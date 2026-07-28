# Production-grade Square-Root Extended Kalman Filter (SR-EKF) v3

Single-file TypeScript library (no runtime dependencies) for fusing IMU + GPS data using a 2D CTRA (Constant Turn Rate and Acceleration) motion model. Optimized for mobile — minimal allocations, fixed-size matrices, Float64Array-backed.

## State Vector (8-state)

```
[x, y, v, ψ, β, a_bias_x, g_bias_z, magDeclination]
```

| Index | Symbol | Description | Unit |
|-------|--------|-------------|------|
| 0 | `x` | East position | m |
| 1 | `y` | North position | m |
| 2 | `v` | Forward speed | m/s |
| 3 | `ψ` | Heading (yaw) | rad |
| 4 | `β` | Sideslip angle (ψ+β = velocity direction) | rad |
| 5 | `a_bias_x` | Accelerometer x-bias (forward) | m/s² |
| 6 | `g_bias_z` | Gyroscope z-bias (yaw rate) | rad/s |
| 7 | `magDeclination` | Magnetic declination (true north − mag north) | rad |

Removed `a_bias_y` — unobservable from forward acceleration alone; removing improved observability of `a_bias_x`.

## Sensor Fusion Architecture

| Rate | Sensor | Mode |
|------|--------|------|
| High (50–400 Hz) | IMU (ax, ay, az, gx, gy, gz) | Predict — bias-corrected readings drive CTRA kinematics |
| Low (1–10 Hz) | GPS (x, y, vx, vy, accuracy) | Update — Mahalanobis-gated correction |
| Medium (10–100 Hz) | Magnetometer (bearing) | Update — heading observation + declination auto-calibration |
| On change | Device orientation (azimuth, pitch, roll) | Frame alignment — rotates IMU from device to ENU |

IMU readings are used as **control inputs** with on-line bias compensation:
- ω = gz − `g_bias_z` (turn rate)
- a_forward = ax − `a_bias_x` (forward acceleration)

## Core Algorithm

### Square-Root EKF (QR-based)

- **Covariance stored as lower-triangular Cholesky factor** `S` such that `P = S·Sᵀ`
- **Prediction**: QR decomposition of `[F·S ; √Q]ᵀ` gives the propagated square root `S⁻`
- **Update**: Pre-array QR `[R_c, 0; H·S⁻, S⁻]` → `Aᵀ = Q·R` → `S⁺ = R[m:m+n, m:m+n]ᵀ`
- **State update**: Standard EKF gain `K = P⁻·Hᵀ·(H·P⁻·Hᵀ + R)⁻¹` — computed via Cholesky solve (chol4x4 + cholSolve4) instead of explicit 4×4 Gaussian elimination inverse
- **Angle normalization**: Heading ψ wrapped to `(−π, π]` after every step
- **QR diagonal sign fix**: After every QR-to-S copy, columns with negative diagonal are negated (`ensurePosDiag`) — QR can produce negative diagonals in R, which becomes invalid Cholesky factors; column negation preserves `P = S·Sᵀ`
- **safeguardState()**: NaN/Inf detection on state and Cholesky factors after prediction; clamps invalid diagonals to 1e-6, sets coasting flag if corrupted

### CTRA Prediction Model

Velocity direction = ψ + β (heading + sideslip)

```
α = ψ + β  (velocity direction)
v̅ = v + ½·a·dt  (average velocity over the step — CTRA integral fix)

if |ω| > ε:
  x' = x + (v̅/ω)·( sin(α + ω·dt) − sin(α) )
  y' = y + (v̅/ω)·( −cos(α + ω·dt) + cos(α) )
else:
  x' = x + v̅·cos(α)·dt
  y' = y + v̅·sin(α)·dt

v'  = v + a_forward·dt
ψ'  = ψ + ω·dt·omegaScale, where omegaScale = (v==0) ? 1 : max(min(|v|/0.5, 1), stillness)
β'  = β · exp(−dt/τ)   (mean-reversion toward 0)
τ = base × (1 − 0.6 × angAccelNorm), where base = |v| < 0.3 ? 0.1 : (|ω| > ε ? 1.5 : max(0.5, 1.5 − (|v|−1.5) / 3.5)), angAccelNorm = min(|dω/dt| / 3.0, 1)
a_bias_x' = a_bias_x · exp(−dt/50)   (mean-reversion when accelEnergy < 0.05)
g_bias_z' = g_bias_z   (random walk, corrected by ZARU during stationary)
magDeclination' = magDeclination   (random walk)
```

- **v̅ = v + ½·a·dt**: Uses the average velocity over the timestep (instead of start-of-step v) in the CTRA position deltas. For straight-line motion this makes the position integration exact for constant acceleration. Jacobian entries for position derivatives (`∂x/∂ψ`, `∂x/∂β`, `∂x/∂g_bias_z`, and y-sym) also use `v̅` to match the forward model; `∂x/∂a_bias_x = -½·dt·∂x/∂v` and `∂y/∂a_bias_x = -½·dt·∂y/∂v` remain correct because `∂x/∂v = ∂x/∂v̅`.
 - **v is no longer hard-clamped to ≥ 0 in `predict`**: Removed `Math.max(0, v)` from predict, lateral accel/ZUPT scalar QR, and ZARU so the filter can express genuine reverse motion. However, the GPS update maintains the invariant **v ≥ 0** with **ψ = direction of motion**: (1) the 180° flip-recovery fires on anti-parallel GPS velocity whenever `v > 0.8` (no position-χ² gate needed), flipping ψ by π before v can go negative; and (2) a post-update backstop reparameterizes any residual `v < 0` into the equivalent `(−v, ψ+π, β+π)` branch (velocity vector unchanged), with the V row/column of the Cholesky factor sign-flipped. This eliminates the reported "U-turn → negative velocity / heading snaps to nose" symptom. All speed-dependent calculations consistently use `Math.abs(v)`: `speedScale`, `maxPlausibleSpeed`, β time constant `|v|<0.3`, lateral accel gate `|v|>0.2`, outlier guard, `speedRamp`, `velFactor`, etc.
- **a_bias_x mean-reversion**: When `accelEnergy < 0.05` (near-zero dynamic acceleration), the bias drifts toward zero with a 50s time constant. This prevents unbounded bias drift during extended constant-speed travel, which would otherwise cause velocity overshoot during GPS outages. Jacobian `F[A_BIAS_X][A_BIAS_X] = exp(−dt/50)` mirrors the state decay.

### GPS Measurement Model

```
z = [x_gps, y_gps, vx_gps, vy_gps]
h(x) = [x, y, v·cos(ψ+β), v·sin(ψ+β)]
```

Jacobian H (4×8):
```
[[1, 0, 0,        0,           0,           0, 0, 0],
 [0, 1, 0,        0,           0,           0, 0, 0],
 [0, 0, cos(ψ+β), -v·sin(ψ+β), -v·sin(ψ+β), 0, 0, 0],
 [0, 0, sin(ψ+β),  v·cos(ψ+β),  v·cos(ψ+β), 0, 0, 0]]
```

### Direction-Aware Outlier Rejection

Decomposes the GPS position innovation into along-track (forward) and cross-track components relative to the velocity direction. Each component has its own plausibility threshold and inflation cap:
- **Forward**: `maxPlausibleSpeed = max(2v, 1) + 2`, threshold = `maxPlausibleSpeed × dtBase × 2`, cap at 5× — trusts forward motion more
- **Cross-track**: threshold = `maxPlausibleSpeed × dtBase × 0.5`, cap at 10× — cross-track jumps (typical GPS multipath) are penalized more aggressively

The guard inflates `posR` multiplicatively, widening the measurement uncertainty rather than rejecting the measurement outright.

### Separate Position/Velocity χ² Gating

The 4-DOF joint Mahalanobis gate (`chiSq < gateThreshold`) is augmented with a 2-DOF position-only χ²: if position alone passes but velocity fails, the update is still accepted. This prevents noisy GPS velocity from blocking valid position corrections at low speeds where velocity error dominates. The position-only χ² uses the same anisotropic noise model (along-track 0.5·posR, cross-track 1.33·posR, heading-rotated) as the main update for gating consistency. **Note**: this separate gating is only active when robust M-estimation is disabled (`!rw.enabled`); when robust weighting is enabled, the re-weighted chiSq is always accepted and the separate gate is bypassed.

### Anisotropic GPS Covariance

Position measurement noise is modeled as a 2×2 covariance rotated by heading:
- Along-track σ = `0.5 × posR` (trust forward motion)
- Cross-track σ = `1.33 × posR` (higher lateral uncertainty)

The Cholesky factor of this rotated covariance is computed inline and placed in the QR pre-array upper-left 2×2 block:

```
R_chol[0][0] = ra (cos-rotated)
R_chol[1][0] = 0
R_chol[0][1] = -sin-rotated * rc
R_chol[1][1] = rc (cos-rotated)
```

### Anisotropic Process Noise

Cross-track process noise is scaled to 0.3× the forward position Q. The 2×2 Cholesky factor is rotated by heading and placed in the prediction QR pre-array, preventing excessive cross-track covariance growth during straight-line motion.

### Nonholonomic Constraint

During straight-line driving (`|ω| < 0.1`, `|v| > 0.15`), a scalar QR pseudo-measurement injects `β ≈ 0` with noise `R = 0.1 rad`. This prevents sideslip drift on long straight segments where β is otherwise unobservable.

### GPS Latency Compensation (Rewind)

When GPS timestamps lag behind IMU time, the filter rewinds state to the GPS timestamp (via `rewindTo()`), processes the GPS update, then replays IMU predictions forward. This avoids temporal misalignment between IMU and GPS.

**Buffer**: Circular buffer of predict states (x, S, IMU inputs, dt, orientation angles). Saved on every `predict()` call, up to 256 entries. Orientation angles (azimuth, pitch, roll) are stored so that `replayFromScratch()` can reconstruct the historical rotation matrix for each replayed step, preventing state corruption if the device orientation changed during the latency window. On rewind, the buffer preserves prior history (`bufTail` unchanged, `bufLen` truncated to the restored index) rather than clearing entirely, enabling subsequent rewind operations for later delayed GPS measurements.

**Critical safeguard**: After rewind, if the position discrepancy between the restored state and the incoming GPS exceeds 10m (>100 m²), the GPS position is applied directly. This handles the case where the buffer only contains pre-GPS-init states, preventing the filter from being stuck far from the true position with an inflated outlier guard.

### Small-ω Jacobian Continuity Fix

When `|ω| ≤ EPS`, the small-ω Jacobian now includes gyro bias derivatives so that position innovation can drive `g_bias_z` correction even during straight-line motion:

```
∂x'/∂g_bias_z = 0.5·v̅·sin(α)·dt²
∂y'/∂g_bias_z = -0.5·v̅·cos(α)·dt²
```

All position derivatives (`∂x/∂ψ`, `∂x/∂β`, `∂x/∂g_bias_z`, and y-sym) use `v̅ = v + ½·a·dt` in both the big-ω and small-ω branches, matching the CTRA position kinematics. Without these terms, `g_bias_z` was only correctable during turns, causing observability gaps during long straight segments.

## Advanced Features

### Adaptive Process Noise Scaling (replaces walking/driving auto-detection)

EMA-tracked IMU energy metrics dynamically scale process noise:
- `accelEnergy` — 0.9 EMA of `max(sqrt(max(varAx+varAy, 0)) / 5, |a_forward| / 5)`, clamped to [0, 5]  
  The `|a_forward|/5` term catches constant braking/deceleration where variance is zero but sustained acceleration is present
- `gyroEnergy` — 0.9 EMA of `sqrt(max(varGz, 0)) / 0.5`, clamped to [0, 5]
- `stepEnergy` — `min(stepFreq / 3.0, 1.0)` from step detection

Scales applied multiplicatively to process noise diagonals:
- Position: `1 + positionAccel × accelEnergy`
- Velocity: `1 + velocityAccel × accelEnergy + velocityStep × stepEnergy`
- Heading: `1 + headingGyro × gyroEnergy + headingStep × stepEnergy + 0.3 × accelEnergy + 1.5 × angAccelBoost` (accelEnergy captures body sway; angAccelBoost inflates Q during corner entry/exit transients)
- Sideslip: `1 + sideslipGyro × gyroEnergy + sideslipStep × stepEnergy + 2.0 × angAccelBoost`
- AccelBias: no scaling (pure random walk)
- GyroBias: `× (1 + 0.3 × gyroEnergy + 0.5 × angAccelBoost)` (gyro energy increases yaw-rate uncertainty; angAccelBoost allows faster bias correction during transients)

Where `angAccelBoost = min(|dω/dt| / 2.0, 1)` — ramps 0→1 for angular acceleration 0→2 rad/s². Zero on straights, full during corner entry/exit.

Configurable via `adaptiveScaling` (defaults tuned for handheld/wearable):
```ts
{ positionAccel: 2.0, velocityAccel: 1.0, velocityStep: 2.5,
  headingGyro: 0.5, headingStep: 0.5, sideslipGyro: 0.5, sideslipStep: 1.8 }
```

### Zero-Velocity Update (ZUPT)
- Stationary detection via exponential `getStillness()` metric derived from 3D IMU variance (1s time-based window, not frame-count based)
- Hysteresis engagement (ZUPT_ON=0.15, ZUPT_OFF=0.05) prevents toggling
- Scalar QR update injects v=0 with 0.01 m/s noise → corrects biases through cross-covariance
- **Smooth engagement**: ZUPT weight = `stillness × speedGate(|v|)`, where `stillness = exp(−(accel3DVar/0.09 + gyro3DEnergy/0.0025))` (1 at rest, → 0 during motion; uses 3D ax+ay+az variance and 3D gx+gy+gz energy) and `speedGate(|v|) = clamp((8.0 − |v|) / 4.0, 0, 1)` (full below 4 m/s, zero above 8 m/s). At rest, stillness ≈ 1 and speedGate ≈ 1, so ZUPT engages whenever the IMU variance is low (constant readings). A genuine (varying) acceleration has high 3D variance → stillness → 0, correctly excluding ZUPT. The speed-gate cutoff is raised to 8 m/s so ZUPT can RECOVER a corrupted/large `v` left over when the phone is placed down (e.g. a stale `v≈6 m/s`). Chi-square innovation gate (threshold 9.0 = 3σ) rejects velocity innovations inconsistent with zero-velocity hypothesis. Continuous on/off transitions; measurement noise R = 0.01 / weight so ZUPT is precise when confident, loose when uncertain.
- **GPS-gated ZUPT**: when GPS is present and reports speed `> 2.0 m/s` (`lastGpsSpeed`), the device is genuinely moving, so ZUPT is suppressed even if `ax` is constant — constant acceleration cruising must integrate normally. At rest GPS speed is ~0 (or GPS absent → `lastGpsSpeed` stays 0), so ZUPT engages and learns the bias. This resolves the IMU-only ambiguity between "constant accel bias at rest" and "genuine constant acceleration".
- **Continuous velR inflation**: Two-stage GPS velocity noise inflation. First, a speed ramp: `velR *= (1 + 2 × speedRamp)` where `speedRamp = max(0, 1 − |v|/10)` (1 at v=0, 0 at v≥10), preventing velocity jumps at moderate city speeds. Then, stationary-weight inflation: `velR *= (1 + min(4, 9 × stationaryWeight))` where `stationaryWeight = max(0, (1.0 − smoothedSpeed) / 1.0)` and `smoothedSpeed` is a 3s EMA of `max(lastGpsSpeed, |v| × 0.3)` (hybrid GPS+EKF speed). At rest (smoothedSpeed≈0), velR is inflated up to 5×; at smoothedSpeed≥1.0, no inflation. A further ×4 inflation is applied when GPS Doppler opposes the IMU heading direction at low speed (`|v| < 2.5` and dot-product < 0), distrusts multipath-corrupted velocity directions in urban canyons
- **No velocity clamping**: `applyZupt()` does NOT force `v=0` — the Kalman update pulls v toward 0 naturally via `rVel = 0.01/weight`. This allows IMU acceleration to gradually build v during ZUPT, preventing the "stuck at red light" problem.
- **Disengagement uses hysteresis**: ZUPT disengages when `gpsMoving` (GPS speed > 2.0 m/s) or `zuptWeight < ZUPT_OFF` (0.05). Since `zuptWeight = stillness × speedGate`, disengagement occurs when stillness drops (motion begins) or speed rises above the gate. On the transition, V/PSI Cholesky rows are inflated (V=1.5, X/Y=3.0, PSI=0.3) to prevent first GPS after stop from overshooting.
- **`getStillness()` is variance-based ONLY**: `stillness = exp(−(accel3DVar / T_a² + gyro3DEnergy / T_w²))` where `T_a = 0.3` and `T_w = 0.05`, with `accel3DVar` = 3D (ax+ay+az) window variance and `gyro3DEnergy` = 3D (gx+gy+gz) mean energy. Returns 0.5 (neutral) when < 5 samples. A stationary device has near-zero 3D variance → stillness ≈ 1; motion drives variance up → stillness → 0. The magnitude term (`|a|/5`) is still used for adaptive process-noise scaling and `accelEnergy`-based ZUPT GPS-gating, but NOT for stillness. Note ZARU is **no longer** gated on `getStillness()` (see ZARU bullet) — it uses the filter's own `|ω|`.
- **Heading covariance floor**: After ZUPT QR update, `S[PSI][PSI]` is clamped to ≥ 0.05 rad — prevents heading from locking up during long stops
  - **ZARU (Zero Angular Rate Update)**: A scalar QR update injects `omega = 0` with noise `R = 0.01 / zaruGate` rad/s, making gyro bias (`g_bias_z`) observable when the device is not rotating — the filter drives `gBiasZ → gz` when the true angular rate is known to be zero, preventing gyro bias drift through stop-start cycles. ZARU is a **standalone method called directly from `predict()`**, gated on the filter's own bias-corrected rate `zaruGate = 1 − min(|ω|/0.2, 1)`, firing when `zaruGate > 0.5` (i.e. `|ω| < 0.1 rad/s`) — NOT on the raw-IMU `getStillness()` and NOT nested inside ZUPT. Being standalone means ZARU runs even when GPS confirms motion and ZUPT is disabled, so gyro bias is always corrected when the device isn't rotating. Rationale: on a real phone the gyro has high *variance* even at rest, so the variance-based `stillness` saturates to 0 → an IMU-stillness-gated ZARU never fires → gyro bias stays uncorrected → ψ integrates the residual rate → the **"large circle" drift** at rest. Gating on the bias-corrected rate `|ω|` (the exact quantity ψ integrates) is noise-amplitude-independent: when the device is truly still, `|ω|` is small whenever `gBiasZ` is even roughly right, and the update then pulls `gBiasZ` the rest of the way to `gz`. A lenient 0.1 rad/s threshold lets it bootstrap from a wrong initial bias (e.g. `gBiasZ=0`, true bias 0.0583 → `|ω|≈0.0583 < 0.1` → ZARU fires → bias learned).

### Magnetometer Integration
- `updateMag(bearing, timestampMs)` applies heading observation via scalar QR with observation model `h(x) = ψ − magDeclination`, so `H[PSI]=1, H[MAG_DECL]=-1`
- Mag declination is the 8th EKF state, calibrated via cross-covariance with heading. When GPS velocity direction and mag bearing disagree, the filter distributes the correction between ψ and magDeclination based on their relative uncertainties.
- Initial covariance default: `0.25 rad²` (~30° std), configurable via `initialCovariance.magDeclination`
- Process noise: random walk with `processNoise.magDeclination` (default `1e-4`, was `1e-6` — raised so declination is actually observable and absorbs bearing offsets via cross-covariance instead of freezing and yanking ψ)
- Heading auto-initialized from compass on first `updateMag()` call (before GPS init)
- Configurable `measurementNoise.heading` (default 0.1 rad)
- **Adaptive measurement noise**: For large innovations, `r = rBase × stillnessFactor × (1 + min(|innov|/0.175, 1))` where `stillnessFactor = max(0.5, 1 − stillness × 0.5)` — at rest (`stillness→1`), `r` drops to 0.5× the baseline so heading converges to the compass faster. A 10° innovation doubles `r` (reducing Kalman gain for smooth convergence); 30° caps at 4× baseline. Prevents heading jumps from abrupt bearing corrections

### Magnetometer ↔ GPS Heading Authority (speed-gated)
- The magnetometer has **no β term** in its observation model (`H[PSI]=1, H[MAG_DECL]=-1` only), so it observes the *body/nose* heading `ψ = bearing + magDeclination`. GPS observes the *velocity* direction `ψ+β` via `H[2:3][PSI]=H[2:3][BETA]`. When sideslip β ≠ 0 the two references disagree; if both updated ψ at speed they fight and ψ oscillates (the reported red-light / >10 km/h oscillation bug).
- `magUpdateSingle` computes `magHeadingTrust = max(0, 1 − |v|/1.5)` (1 at rest, 0 at |v| ≥ 1.5 m/s ≈ 5.4 km/h). When `magHeadingTrust < 0.2` the **entire mag update is skipped** (including the init snap), so at speed GPS owns live heading ψ and mag cannot oscillate it. The `gpsStale` term is used only to tighten the rejection gate, not to decide authority.
- Policy: **GPS owns ψ at speed; mag owns ψ at rest/low-speed and calibrates `magDeclination`**. Because mag is skipped above 5.4 km/h, `magDeclination` is calibrated only at low speed — it is a slow random walk, so this is sufficient. ψ correctly represents body heading at speed; the *moving direction* is `ψ+β` (read both states).
- The init snap (`ψ = wrapAngle(bearing + magDeclination)` when heading/declination covariance is still near its initial value) is likewise gated by `magHeadingTrust < 0.2` so it never rotates ψ at speed.
- **Rest-rate cross-check (drifting-compass guard)**: `magUpdateSingle` runs an angular-rate cross-check at the very top (before the init snap / trust gate). It compares the compass's angular RATE between consecutive `updateMag()` calls against the bias-corrected gyro rate (`lastOmega`). At rest the gyro rate ≈ 0, so any compass rotation above a small margin is drift — uncalibrated mag, hard-iron, or table interference — NOT real motion. The margin is speed-dependent: at `|v| < 0.1` m/s the gyro yaw rate may be unavailable on some Android devices, so the margin widens to `gyroRate + 2.0` rad/s; at speed the standard tight margin is `gyroRate + 0.05` rad/s. Genuine turns pass because there the gyro rate is high and the compass rotation matches it (`magRate ≈ gyroRate`, not `> gyroRate + margin`). The check requires a valid prior call and a plausible inter-call `dtCall` (1ms–60s).

### Sideslip Angle (β)
- 6th state captures the difference between vehicle heading (ψ) and velocity direction (ψ+β)
- Arises from tire slip during turns; GPS velocity direction vs gyro-integrated heading gives observability
- Mean-reversion toward 0 with speed/yaw-rate-dependent time constant: base `τ = |v| < 0.3 ? 0.1 : (|ω| > ε ? 1.5 : max(0.5, 1.5 − (|v|−1.5) / 3.5))` s, then multiplied by `(1 − 0.6 × angAccelNorm)` where `angAccelNorm = min(|dω/dt| / 3.0, 1)` — at max angular acceleration (3 rad/s²), τ drops to 40% of its base value (e.g. 1.5s → 0.6s during turns). This allows β to track the rapidly changing slip angle at corner entry/exit instead of lagging and corrupting heading. Jacobian `computeJacobian()` uses the same modulated τ.
- Jacobian ∂/∂β = ∂/∂ψ (ψ and β appear symmetrically in position/velocity kinematics)
- Process noise scaled by IMU energy and angular acceleration (via adaptive scaling)

### Angular Acceleration Tracking (corner entry/exit transient detection)

Tracks `|dω/dt|` (angular acceleration magnitude) by differencing consecutive bias-corrected gyro rates. Clamped to [0, 10] rad/s². Drives three mechanisms that prevent heading from "getting stuck" at corner transitions:

1. **Q boost** (`angAccelBoost = min(|dω/dt|/2, 1)`): Inflates heading Q by ×(1+1.5·boost), sideslip Q by ×(1+2·boost), and gyro bias Q by ×(1+0.5·boost) during transients. Zero on straights, full at 2 rad/s². Allows faster state correction when the CTRA model is least accurate.
2. **Faster β adaptation**: Modulates sideslip time constant by `(1 − 0.6·angAccelNorm)` where `angAccelNorm = min(|dω/dt|/3, 1)`. At max angular acceleration, β converges 2.5× faster (τ=0.6s vs 1.5s), preventing sideslip lag from corrupting heading during corner entry.
3. **Gyro bias Q boost**: Extra process noise on gyro bias during transients allows GPS velocity direction to correct accumulated bias errors that would otherwise cause heading overshoot at corner exit.

### Speed-Scaled Process Noise
- Position and velocity Q-diagonal entries auto-scale with speed:
  `speedScale = min(sqrt(max(|v|, 0.05) / 5.0), 2)`
- At 5 m/s: scale = 1.0 (baseline). Capped at 2.0 for |v| ≥ 20 m/s to prevent excessive covariance growth at highway speeds.
- Heading, bias, and sideslip noise are not speed-scaled

### Frame Alignment
- `setOrientation(azimuth, pitch, roll)` sets device-to-ENU rotation matrix using Z-X-Y Euler sequence (matches W3C DeviceOrientation convention): pitch around X-axis (front-to-back), roll around Y-axis (left-to-right), azimuth clockwise-positive around Z (compass heading)
- When orientation is set, a **6-axis IMU (including az) is required** — omitting `az` causes NaN propagation as a guard against gravity leakage
- Bias reset (a_bias_x, g_bias_z → 0 with fresh covariance) only triggers when the device orientation changes **relative to the vehicle** (>5° over 1s in device-to-vehicle frame, factoring out heading changes). This prevents vehicle turns from destroying learned bias calibration
- Orientation angles are stored alongside raw IMU in the latency compensation buffer, so rewind/replay correctly reconstructs the historical rotation matrix for each buffered step

### Low-Speed Jump Protection (Urban Multipath)

Two defenses prevent position/heading jumps from noisy GPS in low-speed city environments:

1. **Tighter position outlier guard**: Uses `maxPlausibleSpeed = max(2·|v|, 1) + 2` (instead of `max(v, 5) + 20`) with `×2` multiplier (instead of `×5`), capping inflation at 5× forward and 10× cross-track. At walking speed and 1 Hz GPS, a 5m jump inflates posR by 1.7×; a 20m jump inflates by 8.3×
2. **Smooth heading gating**: GPS velocity heading correction ramps from 0 at v=0.5 to full via `headingGain = v < 0.5 ? 0 : min((v − 0.5 + gyroEnergy × 0.5) / 3.5, 1)`, applied to `H[2:3, PSI/BETA]`. The gyroEnergy term allows faster convergence when the device is actively turning (gyro confirms real motion). During the first 30 seconds after GPS initialization, an adaptive init boost further increases headingGain by up to 2× when heading uncertainty is high (`psiStd > 0.3 rad`), accelerating convergence from a cold start. Prevents hard-threshold discontinuity while still avoiding heading corruption from noisy GPS direction at near-zero speed. Position innovation still corrects heading through position-heading cross-covariance at all speeds
3. **GPS heading init gate lowered** (line 357): `spd > 0.1` — catches even very slow movement for heading initialization from GPS velocity direction
 4. **180° flip recovery gate relaxed** (src/sr-ekf.ts `gpsUpdateSingle`): now fires on anti-parallel GPS velocity (`dot < -0.5·v²`) whenever `v > 0.8`, with **no position-χ² gate requirement** (previously required `chiSq > gateThreshold`). The strict anti-parallel dot-product check still prevents misfire on 90° pedestrian turns (which only trigger at >120° separation). This catches pure velocity-direction reversals — e.g. urban Doppler multipath or a U-turn where position stays consistent but velocity flips — which previously were not corrected and resolved into a negative-v state. A post-update backstop then reparameterizes any residual `v < 0` into the `(−v, ψ+π, β+π)` branch so `v ≥ 0` always holds with `ψ` = direction of motion. The backstop only reparameterizes when `|v| ≥ 0.5 m/s`; at negligible speed (stationary device) it merely clamps `v` to 0 without flipping ψ. Otherwise, GPS velocity *noise* at rest drives `v` to tiny negative values and the π-flip would make the heading shake 180° every step.
5. **Rest-weighted velocity H blocking**: When IMU energy is near zero (`accelEnergy + gyroEnergy < 0.05`), the GPS velocity heading columns (`H[2:3][PSI/BETA]`) are smoothly reduced by a speed-dependent factor `(1 − restW)` where `restW = clamp((1.5 − |v|) / 1.0, 0, 1)`. At rest this completely blocks GPS velocity from rotating ψ, handing authority to the magnetometer. At speed the blocking fades, allowing GPS velocity direction to own heading.

### Lateral Acceleration Constraint
- During turns (`|ω| > 0.1`, `|v| > 0.2`, no device-to-ENU orientation), uses `ay ≈ v·ω` as a pseudo-measurement to estimate sideslip from centripetal acceleration mismatch — higher omega gate prevents body sway (0.3–0.8 m/s²) from corrupting heading during straight-line walking
- Enabled by default (`useLateralAccel: true`) — safe to enable because it's gated on `deviceToEnu === null` (IMU in vehicle frame)
- Measurement noise `r = 1.0 m/s²`

### Robust M-Estimation (Cauchy/Huber)
- Replaces binary Mahalanobis gate with a smooth weight function when enabled
- **Cauchy**: `w = 1 / (1 + (χ²/threshold)²)` — heavy-tailed
- **Huber**: `w = threshold / max(χ², threshold)` — linear penalty
- Config: `robustWeight.enabled`, `.type` ('cauchy'|'huber'), `.threshold` (default 9.488)
- **Per-component velocity robustness** (noisy city GPS heading/speed): GPS heading is inferred from Doppler velocity via `H[2:3][ψ,β]`. In urban canyons multipath corrupts velocity while position stays usable. The joint 4-DOF weight above conflates the two — a bad velocity would down-weight the good position fix, and moderate velocity noise below the joint gate passes through at full strength and jitters heading. A **velocity-only** weight `velRobustW` is derived from the 2-DOF velocity χ² (`chiSqVel`, threshold scaled 4-DOF→2-DOF by `0.63 = χ²₀.₉₅,₂/χ²₀.₉₅,₄`) and applied to `velR` alone (`velR /= √(totalWeight·velRobustW)`), so noisy Doppler is distrusted for heading without sacrificing position. No-op on clean GPS (`velRobustW = 1`); `robustWeight` diagnostic still reports the joint weight. Reduces heading RMS jitter ~40% under ±1.8 m/s cross-velocity multipath with position tracking unchanged.

### Adaptive Noise Scaling
- EMA tracks `innov² / S_innov` ratio across GPS updates
- When ratio > 1, inflates R by adaptive scale factor
- Config: `adaptiveNoise.enabled`, `.smoothing` (0–1, default 0.1), `.maxScale` (default 3)

### GPS Initialization

On the first `updateGps()` call after `reset()`, the filter snaps position directly to the GPS fix (no Kalman update). Heading and velocity are initialized from GPS velocity direction when `|v| > 0.1` m/s (low threshold catches slow-walking pedestrians). **Covariance is tightened** to reflect GPS-derived knowledge — position σ set to `min(current, accuracyMeters/3)` (default 3m if no accuracy reported), velocity σ capped at 0.5 m/s, heading σ capped at 0.3 rad (~17°) when velocity-initialized. Without this, the default `initialCovariance.position=100` (σ=10m) persists after init, letting IMU drift accumulate rapidly between GPS fixes and causing the direction-aware guard to reject valid subsequent fixes.

### Coast Recovery
- When coasting and GPS gate rejects, `resetFromGps()` resets state + Cholesky from GPS fix directly; position σ capped at 5m via `min(current, 5)`
- **Biases preserved across GPS re-acquisition**: `resetFromGps()` keeps the learned `aBiasX`, `gBiasZ`, and `magDeclination` values and their covariances — only kinematic states (x, y, v, ψ, β) are reset. This prevents tunnel exit from discarding the bias calibration learned during coasting via ZUPT/ZARU, eliminating the "fresh drift" symptom after GPS recovery.
- Deduplicated from 3 copies into single private method
- **Guard-inflation-masking reset**: If the outlier guard inflates `posR` by >3× (`posR / preGuardPosR > 3`) and the position jump is physically plausible (`dxNorm < dtSinceLastGps × 50 + 2`) and the device is genuinely moving (`|v| > 2.0 m/s`), the GPS fix is accepted directly via `resetFromGps()` — the inflated R would otherwise reduce Kalman gain to near-zero, causing minute-long convergence. At stationary, GPS position noise routinely exceeds the guard thresholds and resetFromGps would corrupt heading by setting psi = atan2(noisy v).
- **Auto-divergence detection**: When NOT coasting and the previous GPS gate failed, if position error exceeds 10m (100 m²), the filter force-enters coasting to trigger a full reset on the next GPS fix — catches cases where a brief glitchy GPS kept lastGpsTimeMs alive but the state has clearly diverged (e.g. basement exit).
- **Post-update large position reset**: After a successful GPS update, if the resulting position error exceeds 10km² (1e8 m²), a hard reset from GPS is applied — catches catastrophic Cholesky corruption

### Tunnel Navigation (GPS-Denied Dead Reckoning)

When GPS is lost (tunnel, parking garage, urban canyon), the filter transitions to pure IMU prediction with these aids:

- **ZARU** continues firing during coasting — corrects gyro bias when the car is stopped or driving straight (`|ω| < 0.1`), maintaining heading accuracy through long tunnels
- **Nonholonomic constraint** constrains β→0 during straight driving; **lateral acceleration constraint** constrains velocity during turns — both fire inside `predict()` regardless of GPS state
- **Magnetometer** provides heading corrections at low speed (`|v| < 1.5` m/s) when stopped in a tunnel — `updateMag()` is independent of GPS
- **Coasting Q reduction**: During coasting, position/velocity process noise is scaled to 30% (`COAST_Q_FACTOR = 0.3`) — without GPS corrections, bias-driven position error is already captured in P via cross-covariance; adding full Q over-inflates uncertainty. Heading/bias Q remain at full strength since ZARU/Mag still provide corrections
- **Coast covariance cap**: Position σ capped at 500m, velocity σ at 50 m/s during coasting — prevents unbounded growth in very long tunnels (2+ minutes)
- **Stale GPS speed decay**: After 10s of coasting, `lastGpsSpeed` decays toward 0 with 5s time constant — enables ZUPT to engage if the car stops in a tunnel (previously blocked permanently when GPS was lost at >2 m/s)

## Configuration

```ts
interface EkfConfig {
  dt?: number
  processNoise?: {
    position?: number           // σ per √s (default 5.0)
    velocity?: number           // σ per √s (default 0.5)
    heading?: number            // σ per √s (default 0.10)
    sideslip?: number           // σ per √s (default 0.1)
    accelBias?: number          // σ per √s (default 1e-4)
    gyroBias?: number           // σ per √s (default 5e-5)
    magDeclination?: number     // σ per √s (default 1e-4)
  }
  measurementNoise?: {
    position?: number           // σ for GPS position (default 3.0)
    velocity?: number           // σ for GPS velocity (default 0.5)
    heading?: number            // σ for magnetometer heading (default 0.1)
  }
  magneticDeclination?: number  // initial mag-to-true north offset (default 0)
  initialCovariance?: {
    position?: number           // σ² (default 100)
    velocity?: number           // σ² (default 10)
    heading?: number            // σ² (default π²)
    sideslip?: number           // σ² (default 0.25)
    accelBias?: number          // σ² (default 0.1)
    gyroBias?: number           // σ² (default 0.01)
    magDeclination?: number     // σ² (default 0.25)
  }
  gateThreshold?: number        // Mahalanobis chi-square threshold (default 9.488)
  coastTimeoutMs?: number       // GPS timeout before divergence flag (default 5000)
  gpsTimeOffsetMs?: number      // GNSS→local clock offset (default 0)
  useLateralAccel?: boolean     // enable centripetal sideslip constraint (default true)
  robustWeight?: {
    enabled?: boolean
    type?: 'cauchy' | 'huber'
    threshold?: number
  }
  adaptiveNoise?: {
    enabled?: boolean
    smoothing?: number
    maxScale?: number
  }
  adaptiveScaling?: {
    positionAccel?: number
    velocityAccel?: number
    velocityStep?: number
    headingGyro?: number
    headingStep?: number
    sideslipGyro?: number
    sideslipStep?: number
  }
}
```

## API Surface

```ts
class SrEkf {
  constructor(config?: EkfConfig)

  predict(ax: number, ay: number, gz: number, dt: number,
          timestampMs: number, az?: number, gx?: number, gy?: number): void

  updateGps(x: number, y: number, vx: number, vy: number,
            timestampMs: number, accuracyMeters?: number): boolean

  updateMag(bearing: number, timestampMs: number): void

  setOrientation(azimuth: number, pitch: number, roll: number): void

  coast(timeoutMs: number, currentTimeMs: number): boolean

  getState(): NavigationSolution

  getDiagnostics(): EkfDiagnostics

  reset(x: number, y: number, v: number, psi: number): void
}

interface NavigationSolution {
  x: number; y: number;
  v: number; psi: number;
  beta: number;
  aBiasX: number;
  gBiasZ: number;
  magDeclination: number;
  p: Float64Array[];             // full covariance (8×8)
}

interface EkfDiagnostics {
  trace: number;
  gpsInnovation: number[];
  gpsChiSq: number;
  gatePassed: boolean;
  coasting: boolean;
  lastGpsTimeMs: number;
  lastImuTimeMs: number;
  stationary: boolean;   // true iff stillness > 0.7 AND |v| < 3.0 (matches ZUPT engagement)
  magDeclination: number;
  robustWeight: number;
  adaNoiseScale: number;
}
```

## Determinism

Fully deterministic — no `Date.now()` calls. All timestamps are explicit parameters.

## Mobile Performance

Designed for 50–400 Hz IMU and 1–10 Hz GPS on resource-constrained devices:

- **Zero allocations on the predict() hot path** — all working matrices (`tmpF`, `tmpQR`, `tmpPreA`, etc.) are preallocated `Float64Array`s created once in the constructor. `predict()`/`applyZupt()`/`magUpdateSingle()` never `new` anything.
- **Zero-allocation state read**: use `getStateInto(out)` with a caller-owned `NavigationSolution` buffer on the render loop instead of `getState()` (which allocates a fresh 8×8 covariance matrix each call). `getDiagnostics()` is allocation-free (reuses an internal `Float64Array(4)` for `gpsInnovation`).
- **Latency compensation buffer** is a flat `Float64Array` (stride `N` for state, `N·(N+1)/2` for the lower-triangular Cholesky factor) — no object-per-entry garbage. On rewind it preserves **all 8 states** including `magDeclination`.
- **Fixed-size matrix math**: N=8 / M=4. The prediction QR is 16×8, GPS update QR is 12×12, scalar updates 9×9 — all O(N³) with no dynamic sizing.
- **In-place Householder QR** (`qrInPlace`) operates directly on the pre-array, avoiding copies.

## Validation

```bash
npm test                 # 70 tests (9 QR verification + 61 unit)
npm run test:watch       # watch mode
```
