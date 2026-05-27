# Production-grade Square-Root Extended Kalman Filter (SR-EKF) v2

Single-file TypeScript library (no runtime dependencies) for fusing IMU + GPS data using a 2D CTRA (Constant Turn Rate and Acceleration) motion model. Optimized for mobile — minimal allocations, fixed-size matrices, Float64Array-backed.

## State Vector (8-state)

```
[x, y, v, ψ, β, a_bias_x, a_bias_y, g_bias_z]
```

| Index | Symbol | Description | Unit |
|-------|--------|-------------|------|
| 0 | `x` | East position | m |
| 1 | `y` | North position | m |
| 2 | `v` | Forward speed | m/s |
| 3 | `ψ` | Heading (yaw) | rad |
| 4 | `β` | Sideslip angle (ψ+β = velocity direction) | rad |
| 5 | `a_bias_x` | Accelerometer x-bias (forward) | m/s² |
| 6 | `a_bias_y` | Accelerometer y-bias (lateral) | m/s² |
| 7 | `g_bias_z` | Gyroscope z-bias (yaw rate) | rad/s |

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
- **State update**: Standard EKF gain `K = P⁻·Hᵀ·(H·P⁻·Hᵀ + R)⁻¹`
- **Angle normalization**: Heading ψ wrapped to `(−π, π]` after every step

### CTRA Prediction Model

Velocity direction = ψ + β (heading + sideslip)

```
α = ψ + β  (velocity direction)

if |ω| > ε:
  x' = x + (v/ω)·( sin(α + ω·dt) − sin(α) )
  y' = y + (v/ω)·( −cos(α + ω·dt) + cos(α) )
else:
  x' = x + v·cos(α)·dt
  y' = y + v·sin(α)·dt

v'  = v + a_forward·dt
ψ'  = ψ + ω·dt
β'  = β · 0.98   (mean-reversion toward 0, τ ≈ 0.5 s)
a_bias_x' = a_bias_x   (random walk)
a_bias_y' = a_bias_y   (random walk)
g_bias_z' = g_bias_z   (random walk)
```

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

## Advanced Features

### Walking / Driving Auto-Detection
- Configurable `mode: 'walk' | 'drive' | 'auto'` (default `'auto'`)
- Heuristic detection: step frequency (ax zero-crossings), speed threshold, gyro energy
- Walking uses higher process noise (less predictable foot motion)
- `auto` blends both process noise models via likelihood

### Zero-Velocity Update (ZUPT)
- Stationary detection via rolling IMU variance (adaptive noise floor, no tuning)
- 10-frame debounce prevents toggling
- Scalar QR update injects v=0 with 1 cm/s noise → corrects biases through cross-covariance
- On stationary→moving transition: V/PSI covariance inflated 5× to prevent GPS gate rejection

### Magnetometer Integration
- `updateMag(bearing, timestampMs)` applies heading observation via scalar QR
- Auto-calibrates `magDeclination` online: compares GPS velocity heading against mag bearing when speed > 3 m/s
- Heading auto-initialized from compass on first `updateMag()` call (before GPS init)
- Configurable `measurementNoise.heading` (default 0.1 rad)

### Sideslip Angle (β)
- 8th state captures the difference between vehicle heading (ψ) and velocity direction (ψ+β)
- Arises from tire slip during turns; GPS velocity direction vs gyro-integrated heading gives observability
- Mean-reversion toward 0 (decay 0.98/step, τ ≈ 0.5 s) — returns to near-zero during straight driving
- Jacobian ∂/∂β = ∂/∂ψ (ψ and β appear symmetrically in position/velocity kinematics)
- `NavigationSolution.beta` exposes the estimate
- Process noise: 0.1 rad/√s driving, 0.3 rad/√s walking (via blend)

### Speed-Scaled Process Noise
- Position and velocity Q-diagonal entries auto-scale with speed:
  `speedScale = sqrt(max(v, 0.5) / 5.0)`
- At 5 m/s: scale = 1.0 (baseline). At walking (1.4 m/s): scale ≈ 0.53. At highway (33 m/s): scale ≈ 2.56
- Produces tighter filtering at low speeds (dead-reckoning is more accurate), looser at high speed (IMU errors → larger position uncertainty)
- Heading, bias, and sideslip noise are not speed-scaled (they depend on vehicle dynamics, not speed)

### Stop-and-Go Position Hold
- When ZUPT triggers, a 3-measurement QR update pins v=0 AND (x,y) at current estimate:
  - v: R=0.01 m/s (existing velocity ZUPT)
  - x: R=0.01 m (new — position hold)
  - y: R=0.01 m (new — position hold)
- Prevents GPS position drift from pulling the filter during stops (e.g., at traffic lights)
- Single fused QR decomposition avoids redundant Householder reflections

### Frame Alignment
- `setOrientation(azimuth, pitch, roll)` sets device-to-ENU rotation matrix
- When orientation is set and 6-axis IMU provided, readings rotate before bias compensation
- Call from `RotationVector` or `getRotationMatrix()` sensor each frame

## Configuration

```ts
interface EkfConfig {
  dt?: number                     // default prediction dt (s)
  mode?: 'walk' | 'drive' | 'auto'
  processNoise?: {
    position?: number             // σ per √s for position (default 1.0)
    velocity?: number             // σ per √s for velocity (default 0.5)
    heading?: number              // σ per √s for heading (default 0.05)
    sideslip?: number             // σ per √s for sideslip (default 0.1)
    accelBias?: number            // σ per √s for accel bias walk (default 1e-4)
    gyroBias?: number             // σ per √s for gyro bias walk (default 1e-5)
  }
  walkingProcessNoise?: {         // same fields, higher defaults
    position?: number             // (default 2.0)
    velocity?: number             // (default 2.0)
    heading?: number              // (default 0.3)
    sideslip?: number             // (default 0.3)
    accelBias?: number            // (default 1e-3)
    gyroBias?: number             // (default 1e-4)
  }
  measurementNoise?: {
    position?: number             // σ for GPS position (default 3.0)
    velocity?: number             // σ for GPS velocity (default 0.5)
    heading?: number              // σ for magnetometer heading (default 0.1 rad)
  }
  magneticDeclination?: number    // initial mag-to-true north offset (default 0)
  initialCovariance?: {
    position?: number             // σ² initial position uncertainty (default 100)
    velocity?: number             // σ² initial velocity uncertainty (default 10)
    heading?: number              // σ² initial heading uncertainty (default π²)
    sideslip?: number             // σ² initial sideslip uncertainty (default 0.25)
    accelBias?: number            // σ² initial accel bias uncertainty (default 0.1)
    gyroBias?: number             // σ² initial gyro bias uncertainty (default 0.01)
  }
  gateThreshold?: number          // Mahalanobis chi-square threshold (default 9.488)
  coastTimeoutMs?: number         // max IMU-only time before flagging divergence (default 5000)
  gpsTimeOffsetMs?: number        // GNSS→local clock offset (default 0)
}
```

## API Surface

```ts
class SrEkf {
  constructor(config?: EkfConfig)

  // High-rate: propagate state using IMU readings
  predict(ax: number, ay: number, gz: number, dt: number,
          timestampMs: number, az?: number, gx?: number, gy?: number): void

  // Low-rate: correct state with GPS fix
  updateGps(x: number, y: number, vx: number, vy: number,
            timestampMs: number, accuracyMeters?: number): boolean

  // Magnetometer heading observation
  updateMag(bearing: number, timestampMs: number): void

  // Set device orientation for IMU frame alignment
  setOrientation(azimuth: number, pitch: number, roll: number): void

  // Coast detection (returns false if covariance has diverged)
  coast(timeoutMs: number, currentTimeMs: number): boolean

  // Retrieve current navigation solution
  getState(): NavigationSolution

  // Retrieve filter diagnostics
  getDiagnostics(): EkfDiagnostics

  // Reset to initial conditions
  reset(x: number, y: number, v: number, psi: number): void
}

interface NavigationSolution {
  x: number; y: number;          // position (m)
  v: number;                     // speed (m/s)
  psi: number;                   // heading (rad)
  beta: number;                  // sideslip angle (rad)
  aBiasX: number; aBiasY: number;// accel biases
  gBiasZ: number;                // gyro bias
  p: Float64Array[];             // full covariance (8×8)
  mode: 'walk' | 'drive';       // active mode
}

interface EkfDiagnostics {
  trace: number;                 // covariance trace
  gpsInnovation: number[];       // 4-element innovation vector (last GPS update)
  gpsChiSq: number;              // chi-square statistic (last GPS update)
  gatePassed: boolean;           // last GPS update passed gating
  coasting: boolean;             // currently in coast mode
  lastGpsTimeMs: number;         // timestamp of last successful GPS update
  lastImuTimeMs: number;         // timestamp of last IMU prediction
  mode: 'walk' | 'drive';       // actual active mode
  walkLikelihood: number;        // 0–1 likelihood of walking
  stationary: boolean;           // ZUPT stationary detection flag
  magDeclination: number;        // current magnetic declination estimate
}
```

## Determinism

The filter is fully deterministic — no `Date.now()` calls. All timestamps are explicit parameters:
- `predict(timestampMs)` — IMU system time
- `updateGps(timestampMs)` — GPS fix time
- `coast(currentTimeMs)` — wall-clock time for timeout comparison

For log replay: pass the same timestamps in sequence and get identical outputs.

## Validation

```bash
npm test                 # 17 unit tests
npm run test:watch       # watch mode
```
