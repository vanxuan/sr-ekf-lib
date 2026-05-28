# Production-grade Square-Root Extended Kalman Filter (SR-EKF) v3

Single-file TypeScript library (no runtime dependencies) for fusing IMU + GPS data using a 2D CTRA (Constant Turn Rate and Acceleration) motion model. Optimized for mobile — minimal allocations, fixed-size matrices, Float64Array-backed.

## State Vector (7-state)

```
[x, y, v, ψ, β, a_bias_x, g_bias_z]
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
g_bias_z' = g_bias_z   (random walk)
```

### GPS Measurement Model

```
z = [x_gps, y_gps, vx_gps, vy_gps]
h(x) = [x, y, v·cos(ψ+β), v·sin(ψ+β)]
```

Jacobian H (4×7):
```
[[1, 0, 0,        0,           0,           0, 0],
 [0, 1, 0,        0,           0,           0, 0],
 [0, 0, cos(ψ+β), -v·sin(ψ+β), -v·sin(ψ+β), 0, 0],
 [0, 0, sin(ψ+β),  v·cos(ψ+β),  v·cos(ψ+β), 0, 0]]
```

### Small-ω Jacobian Continuity Fix

When `|ω| ≤ EPS`, the small-ω Jacobian now includes gyro bias derivatives so that position innovation can drive `g_bias_z` correction even during straight-line motion:

```
∂x'/∂g_bias_z = 0.5·v·sin(α)·dt²
∂y'/∂g_bias_z = -0.5·v·cos(α)·dt²
```

Without these terms, `g_bias_z` was only correctable during turns, causing observability gaps during long straight segments.

## Advanced Features

### Adaptive Process Noise Scaling (replaces walking/driving auto-detection)

EMA-tracked IMU energy metrics dynamically scale process noise:
- `accelEnergy` — 0.9 EMA of `sqrt(max(varAx+varAy, 0)) / 5.0`, clamped to [0, 5]
- `gyroEnergy` — 0.9 EMA of `sqrt(max(varGz, 0)) / 0.5`, clamped to [0, 5]
- `stepEnergy` — `min(stepFreq / 3.0, 1.0)` from step detection

Scales applied multiplicatively to process noise diagonals:
- Position: `1 + positionAccel × accelEnergy`
- Velocity: `1 + velocityAccel × accelEnergy + velocityStep × stepEnergy`
- Heading: `1 + headingGyro × gyroEnergy + headingStep × stepEnergy`
- Sideslip: `1 + sideslipGyro × gyroEnergy + sideslipStep × stepEnergy`
- Biases: no scaling (pure random walk)

Configurable via `adaptiveScaling` (defaults tuned for handheld/wearable):
```ts
{ positionAccel: 2.0, velocityAccel: 1.0, velocityStep: 2.5,
  headingGyro: 0.5, headingStep: 4.8, sideslipGyro: 0.5, sideslipStep: 1.8 }
```

### Zero-Velocity Update (ZUPT)
- Stationary detection via rolling IMU variance (10-frame window, fixed thresholds — no adaptive noise floor tuning needed)
- 10-frame debounce prevents toggling
- Scalar QR update injects v=0 with 0.01 m/s noise → corrects biases through cross-covariance
- On stationary→moving transition: V/PSI Cholesky rows set to fixed values (V=5.0, X/Y=10.0, PSI×3) to prevent GPS gate rejection
- Position hold (x,y) constraints removed in favor of clean reset on transition

### Magnetometer Integration
- `updateMag(bearing, timestampMs)` applies heading observation via scalar QR
- Auto-calibrates `magDeclination` online: compares GPS velocity heading against mag bearing when speed > 3 m/s
- Heading auto-initialized from compass on first `updateMag()` call (before GPS init)
- Configurable `measurementNoise.heading` (default 0.1 rad)

### Sideslip Angle (β)
- 6th state captures the difference between vehicle heading (ψ) and velocity direction (ψ+β)
- Arises from tire slip during turns; GPS velocity direction vs gyro-integrated heading gives observability
- Mean-reversion toward 0 (decay 0.98/step, τ ≈ 0.5 s)
- Jacobian ∂/∂β = ∂/∂ψ (ψ and β appear symmetrically in position/velocity kinematics)
- Process noise scaled by IMU energy (via adaptive scaling)

### Speed-Scaled Process Noise
- Position and velocity Q-diagonal entries auto-scale with speed:
  `speedScale = sqrt(max(v, 0.5) / 5.0)`
- At 5 m/s: scale = 1.0 (baseline). At walking (1.4 m/s): scale ≈ 0.53. At highway (33 m/s): scale ≈ 2.56
- Heading, bias, and sideslip noise are not speed-scaled

### Frame Alignment
- `setOrientation(azimuth, pitch, roll)` sets device-to-ENU rotation matrix
- When orientation is set and 6-axis IMU provided, readings rotate before bias compensation

### Robust M-Estimation (Cauchy/Huber)
- Replaces binary Mahalanobis gate with a smooth weight function when enabled
- **Cauchy**: `w = 1 / (1 + (χ²/threshold)²)` — heavy-tailed
- **Huber**: `w = threshold / max(χ², threshold)` — linear penalty
- Config: `robustWeight.enabled`, `.type` ('cauchy'|'huber'), `.threshold` (default 9.488)

### Adaptive Noise Scaling
- EMA tracks `innov² / S_innov` ratio across GPS updates
- When ratio > 1, inflates R by adaptive scale factor
- Config: `adaptiveNoise.enabled`, `.smoothing` (0–1, default 0.3), `.maxScale` (default 10)

### Coast Recovery
- When coasting and GPS gate rejects, `resetFromGps()` resets state + Cholesky from GPS fix directly
- Deduplicated from 3 copies into single private method

## Configuration

```ts
interface EkfConfig {
  dt?: number
  processNoise?: {
    position?: number           // σ per √s (default 1.0)
    velocity?: number           // σ per √s (default 0.5)
    heading?: number            // σ per √s (default 0.05)
    sideslip?: number           // σ per √s (default 0.1)
    accelBias?: number          // σ per √s (default 1e-4)
    gyroBias?: number           // σ per √s (default 1e-5)
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
  }
  gateThreshold?: number        // Mahalanobis chi-square threshold (default 9.488)
  coastTimeoutMs?: number       // GPS timeout before divergence flag (default 5000)
  gpsTimeOffsetMs?: number      // GNSS→local clock offset (default 0)
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
  p: Float64Array[];             // full covariance (7×7)
}

interface EkfDiagnostics {
  trace: number;
  gpsInnovation: number[];
  gpsChiSq: number;
  gatePassed: boolean;
  coasting: boolean;
  lastGpsTimeMs: number;
  lastImuTimeMs: number;
  stationary: boolean;
  magDeclination: number;
  robustWeight: number;
  adaNoiseScale: number;
}
```

## Determinism

Fully deterministic — no `Date.now()` calls. All timestamps are explicit parameters.

## Validation

```bash
npm test                 # 37 unit tests
npm run test:watch       # watch mode
```
