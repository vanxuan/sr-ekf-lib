# sr-ekf

Production-grade Square-Root Extended Kalman Filter for real-time IMU + GPS sensor fusion on mobile devices. Fuses accelerometer, gyroscope, GPS, and magnetometer data into a smooth, drift-resistant navigation solution using a 2D CTRA (Constant Turn Rate and Acceleration) motion model.

**Zero runtime dependencies.** Single TypeScript file — Float64Array-backed, zero-allocation hot path, fixed-size matrices.

## Install

```bash
npm install sr-ekf
```

## Quick Start

```ts
import { SrEkf } from 'sr-ekf'

const ekf = new SrEkf()

// On first GPS fix — snaps position and initializes heading from velocity
ekf.updateGps(longitude, latitude, vx, vy, timestampMs, accuracyMeters)

// Every IMU sample (50–400 Hz) — bias-corrected CTRA prediction
ekf.predict(ax, ay, gz, dt, timestampMs, az, gx, gy)

// Every GPS fix (1–10 Hz) — Mahalanobis-gated Kalman update
ekf.updateGps(x, y, vx, vy, timestampMs, accuracyMeters)

// Optional: magnetometer for heading at rest / low-speed
ekf.updateMag(bearing, timestampMs)

// Read state
const nav = ekf.getState()
console.log(nav.x, nav.y, nav.v, nav.psi)  // position, speed, heading
```

## State Vector (8-state)

| Index | Symbol | Description | Unit |
|-------|--------|-------------|------|
| 0 | `x` | East position | m |
| 1 | `y` | North position | m |
| 2 | `v` | Forward speed (always ≥ 0) | m/s |
| 3 | `ψ` | Heading (yaw) | rad |
| 4 | `β` | Sideslip angle | rad |
| 5 | `a_bias_x` | Accelerometer x-bias | m/s² |
| 6 | `g_bias_z` | Gyroscope z-bias | rad/s |
| 7 | `magDeclination` | Magnetic declination | rad |

## Key Features

- **Square-Root EKF** — Cholesky-factored covariance for numerical stability
- **CTRA motion model** — constant turn rate and acceleration kinematics
- **GPS latency compensation** — rewinds state to GPS timestamp and replays
- **Zero-Velocity Update (ZUPT)** — variance-based stillness detection, smooth engagement
- **Zero Angular Rate Update (ZARU)** — gyro bias calibration when not rotating
- **Magnetometer fusion** — auto-calibrates magnetic declination as 8th state
- **Robust M-estimation** — Cauchy/Huber weighting against GPS outliers
- **Anisotropic noise** — along-track vs cross-track GPS uncertainty
- **Direction-aware outlier rejection** — decomposes innovation into forward/cross-track
- **Adaptive process noise** — scales with IMU energy metrics
- **Adaptive noise scaling** — inflates measurement noise during sustained innovation
- **Nonholonomic constraint** — β ≈ 0 during straight-line driving
- **Lateral acceleration constraint** — centripetal sideslip estimation
- **Sideslip state** — separates vehicle heading from velocity direction
- **Frame alignment** — rotates IMU from device frame to ENU via device orientation
- **Low-speed jump protection** — 5 defenses against urban GPS multipath
- **Fast GPS re-acquisition** — snaps to the first valid fix after GPS loss (basement exit, tunnel exit) instead of slowly crawling; divergence auto-detection forces a full reset when the state has clearly drifted
- **Barometric speed estimation** — `updateBaro()` constrains forward speed on ramps via `vz = v × sin(pitch)` during GPS outages
- **Robustness in coasting** — GPS fixes are force-accepted (with reset) when coasting, even under Huber down-weighting, so recovery never stalls

## API

```ts
class SrEkf {
  constructor(config?: EkfConfig)
  predict(ax, ay, gz, dt, timestampMs, az?, gx?, gy?): void
  updateGps(x, y, vx, vy, timestampMs, accuracyMeters?): boolean
  updateMag(bearing, timestampMs): void
  updateBaro(altitude, timestampMs): void      // barometric speed on ramps during GPS outage
  setOrientation(azimuth, pitch, roll): void
  coast(timeoutMs, currentTimeMs): boolean
  getState(): NavigationSolution
  getStateInto(out: NavigationSolution): void   // zero-allocation path
  getImuStats(): { n, meanAxRel, stdAx, meanGzRel, stdGz, lastOmega }
  getStillness(): number                        // 0 = in motion, 1 = at rest (IMU variance-based)
  getDiagnostics(): EkfDiagnostics
  reset(x, y, v, psi): void
  resetBiases(aBiasX?, gBiasZ?): void
  inflateCovariance(params): void
}
```

### `getState()` / `getStateInto()`

```ts
interface NavigationSolution {
  x: number; y: number;     // position (m)
  v: number;                 // speed (m/s)
  psi: number;               // heading (rad)
  beta: number;              // sideslip angle (rad)
  aBiasX: number;            // accelerometer bias (m/s²)
  gBiasZ: number;            // gyroscope bias (rad/s)
  magDeclination: number;    // magnetic declination (rad)
  p: Float64Array[];         // 8×8 covariance matrix
}
```

Use `getStateInto()` on the render loop to avoid allocations.

### `getDiagnostics()`

```ts
interface EkfDiagnostics {
  trace: number;             // covariance trace
  gpsInnovation: Float64Array; // [dx, dy, dvx, dvy]
  gpsChiSq: number;          // last GPS χ²
  gatePassed: boolean;       // last gate result
  coasting: boolean;         // true during GPS outage
  lastGpsTimeMs: number;
  lastImuTimeMs: number;
  stationary: boolean;       // motionStillness > 0.7 AND |v| < 3
  motionStillness: number;   // fused metric (0=moving, 1=stopped)
  magDeclination: number;
  robustWeight: number;      // M-estimator weight (1 = no downweighting)
  adaNoiseScale: number;     // adaptive noise inflation factor
}
```

## Configuration

```ts
new SrEkf({
  processNoise: {
    position: 5.0,           // σ per √s
    velocity: 0.5,
    heading: 0.10,
    sideslip: 0.1,
    accelBias: 1e-4,
    gyroBias: 5e-5,
    magDeclination: 1e-4,
  },
  measurementNoise: {
    position: 3.0,           // σ for GPS position
    velocity: 0.5,           // σ for GPS velocity
    heading: 0.1,            // σ for magnetometer
  },
  initialCovariance: {
    position: 100,           // σ²
    velocity: 10,
    heading: 9.87,           // ≈ π²
    sideslip: 0.25,
    accelBias: 0.1,
    gyroBias: 0.01,
    magDeclination: 0.25,
  },
  magneticDeclination: 0,    // initial mag declination (rad)
  gateThreshold: 9.488,      // Mahalanobis χ² gate (95% for 4-DOF)
  coastTimeoutMs: 5000,      // GPS timeout before coasting
  gpsTimeOffsetMs: 0,        // GNSS→local clock offset
  useLateralAccel: true,     // centripetal sideslip constraint
  robustWeight: {
    enabled: true,           // robust M-estimation
    type: 'huber',           // 'cauchy' | 'huber'
    threshold: 9.488,
  },
  adaptiveNoise: {
    enabled: false,
    smoothing: 0.1,          // 0–1
    maxScale: 3,
  },
})
```

## How It Works

1. **IMU predict** (high rate): bias-corrected readings drive CTRA kinematics. Gyro bias is corrected by ZARU when not rotating. Covariance propagated via QR decomposition.

2. **GPS update** (low rate): 4-DOF (x, y, vx, vy) Kalman update with Mahalanobis gating. Position uses anisotropic noise model. Velocity direction corrects heading when speed > 0.5 m/s.

3. **Magnetometer update** (medium rate): heading observation with auto-calibrating declination state. Only trusted at low speed (skipped when `|v| > 1.5 m/s`) to avoid fighting GPS velocity direction.

4. **ZUPT** (on IMU): when the device is stationary (variance-based detection), a zero-velocity pseudo-measurement corrects biases through cross-covariance.

5. **GPS re-acquisition** (after loss): when coasting and a GPS fix's innovation is statistically implausible, the filter snaps directly to the fix (`resetFromGps`) — biases learned during coasting are preserved. The outlier guard no longer masks large position jumps: if it inflates measurement noise >3× and the jump is physically plausible, the fix is force-accepted (when moving >0.5 m/s or the jump exceeds 15m) so recovery is immediate instead of minute-long.

## Mobile Performance

- **Zero allocations** on `predict()` and `getStateInto()` hot paths
- All working matrices preallocated as `Float64Array` in constructor
- Fixed-size: N=8 states, M=4 GPS observations, QR pre-arrays 16×8 / 12×12 / 9×9
- Designed for 50–400 Hz IMU and 1–10 Hz GPS on resource-constrained devices

## Validation

```bash
npm test                 # 90 tests (9 QR verification + 81 unit)
npm run build            # TypeScript → dist/
```

## License

MIT

## Contributing

Contributions welcome. Run `npm test` before submitting PRs. The algorithm is fully documented in [AGENTS.md](AGENTS.md).
