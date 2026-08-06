# Track: Separate device stillness from motion stillness

## Goal

The filter currently conflates two physically distinct quantities under one "stillness" umbrella:

- **Device stillness** (`getStillness()`, IMU 3D-variance) — "is the device physically still?" A mounted phone cruising on a smooth highway is device-still ≈ 1; a hand-held phone at a red light is device-still ≈ 0.
- **Motion stillness** (velocity-domain) — "is the vehicle/device velocity ≈ 0?" Evidence: GPS speed, filter `v`, and device stillness as a proxy when GPS is lost.

This track introduces a fused `motionStillness ∈ [0,1]` metric, maintains it across `predict()`/`updateGps()`, and routes velocity-domain consumers to it while keeping device stillness for device-domain logic. This makes velocity constraints (ZUPT, coasting velocity damping, GPS `stationaryWeight`, `stationary` diagnostics) fire whenever the *vehicle* is stopped — regardless of how the device is held — and keeps raw-IMU logic (mag adaptive noise, rest-rotation `omegaScale`) tied to the device.

## Current state (baseline)

| Consumer | Current metric | File/line |
|---|---|---|
| ZUPT (`v=0` constraint) | `stillness × speedGate`, hard-overridden by `gpsMoving = lastGpsSpeed > 2.0` | `predict()` — `zuptWeight`, `gpsMoving` |
| Coasting velocity damping (`v≈0`) | `accelEnergy + gyroEnergy < 0.1` | `predict()` coasting block |
| GPS `stationaryWeight` (velocity z-blend) | `smoothedSpeed` EMA (motion-domain, separate mechanism) | `updateGps()` |
| Mag adaptive noise / gate skip | `getStillness()` | `magUpdateSingle()` |
| Heading `omegaScale` (rest-rotation) | `vehMoving × stillness` | `predict()` |
| Diagnostics `stationary` | `getStillness() > 0.7 && |v| < 3.0` | `getDiagnostics()` |
| Dead scaffolding | `varAccelEnergy` / `varGyroEnergy` computed but never read | `computeAdaptiveQ()` |

## Design

### 1. Fused `motionStillness` metric

Maintained every step, combining velocity evidence with a smooth blend:

```
speedEvidence = gpsFresh ? smoothedSpeed
             : max(|v| · k, deviceStillnessBlend)   // GPS stale: filter v, then IMU proxy
motionStillness = clamp(1 − speedEvidence / vCut, 0, 1)
```

- **GPS fresh:** driven by `smoothedSpeed` (existing 3s EMA hybrid) — authoritative when GPS alive.
- **GPS stale (coasting):** falls back to filter `v`; when `v` is unreliable (ZUPT-held), falls back to device stillness as a proxy.
- `vCut` ≈ 1.0 m/s reference speed (smooth ramp, replaces the hard `gpsMoving > 2.0` binary).

### 2. Consumer routing

- **Motion-domain** (velocity constraints) → `motionStillness`:
  - ZUPT engagement weight (`zuptWeight`)
  - Coasting velocity damping activation
  - GPS `stationaryWeight` (unify onto the shared metric)
  - Diagnostics `stationary`
- **Device-domain** (device-motion phenomena) → keep `getStillness()`:
  - Mag adaptive measurement noise (`stillnessFactor`)
  - Mag gate-skip / blend logic
  - Heading `omegaScale` rest-rotation tracking
  - ZARU (gated on filter `ω`)

### 3. Behavior notes

- Bias learning at rest is safe under hand tremor: ZUPT injects `v=0`; `a_bias_x` is corrected via `P[V][A_BIAS_X]` cross-covariance, not raw accel variance.
- The `gpsMoving` override is removed in favor of the continuous `motionStillness` ramp; GPS-confirmed cruise must still suppress ZUPT (velocity high → motionStillness low).
- Coasting behavior is preserved: with GPS lost, `motionStillness` degrades to the device-stillness proxy + filter `v`.

## Acceptance criteria

- [ ] `motionStillness` exposed and documented; `stationary` diagnostic uses it.
- [ ] Hand-held phone at a stop (high device variance, GPS speed ≈ 0) engages ZUPT and damps velocity during coasting.
- [ ] Mounted phone cruising on a smooth road (device-still ≈ 1, GPS speed high) does NOT engage ZUPT.
- [ ] Device-domain logic (mag adaptive noise, `omegaScale`, ZARU) behavior unchanged — existing rest/rotation tests still pass.
- [ ] All 78+ existing tests pass; new regression tests added; `npm run build` clean; AGENTS.md/README updated.
