# Code Restructure — Multi-Module Extraction

## Overview

Break up the monolithic `src/sr-ekf.ts` (~2337 lines) into a focused directory of TypeScript modules organized by functional concern. `sr-ekf.ts` remains the single public entry point (exporting `SrEkf` class, `NavigationSolution`, `EkfDiagnostics`, `EkfConfig` interfaces). Zero runtime dependencies preserved.

## Module Plan (by concern)

| Module | Contents | Rationale |
|---|---|---|
| `src/math.ts` | `qrInPlace`, `cholSolve4`, `ensurePosDiag`, `wrapAngle`, 4×4/8×8 matrix helpers, `mulScaQR`, `safeLog` | Pure functions — no instance state. Easiest extraction with bit-identical output. |
| `src/ring-buf.ts` | `RingBuf` class (push, shift, get, length) | Standalone data structure; no coupling to EKF logic. |
| `src/config.ts` | `EkfConfig`, default config values, `I` index constants, module-level constants (`EPS`, `GPS_REST_NOISE`, etc.) | Type definitions and constants — pure, extractable. |
| `src/ctra.ts` | `ctraDelta()`, `computeJacobian()`, `processNoise()`, adaptive-Q helpers | Kinematic pure functions; inputs are state + control vector, outputs are delta + Jacobian. |
| `src/diagnostics.ts` | `getDiagnostics()`, `getDebug()`, `getImuStats()` readout helpers | Pure read functions on state snapshot. |
| `src/sr-ekf.ts` (trimmed) | `SrEkf` class: `predict()`, `updateGps()`, `updateMag()`, `updateBaro()`, `setOrientation()`, `coast()`, `reset()`, plus private orchestration (`gpsUpdateSingle`, `applyZupt`, `applyZaru`, `magUpdateSingle`, `updateMotionStillness`, latency buffer management) | Entry point — wires modules together. Public API unchanged. |

## Test Strategy

- **Module-level tests**: `tests/math.test.ts` (QR decomposition, Cholesky solve, angle wrapping), `tests/ctra.test.ts` (kinematics, Jacobian continuity)
- **Existing tests preserved**: `tests/sr-ekf.test.ts` (81 tests) and `tests/qr-verification.test.ts` (9 tests) must pass with bit-identical output.
- **No new integration tests**: existing tests cover all behaviors through the public API.

## Constraints

1. **Exact numerical preservation**: Every extracted computation must produce bit-identical results. Floating-point order of operations must not change.
2. **Zero dependencies**: No runtime dependencies added.
3. **Public API unchanged**: `import { SrEkf } from 'sr-ekf'` works exactly as before. All exported types (`NavigationSolution`, `EkfDiagnostics`, `EkfConfig`) remain available from the same entry point.
4. **Build output unchanged**: `npm run build` produces `dist/sr-ekf.js` + `dist/sr-ekf.d.ts` with identical behavior. Internal modules may also appear in `dist/` but are not part of the public API.
5. **AGENTS.md stays the source of truth**: Algorithm documentation is not split; AGENTS.md remains the authoritative reference.

## Out of Scope

- Algorithm changes, new features, or behavior modifications.
- Changes to the CI pipeline or package distribution.
- Splitting `predict()` or `updateGps()` into separate files (these remain in `sr-ekf.ts` due to deep state coupling — only their pure sub-computations are extracted).

## Acceptance Criteria

1. `src/math.ts` exported — all matrix/QR/angle functions extracted, existing `qr-verification.test.ts` + new `math.test.ts` pass.
2. `src/ring-buf.ts` exported — `RingBuf` extracted, all window-dependent tests pass.
3. `src/config.ts` exported — `EkfConfig`, defaults, constants, index moved.
4. `src/ctra.ts` exported — CTRA kinematics + Jacobian extracted.
5. `src/diagnostics.ts` exported — diagnostic readouts extracted.
6. `src/sr-ekf.ts` trimmed — imports modules above, exports same public API.
7. `npm test` → 90+ tests pass (81 integration + 9 QR + new module tests).
8. `npx tsc` → BUILD_OK.
9. `git diff` shows no numerical tolerance changes in existing tests.
