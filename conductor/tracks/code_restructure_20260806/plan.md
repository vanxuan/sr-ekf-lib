# Implementation Plan — code_restructure_20260806

> Workflow: Test-Driven Development. Every feature task is split into a "Write Tests" (Red) sub-task followed by an "Implement Feature" (Green) sub-task. Status markers: `[ ]` = pending, `[~]` = in progress, `[x]` = done (with commit SHA).

## Phase 1: Extract pure infrastructure (math, RingBuf, config) `[checkpoint: e915fb5]`

- [x] Task: Extract matrix/QR math into `src/math.ts`
    - [x] Write tests: `tests/math.test.ts` — 6 unit tests for `qrInPlace`, `wrapAngle`, `ensureDiag`, `chol4x4`+`cholSolve4`, `matLowerToFull`, `traceOfP`.
    - [x] Implement: Created `src/math.ts` consolidating all functions from `src/matrix.ts` + `qrInPlace`, `wrapAngle`, `copySfromQR` from `sr-ekf.ts`. Updated import in `sr-ekf.ts`; replaced private method bodies with thin wrappers. Deleted `src/matrix.ts`. `qr-verification.test.ts` and full suite pass unchanged (96 tests).
- [x] Task: Extract RingBuf into `src/ring-buf.ts`
    - [x] Write tests: `tests/ring-buf.test.ts` — 7 unit tests for push, shift, get, wrap-around, clear, capacity.
    - [x] Implement: Created `src/ring-buf.ts` with `RingBuf` class; imported in `sr-ekf.ts`; removed inline definition. Fixed latent shift() bug (computed index before decrement). All window-dependent tests pass unchanged (103 tests).
- [x] Task: Extract config types, defaults, and constants into `src/config.ts`
    - [x] Write tests: n/a (types/constants — implicitly tested by existing suite).
    - [x] Implement: Created `src/config.ts` with `EkfConfig`, defaults, `I` (as const object), `N`/`M`/`PRE`/`MAG_PRE`/`NTRI`, `DEFAULTS`, `EPS`, all motion constants. Re-exported types from `sr-ekf.ts` for downstream consumers.
- [x] Task: Conductor - User Manual Verification 'Phase 1' (Protocol in workflow.md): user confirmed yes; verification report attached to checkpoint `e915fb5`.

## Phase 2: Extract pure-domain logic (CTRA kinematics, diagnostics) `[checkpoint: acb53ad]`

- [x] Task: Extract CTRA kinematics into `src/ctra.ts`
    - [x] Write tests: `tests/ctra.test.ts` — 7 unit tests for `ctraDelta()` (big-ω, small-ω, zero-psi branches) and `computeJacobian()` (diagonal, position-velocity, gyro-bias derivatives, beta decay).
    - [x] Implement: Created `src/ctra.ts` with exported `ctraDelta()` (pure) and `computeJacobian()` (takes state + `F` output array). `computeAdaptiveQ()` kept in `sr-ekf.ts` due to deep state coupling. Updated call sites; removed private methods. All predict-dependent tests pass unchanged (110 tests).
- [x] Task: Extract diagnostic readouts into `src/diagnostics.ts`
    - [x] Write tests: `tests/diagnostics.test.ts` — 8 unit tests for `wmean`, `wstd`, `buildDiagnostics`, `buildDebug`, `buildImuStats`.
    - [x] Implement: Created `src/diagnostics.ts` with pure snapshot-formatting functions; updated `sr-ekf.ts` `getDiagnostics()`/`getDebug()`/`getImuStats()` to delegate. All diagnostic tests pass unchanged (118 tests).
- [x] Task: Conductor - User Manual Verification 'Phase 2' (Protocol in workflow.md): user confirmed yes; verification report attached to checkpoint `acb53ad`.

## Phase 3: Trim sr-ekf.ts entry point and final cleanup

- [ ] Task: Update `sr-ekf.ts` to cleanly import all modules, re-export public types
    - [ ] Write tests: n/a — existing 90+ tests are the regression gate.
    - [ ] Implement: Remove now-duplicate imports/inline definitions; verify `SrEkf` class methods reference imported functions; re-export `NavigationSolution`, `EkfDiagnostics`, `EkfConfig` from the entry point for downstream consumers.
- [ ] Task: Verify build and publish artifacts
    - [ ] Write tests: n/a — `npx tsc` is the gate.
    - [ ] Implement: Run `npm run build`; verify `dist/sr-ekf.js` + `dist/sr-ekf.d.ts` produced; verify no new exports leak. Update `package.json` `files` whitelist if needed; update imports in test imports if `tsconfig.json` path resolution changed.
- [ ] Task: Conductor - User Manual Verification 'Phase 3' (Protocol in workflow.md)

## Definition of Done

- [ ] All tests pass (`npm test`), coverage > 80%
- [ ] `npm run build` clean (tsc)
- [ ] Public API unchanged — existing consumers need no code changes
- [ ] AGENTS.md/README updated if file structure/import paths changed
- [ ] Per-task commits with git notes; phase checkpoints per workflow.md
