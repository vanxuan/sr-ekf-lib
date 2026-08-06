# Implementation Plan — code_restructure_20260806

> Workflow: Test-Driven Development. Every feature task is split into a "Write Tests" (Red) sub-task followed by an "Implement Feature" (Green) sub-task. Status markers: `[ ]` = pending, `[~]` = in progress, `[x]` = done (with commit SHA).

## Phase 1: Extract pure infrastructure (math, RingBuf, config)

- [x] Task: Extract matrix/QR math into `src/math.ts`
    - [x] Write tests: `tests/math.test.ts` — 6 unit tests for `qrInPlace`, `wrapAngle`, `ensureDiag`, `chol4x4`+`cholSolve4`, `matLowerToFull`, `traceOfP`.
    - [x] Implement: Created `src/math.ts` consolidating all functions from `src/matrix.ts` + `qrInPlace`, `wrapAngle`, `copySfromQR` from `sr-ekf.ts`. Updated import in `sr-ekf.ts`; replaced private method bodies with thin wrappers. Deleted `src/matrix.ts`. `qr-verification.test.ts` and full suite pass unchanged (96 tests).
- [ ] Task: Extract RingBuf into `src/ring-buf.ts`
    - [ ] Write tests: `tests/ring-buf.test.ts` — push, shift, get, length, wrap-around behavior, empty-buffer edge cases.
    - [ ] Implement: Create `src/ring-buf.ts` with `RingBuf` class; import in `sr-ekf.ts`; remove inline definition. All window-dependent tests pass unchanged.
- [ ] Task: Extract config types, defaults, and constants into `src/config.ts`
    - [ ] Write tests: n/a (types/constants — implicitly tested by existing suite).
    - [ ] Implement: Create `src/config.ts` with `EkfConfig`, defaults, `I` index constants, module-level constants (EPS, GPS_REST_NOISE, etc.); import in `sr-ekf.ts`.
- [ ] Task: Conductor - User Manual Verification 'Phase 1' (Protocol in workflow.md)

## Phase 2: Extract pure-domain logic (CTRA kinematics, diagnostics)

- [ ] Task: Extract CTRA kinematics into `src/ctra.ts`
    - [ ] Write tests: `tests/ctra.test.ts` — unit tests for `ctraDelta()` (big-ω and small-ω branches, zero-ω edge), Jacobian computation (verify Jacobian entries against finite differences), adaptive-Q formula.
    - [ ] Implement: Create `src/ctra.ts` with `ctraDelta()`, `computeJacobian()`, `computeAdaptiveQ()`, process-noise helpers; import in `sr-ekf.ts`; remove inline definitions. All predict-dependent tests pass unchanged.
- [ ] Task: Extract diagnostic readouts into `src/diagnostics.ts`
    - [ ] Write tests: `tests/diagnostics.test.ts` — unit tests for diagnostic formatting, debug-snapshot correctness.
    - [ ] Implement: Create `src/diagnostics.ts` with `buildDiagnostics()`, `buildDebug()`, `buildImuStats()`; import in `sr-ekf.ts`; remove inline definitions. All diagnostic tests pass unchanged.
- [ ] Task: Conductor - User Manual Verification 'Phase 2' (Protocol in workflow.md)

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
