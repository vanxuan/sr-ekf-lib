# Implementation Plan — velocity_prior_konsist_20260806

> Workflow: Single-task fix. Status markers: `[ ]` = pending, `[~]` = in progress, `[x]` = done (with commit SHA).

## Phase 1: Unify velocity prior gate onto motionStillness

- [ ] Task: Update velocity prior stationarity gate in `predict()`
    - [ ] Write tests: Update 'should preserve genuine motion during coasting' to explicitly assert the prior fires (coastSpeed stays ≈ 8); no new tests needed — existing coasting tests are the regression gate.
    - [ ] Implement: In `predict()`, change gate from `(this.accelEnergy + this.gyroEnergy > 0.05 || Math.abs(this.x[I.V]) > 0.5)` to `(this.motionStillness < COAST_DAMP_STILL)`. Add `&& !this._zuptEngaged`. Remove the now-unused `stationarity` variable.
- [ ] Task: Conductor - User Manual Verification 'Phase 1' (Protocol in workflow.md)

## Definition of Done

- [ ] All tests pass (118+ test suite green)
- [ ] `npx tsc` BUILD_OK
- [ ] Committed with git note
