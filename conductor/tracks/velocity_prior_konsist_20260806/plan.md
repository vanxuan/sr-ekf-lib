# Implementation Plan — velocity_prior_konsist_20260806

> Workflow: Single-task fix. Status markers: `[ ]` = pending, `[~]` = in progress, `[x]` = done (with commit SHA).

## Phase 1: Unify last two `accelEnergy + gyroEnergy` holdouts

- [ ] Task: Update velocity prior and rest-weighted H block gates
    - [ ] Write tests: Update 'should preserve genuine motion during coasting' to explicitly assert the prior fires with the new gate. No new tests needed — existing coasting/heading tests are the regression gate.
    - [ ] Implement:
        1. In `predict()`, velocity prior gate: change `(this.accelEnergy + this.gyroEnergy > 0.05 || Math.abs(this.x[I.V]) > 0.5)` to `(this.motionStillness < COAST_DAMP_STILL)`. Add `!this._zuptEngaged` guard. Remove now-unused `stationarity` variable.
        2. In `gpsUpdateSingle()`, rest-weighted H gate: change `this.accelEnergy + this.gyroEnergy < 0.05` to `this.motionStillness > COAST_DAMP_STILL`.
- [ ] Task: Conductor - User Manual Verification 'Phase 1' (Protocol in workflow.md)

## Definition of Done

- [ ] All tests pass (118+ test suite green)
- [ ] `npx tsc` BUILD_OK
- [ ] Committed with git note
