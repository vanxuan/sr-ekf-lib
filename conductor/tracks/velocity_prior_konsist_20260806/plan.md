# Implementation Plan — velocity_prior_konsist_20260806

> Workflow: Single-task fix. Status markers: `[ ]` = pending, `[~]` = in progress, `[x]` = done (with commit SHA).

## Phase 1: Unify last two `accelEnergy + gyroEnergy` holdouts

- [x] Task: Update velocity prior and rest-weighted H block gates
    - [x] Write tests: Updated test comment; existing coasting/heading tests (118) pass unchanged — regression gate validates no behavioral drift.
    - [x] Implement:
        1. In `predict()`, velocity prior gate: changed to `(this.motionStillness < COAST_DAMP_STILL)`. Added `!this._zuptEngaged`. Removed unused `stationarity`.
        2. In `gpsUpdateSingle()`, rest-weighted H gate: changed to `this.motionStillness > COAST_DAMP_STILL`.
- [~] Task: Conductor - User Manual Verification 'Phase 1' (Protocol in workflow.md)

## Definition of Done

- [ ] All tests pass (118+ test suite green)
- [ ] `npx tsc` BUILD_OK
- [ ] Committed with git note
