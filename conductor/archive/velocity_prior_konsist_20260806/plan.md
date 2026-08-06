# Implementation Plan — velocity_prior_konsist_20260806

> Workflow: Single-task fix. Status markers: `[ ]` = pending, `[~]` = in progress, `[x]` = done (with commit SHA).

## Phase 1: Unify last two `accelEnergy + gyroEnergy` holdouts `[checkpoint: bab2824]`

- [x] Task: Update velocity prior and rest-weighted H block gates
    - [x] Write tests: Updated test comment; existing coasting/heading tests (118) pass unchanged — regression gate validates no behavioral drift.
    - [x] Implement:
        1. In `predict()`, velocity prior gate: changed to `(this.motionStillness < COAST_DAMP_STILL)`. Added `!this._zuptEngaged`. Removed unused `stationarity`.
        2. In `gpsUpdateSingle()`, rest-weighted H gate: changed to `this.motionStillness > COAST_DAMP_STILL`.
- [x] Task: Conductor - User Manual Verification 'Phase 1' (Protocol in workflow.md): user confirmed yes; verification report attached to checkpoint `bab2824`.

## Definition of Done

- [x] All tests pass — 118/118 green
- [x] `npx tsc` BUILD_OK
- [x] Committed with git note
