# Implementation Plan — motion_stillness_20260806

> Workflow: Test-Driven Development. Every feature task is split into a "Write Tests" (Red) sub-task followed by an "Implement Feature" (Green) sub-task. Status markers: `[ ]` = pending, `[~]` = in progress, `[x]` = done (with commit SHA).

## Phase 1: Fused motionStillness metric

- [x] Task: Add `motionStillness` metric and maintenance `2e35848`
    - [x] Write tests: `tests/sr-ekf.test.ts` — motionStillness driven by `smoothedSpeed` when GPS fresh (cruise ⇒ low, stop ⇒ high); falls back to filter `v`/device-stillness proxy when GPS stale (coasting); smooth ramp, no hard binary. `2e35848`
    - [x] Implement: add private `motionStillness` field, compute in `predict()`/`updateGps()` using the fused formula (GPS-fresh `smoothedSpeed`, GPS-stale filter `v` + device-stillness proxy, `vCut` ramp). `2e35848`
- [x] Task: Expose `motionStillness` and update `stationary` diagnostic `2e35848`
    - [x] Write tests: `getDiagnostics().motionStillness` present and consistent; `stationary` follows motionStillness (hand-held-at-stop ⇒ true; mounted-cruise ⇒ false). `2e35848`
    - [x] Implement: add `motionStillness` to `EkfDiagnostics`; rewire `stationary = motionStillness > 0.7`. `2e35848`
- [ ] Task: Conductor - User Manual Verification 'Phase 1' (Protocol in workflow.md)

## Phase 2: Route velocity-domain consumers to motionStillness

- [ ] Task: ZUPT gating keys on motionStillness
    - [ ] Write tests: hand-held phone at a stop (high device variance, GPS speed ≈ 0) engages ZUPT and learns bias; mounted phone cruising (device-still ≈ 1, GPS speed high) does NOT engage ZUPT; removes reliance on hard `gpsMoving > 2.0`.
    - [ ] Implement: `zuptWeight = motionStillness × speedGate` in `predict()`; drop the `gpsMoving` binary override in favor of the continuous ramp.
- [ ] Task: Coasting velocity damping keys on motionStillness
    - [ ] Write tests: hand-held stop during coasting (device variance high, filter `v` → 0) still damps velocity to 0; genuine motion during coasting (motionStillness low) preserved.
    - [ ] Implement: coasting `stationarity` gate uses `motionStillness` instead of `accelEnergy + gyroEnergy`.
- [ ] Task: Unify GPS `stationaryWeight` onto motionStillness
    - [ ] Write tests: GPS velocity z-blend / velR inflation consistent with motionStillness at rest-exit and cruise (no regression in stationary-weight behavior).
    - [ ] Implement: derive `stationaryWeight` from the shared `motionStillness` metric instead of the inline `smoothedSpeed` EMA (keep smoothing).
- [ ] Task: Diagnostics `stationary` uses motionStillness (finalize)
    - [ ] Write tests: `stationary` true for a hand-held stop, false for mounted cruise, matches ZUPT engagement.
    - [ ] Implement: verify/refine `stationary` wiring; ensure `getDebug()`/`getImuStats()` unchanged.
- [ ] Task: Conductor - User Manual Verification 'Phase 2' (Protocol in workflow.md)

## Phase 3: Device-domain regression hardening and docs

- [ ] Task: Verify device-domain consumers unchanged
    - [ ] Write tests: existing table-rotation, rest-mag blend, ZARU, and omegaScale tests pass unchanged with the new metric (device stillness still drives mag adaptive noise and rest-rotation).
    - [ ] Implement: confirm no behavioral drift; adjust only if a regression is caught by the new tests.
- [ ] Task: Update AGENTS.md and README
    - [ ] Write tests: n/a (docs) — full suite green gate.
    - [ ] Implement: document `motionStillness` (definition, fusion formula, consumer routing table, coasting fallback, hand-held-at-stop behavior) in AGENTS.md; sync README validation count; update test count if changed.
- [ ] Task: Conductor - User Manual Verification 'Phase 3' (Protocol in workflow.md)

## Definition of Done

- [ ] All tests pass (`CI=true npm test`), coverage > 80%
- [ ] `npm run build` clean (tsc)
- [ ] AGENTS.md/README synced with behavior changes
- [ ] Per-task commits with git notes; phase checkpoints per workflow.md
