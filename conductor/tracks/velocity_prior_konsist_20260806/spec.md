# Unify Velocity Prior Gate onto motionStillness

## Overview

The velocity prior (a scalar QR pseudo-measurement that pulls `v` toward `coastSpeed` during coasting) uses the **old** `accelEnergy + gyroEnergy` gate to decide whether the car is moving. This was missed during the `motion_stillness` track — all other velocity-domain consumers (ZUPT, coasting damping, GPS stationaryWeight, stationary diagnostic) were migrated to the fused `motionStillness` metric. The velocity prior is the last holdout.

## Change

1. **Gate**: Replace `this.accelEnergy + this.gyroEnergy > 0.05 || Math.abs(this.x[I.V]) > 0.5` with `this.motionStillness < COAST_DAMP_STILL` (0.5) — consistent with the coasting damping gate from the `motion_stillness` track.
2. **ZUPT guard**: Add `!this._zuptEngaged` to prevent the prior from competing with ZUPT during a stop.

## Test Strategy

- Update the existing 'should preserve genuine motion during coasting' test to verify the prior still fires correctly with the new gate.
- Existing coasting tests (damping, tunnel stop, fallback) must pass unchanged.

## Constraints

1. Exact numerical preservation for all existing tests — the gate change may alter behavior during coasting near the 0.5 threshold, which must be validated.
2. Zero behavioral change when GPS is fresh (prior only fires during coasting).

## Acceptance Criteria

1. Velocity prior gate uses `motionStillness < COAST_DAMP_STILL`.
2. `!_zuptEngaged` guard added.
3. All 118 existing tests pass.
4. The 'preserve genuine motion during coasting' test still validates the prior holds v > 7 after 12s tunnel cruise.
