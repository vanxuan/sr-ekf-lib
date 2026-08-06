# Unify Last Two `accelEnergy + gyroEnergy` Holdouts onto motionStillness

## Overview

Two velocity-domain consumers still use the old `accelEnergy + gyroEnergy` gate instead of the fused `motionStillness` metric, missed during the `motion_stillness` track:

1. **Velocity prior** (coasting) — decides whether the car is moving to anchor `v` at `coastSpeed`. Uses `accelEnergy + gyroEnergy > 0.05 || |v| > 0.5`.
2. **Rest-weighted velocity H blocking** (updateGps) — blocks GPS velocity heading columns from rotating ψ at rest. Uses `accelEnergy + gyroEnergy < 0.05`.

Both should use `motionStillness` for consistency with all other velocity-domain consumers (ZUPT, coasting damping, GPS stationaryWeight, stationary diagnostic).

## Changes

1. **Velocity prior gate**: Replace `(this.accelEnergy + this.gyroEnergy > 0.05 || Math.abs(this.x[I.V]) > 0.5)` with `(this.motionStillness < COAST_DAMP_STILL)`. Add `!this._zuptEngaged` guard. Remove the now-unused `stationarity` variable.
2. **Rest-weighted H gate**: Replace `this.accelEnergy + this.gyroEnergy < 0.05` with `this.motionStillness > COAST_DAMP_STILL`. Keep `lastGpsSpeed < 0.3` guard.

## Test Strategy

- Update the existing 'should preserve genuine motion during coasting' test to explicitly assert the prior fires with the new gate.
- Existing coasting/heading tests must pass unchanged.

## Constraints

1. Exact numerical preservation for all existing tests — gate changes may alter behavior near thresholds, must be validated.
2. Zero behavioral change when GPS is fresh (prior only fires during coasting; H block only fires when `lastGpsSpeed < 0.3`).

## Acceptance Criteria

1. Velocity prior gate uses `motionStillness < COAST_DAMP_STILL` with `!_zuptEngaged` guard.
2. Rest-weighted H block gate uses `motionStillness > COAST_DAMP_STILL`.
3. All 118 existing tests pass.
4. The 'preserve genuine motion during coasting' test still validates the prior holds v > 7 after 12s tunnel cruise.
