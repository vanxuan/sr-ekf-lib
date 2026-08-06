import { describe, it, expect } from 'vitest'
import { ctraDelta, computeJacobian } from '../src/ctra'
import { EPS } from '../src/config'

describe('CTRA kinematics', () => {
  // ─── ctraDelta ───────────────────────────────────────────────
  it('ctraDelta with nonzero omega computes position delta', () => {
    const [dx, dy] = ctraDelta(0.5, 5, 0.2, 0.1)
    // Approximate check: moving at 5 m/s for 0.1s with slight turn
    expect(dx).toBeGreaterThan(0)
    expect(dy).toBeGreaterThan(0)
    expect(Math.sqrt(dx * dx + dy * dy)).toBeCloseTo(0.5, 1)
  })

  it('ctraDelta with near-zero omega uses small-angle branch', () => {
    const [dx, dy] = ctraDelta(0.5, 5, 1e-6, 0.1)
    // Straight line: dx = v·cos(ψ)·dt, dy = v·sin(ψ)·dt
    expect(dx).toBeCloseTo(5 * Math.cos(0.5) * 0.1, 6)
    expect(dy).toBeCloseTo(5 * Math.sin(0.5) * 0.1, 6)
  })

  it('ctraDelta with zero psi gives x-only displacement', () => {
    const [dx, dy] = ctraDelta(0, 10, 0.3, 0.05)
    // Turning right from heading 0
    expect(dx).toBeGreaterThan(0)
    expect(Math.abs(dy)).toBeLessThan(0.2)
  })

  // ─── computeJacobian ─────────────────────────────────────────
  it('computeJacobian has non-decay diagonals as 1', () => {
    const F: Float64Array[] = [new Float64Array(8), new Float64Array(8), new Float64Array(8), new Float64Array(8), new Float64Array(8), new Float64Array(8), new Float64Array(8), new Float64Array(8)]
    computeJacobian(0.5, 5, 0, 1e-6, 0.1, 1.5, 0.98, EPS, F)
    // Decay states (beta and accelBias) have non-1 diagonal
    for (const i of [0, 1, 2, 3, 6, 7]) expect(F[i][i]).toBeCloseTo(1, 6)
    expect(F[4][4]).toBeCloseTo(Math.exp(-0.1 / 1.5), 6)
    expect(F[5][5]).toBeCloseTo(0.98, 6)
  })

  it('computeJacobian small-omega branch position-velocity coupling', () => {
    const F: Float64Array[] = [new Float64Array(8), new Float64Array(8), new Float64Array(8), new Float64Array(8), new Float64Array(8), new Float64Array(8), new Float64Array(8), new Float64Array(8)]
    computeJacobian(0, 5, 0, 1e-6, 0.1, 0.8, 0.98, EPS, F)
    // Straight line: dx/dv = cos(ψ)·dt, dy/dv = sin(ψ)·dt
    expect(F[0][2]).toBeCloseTo(1 * 0.1, 6)
    expect(F[1][2]).toBeCloseTo(0 * 0.1, 6)
  })

  it('computeJacobian small-omega branch has gyro bias derivatives', () => {
    const F: Float64Array[] = [new Float64Array(8), new Float64Array(8), new Float64Array(8), new Float64Array(8), new Float64Array(8), new Float64Array(8), new Float64Array(8), new Float64Array(8)]
    const psi = 1.2, v = 4, a = 0.5, dt = 0.1
    computeJacobian(psi, v, a, 1e-6, dt, 0.8, 0.98, EPS, F)
    const vAvg = v + 0.5 * a * dt
    const sp = Math.sin(psi), cp = Math.cos(psi)
    // ∂x/∂g_bias_z = -0.5·v̅·sin(ψ)·dt²
    expect(F[0][6]).toBeCloseTo(-0.5 * vAvg * sp * dt * dt, 6)
    // ∂y/∂g_bias_z = 0.5·v̅·cos(ψ)·dt²
    expect(F[1][6]).toBeCloseTo(0.5 * vAvg * cp * dt * dt, 6)
  })

  it('computeJacobian beta mean-reversion term', () => {
    const F: Float64Array[] = [new Float64Array(8), new Float64Array(8), new Float64Array(8), new Float64Array(8), new Float64Array(8), new Float64Array(8), new Float64Array(8), new Float64Array(8)]
    const betaTau = 1.5
    computeJacobian(0, 5, 0, 0.05, 0.1, betaTau, 0.98, EPS, F)
    expect(F[4][4]).toBeCloseTo(Math.exp(-0.1 / betaTau), 6)
  })
})
