import { describe, it, expect } from 'vitest'
import { qrInPlace, wrapAngle, ensureDiag, chol4x4, cholSolve4, matCreate, matLowerToFull, traceOfP } from '../src/math'

describe('math utilities', () => {
  // ─── qrInPlace ───────────────────────────────────────────────
  it('qrInPlace produces upper-triangular R', () => {
    const A: Float64Array[] = matCreate(4, 4)
    for (let i = 0; i < 4; i++)
      for (let j = 0; j < 4; j++)
        A[i][j] = (i + 1) * (j + 1) * 2 + i - j
    const vBuf = new Float64Array(4)
    qrInPlace(A, 4, 4, vBuf)
    for (let i = 1; i < 4; i++)
      for (let j = 0; j < i; j++)
        expect(Math.abs(A[i][j])).toBeLessThan(1e-10)
  })

  // ─── wrapAngle ───────────────────────────────────────────────
  it('wrapAngle normalizes to (-π, π]', () => {
    expect(wrapAngle(0)).toBeCloseTo(0)
    expect(wrapAngle(Math.PI)).toBeCloseTo(Math.PI)
    expect(wrapAngle(-Math.PI)).toBeCloseTo(Math.PI)
    expect(wrapAngle(3 * Math.PI)).toBeCloseTo(Math.PI)
    expect(wrapAngle(-3 * Math.PI)).toBeCloseTo(Math.PI)
    expect(wrapAngle(1.5)).toBeCloseTo(1.5)
    expect(wrapAngle(-1.5)).toBeCloseTo(-1.5)
    expect(wrapAngle(10)).toBeGreaterThan(-Math.PI)
    expect(wrapAngle(10)).toBeLessThanOrEqual(Math.PI)
  })

  // ─── ensureDiag ──────────────────────────────────────────────
  it('ensureDiag flips columns with negative diagonal', () => {
    const S: Float64Array[] = matCreate(3, 3)
    S[0][0] = 2; S[1][0] = 1; S[1][1] = -3; S[2][0] = 0.5; S[2][1] = 0.2; S[2][2] = 1
    ensureDiag(S)
    for (let i = 0; i < 3; i++) expect(S[i][i]).toBeGreaterThanOrEqual(0)
    // Column 0 (S[0][0]=2 >= 0) → not flipped
    expect(S[1][0]).toBe(1)
    expect(S[2][0]).toBe(0.5)
    // Column 1 (S[1][1]=-3 < 0) → negated: S[1][1]=3, S[2][1]=-0.2
    expect(S[1][1]).toBe(3)
    expect(S[2][1]).toBe(-0.2)
  })

  // ─── chol4x4 + cholSolve4 ────────────────────────────────────
  it('chol4x4 and cholSolve4 solve a 4×4 linear system', () => {
    const L = matCreate(4, 4)
    const A: Float64Array[] = matCreate(4, 4)
    for (let i = 0; i < 4; i++)
      for (let j = 0; j < 4; j++)
        A[i][j] = (i + 1) * (j + 1) * 3 + (i === j ? 10 : 0)
    const ok = chol4x4(L, A)
    expect(ok).toBe(true)
    const b = new Float64Array([1, 2, 3, 4])
    cholSolve4(L, b)
    const x = [b[0], b[1], b[2], b[3]]
    const bOrig = [1, 2, 3, 4]
    for (let i = 0; i < 4; i++) {
      let s = 0
      for (let j = 0; j < 4; j++) s += A[i][j] * x[j]
      expect(s).toBeCloseTo(bOrig[i], 2)
    }
  })

  // ─── matLowerToFull ──────────────────────────────────────────
  it('matLowerToFull reconstructs P = S·Sᵀ', () => {
    const S: Float64Array[] = matCreate(3, 3)
    S[0][0] = 2; S[1][0] = 1; S[1][1] = 3; S[2][0] = 0.5; S[2][1] = 0.2; S[2][2] = 1.5
    const P = matLowerToFull(S)
    for (let i = 0; i < 3; i++)
      for (let j = 0; j < 3; j++)
        expect(P[i][j]).toBeCloseTo(P[j][i], 10)
    for (let i = 0; i < 3; i++)
      for (let j = 0; j < 3; j++) {
        let s = 0
        for (let k = 0; k <= Math.min(i, j); k++) s += S[i][k] * S[j][k]
        expect(P[i][j]).toBeCloseTo(s, 10)
      }
  })

  // ─── traceOfP ────────────────────────────────────────────────
  it('traceOfP computes trace from Cholesky factor', () => {
    const S: Float64Array[] = matCreate(3, 3)
    S[0][0] = 2; S[1][0] = 0; S[1][1] = 3; S[2][0] = 0; S[2][1] = 0; S[2][2] = 4
    expect(traceOfP(S)).toBeCloseTo(4 + 9 + 16)
  })
})
