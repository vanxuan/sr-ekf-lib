import { describe, it, expect } from 'vitest'
import { buildDiagnostics, buildDebug, buildImuStats, wmean, wstd } from '../src/diagnostics'

describe('diagnostics', () => {
  // ─── wmean / wstd ────────────────────────────────────────────
  it('wmean computes mean of RingBuf values', () => {
    // Mock RingBuf for testing
    const buf = { length: 5, get: (i: number) => [10, 20, 30, 40, 50][i] } as any
    expect(wmean(buf)).toBe(30)
  })

  it('wmean returns 0 for empty buffer', () => {
    const buf = { length: 0, get: () => 0 } as any
    expect(wmean(buf)).toBe(0)
  })

  it('wstd computes std of RingBuf values', () => {
    const buf = { length: 5, get: (i: number) => [10, 20, 30, 40, 50][i] } as any
    // mean=30, variance = (400+100+0+100+400)/5=200, std = sqrt(200) ≈ 14.14
    expect(wstd(buf, 30)).toBeCloseTo(Math.sqrt(200), 6)
  })

  it('wstd returns 0 for fewer than 2 samples', () => {
    const buf = { length: 1, get: () => 5 } as any
    expect(wstd(buf, 5)).toBe(0)
  })

  // ─── buildDiagnostics ────────────────────────────────────────
  it('buildDiagnostics returns formatted EkfDiagnostics', () => {
    const innov = new Float64Array([0.1, 0.2, 0.3, 0.4])
    const d = buildDiagnostics(42, innov, 8.5, true, false, 1000, 2000, 0, 0.9, 0.1, 1.0, 1.2)
    expect(d.trace).toBe(42)
    expect(d.gpsChiSq).toBe(8.5)
    expect(d.gatePassed).toBe(true)
    expect(d.coasting).toBe(false)
    expect(d.motionStillness).toBe(0.9)
    expect(d.stationary).toBe(true) // ms=0.9 > 0.7 && |v|=0 < 3
    expect(d.robustWeight).toBe(1.0)
  })

  it('buildDiagnostics stationary=false when v >= 3', () => {
    const innov = new Float64Array(4)
    const d = buildDiagnostics(0, innov, 0, true, false, 0, 0, 3.5, 1.0, 0, 1.0, 1.0)
    expect(d.stationary).toBe(false)
  })

  // ─── buildDebug ──────────────────────────────────────────────
  it('buildDebug returns formatted debug snapshot', () => {
    const dbg = buildDebug(0.8, 0.05, 0.01, 5, 90, 0.7, 0.5, 0, 100, 0, 0.9, 10, 15, 0.3)
    expect(dbg.stillness).toBe(0.8)
    expect(dbg.aBiasX).toBe(0.05)
    expect(dbg.v).toBe(5)
    expect(dbg.psiDeg).toBe(90)
    expect(dbg.zuptWeight).toBe(0.7)
    expect(dbg.n).toBe(100)
  })

  // ─── buildImuStats ───────────────────────────────────────────
  it('buildImuStats returns formatted IMU stats', () => {
    const buf = { length: 5, get: (i: number) => [10, 20, 30, 40, 50][i] } as any
    const stats = buildImuStats(buf, buf, 0.05, 0.01, 0.1)
    expect(stats.n).toBe(5)
    expect(stats.meanAxRel).toBeCloseTo(30 - 0.05, 6) // mean - aBiasX
    expect(stats.lastOmega).toBe(0.1)
  })
})
