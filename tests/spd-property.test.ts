import { describe, it, expect } from 'vitest'
import { SrEkf } from '../src/sr-ekf'
import { traceOfP } from '../src/math'

function isSPD(P: Float64Array[]): boolean {
  const n = P.length
  const L: Float64Array[] = Array.from({ length: n }, () => new Float64Array(n))
  for (let i = 0; i < n; i++) {
    for (let j = 0; j <= i; j++) {
      let s = P[i][j]
      for (let k = 0; k < j; k++) s -= L[i][k] * L[j][k]
      if (i === j) {
        if (s <= 1e-12) return false
        L[i][j] = Math.sqrt(s)
      } else {
        if (Math.abs(L[j][j]) < 1e-12) return false
        L[i][j] = s / L[j][j]
      }
      if (!isFinite(L[i][j])) return false
    }
  }
  return true
}

function checkState(ekf: SrEkf, label: string): void {
  const st = ekf.getState()
  // No NaN/Inf in state vector
  for (const k of ['x', 'y', 'v', 'psi', 'beta', 'aBiasX', 'gBiasZ', 'magDeclination'] as const) {
    expect(isFinite(st[k]), `${label}: ${k} isFinite`).toBe(true)
  }
  // P symmetric and diagonal positive
  const P = st.p
  const n = P.length
  for (let i = 0; i < n; i++) {
    expect(P[i][i], `${label}: P[${i}][${i}] >0`).toBeGreaterThan(0)
    for (let j = 0; j < n; j++) {
      expect(P[i][j], `${label}: P symmetric`).toBeCloseTo(P[j][i], 10)
      expect(isFinite(P[i][j]), `${label}: P[${i}][${j}] finite`).toBe(true)
    }
  }
  // traceOfP matches sum diag
  // trace via diagnostics should be sum of squares of S, which equals trace(P)
  const trace = traceOfP(P as any) // not needed; we just check P SPD
  // SPD check via Cholesky
  expect(isSPD(P), `${label}: P is SPD`).toBe(true)
  // v^T P v >0 for random probes
  for (let t = 0; t < 3; t++) {
    const v = new Float64Array(n)
    for (let i = 0; i < n; i++) v[i] = Math.sin(i * 1.3 + t * 2.7)
    let q = 0
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) q += v[i] * P[i][j] * v[j]
    expect(q, `${label}: vTPv >0 probe ${t}`).toBeGreaterThan(0)
  }
}

describe('SPD property after pseudo-measurements', () => {
  it('P remains SPD through predict + ZUPT + ZARU cycle', () => {
    const ekf = new SrEkf()
    ekf.updateGps(0, 0, 0, 0, 1000, 3)
    checkState(ekf, 'after gps init')
    // stationary predicts should trigger ZUPT/ZARU
    for (let i = 0; i < 50; i++) {
      ekf.predict(0, 0, 0, 0.02, 1100 + i * 20, 0, 0, 0)
      if (i % 10 === 0) checkState(ekf, `stationary predict ${i}`)
    }
    checkState(ekf, 'after stationary predicts')
  })

  it('P remains SPD through lateral-accel and nonholonomic branches', () => {
    const ekf = new SrEkf({ useLateralAccel: true })
    ekf.updateGps(0, 0, 5, 0, 1000, 3)
    // drive straight: nonholonomic should fire (|ω|<0.1, |v|>0.15)
    for (let i = 0; i < 20; i++) {
      ekf.predict(0.1, 0.02, 0.02, 0.02, 2000 + i * 20, 0, 0, 0)
    }
    checkState(ekf, 'after nonholonomic')
    // turn: lateral accel (|ω|>0.1, |v|>0.2) should fire
    for (let i = 0; i < 30; i++) {
      ekf.predict(0.1, 0.6, 0.4, 0.02, 3000 + i * 20, 0, 0, 0)
    }
    checkState(ekf, 'after lateral accel')
  })

  it('P remains SPD through GPS + mag + baro updates', () => {
    const ekf = new SrEkf()
    ekf.setOrientation(0.1, 0.05, -0.02)
    ekf.updateGps(0, 0, 3, 0, 1000, 3)
    checkState(ekf, 'init')
    for (let i = 0; i < 10; i++) {
      // coasting-like predicts to grow coastTime then baro
      ekf.predict(0.2, 0.05, 0.01, 0.02, 1100 + i * 20, 0, 0, 0)
    }
    // mag update at low speed
    ekf.updateMag(0.1, 1500)
    checkState(ekf, 'after mag')
    // baro during coasting (force coasting)
    ekf.coast(0, 5000)
    ekf.updateBaro(100, 1600)
    ekf.updateBaro(100.5, 2600) // dt 1s, vz=0.5, pitch ~0.05 => should apply
    checkState(ekf, 'after baro coasting')
    // baro outside coasting (opportunistic) should also stay SPD
    ekf.updateGps(10, 0, 3, 0, 3000, 3)
    ekf.updateBaro(101, 3600)
    ekf.updateBaro(101.6, 4600)
    checkState(ekf, 'after baro non-coasting')
    // GPS update with velocity
    ekf.updateGps(12, 0.2, 3.1, 0.1, 5000, 3)
    checkState(ekf, 'after gps update')
  })

  it('P remains SPD through coasting damping, velocity prior, and resetFromGps', () => {
    const ekf = new SrEkf()
    ekf.updateGps(0, 0, 8, 0, 1000, 3)
    // simulate tunnel: coast + predicts with bias drift
    ekf.coast(0, 8000) // force coasting
    for (let i = 0; i < 100; i++) {
      ekf.predict(0.05, 0.02, 0.01, 0.02, 9000 + i * 20, 0, 0, 0)
      if (i === 50) ekf.coast(0, 10000) // ensure still coasting
    }
    checkState(ekf, 'after long coasting with damping/prior')
    // large GPS jump should trigger resetFromGps and keep SPD
    ekf.updateGps(500, 0, 8, 0, 12000, 3)
    checkState(ekf, 'after resetFromGps')
  })

  it('P remains SPD under rapid rewind/replay with orientation changes', () => {
    const ekf = new SrEkf()
    ekf.setOrientation(0, 0.02, 0)
    ekf.updateGps(0, 0, 2, 0, 1000, 3)
    for (let i = 0; i < 20; i++) ekf.predict(0.1, 0, 0.05, 0.02, 1100 + i * 20, 0, 0, 0)
    ekf.setOrientation(0.2, 0.02, 0)
    for (let i = 20; i < 40; i++) ekf.predict(0.1, 0, 0.05, 0.02, 1100 + i * 20, 0, 0, 0)
    // delayed GPS triggers rewind
    ekf.updateGps(5, 0.1, 2, 0.1, 1300, 3) // timestamp behind last IMU
    checkState(ekf, 'after rewind/replay')
  })
})
