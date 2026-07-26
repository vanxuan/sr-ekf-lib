import { describe, it, expect } from 'vitest'
import { SrEkf } from '../src/sr-ekf'

function maxDiff(A: Float64Array[], B: Float64Array[]): number {
  let d = 0
  for (let i = 0; i < A.length; i++)
    for (let j = 0; j < A[0].length; j++)
      d = Math.max(d, Math.abs(A[i][j] - B[i][j]))
  return d
}

function computeActualPosR(basePosR: number, v: number, psi: number, beta: number,
                           gpsX: number, gpsY: number, predX: number, predY: number,
                           dtSinceLastGpsSec: number): number {
  const maxPlausibleSpeed = Math.max(v * 2, 1.0) + 2;
  const dx = gpsX - predX, dy = gpsY - predY;
  const psiBeta = psi + beta;
  const cp = Math.cos(psiBeta), sp = Math.sin(psiBeta);
  const forward = dx * cp + dy * sp;
  const cross = -dx * sp + dy * cp;
  const dtBase = Math.max(dtSinceLastGpsSec, 0.1);
  const maxForward = maxPlausibleSpeed * dtBase * 2;
  const maxCross = maxPlausibleSpeed * dtBase * 0.5;
  if (Math.abs(forward) > maxForward) basePosR *= Math.min(Math.abs(forward) / maxForward, 20);
  if (Math.abs(cross) > maxCross) basePosR *= Math.min(Math.abs(cross) / maxCross, 40);
  return basePosR;
}

// Joseph form: P⁺ = (I - K·H)·P⁻·(I - K·H)ᵀ + K·R·Kᵀ
// R is a 4×4 full covariance matrix (Float64Array[]) matching the anisotropic QR pre-array
function josephForm(P: Float64Array[], H: Float64Array[], R: Float64Array[]): Float64Array[] {
  const n = P.length, m = 4
  const HP = new Array(m)
  for (let i = 0; i < m; i++) {
    HP[i] = new Float64Array(n)
    for (let j = 0; j < n; j++) {
      let s = 0
      for (let k = 0; k < n; k++) s += H[i][k] * P[k][j]
      HP[i][j] = s
    }
  }
  const S_inn = new Array(m)
  for (let i = 0; i < m; i++) {
    S_inn[i] = new Float64Array(m)
    for (let j = 0; j <= i; j++) {
      let s = 0
      for (let k = 0; k < n; k++) s += HP[i][k] * H[j][k]
      S_inn[i][j] = s + R[i][j]
      S_inn[j][i] = S_inn[i][j]
    }
  }
  const L = new Array(m)
  for (let i = 0; i < m; i++) {
    L[i] = new Float64Array(m)
    for (let j = 0; j <= i; j++) {
      let s = 0
      for (let k = 0; k < j; k++) s += L[i][k] * L[j][k]
      if (i === j) L[i][i] = Math.sqrt(S_inn[i][i] - s)
      else L[i][j] = (S_inn[i][j] - s) / L[j][j]
    }
  }
  const K = new Array(n)
  for (let i = 0; i < n; i++) {
    K[i] = new Float64Array(m)
    const y = new Float64Array(m)
    for (let j = 0; j < m; j++) {
      let s = 0
      for (let k = 0; k < n; k++) s += P[i][k] * H[j][k]
      for (let k = 0; k < j; k++) s -= L[j][k] * y[k]
      y[j] = s / L[j][j]
    }
    for (let j = m - 1; j >= 0; j--) {
      let s = y[j]
      for (let k = j + 1; k < m; k++) s -= L[k][j] * K[i][k]
      K[i][j] = s / L[j][j]
    }
  }
  const IKH = new Array(n)
  for (let i = 0; i < n; i++) {
    IKH[i] = new Float64Array(n)
    for (let j = 0; j < n; j++) {
      let s = i === j ? 1 : 0
      for (let k = 0; k < m; k++) s -= K[i][k] * H[k][j]
      IKH[i][j] = s
    }
  }
  const Pplus = new Array(n)
  for (let i = 0; i < n; i++) {
    Pplus[i] = new Float64Array(n)
    for (let j = 0; j < n; j++) {
      let s = 0
      for (let a = 0; a < n; a++) {
        let t1 = 0
        for (let b = 0; b < n; b++) t1 += IKH[i][b] * P[b][a]
        s += t1 * IKH[j][a]
      }
      for (let a = 0; a < m; a++)
        for (let b = 0; b < m; b++)
          s += K[i][a] * R[a][b] * K[j][b]
      Pplus[i][j] = s
    }
  }
  return Pplus
}

describe('SR-EKF correctness verification', () => {
  it('P⁺ from QR matches Joseph form', () => {
    const ekf = new SrEkf({
      measurementNoise: { position: 1.0, velocity: 0.5 },
      processNoise: { position: 0.001, velocity: 0.01, heading: 0.001, sideslip: 0.01, accelBias: 0, gyroBias: 0 },
      initialCovariance: { position: 0.01, velocity: 0.01, heading: 0.01, sideslip: 0.01, accelBias: 0.001, gyroBias: 0.001 }
    })
    ekf.reset(10, 20, 5, 0.5)
    ekf.predict(0.1, 0, 0.02, 0.1, 0)
    ekf.updateGps(11, 21, 4.9, 0.3, 1)
    ekf.predict(0.05, 0, 0.01, 0.1, 2)
    const sBefore = ekf.getState()
    const P_pre = sBefore.p
    ekf.updateGps(12, 22, 4.8, 0.4, 3)
    const P_qr = ekf.getState().p

    const psi = sBefore.psi, beta = sBefore.beta, v = sBefore.v
    const psiBeta = psi + beta
    const cp = Math.cos(psiBeta), sp = Math.sin(psiBeta)
    const H: Float64Array[] = new Array(4)
    for (let i = 0; i < 4; i++) H[i] = new Float64Array(8)
    H[0][0] = 1; H[1][1] = 1
    H[2][2] = cp; H[2][3] = -v * sp; H[2][4] = -v * sp
    H[3][2] = sp; H[3][3] = v * cp; H[3][4] = v * cp
    // Match the EKF's speed-dependent velR inflation (ramp to v=10)
    const speedRamp = Math.max(0, 1 - sBefore.v / 10);
    const velR = 0.5 * (1 + 2 * speedRamp);
    // Compute posR after direction-aware guard (replicates updateGps guard logic)
    const actualPosR = computeActualPosR(1.0, v, psi, beta, 12, 22, sBefore.x, sBefore.y, 0.002);
    const ra = actualPosR * 0.5, rc = actualPosR * 1.33;
    const ra2 = ra * ra, rc2 = rc * rc;
    const R: Float64Array[] = new Array(4)
    for (let i = 0; i < 4; i++) R[i] = new Float64Array(4)
    R[0][0] = cp * cp * ra2 + sp * sp * rc2;
    R[0][1] = cp * sp * (ra2 - rc2); R[1][0] = R[0][1];
    R[1][1] = sp * sp * ra2 + cp * cp * rc2;
    R[2][2] = velR * velR; R[3][3] = velR * velR;

    expect(maxDiff(P_qr, josephForm(P_pre, H, R))).toBeLessThan(1e-4)
  })

  it('P⁺ from QR matches Joseph form with accuracy scaling', () => {
    const ekf = new SrEkf({
      measurementNoise: { position: 3.0, velocity: 0.5 },
      processNoise: { position: 0.001, velocity: 0.01, heading: 0.001, sideslip: 0.01, accelBias: 0, gyroBias: 0 },
      initialCovariance: { position: 0.01, velocity: 0.01, heading: 0.01, sideslip: 0.01, accelBias: 0.001, gyroBias: 0.001 }
    })
    ekf.reset(0, 0, 10, 0)
    ekf.predict(0, 0, 0, 0.1, 0)
    ekf.updateGps(0, 0, 10, 0, 1)
    ekf.predict(0.1, 0, 0.01, 0.1, 2)
    const sPredict = ekf.getState()
    const P_pre = sPredict.p
    ekf.updateGps(1.1, 0.1, 10.1, 0.1, 3, 6)
    const P_qr = ekf.getState().p

    const ekf2 = new SrEkf({
      measurementNoise: { position: 3.0, velocity: 0.5 },
      processNoise: { position: 0.001, velocity: 0.01, heading: 0.001, sideslip: 0.01, accelBias: 0, gyroBias: 0 },
      initialCovariance: { position: 0.01, velocity: 0.01, heading: 0.01, sideslip: 0.01, accelBias: 0.001, gyroBias: 0.001 }
    })
    ekf2.reset(0, 0, 10, 0)
    ekf2.predict(0, 0, 0, 0.1, 0)
    ekf2.updateGps(0, 0, 10, 0, 1)
    ekf2.predict(0.1, 0, 0.01, 0.1, 2)
    const s2 = ekf2.getState()

    const posR = 3 * Math.max(6 / 3, 0.1)
    const baseVelR = 0.5 * Math.max(6 / 3, 0.1)
    // Apply speed ramp as code does (v=10 → speedRamp=0)
    const velR = baseVelR * (1 + 4 * Math.max(0, 1 - s2.v / 10))
    // Compute posR after direction-aware guard
    const actualPosR = computeActualPosR(posR, s2.v, s2.psi, s2.beta, 1.1, 0.1, sPredict.x, sPredict.y, 0.002);
    const ra = actualPosR * 0.5, rc = actualPosR * 1.33;
    const ra2 = ra * ra, rc2 = rc * rc;
    const psi = s2.psi, beta = s2.beta, v = s2.v
    const psiBeta = psi + beta
    const cp = Math.cos(psiBeta), sp = Math.sin(psiBeta)
    const H: Float64Array[] = new Array(4)
    for (let i = 0; i < 4; i++) H[i] = new Float64Array(8)
    H[0][0] = 1; H[1][1] = 1
    H[2][2] = cp; H[2][3] = -v * sp; H[2][4] = -v * sp
    H[3][2] = sp; H[3][3] = v * cp; H[3][4] = v * cp
    const R2: Float64Array[] = new Array(4)
    for (let i = 0; i < 4; i++) R2[i] = new Float64Array(4)
    R2[0][0] = cp * cp * ra2 + sp * sp * rc2;
    R2[0][1] = cp * sp * (ra2 - rc2); R2[1][0] = R2[0][1];
    R2[1][1] = sp * sp * ra2 + cp * cp * rc2;
    R2[2][2] = velR * velR; R2[3][3] = velR * velR;

    expect(maxDiff(P_qr, josephForm(P_pre, H, R2))).toBeLessThan(1e-4)
  })

  // ─── Stress / Edge-case tests ──────────────────────────────

  it('handles varying dt without destabilizing', () => {
    const ekf = new SrEkf({
      measurementNoise: { position: 3.0, velocity: 1.0 },
      processNoise: { position: 0.1, velocity: 0.5, heading: 0.05, sideslip: 0.1, accelBias: 1e-5, gyroBias: 1e-6 },
      gateThreshold: 25
    })
    ekf.reset(0, 0, 5, 0)
    ekf.updateGps(0, 0, 5, 0, 0)
    const dts = [0.05, 0.01, 0.1, 0.02, 0.005]
    let clock = 1, trueX = 0
    for (let cycle = 0; cycle < 200; cycle++) {
      const dt = dts[cycle % dts.length]
      ekf.predict(0, 0, 0, dt, clock)
      trueX += 5 * dt
      clock += dt
      if (cycle % 20 === 19) {
        const ok = ekf.updateGps(trueX, 0, 5, 0, clock)
        expect(ok).toBe(true)
      }
    }
    const s = ekf.getState()
    expect(isFinite(s.x)).toBe(true)
    expect(isFinite(s.v)).toBe(true)
    expect(s.v).toBeGreaterThan(0)
  })

  it('handles reverse timestamps (out-of-order reject)', () => {
    const ekf = new SrEkf()
    ekf.reset(0, 0, 10, 0)
    ekf.predict(0, 0, 0, 0.1, 100)
    const s1 = ekf.getState()
    ekf.predict(0, 0, 0, 0.1, 50) // earlier timestamp → should be rejected
    const s2 = ekf.getState()
    expect(s2.x).toBe(s1.x)
    expect(s2.y).toBe(s1.y)
    expect(s2.v).toBe(s1.v)
  })

  it('handles dt=0 gracefully', () => {
    const ekf = new SrEkf()
    ekf.reset(0, 0, 10, 0)
    ekf.predict(0, 0, 0, 0, 100) // dt=0 → no-op
    ekf.predict(0, 0, 0, 0.1, 101) // normal
    const s = ekf.getState()
    expect(s.v).toBe(10)
  })

  it('handles very small GPS accuracy (0.01m)', () => {
    const ekf = new SrEkf({
      measurementNoise: { position: 3.0, velocity: 0.5 }
    })
    ekf.reset(0, 0, 5, 0)
    ekf.updateGps(0, 0, 5, 0, 0)
    ekf.predict(0, 0, 0, 0.1, 1)
    const ok = ekf.updateGps(0.6, 0, 5.1, 0, 2, 0.01)
    expect(ok).toBe(true)
    const s = ekf.getState()
    expect(isFinite(s.x)).toBe(true)
    expect(s.x).toBeCloseTo(0.6, 1)
  })

  it('handles very large GPS accuracy (1000m)', () => {
    const ekf = new SrEkf({
      measurementNoise: { position: 3.0, velocity: 0.5 }
    })
    ekf.reset(0, 0, 5, 0)
    ekf.updateGps(0, 0, 5, 0, 0)
    ekf.predict(0, 0, 0, 0.1, 1)
    const ok = ekf.updateGps(1000, 1000, 5, 0, 2, 1000)
    expect(ok).toBe(true)
    const s = ekf.getState()
    expect(isFinite(s.x)).toBe(true)
    // Large accuracy → measurement inflates R → small correction
    expect(Math.abs(s.x)).toBeLessThan(200)
  })

  it('recovers from successive reset() calls', () => {
    const ekf = new SrEkf()
    for (let r = 0; r < 10; r++) {
      ekf.reset(r * 10, r * 10, 5, 0)
      for (let i = 0; i < 10; i++) ekf.predict(0.1, 0, 0.01, 0.1, i)
      ekf.updateGps(r * 10 + 5, r * 10 + 5, 5, 0, 10)
    }
    const s = ekf.getState()
    expect(isFinite(s.x)).toBe(true)
    expect(isFinite(s.y)).toBe(true)
  })

  it('handles consecutive GPS updates without intermediate predict', () => {
    const ekf = new SrEkf({
      measurementNoise: { position: 2.0, velocity: 0.5 }
    })
    ekf.reset(0, 0, 5, 0)
    ekf.updateGps(0, 0, 5, 0, 0)
    ekf.updateGps(1, 0, 5, 0, 1)  // no predict in between
    const s = ekf.getState()
    expect(isFinite(s.x)).toBe(true)
  })
})
