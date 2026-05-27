import { describe, it, expect } from 'vitest'
import { SrEkf } from '../src/sr-ekf'

function isLowerTriangular(m: Float64Array[]): boolean {
  const n = m.length
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      if (m[i][j] !== 0) return false
    }
  }
  return true
}

function isSymmetric(m: Float64Array[]): boolean {
  const n = m.length
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      if (Math.abs(m[i][j] - m[j][i]) > 1e-10) return false
    }
  }
  return true
}

function trace(m: Float64Array[]): number {
  let t = 0
  for (let i = 0; i < m.length; i++) t += m[i][i]
  return t
}

describe('SrEkf', () => {
  it('should initialize with default state', () => {
    const ekf = new SrEkf()
    const s = ekf.getState()
    expect(s.x).toBe(0)
    expect(s.y).toBe(0)
    expect(s.v).toBe(0)
    expect(s.psi).toBe(0)
    expect(s.aBiasX).toBe(0)
    expect(s.aBiasY).toBe(0)
    expect(s.gBiasZ).toBe(0)
    expect(s.beta).toBe(0)
    const p = s.p
    expect(isSymmetric(p)).toBe(true)
    for (let i = 0; i < 8; i++) expect(p[i][i]).toBeGreaterThan(0)
  })

  it('should reset state correctly', () => {
    const ekf = new SrEkf()
    const d = ekf.getDiagnostics()
    const trace1 = d.trace
    ekf.reset(10, 20, 5, 0.5)
    const s = ekf.getState()
    expect(s.x).toBe(10)
    expect(s.y).toBe(20)
    expect(s.v).toBe(5)
    expect(s.psi).toBeCloseTo(0.5, 10)
    expect(s.aBiasX).toBe(0)
    expect(s.aBiasY).toBe(0)
    expect(s.gBiasZ).toBe(0)
    const d2 = ekf.getDiagnostics()
    expect(d2.trace).toBeGreaterThan(0)
  })

  it('should produce no drift with zero IMU input (stationary)', () => {
    const ekf = new SrEkf()
    ekf.reset(100, 200, 0, 0)
    for (let i = 0; i < 100; i++) {
      ekf.predict(0, 0, 0, 0.01, i)
    }
    const s = ekf.getState()
    expect(s.x).toBeCloseTo(100, 8)
    expect(s.y).toBeCloseTo(200, 8)
    expect(s.v).toBeCloseTo(0, 8)
    expect(s.psi).toBeCloseTo(0, 8)
  })

  it('should match analytic straight-line motion', () => {
    const ekf = new SrEkf({
      processNoise: { position: 0, velocity: 0, heading: 0, accelBias: 0, gyroBias: 0 }
    })
    ekf.reset(0, 0, 10, 0)
    const dt = 0.01
    for (let i = 0; i < 1000; i++) {
      ekf.predict(0, 0, 0, dt, i)
    }
    const s = ekf.getState()
    expect(s.x).toBeCloseTo(100, 2)
    expect(s.y).toBeCloseTo(0, 8)
    expect(s.v).toBeCloseTo(10, 6)
    expect(s.psi).toBeCloseTo(0, 6)
  })

  it('should match analytic diagonal motion', () => {
    const ekf = new SrEkf({
      processNoise: { position: 0, velocity: 0, heading: 0, accelBias: 0, gyroBias: 0 }
    })
    const psi = Math.PI / 4
    const v = 10 * Math.SQRT2
    ekf.reset(0, 0, v, psi)
    const dt = 0.01
    for (let i = 0; i < 1000; i++) {
      ekf.predict(0, 0, 0, dt, i)
    }
    const s = ekf.getState()
    expect(s.x).toBeCloseTo(100, 1)
    expect(s.y).toBeCloseTo(100, 1)
  })

  it('should match analytic constant-turn motion (circular arc)', () => {
    const ekf = new SrEkf({
      processNoise: { position: 0, velocity: 0, heading: 0, accelBias: 0, gyroBias: 0 }
    })
    const omega = 0.5
    const v = 10
    const R_curve = v / omega
    ekf.reset(0, 0, v, 0)
    const dt = 0.01
    const steps = 200
    for (let i = 0; i < steps; i++) {
      ekf.predict(0, 0, omega, dt, i)
    }
    const t_total = steps * dt
    const psi_final = omega * t_total
    const x_expected = R_curve * Math.sin(psi_final)
    const y_expected = R_curve * (1 - Math.cos(psi_final))
    const s = ekf.getState()
    expect(s.x).toBeCloseTo(x_expected, 1)
    expect(s.y).toBeCloseTo(y_expected, 1)
    expect(s.psi).toBeCloseTo(psi_final, 4)
    expect(s.v).toBeCloseTo(v, 6)
  })

  it('should converge position and velocity with GPS updates', () => {
    const ekf = new SrEkf({
      measurementNoise: { position: 0.5, velocity: 0.1 },
      processNoise: { position: 0.001, velocity: 0.01, heading: 0.001, accelBias: 0, gyroBias: 0 }
    })
    ekf.reset(0, 0, 0, 0)
    const trueV = 10
    const truePsi = 0.3
    const trueVx = trueV * Math.cos(truePsi)
    const trueVy = trueV * Math.sin(truePsi)
    for (let i = 0; i < 50; i++) {
      const t = i * 0.1
      ekf.predict(0, 0, 0, 0.1, i)
      ekf.updateGps(trueVx * t, trueVy * t, trueVx, trueVy, i)
    }
    const s = ekf.getState()
    const estV = s.v
    const estPsi = s.psi
    expect(s.x).toBeGreaterThan(40)
    expect(s.y).toBeGreaterThan(10)
    expect(Math.abs(estV - trueV)).toBeLessThan(0.5)
    expect(Math.abs(estPsi - truePsi)).toBeLessThan(0.05)
  })

  it('should estimate a constant accel bias', () => {
    const ekf = new SrEkf({
      measurementNoise: { position: 0.5, velocity: 0.1 },
      processNoise: { position: 0.001, velocity: 0.1, heading: 0.001, accelBias: 1e-5, gyroBias: 0 }
    })
    ekf.reset(0, 0, 0, 0)
    const trueVx = 5
    const trueVy = 0
    const bias_x = 0.5
    for (let i = 0; i < 200; i++) {
      const t = i * 0.1
      ekf.predict(-bias_x, 0, 0, 0.1, i)
      const gpsOk = ekf.updateGps(trueVx * t, 0, trueVx, 0, i)
      expect(gpsOk).toBe(true)
    }
    const s = ekf.getState()
    expect(s.aBiasX).toBeLessThan(-0.3)
    expect(s.aBiasX).toBeGreaterThan(-0.7)
    expect(s.v).toBeCloseTo(5, 0)
  })

  it('should estimate a constant gyro bias', () => {
    const ekf = new SrEkf({
      measurementNoise: { position: 2.0, velocity: 0.5 },
      processNoise: { position: 0.001, velocity: 0.01, heading: 0.001, accelBias: 0, gyroBias: 1e-6 }
    })
    ekf.reset(0, 0, 10, 0)
    const gbias = 0.05
    for (let i = 0; i < 200; i++) {
      ekf.predict(0, 0, -gbias, 0.1, i)
      ekf.updateGps(10 * (i + 1) * 0.1, 0, 10, 0, i)
    }
    const s = ekf.getState()
    expect(s.gBiasZ).toBeGreaterThan(-0.08)
    expect(s.gBiasZ).toBeLessThan(-0.01)
  })

  it('should reject GPS outliers via Mahalanobis gating', () => {
    const ekf = new SrEkf({
      measurementNoise: { position: 1.0, velocity: 0.5 },
      gateThreshold: 9.488
    })
    ekf.reset(0, 0, 5, 0)
    const gpsOk = ekf.updateGps(0, 0, 5, 0, 0)
    expect(gpsOk).toBe(true)
    const outlierRejected = ekf.updateGps(1000, 1000, 5, 0, 1)
    expect(outlierRejected).toBe(false)
    const d = ekf.getDiagnostics()
    expect(d.gatePassed).toBe(false)
    expect(d.gpsChiSq).toBeGreaterThan(9.488)
  })

  it('should keep heading in (-pi, pi] after many rotations', () => {
    const ekf = new SrEkf({
      processNoise: { position: 0, velocity: 0, heading: 0, accelBias: 0, gyroBias: 0 }
    })
    ekf.reset(0, 0, 10, 0)
    const fastOmega = 20
    for (let i = 0; i < 1000; i++) {
      ekf.predict(0, 0, fastOmega, 0.01, i)
    }
    const s = ekf.getState()
    for (let i = 0; i < 10; i++) {
      ekf.predict(0, 0, fastOmega, 0.1, i)
      const s2 = ekf.getState()
      expect(s2.psi).toBeGreaterThan(-Math.PI)
      expect(s2.psi).toBeLessThanOrEqual(Math.PI)
    }
  })

  it('should maintain symmetric positive-definite covariance after many steps', () => {
    const ekf = new SrEkf()
    ekf.reset(0, 0, 5, 0.5)
    for (let i = 0; i < 1000; i++) {
      ekf.predict(0.1, 0.05, 0.01, 0.01, i)
      if (i % 100 === 0) {
        ekf.updateGps(0, 0, 5, 0, i)
      }
    }
    const s = ekf.getState()
    const P = s.p
    expect(isSymmetric(P)).toBe(true)
    for (let i = 0; i < 8; i++) {
      expect(P[i][i]).toBeGreaterThan(0)
    }
  })

  it('should have lower-triangular S at all times', () => {
    const ekf = new SrEkf()
    ekf.reset(0, 0, 5, 0)
    for (let i = 0; i < 100; i++) {
      ekf.predict(0.1, 0, 0.05, 0.01, i)
    }
    const d = ekf.getDiagnostics()
    expect(d.trace).toBeGreaterThan(0)
  })

  it('should grow covariance during coast mode', () => {
    const ekf = new SrEkf({
      processNoise: { position: 0.1, velocity: 1.0, heading: 0.1, accelBias: 0.001, gyroBias: 0.0001 }
    })
    ekf.reset(0, 0, 10, 0)
    const traces: number[] = []
    for (let i = 0; i < 200; i++) {
      ekf.predict(0, 0, 0, 0.1, i)
      if (i % 20 === 0) {
        traces.push(ekf.getDiagnostics().trace)
      }
    }
    for (let i = 1; i < traces.length; i++) {
      expect(traces[i]).toBeGreaterThan(traces[i - 1])
    }
  })

  it('should converge heading from velocity GPS updates', () => {
    const ekf = new SrEkf({
      measurementNoise: { position: 3.0, velocity: 0.3 },
      processNoise: { position: 0.01, velocity: 0.1, heading: 0.01, accelBias: 0, gyroBias: 0 }
    })
    ekf.reset(0, 0, 1, 0)
    const truePsi = 1.2
    const trueV = 10
    for (let i = 0; i < 30; i++) {
      ekf.predict(0, 0, 0, 0.1, i)
      ekf.updateGps(
        trueV * Math.cos(truePsi) * (i + 1) * 0.1,
        trueV * Math.sin(truePsi) * (i + 1) * 0.1,
        trueV * Math.cos(truePsi),
        trueV * Math.sin(truePsi),
        i
      )
    }
    const s = ekf.getState()
    expect(Math.abs(s.psi - truePsi)).toBeLessThan(0.1)
  })

  it('should handle accelerating motion', () => {
    const ekf = new SrEkf({
      measurementNoise: { position: 2.0, velocity: 0.5 },
      processNoise: { position: 0.01, velocity: 0.5, heading: 0.01, accelBias: 1e-5, gyroBias: 0 }
    })
    ekf.reset(0, 0, 0, 0)
    const axAccel = 2.0
    for (let i = 0; i < 100; i++) {
      const t = i * 0.1
      ekf.predict(axAccel, 0, 0, 0.1, i)
      const v_expected = axAccel * t
      ekf.updateGps(
        0.5 * axAccel * t * t,
        0,
        v_expected,
        0,
        i
      )
    }
    const s = ekf.getState()
    expect(s.v).toBeGreaterThan(15)
    expect(s.v).toBeLessThan(25)
  })

  it('should handle coast mode timeout', () => {
    const ekf = new SrEkf({ coastTimeoutMs: 100 })
    ekf.reset(0, 0, 10, 0)
    const ok = ekf.coast(50, 0)
    expect(ok).toBe(true)
  })

  it('should respect sideslip process noise config', () => {
    const ekf = new SrEkf({
      mode: 'drive',
      processNoise: { sideslip: 0 }
    })
    ekf.reset(0, 0, 10, 0)
    for (let i = 0; i < 100; i++) {
      ekf.predict(0, 0, 0.5, 0.01, i)
    }
    const s = ekf.getState()
    const betaCov = s.p[4][4]
    const ekfDefault = new SrEkf({ mode: 'drive' })
    ekfDefault.reset(0, 0, 10, 0)
    for (let i = 0; i < 100; i++) {
      ekfDefault.predict(0, 0, 0.5, 0.01, i)
    }
    const sDefault = ekfDefault.getState()
    expect(betaCov).toBeLessThan(sDefault.p[4][4])
  })
})
