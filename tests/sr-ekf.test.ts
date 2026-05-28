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
    expect(s.gBiasZ).toBe(0)
    expect(s.beta).toBe(0)
    const p = s.p
    expect(isSymmetric(p)).toBe(true)
    for (let i = 0; i < 7; i++) expect(p[i][i]).toBeGreaterThan(0)
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
    for (let i = 0; i < 7; i++) {
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

  it('should recover heading from 180° error using GPS velocity', () => {
    const ekf = new SrEkf({
      measurementNoise: { position: 3.0, velocity: 0.5 }
    })
    ekf.reset(0, 0, 5, Math.PI)
    ekf.updateGps(0, 0, -5, 0, 0)
    ekf.updateGps(10, 0, 5, 0, 1)
    const s = ekf.getState()
    expect(Math.abs(s.psi)).toBeLessThan(0.1)
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
      processNoise: { sideslip: 0 }
    })
    ekf.reset(0, 0, 10, 0)
    for (let i = 0; i < 100; i++) {
      ekf.predict(0, 0, 0.5, 0.01, i)
    }
    const s = ekf.getState()
    const betaCov = s.p[4][4]
    const ekfDefault = new SrEkf()
    ekfDefault.reset(0, 0, 10, 0)
    for (let i = 0; i < 100; i++) {
      ekfDefault.predict(0, 0, 0.5, 0.01, i)
    }
    const sDefault = ekfDefault.getState()
    expect(betaCov).toBeLessThan(sDefault.p[4][4])
  })

  it('should detect stationary via IMU variance and trigger ZUPT', () => {
    const ekf = new SrEkf()
    ekf.reset(100, 200, 0, 0)
    for (let i = 0; i < 50; i++) {
      ekf.predict(0, 0, 0, 0.01, i)
    }
    const d = ekf.getDiagnostics()
    expect(d.stationary).toBe(true)
    const s = ekf.getState()
    expect(s.v).toBeCloseTo(0, 6)
    const posCov = s.p[0][0]
    expect(posCov).toBeLessThan(2)
  })

  it('should restore motion after stationary ends', () => {
    const ekf = new SrEkf()
    ekf.reset(100, 200, 0, 0)
    for (let i = 0; i < 50; i++) {
      ekf.predict(0, 0, 0, 0.01, i)
    }
    expect(ekf.getDiagnostics().stationary).toBe(true)
    for (let i = 0; i < 30; i++) {
      ekf.predict(1, 0, 0, 0.01, i + 50)
    }
    const d = ekf.getDiagnostics()
    expect(d.stationary).toBe(false)
  })

  it('should converge heading from magnetometer updates', () => {
    const ekf = new SrEkf({ measurementNoise: { heading: 0.1 } })
    ekf.reset(0, 0, 5, Math.PI)
    for (let i = 0; i < 50; i++) {
      ekf.updateMag(0, i)
    }
    const s = ekf.getState()
    expect(Math.abs(s.psi)).toBeLessThan(0.5)
  })

  it('should auto-initialize heading from first mag update (no GPS)', () => {
    const ekf = new SrEkf()
    ekf.reset(0, 0, 0, 0)
    expect(ekf.getState().psi).toBe(0)
    ekf.updateMag(1.5, 0)
    const s = ekf.getState()
    expect(s.psi).toBeCloseTo(1.5, 5)
  })

  it('should rotate IMU readings with non-identity orientation', () => {
    const ekf = new SrEkf({
      processNoise: { position: 0, velocity: 0, heading: 0, accelBias: 0, gyroBias: 0 }
    })
    ekf.reset(0, 0, 0, 0)
    ekf.setOrientation(Math.PI / 2, 0, 0)
    ekf.predict(1, 0, 0, 0.1, 0, 0, 0, 0)
    const s = ekf.getState()
    expect(s.v).toBeCloseTo(0, 4)
  })

  it('should rotate device-y acceleration to east with azimuth=π/2', () => {
    const ekf = new SrEkf({
      processNoise: { position: 0, velocity: 0, heading: 0, accelBias: 0, gyroBias: 0 }
    })
    ekf.reset(0, 0, 0, 0)
    ekf.setOrientation(Math.PI / 2, 0, 0)
    ekf.predict(0, 1, 0, 0.1, 0)
    ekf.predict(0, 1, 0, 0.1, 1)
    const s = ekf.getState()
    expect(s.x).toBeGreaterThan(0.005)
  })

  it('should preserve behavior with identity orientation', () => {
    const ekf = new SrEkf({
      processNoise: { position: 0, velocity: 0, heading: 0, accelBias: 0, gyroBias: 0 }
    })
    ekf.reset(0, 0, 0, 0)

    const ekfRef = new SrEkf({
      processNoise: { position: 0, velocity: 0, heading: 0, accelBias: 0, gyroBias: 0 }
    })
    ekfRef.reset(0, 0, 0, 0)

    ekf.setOrientation(0, 0, 0)
    ekf.predict(2, 0, 0.5, 0.1, 0, 0, 0, 0)
    ekfRef.predict(2, 0, 0.5, 0.1, 0)

    const s1 = ekf.getState()
    const s2 = ekfRef.getState()
    expect(s1.x).toBeCloseTo(s2.x, 10)
    expect(s1.y).toBeCloseTo(s2.y, 10)
    expect(s1.v).toBeCloseTo(s2.v, 10)
    expect(s1.psi).toBeCloseTo(s2.psi, 10)
  })

  it('should separate sideslip from heading using mag + GPS', () => {
    const ekf = new SrEkf({
      measurementNoise: { position: 10.0, velocity: 0.3, heading: 0.05 },
      processNoise: { position: 0.1, velocity: 0.1, heading: 0.001, sideslip: 0.2, accelBias: 0, gyroBias: 0 }
    })
    ekf.reset(0, 0, 10, 0)
    const omega = 0.3
    const trueBeta = 0.08
    for (let i = 0; i < 200; i++) {
      const t = i * 0.1
      ekf.predict(0, 0, omega, 0.1, i)
      const psi_t = omega * t
      const alpha = psi_t + trueBeta
      const x = (10 / omega) * (Math.sin(omega * t + trueBeta) - Math.sin(trueBeta))
      const y = (10 / omega) * (-Math.cos(omega * t + trueBeta) + Math.cos(trueBeta))
      ekf.updateMag(psi_t, i)
      ekf.updateGps(x, y, 10 * Math.cos(alpha), 10 * Math.sin(alpha), i, 15)
    }
    const s = ekf.getState()
    expect(s.beta).toBeGreaterThan(0.03)
  })

  it('should report coasting when timeout exceeded', () => {
    const ekf = new SrEkf({
      processNoise: { position: 10, velocity: 10, heading: 1, accelBias: 0.1, gyroBias: 0.01 }
    })
    ekf.reset(0, 0, 100, 0)
    ekf.updateGps(0, 0, 100, 0, 0)
    for (let i = 1; i < 5000; i++) {
      ekf.predict(0, 0, 0, 0.01, i)
    }
    const ok = ekf.coast(100, 5000)
    expect(ok).toBe(true)
    const d = ekf.getDiagnostics()
    expect(d.coasting).toBe(true)
  })

  it('should accept larger innovations with high accuracyMeters', () => {
    const ekf = new SrEkf({
      measurementNoise: { position: 3.0, velocity: 0.5 },
      gateThreshold: 9.488
    })
    ekf.reset(0, 0, 5, 0)
    ekf.updateGps(0, 0, 5, 0, 0, 3)
    const result = ekf.updateGps(1000, 1000, 5, 0, 1, 500)
    expect(result).toBe(true)
    expect(ekf.getDiagnostics().gatePassed).toBe(true)

    const ekf2 = new SrEkf({
      measurementNoise: { position: 3.0, velocity: 0.5 },
      gateThreshold: 9.488
    })
    ekf2.reset(0, 0, 5, 0)
    ekf2.updateGps(0, 0, 5, 0, 0)
    const result2 = ekf2.updateGps(1000, 1000, 5, 0, 1)
    expect(result2).toBe(false)
  })

  it('should calibrate magnetic declination from GPS heading', () => {
    const ekf = new SrEkf({
      measurementNoise: { position: 2.0, velocity: 0.3 }
    })
    ekf.reset(0, 0, 10, 0)
    const trueHdg = 0.5
    const magDeclination = 0.2
    for (let i = 0; i < 10; i++) {
      const t = i * 0.1
      ekf.predict(0, 0, 0, 0.1, i)
      ekf.updateGps(
        10 * Math.cos(trueHdg) * t,
        10 * Math.sin(trueHdg) * t,
        10 * Math.cos(trueHdg),
        10 * Math.sin(trueHdg),
        i
      )
    }
    ekf.updateMag(trueHdg - magDeclination, 100)
    expect(ekf.getDiagnostics().magDeclination).toBeCloseTo(0, 2)
    for (let i = 0; i < 40; i++) {
      const t = i * 0.1 + 1.0
      ekf.predict(0, 0, 0, 0.1, i + 100)
      ekf.updateGps(
        10 * Math.cos(trueHdg) * t,
        10 * Math.sin(trueHdg) * t,
        10 * Math.cos(trueHdg),
        10 * Math.sin(trueHdg),
        i + 100,
        5
      )
      ekf.updateMag(trueHdg - magDeclination, i + 100)
    }
    const d = ekf.getDiagnostics()
    expect(d.magDeclination).toBeGreaterThan(0.05)
  })

  it('should return NavigationSolution with correct shape', () => {
    const ekf = new SrEkf()
    ekf.reset(10, 20, 5, 0.5)
    const s = ekf.getState()
    expect(s).toHaveProperty('x')
    expect(s).toHaveProperty('y')
    expect(s).toHaveProperty('v')
    expect(s).toHaveProperty('psi')
    expect(s).toHaveProperty('beta')
    expect(s).toHaveProperty('aBiasX')
    expect(s).toHaveProperty('gBiasZ')
    expect(s).toHaveProperty('p')
    expect(s.p.length).toBe(7)
    expect(s.p[0].length).toBe(7)
    expect(isSymmetric(s.p)).toBe(true)
  })

  it('should return EkfDiagnostics with correct shape', () => {
    const ekf = new SrEkf()
    ekf.reset(0, 0, 5, 0)
    ekf.predict(0, 0, 0, 0.01, 0)
    ekf.updateGps(0, 0, 5, 0, 1)
    const d = ekf.getDiagnostics()
    expect(d).toHaveProperty('trace')
    expect(d).toHaveProperty('gpsInnovation')
    expect(d.gpsInnovation).toHaveLength(4)
    expect(d).toHaveProperty('gpsChiSq')
    expect(d).toHaveProperty('gatePassed')
    expect(d).toHaveProperty('coasting')
    expect(d).toHaveProperty('lastGpsTimeMs')
    expect(d).toHaveProperty('lastImuTimeMs')
    expect(d).toHaveProperty('stationary')
    expect(d).toHaveProperty('magDeclination')
    expect(d).toHaveProperty('robustWeight')
    expect(d).toHaveProperty('adaNoiseScale')
  })

  it('should downweight large innovations with robust M-estimation (Cauchy)', () => {
    const ekf = new SrEkf({
      measurementNoise: { position: 1.0, velocity: 0.5 },
      gateThreshold: 9.488,
      robustWeight: { enabled: true, type: 'cauchy', threshold: 9.488 }
    })
    ekf.reset(0, 0, 5, 0)
    ekf.updateGps(0, 0, 5, 0, 0)
    const w1 = ekf.getDiagnostics().robustWeight
    expect(w1).toBeCloseTo(1, 2)

    ekf.updateGps(1000, 1000, 5, 0, 1)
    const w2 = ekf.getDiagnostics().robustWeight
    expect(w2).toBeLessThan(0.5)
    expect(w2).toBeGreaterThan(0)
    const d = ekf.getDiagnostics()
    expect(d.gatePassed).toBe(true)
  })

  it('should downweight with Huber robust M-estimation', () => {
    const ekf = new SrEkf({
      measurementNoise: { position: 1.0, velocity: 0.5 },
      gateThreshold: 9.488,
      robustWeight: { enabled: true, type: 'huber', threshold: 9.488 }
    })
    ekf.reset(0, 0, 5, 0)
    for (let i = 0; i < 10; i++) {
      ekf.predict(0, 0, 0, 0.1, i + 1)
      ekf.updateGps(0.5 * (i + 1), 0, 5, 0, i + 1)
    }
    ekf.updateGps(1000, 1000, 5, 0, 100)
    const w = ekf.getDiagnostics().robustWeight
    expect(w).toBeLessThan(0.5)
    expect(w).toBeGreaterThan(0)
    const s = ekf.getState()
    expect(Math.abs(s.x - 5)).toBeLessThan(20)
    expect(ekf.getDiagnostics().gatePassed).toBe(true)
  })

  it('should reject extreme outliers via binary gate when robust M is disabled', () => {
    const ekf = new SrEkf({
      measurementNoise: { position: 1.0, velocity: 0.5 },
      gateThreshold: 9.488,
      robustWeight: { enabled: false }
    })
    ekf.reset(0, 0, 5, 0)
    ekf.updateGps(0, 0, 5, 0, 0)
    const result = ekf.updateGps(1e6, 1e6, 5, 0, 1)
    expect(result).toBe(false)
  })

  it('should inflate R via adaptive noise when innovations are large', () => {
    const ekf = new SrEkf({
      measurementNoise: { position: 1.0, velocity: 0.5 },
      adaptiveNoise: { enabled: true, smoothing: 0.3, maxScale: 20 }
    })
    ekf.reset(0, 0, 5, 0)
    ekf.updateGps(0, 0, 5, 0, 0)
    expect(ekf.getDiagnostics().adaNoiseScale).toBe(1)
    for (let i = 0; i < 20; i++) {
      ekf.predict(0, 0, 0, 0.1, i + 1)
      ekf.updateGps(100, 100, 5, 0, i + 1)
    }
    expect(ekf.getDiagnostics().adaNoiseScale).toBeGreaterThan(1.5)
  })

  it('should track IMU-only stationary correctly', () => {
    const ekf = new SrEkf({
      processNoise: { position: 0, velocity: 0, heading: 0, accelBias: 0, gyroBias: 0 }
    })
    ekf.reset(100, 200, 0, 0)
    for (let i = 0; i < 50; i++) {
      ekf.predict(0, 0, 0, 0.01, i)
    }
    const s = ekf.getState()
    expect(s.x).toBeCloseTo(100, 8)
    expect(s.y).toBeCloseTo(200, 8)
    expect(s.v).toBeCloseTo(0, 8)
  })
})
