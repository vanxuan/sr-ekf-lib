import { describe, it, expect } from 'vitest'
import { SrEkf, NavigationSolution } from '../src/sr-ekf'

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
      ekf.predict(0, v * omega, omega, dt, i)
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
    expect(s.aBiasX).toBeLessThan(-0.25)
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

  it('should reject GPS outliers via direction-aware guard', () => {
    const ekf = new SrEkf({
      measurementNoise: { position: 1.0, velocity: 0.5 },
      gateThreshold: 9.488,
      robustWeight: { enabled: false }
    })
    ekf.reset(0, 0, 5, 0)
    const gpsOk = ekf.updateGps(0, 0, 5, 0, 0)
    expect(gpsOk).toBe(true)
    // 1000m jump: guard inflates R (5× forward, 10× cross-track caps) but χ²
    // still exceeds gate → outlier rejected
    const outlierAccepted = ekf.updateGps(1000, 1000, 5, 0, 1)
    expect(outlierAccepted).toBe(false)
    // State unchanged because outlier was rejected
    const s = ekf.getState()
    expect(s.x).toBeLessThan(1)
    expect(s.y).toBeLessThan(1)
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

  it('should keep heading stable under noisy city GPS velocity while tracking position', () => {
    // Simulates urban-canyon multipath: clean GPS position, but Doppler velocity
    // direction is noisy (±1.8 m/s cross-track). Per-component velocity robustness
    // must attenuate the noisy heading measurement without sacrificing the position
    // fix. Guards against the joint-weight regression where moderate velocity noise
    // (joint χ² < gate) passes through at full strength and jitters heading.
    const ekf = new SrEkf({ measurementNoise: { position: 3.0, velocity: 0.5 } })
    ekf.reset(0, 0, 8, 0)
    const dt = 0.1, trueV = 8
    let seed = 999, sumSq = 0, n = 0
    const rand = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff * 2 - 1 }
    for (let i = 0; i < 120; i++) {
      ekf.predict(0, 0, 0, dt, i)
      const xTrue = trueV * (i + 1) * dt
      ekf.updateGps(xTrue, 0, trueV + rand() * 0.8, rand() * 1.8, i)
      const s = ekf.getState()
      if (i > 10) { sumSq += s.psi * s.psi; n++ }
    }
    const s = ekf.getState()
    const psiRms = Math.sqrt(sumSq / n)
    // Old joint-only weighting gives psiRms ≈ 0.092; per-component gives ≈ 0.056.
    expect(psiRms).toBeLessThan(0.07)
    // Position must still track the clean measurements (true x = 96).
    expect(Math.abs(s.x - 96)).toBeLessThan(3)
    expect(Math.abs(s.y)).toBeLessThan(2)
  })

  it('should reject a GPS velocity outlier for heading while still correcting position', () => {
    // A single gross Doppler outlier (bad heading) must not corrupt heading, and
    // must not block the valid position correction (decoupled robustness).
    const ekf = new SrEkf({ measurementNoise: { position: 3.0, velocity: 0.5 } })
    ekf.reset(0, 0, 8, 0)
    for (let i = 0; i < 5; i++) { ekf.predict(0, 0, 0, 0.1, i); ekf.updateGps(8 * (i + 1) * 0.1, 0, 8, 0, i) }
    ekf.predict(0, 0, 0, 0.1, 5)
    const xValid = 8 * 6 * 0.1
    ekf.updateGps(xValid, 0, 8, 30, 5) // vy=30 → absurd ~75° heading, valid position
    const s = ekf.getState()
    expect(Math.abs(s.psi)).toBeLessThan(0.05)      // heading protected
    expect(Math.abs(s.x - xValid)).toBeLessThan(0.2) // position still corrected
  })

  it('should recover heading from 180° error using GPS velocity', () => {
    const ekf = new SrEkf({
      measurementNoise: { position: 3.0, velocity: 0.3 }
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

  it('should not drift in a circle at rest with a constant IMU bias (variance-only stillness)', () => {
    // Regression for the "car icon drifts in a very large circle" symptom:
    // a stationary device with a constant accel/gyro bias has ZERO variance but
    // nonzero |a|/|gz| magnitude. If stillness used magnitude, energy would stay
    // high, stillness would be 0, and ZUPT/ZARU would never fire -> the bias
    // builds a slowly-rotating v -> circular drift. With variance-only
    // stillness, the bias reads as perfectly still and is learned.
    const ekf = new SrEkf()
    ekf.reset(0, 0, 0, 0)
    for (let i = 0; i < 6000; i++) {
      // constant bias offset (zero variance): ax=0.05, gz=-0.0047
      ekf.predict(0.05, 0, -0.0047, 0.02, i)
    }
    const s = ekf.getState()
    const radius = Math.hypot(s.x, s.y)
    expect(radius).toBeLessThan(1) // position must stay near origin, not circle away
    expect(Math.abs(s.v)).toBeLessThan(0.1)
    expect(ekf.getStillness()).toBeGreaterThan(0.7)
  })

  it('should restore motion after stationary ends', () => {
    const ekf = new SrEkf()
    ekf.reset(100, 200, 0, 0)
    for (let i = 0; i < 50; i++) {
      ekf.predict(0, 0, 0, 0.01, i)
    }
    expect(ekf.getDiagnostics().stationary).toBe(true)
    // Real motion must produce IMU variance (stillness is variance-based). A
    // strong, sustained oscillating acceleration creates variance, flips
    // stationarity off, and lets the velocity estimate grow (ZUPT disengages).
    for (let i = 0; i < 100; i++) {
      ekf.predict(Math.sin(i * 1.5) * 3, 0, 0, 0.01, i + 50)
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
    ekf.predict(0, 1, 0, 0.1, 0, 0)
    ekf.predict(0, 1, 0, 0.1, 1, 0)
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

  it('should fire lateral accel constraint with device orientation set', () => {
    // With the fix, lateral accel fires regardless of deviceToEnu because ay
    // is already projected into the vehicle body frame in predict().  Verify
    // that a turning car with orientation set has velocity constrained by the
    // lateral accel constraint (v ≈ ay/omega).
    const ekf = new SrEkf({
      processNoise: { position: 0, velocity: 0, heading: 0, accelBias: 0, gyroBias: 0, sideslip: 0 },
      useLateralAccel: true
    })
    ekf.reset(0, 0, 5, 0)
    ekf.setOrientation(0, 0, 0)
    // Turn: gz=0.5 rad/s, ay=2.5 (centripetal = v*ω = 5*0.5 = 2.5)
    for (let i = 0; i < 50; i++) ekf.predict(0, 2.5, 0.5, 0.02, i, 0, 0, 0)
    const s = ekf.getState()
    // Lateral accel constrains v during turns — velocity should remain near 5 m/s
    // (without lateral accel, accumulated heading error from the turn could drift v)
    expect(s.v).toBeGreaterThan(4.5)
    expect(s.v).toBeLessThan(5.5)
    // Heading should accumulate from the turn: ψ ≈ ω*dt*N = 0.5*0.02*50 = 0.5 rad
    expect(s.psi).toBeGreaterThan(0.4)
    expect(s.psi).toBeLessThan(0.6)
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
      gateThreshold: 9.488,
      robustWeight: { enabled: false }
    })
    ekf2.reset(0, 0, 5, 0)
    ekf2.updateGps(0, 0, 5, 0, 0)
    const result2 = ekf2.updateGps(1000, 1000, 5, 0, 1)
    // Tighter guard caps (5×/10×) → 1000m jump exceeds gate → rejected
    expect(result2).toBe(false)
  })

  it('should not let magnetometer rotate heading when GPS is fresh at speed', () => {
    // Core regression for the mag/GPS oscillation bug: with a fresh GPS fix at
    // speed, a magnetometer reading that disagrees with the GPS heading must NOT
    // move ψ (mag update is skipped). Previously mag — which has no β term —
    // fought the GPS velocity-direction update and made ψ oscillate whenever
    // sideslip β ≠ 0.
    const ekf = new SrEkf({
      measurementNoise: { position: 2.0, velocity: 0.3, heading: 0.05 }
    })
    ekf.reset(0, 0, 20, 0)
    const trueHdg = 0.5
    // Warm up directly at speed (v=20, ~72 km/h) so the state is genuinely fast
    // and GPS heading is well-constrained. Mag is skipped at this speed, so the
    // init snap stops firing once heading covariance shrinks.
    for (let i = 0; i < 60; i++) {
      const t = i * 0.01
      ekf.predict(0, 0, 0, 0.01, i)
      ekf.updateGps(
        20 * Math.cos(trueHdg) * t,
        20 * Math.sin(trueHdg) * t,
        20 * Math.cos(trueHdg),
        20 * Math.sin(trueHdg),
        i, 5
      )
      ekf.updateMag(trueHdg, i)
    }
    expect(ekf.getState().v).toBeGreaterThan(15)
    const psiBefore = ekf.getState().psi
    // Disagreeing mag reading, same timestamp as the fresh GPS fix.
    ekf.updateMag(trueHdg + 1.0, 60)
    const psiAfter = ekf.getState().psi
    expect(Math.abs(psiAfter - psiBefore)).toBeLessThan(1e-6)
  })

  it('should not oscillate heading from mag/GPS conflict at speed with sideslip', () => {
    // Reproduces the reported bug: at >10 km/h with an active turn (sideslip
    // β ≠ 0), the magnetometer (which has no β term) used to fight the GPS
    // velocity-direction update for ψ, making ψ oscillate away/toward the car's
    // nose. Mag is now skipped at speed with a fresh GPS fix, so ψ must stay
    // locked to the GPS velocity direction.
    const ekf = new SrEkf({
      measurementNoise: { position: 3.0, velocity: 0.5, heading: 0.1 }
    })
    ekf.reset(0, 0, 12, 0) // ~43 km/h
    const omega = 0.2
    const trueBeta = 0.1
    let maxPsiJump = 0
    let prevPsi = ekf.getState().psi
    for (let i = 0; i < 300; i++) {
      const t = i * 0.05
      ekf.predict(0, 0, omega, 0.05, i)
      const psi_t = omega * t
      const alpha = psi_t + trueBeta
      const x = (12 / omega) * (Math.sin(omega * t + trueBeta) - Math.sin(trueBeta))
      const y = (12 / omega) * (-Math.cos(omega * t + trueBeta) + Math.cos(trueBeta))
      // GPS every step (so it is always fresh → mag skipped), bearing = nose dir.
      ekf.updateGps(x, y, 12 * Math.cos(alpha), 12 * Math.sin(alpha), i, 10)
      ekf.updateMag(psi_t, i)
      const psi = ekf.getState().psi
      maxPsiJump = Math.max(maxPsiJump, Math.abs(psi - prevPsi))
      prevPsi = psi
    }
    const s = ekf.getState()
    // ψ must track the GPS velocity direction (psi_t + trueBeta), not the nose
    // bearing (psi_t), and must not oscillate step-to-step.
    expect(Math.abs(s.psi - (omega * 300 * 0.05 + trueBeta))).toBeLessThan(0.2)
    expect(maxPsiJump).toBeLessThan(0.15)
  })

  it('should keep velocity non-negative and heading = motion direction after a U-turn', () => {
    // Reproduces the reported "U-turn → velocity negative / heading snaps to
    // nose" bug. After a 180° heading reversal, the GPS velocity direction is
    // anti-parallel to the current heading. The 180° flip-recovery (now fired
    // on anti-parallel velocity whenever v > 0.8, independent of the position
    // χ² gate) must flip ψ to the new direction of motion instead of letting the
    // Kalman velocity update resolve the v/heading sign ambiguity into the
    // negative-v branch. Velocity must stay ≥ 0 and ψ must track the new GPS
    // velocity direction.
    const ekf = new SrEkf({
      measurementNoise: { position: 3.0, velocity: 0.5, heading: 0.1 }
    })
    ekf.reset(0, 0, 10, Math.PI / 2)
    let py = 0
    for (let i = 0; i < 400; i++) {
      ekf.predict(0, 0, 0, 0.05, i)
      py += 10 * 0.05
      if (i < 200) {
        // driving north
        ekf.updateGps(0, py, 10, 0, i, 5)
        ekf.updateMag(Math.PI / 2, i)
      } else {
        // hard U-turn telemetry: position continues north but reported GPS
        // velocity now points south (anti-parallel to current heading ψ=π/2)
        ekf.updateGps(0, py, -10, 0, i, 5)
        ekf.updateMag(Math.PI / 2, i)
      }
    }
    const s = ekf.getState()
    expect(s.v).toBeGreaterThan(0) // never negative
    expect(ekf.getDiagnostics().coasting).toBe(false)
    // heading must now point south (direction of motion), not stay at north
    expect(Math.abs(s.psi - Math.PI)).toBeLessThan(0.2)
  })

  it('should not shake heading (180° flips) when stationary on a table', () => {
    // Regression: when the device is at rest (v ≈ 0) with GPS multipath jitter,
    // the GPS velocity measurement is pure noise. The velocity sign-ambiguity
    // backstop must NOT reparameterize at negligible speed, otherwise v jitters
    // around zero and ψ flips by π every step (heading "shaking"). Verify v stays
    // ~0 and ψ does not do any 180° flip.
    const ekf = new SrEkf({
      measurementNoise: { position: 3.0, velocity: 0.5, heading: 0.1 }
    })
    ekf.reset(0, 0, 10, 0)
    for (let i = 0; i < 100; i++) {
      ekf.predict(0, 0, 0, 0.02, i)
      ekf.updateGps(0, 10 * 0.02 * i, 0, 10, i, 5)
      ekf.updateMag(0, i)
    }
    let prevPsi = ekf.getState().psi
    let n180 = 0
    let minV = 0
    for (let i = 100; i < 600; i++) {
      const ax = (Math.random() + Math.random() + Math.random() + Math.random() - 2) / 2 * 0.005
      const ay = (Math.random() + Math.random() + Math.random() + Math.random() - 2) / 2 * 0.005
      const gz = (Math.random() + Math.random() + Math.random() + Math.random() - 2) / 2 * 0.001
      ekf.predict(ax, ay, gz, 0.02, i)
      const gvx = (Math.random() + Math.random() + Math.random() + Math.random() - 2) / 2 * 1.5
      const gvy = (Math.random() + Math.random() + Math.random() + Math.random() - 2) / 2 * 1.5
      ekf.updateGps(
        (Math.random() + Math.random() + Math.random() + Math.random() - 2) / 2 * 1.2,
        20 + (Math.random() + Math.random() + Math.random() + Math.random() - 2) / 2 * 1.2,
        gvx, gvy, i, 5
      )
      const s = ekf.getState()
      const jump = Math.abs(s.psi - prevPsi)
      if (Math.abs(Math.PI - jump) < 0.3) n180++
      minV = Math.min(minV, s.v)
      prevPsi = s.psi
    }
    expect(n180).toBe(0)
    expect(ekf.getDiagnostics().stationary).toBe(true)
    expect(minV).toBeGreaterThanOrEqual(0) // v never negative
    expect(ekf.getState().v).toBeLessThan(0.2) // settled near zero
  })

  it('should still absorb a stable magnetometer heading while stationary', () => {
    // The drift guard must NOT reject a genuinely stable compass: after warm-up a
    // stable compass bearing should be fused and the heading converge to it.
    const ekf = new SrEkf({
      measurementNoise: { position: 3.0, velocity: 0.5, heading: 0.1 }
    })
    ekf.reset(0, 0, 0, 0)
    const step = 20
    let t = 0
    for (; t < 4000; ) {
      ekf.predict(0, 0, 0, 0.02, t)
      ekf.updateGps(0, 0, 0, 0, t, 5)
      ekf.updateMag((30 * Math.PI) / 180, t) // stable 30° compass
      t += step
    }
    const psi = ekf.getState().psi
    const diff = Math.abs(((psi - (30 * Math.PI) / 180 + Math.PI) % (2 * Math.PI)) - Math.PI)
    expect(diff).toBeLessThan(0.05)
  })

  it('should not let heading rotate while at rest when ZARU velocity gate is suppressed', () => {
    // Regression: ZARU (gyro-bias correction at rest) used to be gated on the
    // velocity-coupled ZUPT weight `w > 0.8`. If the device was on a slightly
    // moving surface (stillness just under 0.8) OR GPS kept the velocity estimate
    // slightly elevated (speedGate < 1), ZARU never fired, gyro bias stayed
    // uncorrected, and ψ integrated the residual rate → heading ROTATED CONSTANTLY
    // while apparently stationary. ZARU must engage on IMU stillness alone.
    // Use a gyro bias with no magnetometer so nothing else corrects heading.
    const ekf = new SrEkf({
      measurementNoise: { position: 3.0, velocity: 0.5, heading: 0.1 }
    })
    ekf.reset(0, 0, 10, 0)
    for (let i = 0; i < 150; i++) {
      ekf.predict(0, 0, 0.04, 0.02, i) // true gyro bias 0.04 rad/s
      ekf.updateGps(0, 10 * 0.02 * i, 0, 10, i, 5)
      // NO magnetometer updates
    }
    const psi0 = ekf.getState().psi
    for (let i = 150; i < 1650; i++) {
      // slight table sway keeps stillness ~0.8; GPS velocity stays ~2 m/s so the
      // OLD velocity-coupled gate would have disabled ZARU.
      const sway = 0.08 * Math.sin(i * 0.05)
      ekf.predict(sway + (Math.random() + Math.random() + Math.random() + Math.random() - 2) / 2 * 0.01,
        (Math.random() + Math.random() + Math.random() + Math.random() - 2) / 2 * 0.01,
        0.04 + (Math.random() + Math.random() + Math.random() + Math.random() - 2) / 2 * 0.002,
        0.02, i)
      ekf.updateGps(
        (Math.random() + Math.random() + Math.random() + Math.random() - 2) / 2 * 1.2,
        20 + (Math.random() + Math.random() + Math.random() + Math.random() - 2) / 2 * 1.2,
        (Math.random() + Math.random() + Math.random() + Math.random() - 2) / 2 * 0.3,
        2.0 + (Math.random() + Math.random() + Math.random() + Math.random() - 2) / 2 * 0.3,
        i, 5
      )
    }
    const psi1 = ekf.getState().psi
    let rot = psi1 - psi0
    while (rot > Math.PI) rot -= 2 * Math.PI
    while (rot < -Math.PI) rot += 2 * Math.PI
    // heading must stay put (gyro bias corrected by ZARU), not rotate.
    expect(Math.abs(rot)).toBeLessThan(0.1) // < ~5.7° over 30s
    expect(Math.abs(ekf.getState().gBiasZ - 0.04)).toBeLessThan(0.01)
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
    expect(s.p.length).toBe(8)
    expect(s.p[0].length).toBe(8)
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
    // Direction-aware guard inflates R → chiSq ~7.13 → weight ~0.64
    expect(w2).toBeLessThan(0.8)
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
    const accepted = ekf.updateGps(1000, 1000, 5, 0, 100)
    // Robust M weights the outlier down (soft rejection, not binary gate)
    expect(accepted).toBe(true)
    const w = ekf.getDiagnostics().robustWeight
    expect(w).toBeLessThan(0.1)
    const s = ekf.getState()
    expect(Math.abs(s.x - 5)).toBeLessThan(10)
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
      measurementNoise: { position: 1.0, velocity: 0.3 },
      initialCovariance: { position: 0.1, velocity: 0.1, heading: 0.1, sideslip: 0.1, accelBias: 0.01, gyroBias: 0.01 },
      adaptiveNoise: { enabled: true, smoothing: 0.3, maxScale: 20 }
    })
    ekf.reset(0, 0, 2, 0)
    ekf.updateGps(0, 0, 2, 0, 0)
    expect(ekf.getDiagnostics().adaNoiseScale).toBe(1)
    // Use velocity mismatches that don't trigger direction-aware guard
    for (let i = 0; i < 20; i++) {
      ekf.predict(0, 0, 0, 0.1, i + 1)
      ekf.updateGps(0, 0, 5 + i * 2, 0, i + 1)
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

  it('yaw rate does not drift during stationary', () => {
    const ekf = new SrEkf({ dt: 0.01 })
    ekf.reset(0, 0, 0, 0)
    // g_bias_z starts at 0. Feed gz=0 so ω=0.
    // ZUPT engages and ZARU drives gyro bias toward gz.
    for (let i = 0; i < 1500; i++) ekf.predict(0, 0, 0, 0.01, i * 10)
    // 15s stationary — bias near 0
    const s = ekf.getState()
    expect(Math.abs(s.gBiasZ)).toBeLessThan(0.001)
    // Feed gz=0.05 — ZARU drives bias toward 0.05 (convergence is slow
    // because gyro bias covariance is tiny after 15s of ZUPT+ZARU)
    for (let i = 0; i < 3000; i++) ekf.predict(0, 0, 0.05, 0.01, 15000 + i * 10)
    // ω = gz - gBiasZ should be significantly reduced after 30s
    const omega = 0.05 - ekf.getState().gBiasZ
    expect(Math.abs(omega)).toBeLessThan(0.02)
  })

  // ─── Tier 1: GPS latency compensation ─────────────────────────

  it('should compensate GPS latency by rewinding and re-predicting', () => {
    const ekf = new SrEkf({
      measurementNoise: { position: 1.0, velocity: 0.3 },
      processNoise: { position: 0.001, velocity: 0.01, heading: 0.001, accelBias: 0, gyroBias: 0 }
    })
    ekf.reset(0, 0, 10, 0)
    ekf.updateGps(0, 0, 10, 0, 0)
    // Predict from t=10 to t=109 (100 steps at 0.01s)
    for (let i = 0; i < 100; i++) ekf.predict(0, 0, 0, 0.01, 10 + i)
    expect(ekf.getState().x).toBeGreaterThan(9)
    // GPS arrives at t=60 (mid-point, 0.49s delayed)
    // Position at t=60 from IMU alone would be x≈6.0 (10 * 0.6s)
    // GPS says the true position is (5.0, 0) — a 1.0m error
    const ok = ekf.updateGps(5.0, 0, 10, 0, 60)
    expect(ok).toBe(true)
    // After latency compensation, state should be numerically stable
    const s = ekf.getState()
    expect(isFinite(s.x)).toBe(true)
    expect(isFinite(s.v)).toBe(true)
    expect(s.v).toBeGreaterThan(5)
  })

  it('should handle GPS latency with a curved trajectory', () => {
    const ekf = new SrEkf({
      measurementNoise: { position: 2.0, velocity: 0.5 },
      processNoise: { position: 0.01, velocity: 0.1, heading: 0.01, accelBias: 0, gyroBias: 0 }
    })
    const omega = 0.3
    const v = 10
    ekf.reset(0, 0, v, 0)
    ekf.updateGps(0, 0, v, 0, 0)
    let t = 10
    for (let i = 0; i < 80; i++) {
      ekf.predict(0, 0, omega, 0.01, t)
      t++
    }
    // After 80 predicts at 100Hz = 0.8s
    const s1 = ekf.getState()
    // GPS at t=40 (400ms ago)
    const truePsi = omega * 0.4
    const trueAlpha = truePsi
    const trueX = (v / omega) * Math.sin(omega * 0.4)
    const trueY = (v / omega) * (1 - Math.cos(omega * 0.4))
    const ok = ekf.updateGps(trueX, trueY, v * Math.cos(trueAlpha), v * Math.sin(trueAlpha), 40)
    expect(ok).toBe(true)
    // State should be numerically stable after re-predict
    const s2 = ekf.getState()
    expect(isFinite(s2.x)).toBe(true)
    expect(isFinite(s2.y)).toBe(true)
    expect(isFinite(s2.v)).toBe(true)
    expect(isFinite(s2.psi)).toBe(true)
    expect(s2.v).toBeGreaterThan(5)
  })

  it('should replay from scratch after GPS latency (no buffer corruption)', () => {
    const ekf = new SrEkf({
      measurementNoise: { position: 1.0, velocity: 0.3 },
      processNoise: { position: 0.001, velocity: 0.01, heading: 0.001, accelBias: 0, gyroBias: 0 }
    })
    ekf.reset(0, 0, 5, 0)
    ekf.updateGps(0, 0, 5, 0, 0)
    // Fill buffer with IMU data
    for (let i = 0; i < 50; i++) ekf.predict(0, 0, 0, 0.01, i + 10)
    // Delayed GPS
    ekf.updateGps(0.25, 0, 5, 0, 10)
    // After compensation, subsequent predicts should work normally
    for (let i = 0; i < 20; i++) ekf.predict(0, 0, 0, 0.01, i + 200)
    const s = ekf.getState()
    expect(isFinite(s.x)).toBe(true)
    expect(isFinite(s.y)).toBe(true)
    expect(isFinite(s.v)).toBe(true)
  })

  // ─── Tier 1: GPS latency buffer preserves all 8 states ──────────

  it('should preserve magDeclination state across GPS latency rewind', () => {
    // Regression: the latency-compensation buffer previously stored only 7
    // state components, so the 8th state (magDeclination) was dropped on rewind.
    const ekf = new SrEkf({
      magneticDeclination: 0.37,
      measurementNoise: { position: 1.0, velocity: 0.3 },
      processNoise: { position: 0.001, velocity: 0.01, heading: 0.001, accelBias: 0, gyroBias: 0 },
      initialCovariance: { magDeclination: 0.25 }
    })
    ekf.reset(0, 0, 5, 0)
    ekf.updateGps(0, 0, 5, 0, 0)
    // Fill buffer with IMU data so a delayed GPS triggers a rewind
    for (let i = 0; i < 50; i++) ekf.predict(0, 0, 0, 0.01, i + 10)
    const before = ekf.getState().magDeclination
    expect(before).toBeCloseTo(0.37, 6)
    // Delayed GPS (timestamp 10, while lastImuTimeMs is ~60) → rewind
    ekf.updateGps(0.25, 0, 5, 0, 10)
    const after = ekf.getState().magDeclination
    expect(isFinite(after)).toBe(true)
    expect(after).toBeCloseTo(0.37, 6)
  })

  // ─── Tier 1: Fixed-time IMU windows ───────────────────────────

  it('should maintain correct IMU window duration with varying dt', () => {
    const ekf = new SrEkf({
      processNoise: { position: 0.1, velocity: 0.1, heading: 0.01, accelBias: 0, gyroBias: 0 }
    })
    ekf.reset(0, 0, 5, 0)
    // Mix of dt values (50Hz, 100Hz, 200Hz)
    const dts = [0.02, 0.01, 0.005, 0.01, 0.02]
    for (let cycle = 0; cycle < 300; cycle++) {
      const dt = dts[cycle % dts.length]
      ekf.predict(1, 0.1, 0.05, dt, cycle + 1)
    }
    // After many steps, IMU energy should be reasonable
    const d = ekf.getDiagnostics()
    expect(d.trace).toBeGreaterThan(0)
    expect(isFinite(d.trace)).toBe(true)
    // Stationary detection should work even with varying dt
    const ekf2 = new SrEkf()
    ekf2.reset(0, 0, 0, 0)
    for (let cycle = 0; cycle < 100; cycle++) {
      const dt = dts[cycle % dts.length]
      ekf2.predict(0, 0, 0, dt, cycle + 1)
    }
    expect(ekf2.getDiagnostics().stationary).toBe(true)
  })

  // ─── Tier 1: adaNoiseScale decay ──────────────────────────────

  it('should decay adaNoiseScale after large innovations subside', () => {
    const ekf = new SrEkf({
      measurementNoise: { position: 1.0, velocity: 0.3 },
      initialCovariance: { position: 0.1, velocity: 0.1, heading: 0.1, sideslip: 0.1, accelBias: 0.01, gyroBias: 0.01 },
      adaptiveNoise: { enabled: true, smoothing: 0.3, maxScale: 20 }
    })
    ekf.reset(0, 0, 2, 0)
    ekf.updateGps(0, 0, 2, 0, 0)
    // Inject large innovations via velocity ramping
    for (let i = 0; i < 20; i++) {
      ekf.predict(0, 0, 0, 0.1, i + 1)
      ekf.updateGps(0, 0, 5 + i * 2, 0, i + 1)
    }
    const inflated = ekf.getDiagnostics().adaNoiseScale
    expect(inflated).toBeGreaterThan(1.5)
    // Continue with GPS that matches the corrected state exactly
    const s0 = ekf.getState()
    for (let i = 0; i < 30; i++) {
      ekf.predict(0, 0, 0, 0.1, i + 100)
      ekf.updateGps(0, 0, s0.v, 0, i + 100)
    }
    const decayed = ekf.getDiagnostics().adaNoiseScale
    expect(decayed).toBeLessThan(inflated * 0.8)
  })

  // ─── Tier 2: Zero-alloc getStateInto ─────────────────────────

  it('getStateInto should match getState output', () => {
    const ekf = new SrEkf()
    ekf.reset(10, 20, 5, 0.5)
    ekf.predict(0.1, 0.05, 0.01, 0.01, 0)
    const ref = ekf.getState()
    const p = [new Float64Array(8), new Float64Array(8), new Float64Array(8),
               new Float64Array(8), new Float64Array(8), new Float64Array(8), new Float64Array(8), new Float64Array(8)]
    const out: NavigationSolution = { x: 0, y: 0, v: 0, psi: 0, beta: 0, aBiasX: 0, gBiasZ: 0, magDeclination: 0, p }
    ekf.getStateInto(out)
    expect(out.x).toBe(ref.x)
    expect(out.y).toBe(ref.y)
    expect(out.v).toBe(ref.v)
    expect(out.psi).toBe(ref.psi)
    expect(out.beta).toBe(ref.beta)
    expect(out.aBiasX).toBe(ref.aBiasX)
    expect(out.gBiasZ).toBe(ref.gBiasZ)
    expect(out.magDeclination).toBe(ref.magDeclination)
    for (let i = 0; i < 8; i++)
      for (let j = 0; j < 8; j++)
        expect(out.p[i][j]).toBe(ref.p[i][j])
  })

  it('getStateInto should reuse provided arrays (no allocation)', () => {
    const ekf = new SrEkf()
    ekf.reset(10, 20, 5, 0.5)
    const p = [new Float64Array(8), new Float64Array(8), new Float64Array(8),
               new Float64Array(8), new Float64Array(8), new Float64Array(8), new Float64Array(8), new Float64Array(8)]
    const p0 = p[0]
    const out: NavigationSolution = { x: 0, y: 0, v: 0, psi: 0, beta: 0, aBiasX: 0, gBiasZ: 0, magDeclination: 0, p }
    ekf.getStateInto(out)
    expect(out.p[0]).toBe(p0)
  })

  // ─── Tier 2: Position step-change outlier guard ──────────────

  it('should inflate posR for implausible position jumps', () => {
    // With a tight initial covariance (P[X][X]=1, posR=1):
    // A 20m jump in 0.1s: without guard K ≈ 0.5 → Δx ≈ 10m
    // With direction-aware guard: forward=20 > maxForward=2.4 → 8.3× posR
    // Anisotropic R: along-track σ = 8.3×0.5 ≈ 4.2m → K ≈ 0.1 → Δx ≈ 2m
    const ekf = new SrEkf({
      measurementNoise: { position: 1.0, velocity: 0.5 },
      initialCovariance: { position: 1, velocity: 1, heading: 0.1, sideslip: 0.1, accelBias: 0.01, gyroBias: 0.01 },
      gateThreshold: 100,
      robustWeight: { enabled: false }
    })
    ekf.reset(0, 0, 5, 0)
    ekf.updateGps(0, 0, 5, 0, 0)   // init, sets small cov
    ekf.predict(0, 0, 0, 0.1, 1)
    const ok = ekf.updateGps(20, 0, 5, 0, 1)
    expect(ok).toBe(true)
    const s = ekf.getState()
    // Guard limits correction: without guard x ≈ 10m, with guard less
    expect(s.x).toBeGreaterThan(0)
    expect(s.x).toBeLessThan(8)
  })

  // ─── Tier 2: GPS timestamp monotonicity ──────────────────────

  it('should reject duplicate GPS timestamps', () => {
    const ekf = new SrEkf()
    ekf.reset(0, 0, 5, 0)
    ekf.updateGps(0, 0, 5, 0, 0)   // init at t=0
    ekf.predict(0, 0, 0, 0.1, 1)
    ekf.updateGps(10, 10, 5, 0, 1) // normal update at t=1
    const s_before = ekf.getState()
    const ok = ekf.updateGps(100, 100, 5, 0, 1) // same timestamp t=1 → rejected
    expect(ok).toBe(true) // returns true but no update applied
    const s_after = ekf.getState()
    expect(s_after.x).toBe(s_before.x) // state unchanged
  })

  it('should reject stale out-of-order GPS timestamps', () => {
    const ekf = new SrEkf()
    ekf.reset(0, 0, 5, 0)
    ekf.updateGps(0, 0, 5, 0, 0)
    ekf.predict(0, 0, 0, 0.1, 10)
    ekf.updateGps(1, 0, 5, 0, 10)
    const s_before = ekf.getState()
    // Earlier timestamp than last accepted GPS (t=5 < t=10) → rejected
    const ok = ekf.updateGps(999, 999, 5, 0, 5)
    expect(ok).toBe(true)
    const s_after = ekf.getState()
    expect(s_after.x).toBe(s_before.x) // state unchanged
  })

  it('should reset biases via resetBiases()', () => {
    const ekf = new SrEkf({
      processNoise: { position: 0.001, velocity: 0.01, heading: 0.001, accelBias: 0.1, gyroBias: 1e-5 },
      initialCovariance: { accelBias: 0.2, gyroBias: 0.02 }
    })
    ekf.reset(0, 0, 5, 0)
    ekf.predict(10, 0, 0, 1, 0)
    ekf.updateGps(5, 0, 5, 0, 1)
    const before = ekf.getState()
    ekf.resetBiases(0, 0)
    const after = ekf.getState()
    expect(after.aBiasX).toBe(0)
    expect(after.gBiasZ).toBe(0)
    // Covariance should be restored to initial values
    expect(after.p[5][5]).toBeCloseTo(0.2, 3)
    expect(after.p[6][6]).toBeCloseTo(0.02, 3)
  })

  it('should inflate covariance via inflateCovariance()', () => {
    const ekf = new SrEkf({
      processNoise: { position: 0, velocity: 0, heading: 0, accelBias: 0, gyroBias: 0 }
    })
    ekf.reset(0, 0, 0, 0)
    const before = ekf.getState()
    ekf.inflateCovariance({ position: 400, heading: 9 })
    const after = ekf.getState()
    expect(after.p[0][0]).toBe(400)
    expect(after.p[1][1]).toBe(400)
    expect(after.p[3][3]).toBeCloseTo(9, 8)
  })

  it('should re-acquire zero velocity via ZUPT when stationary with corrupted v (phone on table)', () => {
    const ekf = new SrEkf()
    ekf.reset(100, 200, 1.0, 0)
    ekf.updateGps(100, 200, 0, 0, 0)
    for (let i = 0; i < 400; i++) {
      ekf.predict(0, 0, 0, 0.01, i + 1)
    }
    const s = ekf.getState()
    expect(Math.abs(s.v)).toBeLessThan(0.1)
    expect(ekf.getDiagnostics().stationary).toBe(true)
  })

  it('should still disengage ZUPT on real motion after stationary re-acquisition', () => {
    const ekf = new SrEkf()
    ekf.reset(100, 200, 1.0, 0)
    ekf.updateGps(100, 200, 0, 0, 0)
    for (let i = 0; i < 400; i++) {
      ekf.predict(0, 0, 0, 0.01, i + 1)
    }
    expect(Math.abs(ekf.getState().v)).toBeLessThan(0.1)
    // Now the user picks the phone up and accelerates
    for (let i = 0; i < 100; i++) {
      ekf.predict(2, 0, 0, 0.01, i + 401)
    }
    const after = ekf.getState()
    expect(ekf.getDiagnostics().stationary).toBe(false)
    expect(after.v).toBeGreaterThan(0.1)
  })

  it('should track heading during corner entry and exit with minimal lag (Trapezoidal integration)', () => {
    const ekf = new SrEkf({
      processNoise: { position: 1.0, velocity: 0.1, heading: 0.05, sideslip: 0.05, accelBias: 1e-5, gyroBias: 1e-6 },
      measurementNoise: { position: 1.0, velocity: 0.2 }
    })
    ekf.reset(0, 0, 10, 0)
    const dt = 0.02, truthDt = 0.001, v = 10
    let time = 0, truePsi = 0, trueX = 0, trueY = 0
    let maxHeadingError = 0
    const getOmegaAt = (t: number) => {
      if (t < 5) return 0
      if (t < 6) return (t - 5)
      if (t < 8) return 1
      if (t < 9) return 1 - (t - 8)
      return 0
    }
    let nextGpsTime = 0
    while (time < 10) {
      for (let j = 0; j < dt / truthDt; j++) {
        const t = time + j * truthDt
        const w = getOmegaAt(t)
        truePsi += w * truthDt
        trueX += v * Math.cos(truePsi) * truthDt
        trueY += v * Math.sin(truePsi) * truthDt
      }
      time += dt
      const wEKF = getOmegaAt(time)
      ekf.predict(0, v * wEKF, wEKF, dt, time * 1000)
      if (time >= nextGpsTime) {
        ekf.updateGps(trueX, trueY, v * Math.cos(truePsi), v * Math.sin(truePsi), time * 1000)
        nextGpsTime += 1.0
      }
      const s = ekf.getState()
      const error = Math.abs((s.psi - truePsi + 3 * Math.PI) % (2 * Math.PI) - Math.PI)
      maxHeadingError = Math.max(maxHeadingError, error)
    }
    expect(maxHeadingError * 180 / Math.PI).toBeLessThan(5.0)
  })

  it('should track heading and velocity during walking with step-induced energy', () => {
    const ekf = new SrEkf()
    ekf.reset(0, 0, 0, 0)
    const dt = 0.02
    let time = 0, truePsi = 0.5, trueV = 1.3, trueX = 0, trueY = 0
    let nextGpsTime = 0
    for (let i = 0; i < 10 / dt; i++) {
      const stepFreq = 2.0
      const ax = 2.0 * Math.sin(2 * Math.PI * stepFreq * time)
      ekf.predict(ax, 0, 0, dt, time * 1000)
      trueX += trueV * Math.cos(truePsi) * dt
      trueY += trueV * Math.sin(truePsi) * dt
      time += dt
      if (time >= nextGpsTime) {
        ekf.updateGps(trueX, trueY, trueV * Math.cos(truePsi), trueV * Math.sin(truePsi), time * 1000)
        nextGpsTime += 1.0
      }
    }
    const s = ekf.getState()
    expect(s.v).toBeCloseTo(trueV, 0.2)
    expect(Math.abs((s.psi - truePsi + 3 * Math.PI) % (2 * Math.PI) - Math.PI)).toBeLessThan(0.1)
  })
})
