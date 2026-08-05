const N = 8, M = 4, PRE = M + N, MAG_PRE = 1 + N;
const NTRI = N * (N + 1) / 2;   // entries in lower-triangular N×N matrix

const enum I { X, Y, V, PSI, BETA, A_BIAS_X, G_BIAS_Z, MAG_DECL }

export interface EkfConfig {
  dt?: number
  processNoise?: {
    position?: number; velocity?: number; heading?: number
    sideslip?: number; accelBias?: number; gyroBias?: number; magDeclination?: number
  }
  measurementNoise?: { position?: number; velocity?: number; heading?: number }
  magneticDeclination?: number
  initialCovariance?: {
    position?: number; velocity?: number; heading?: number
    sideslip?: number; accelBias?: number; gyroBias?: number; magDeclination?: number
  }
  gateThreshold?: number; coastTimeoutMs?: number; gpsTimeOffsetMs?: number
  robustWeight?: { enabled?: boolean; type?: 'huber' | 'cauchy'; threshold?: number }
  adaptiveNoise?: { enabled?: boolean; smoothing?: number; maxScale?: number }
  adaptiveScaling?: {
    positionAccel?: number; velocityAccel?: number; velocityStep?: number
    headingGyro?: number; headingStep?: number; sideslipGyro?: number; sideslipStep?: number
  }
  useLateralAccel?: boolean
}

export interface NavigationSolution {
  x: number; y: number; v: number; psi: number; beta: number
  aBiasX: number; gBiasZ: number; magDeclination: number; p: Float64Array[]
}

export interface EkfDiagnostics {
  trace: number; gpsInnovation: Float64Array; gpsChiSq: number
  gatePassed: boolean; coasting: boolean; lastGpsTimeMs: number
  lastImuTimeMs: number; stationary: boolean; magDeclination: number
  robustWeight: number; adaNoiseScale: number
}

const DEFAULTS = {
  dt: 0.01,
  processNoise: { position: 5.0, velocity: 0.5, heading: 0.10, sideslip: 0.1, accelBias: 1e-4, gyroBias: 5e-5, magDeclination: 1e-4 },
  measurementNoise: { position: 3.0, velocity: 0.5, heading: 0.1 },
  magneticDeclination: 0,
  initialCovariance: { position: 100, velocity: 10, heading: Math.PI * Math.PI, sideslip: 0.25, accelBias: 0.1, gyroBias: 0.01, magDeclination: 0.25 },
  gateThreshold: 9.488, coastTimeoutMs: 5000, gpsTimeOffsetMs: 0,
  robustWeight: { enabled: true, type: 'huber', threshold: 9.488 },
  adaptiveNoise: { enabled: false, smoothing: 0.1, maxScale: 3 },
  adaptiveScaling: { positionAccel: 2.0, velocityAccel: 1.0, velocityStep: 2.5, headingGyro: 0.5, headingStep: 1.5, sideslipGyro: 0.5, sideslipStep: 1.8 },
  useLateralAccel: true
};

const EPS = 1e-4, TWO_PI = 2 * Math.PI;

import { matCreate, matLowerToFull, matLowerToFullInto, chol4x4, cholSolve4, ensureDiag } from './matrix';

class RingBuf {
  private buf: Float64Array;
  private head = 0;
  length = 0;
  private readonly mask: number;

  constructor(capacity: number) {
    this.buf = new Float64Array(capacity);
    this.mask = capacity - 1;
  }

  push(v: number): void {
    this.buf[this.head] = v;
    this.head = (this.head + 1) & this.mask;
    if (this.length < this.buf.length) this.length++;
  }

  shift(): number {
    if (this.length === 0) return 0;
    this.length--;
    const idx = (this.head - this.length + this.buf.length) & this.mask;
    return this.buf[idx];
  }

  get(i: number): number {
    const idx = (this.head - this.length + i + this.buf.length) & this.mask;
    return this.buf[idx];
  }

  clear(): void {
    this.head = 0;
    this.length = 0;
  }
}

export class SrEkf {
  private readonly x = new Float64Array(N);
  private readonly S: Float64Array[];
  private config: Required<EkfConfig>;
  private lastGpsTimeMs = 0;
  private lastGpsSpeed = 0;
  private smoothedSpeed?: number;
  private lastImuTimeMs = -Infinity;
  private lastChiSq = 0;
  private lastGatePassed = false;
  private coasting = false;
  private gpsInitialized = false;
  private accelEnergy = 0;
  private gyroEnergy = 0;
  private varAccelEnergy = 0;
  private varGyroEnergy = 0;
  private stepEnergy = 0;
  private _accel3DVar = 0;
  private _gyro3DEnergy = 0;
  private magRejectCount = 0;
  private _debugMagTrust = 1;
  private _debugInnovDeg = 0;
  private _debugGateThreshDeg = 0;
  private _debugMagAlpha = 0;

  private readonly tmpF = matCreate(N, N);
  private readonly tmpFS = matCreate(N, N);
  private readonly tmpSqrtQ = matCreate(N, N);
  private readonly tmpH = matCreate(M, N);
  private readonly tmpA = matCreate(M, N);
  private readonly tmpInnov = new Float64Array(M);
  private readonly tmpWork4x4 = matCreate(M, M);
  private readonly tmpHPH = matCreate(M, M);
  private readonly tmpQR = matCreate(2 * N, N);
  private readonly tmpHouseV = new Float64Array(2 * N);
  private readonly tmpZ = new Float64Array(M);
  private readonly tmpPreA = matCreate(PRE, PRE);
  private readonly tmpPreAT = matCreate(PRE, PRE);
  private readonly tmpMagHS = new Float64Array(N);
  private readonly tmpMagAT = matCreate(MAG_PRE, MAG_PRE);
  private readonly tmpW = new Float64Array(M);
  private readonly tmpBuf = new Float64Array(N);
  private readonly tmpL4x4 = matCreate(M, M);
  private readonly _tmpR3x3 = matCreate(3, 3);
  // Reusable scratch so getDiagnostics() allocates nothing on the hot path
  // (typically called every render frame on mobile). getState() keeps
  // returning a fresh covariance matrix for caller safety; use getStateInto()
  // with a caller-owned buffer for a zero-allocation fast path.
  private readonly _innovCache = new Float64Array(M);
  private readonly _tmpRdv3x3 = matCreate(3, 3);
  private readonly tmpLatHS = new Float64Array(N);
  private readonly tmpLatPre = matCreate(MAG_PRE, MAG_PRE);

  // ─── GPS latency compensation buffer ───────────────────────────
  private static readonly BUF_CAP = 256;
  private readonly bufT = new Float64Array(SrEkf.BUF_CAP);
  private readonly bufX = new Float64Array(SrEkf.BUF_CAP * N);
  private readonly bufS = new Float64Array(SrEkf.BUF_CAP * NTRI);
  private readonly bufIn = new Float64Array(SrEkf.BUF_CAP * 10);
  private readonly bufReplay = new Float64Array(SrEkf.BUF_CAP * 11);
  private bufTail = 0;
  private bufLen = 0;

  // ─── Time-tracked IMU windows (circular buffers) ────────────────
  private readonly imuTS = new RingBuf(256);
  private readonly stepTS = new RingBuf(128);

  private readonly stepBuffer = new RingBuf(128);
  private stepFreq = 0;
  private readonly axWindow = new RingBuf(256);
  private readonly ayWindow = new RingBuf(256);
  private readonly azWindow = new RingBuf(256);
  private readonly gxWindow = new RingBuf(256);
  private readonly gyWindow = new RingBuf(256);
  private readonly gzWindow = new RingBuf(256);
  private _zuptEngaged = false;
  private _wasZuptEngaged = false;
  private lastMagBearing = 0;
  private lastMagTimeMs = 0;
  private lastOmega = 0;             // bias-corrected gyro rate (rad/s) from the last predict()
  private prevOmega = 0;             // omega from the previous predict() for angular acceleration
  private angAccel = 0;              // |dω/dt| (rad/s²) — angular acceleration magnitude
  private smoothAngAccel = 0;        // EMA of angAccel to sustain boost after transients
  private lastPredictTimeMs = 0;
  private prevCallMagBearing = 0;     // compass bearing at the previous updateMag() call
  private prevCallMagTimeMs = 0;     // timestamp at the previous updateMag() call
  private gpsInitTimeMs = 0;  // Time when GPS first initialized, used for heading init boost
  private deviceToEnu: Float64Array[] | null = null;
  private deviceToVehicle: Float64Array[] | null = null;
  private curAzimuth = NaN;
  private curPitch = NaN;
  private curRoll = NaN;
  // ─── Velocity prior (coasting speed anchor) ───────────────────
  private coastSpeed = 0;            // speed at GPS loss — soft prior target during coasting
  private coastSpeedReady = false;   // true once coastSpeed has been captured
  // ─── Barometric altitude tracking ─────────────────────────────
  private lastBaroAlt = NaN;         // last barometric altitude (m)
  private lastBaroTimeMs = 0;        // last barometer timestamp (ms)
  private adaNoiseScale = 1;
  private adaConvergeCount = 0;
  private _traceCache = 0;
  private _zuptWeight = 0;
  private _speedGate = 0;
  private _accelGate = 0;
  private robustWeight = 1;
  private _varAx = 0;
  private _varAy = 0;
  private _varGz = 0;
  constructor(config?: EkfConfig) {
    this.config = {
      dt: config?.dt ?? DEFAULTS.dt,
      processNoise: { ...DEFAULTS.processNoise, ...config?.processNoise },
      measurementNoise: { ...DEFAULTS.measurementNoise, ...config?.measurementNoise },
      initialCovariance: { ...DEFAULTS.initialCovariance, ...config?.initialCovariance },
      gateThreshold: config?.gateThreshold ?? DEFAULTS.gateThreshold,
      coastTimeoutMs: config?.coastTimeoutMs ?? DEFAULTS.coastTimeoutMs,
      gpsTimeOffsetMs: config?.gpsTimeOffsetMs ?? DEFAULTS.gpsTimeOffsetMs,
      magneticDeclination: config?.magneticDeclination ?? DEFAULTS.magneticDeclination,
      robustWeight: { ...DEFAULTS.robustWeight, ...config?.robustWeight } as Required<NonNullable<EkfConfig['robustWeight']>>,
      adaptiveNoise: { ...DEFAULTS.adaptiveNoise, ...config?.adaptiveNoise },
      adaptiveScaling: { ...DEFAULTS.adaptiveScaling, ...config?.adaptiveScaling },
      useLateralAccel: config?.useLateralAccel ?? DEFAULTS.useLateralAccel
    };

    this.S = matCreate(N, N);
    this.reset(0, 0, 0, 0);
  }

  reset(x: number, y: number, v: number, psi: number): void {
    const ic = this.config.initialCovariance;
    this.x[I.X] = x;
    this.x[I.Y] = y;
    this.x[I.V] = v;
    this.x[I.PSI] = this.wrapAngle(psi);
    this.x[I.BETA] = 0;
    this.x[I.A_BIAS_X] = 0;
    this.x[I.G_BIAS_Z] = 0;
    this.x[I.MAG_DECL] = this.config.magneticDeclination ?? 0;

    for (let i = 0; i < N; i++) this.S[i].fill(0);

    this.S[I.X][I.X] = Math.sqrt(ic.position!);
    this.S[I.Y][I.Y] = Math.sqrt(ic.position!);
    this.S[I.V][I.V] = Math.sqrt(ic.velocity!);
    this.S[I.PSI][I.PSI] = Math.sqrt(ic.heading!);
    this.S[I.BETA][I.BETA] = Math.sqrt(ic.sideslip!);
    this.S[I.A_BIAS_X][I.A_BIAS_X] = Math.sqrt(ic.accelBias!);
    this.S[I.G_BIAS_Z][I.G_BIAS_Z] = Math.sqrt(ic.gyroBias!);
    this.S[I.MAG_DECL][I.MAG_DECL] = Math.sqrt(ic.magDeclination ?? 0.25);

    this._traceCache = 0;
    for (let i = 0; i < N; i++) this._traceCache += this.S[i][i] * this.S[i][i];

    this.lastGpsTimeMs = 0;
    this.gpsInitTimeMs = 0;
    this.lastImuTimeMs = -Infinity;
    this.lastChiSq = 0;
    this.lastGatePassed = false;
    this.coasting = false;
    this.gpsInitialized = false;
    this.accelEnergy = 0;
    this.gyroEnergy = 0;
    this.varAccelEnergy = 0;
    this.varGyroEnergy = 0;
    this.stepBuffer.clear();
    this.stepFreq = 0;
    this.axWindow.clear();
    this.ayWindow.clear();
    this.azWindow.clear();
    this.gzWindow.clear();
    this.gxWindow.clear();
    this.gyWindow.clear();
    this.lastMagBearing = 0;
    this.lastMagTimeMs = 0;
    this.lastOmega = 0;
    this.prevOmega = 0;
    this.angAccel = 0;
    this.smoothAngAccel = 0;
    this.lastPredictTimeMs = 0;
    this.prevCallMagBearing = 0;
    this.prevCallMagTimeMs = 0;
    this.magRejectCount = 0;
    this.deviceToEnu = null;
    this.deviceToVehicle = null;
    this.curAzimuth = NaN;
    this.curPitch = NaN;
    this.curRoll = NaN;
    this.adaNoiseScale = 1;
    this.adaConvergeCount = 0;
    this.robustWeight = 1;
    this._zuptEngaged = false;
    this._wasZuptEngaged = false;
    this.coastSpeed = 0;
    this.coastSpeedReady = false;
    this.lastBaroAlt = NaN;
    this.lastBaroTimeMs = 0;
    this.bufTail = 0;
    this.bufLen = 0;
    this.imuTS.clear();
    this.stepTS.clear();
  }

  resetBiases(aBiasX?: number, gBiasZ?: number): void {
    this.x[I.A_BIAS_X] = aBiasX ?? 0;
    this.x[I.G_BIAS_Z] = gBiasZ ?? 0;
    const ic = this.config.initialCovariance;
    this.S[I.A_BIAS_X][I.A_BIAS_X] = Math.sqrt(ic.accelBias!);
    this.S[I.G_BIAS_Z][I.G_BIAS_Z] = Math.sqrt(ic.gyroBias!);
    for (let i = 0; i < I.A_BIAS_X; i++) {
      this.S[I.A_BIAS_X][i] = 0;
      this.S[I.G_BIAS_Z][i] = 0;
    }
    this.S[I.G_BIAS_Z][I.A_BIAS_X] = 0;
  }

  inflateCovariance(params: { position?: number; heading?: number }): void {
    if (params.position !== undefined) {
      const newStd = Math.sqrt(params.position);
      this.S[I.X][I.X] = newStd;
      this.S[I.Y][I.Y] = newStd;
    }
    if (params.heading !== undefined) {
      this.S[I.PSI][I.PSI] = Math.sqrt(params.heading);
    }
  }

  setOrientation(azimuth: number, pitch: number, roll: number): void {
    // Z-X-Y rotation (W3C DeviceOrientation convention):
    // R = Rz(azimuth) * Rx(pitch) * Ry(roll)
    // azimuth: rotation around Z (clockwise-positive compass heading)
    // pitch:   rotation around X (front-to-back tilt)
    // roll:    rotation around Y (left-to-right tilt)
    const ca = Math.cos(azimuth), sa = -Math.sin(azimuth);
    const cp = Math.cos(pitch), sp = Math.sin(pitch);
    const cr = Math.cos(roll), sr = Math.sin(roll);
    const R = this._tmpR3x3;
    R[0][0] = ca * cr - sa * sp * sr;  R[0][1] = -sa * cp;           R[0][2] = ca * sr + sa * sp * cr;
    R[1][0] = sa * cr + ca * sp * sr;  R[1][1] = ca * cp;            R[1][2] = sa * sr - ca * sp * cr;
    R[2][0] = -cp * sr;                 R[2][1] = sp;                R[2][2] = cp * cr;

    this.curAzimuth = azimuth;
    this.curPitch = pitch;
    this.curRoll = roll;

    // Device-to-vehicle orientation: R_dev_to_veh = Rz(-ψ) * R_dev_to_ENU
    // This stays constant when the vehicle turns. We compare it across calls
    // to detect actual device repositioning (not vehicle turns).
    const psi = this.x[I.PSI];
    const cPsi = Math.cos(psi), sPsi = Math.sin(psi);
    const Rdv = this._tmpRdv3x3;
    for (let col = 0; col < 3; col++) {
      Rdv[0][col] = cPsi * R[0][col] + sPsi * R[1][col];
      Rdv[1][col] = -sPsi * R[0][col] + cPsi * R[1][col];
      Rdv[2][col] = R[2][col];
    }

    if (this.deviceToVehicle) {
      const prev = this.deviceToVehicle;
      const trace = Rdv[0][0]*prev[0][0] + Rdv[0][1]*prev[0][1] + Rdv[0][2]*prev[0][2]
                  + Rdv[1][0]*prev[1][0] + Rdv[1][1]*prev[1][1] + Rdv[1][2]*prev[1][2]
                  + Rdv[2][0]*prev[2][0] + Rdv[2][1]*prev[2][1] + Rdv[2][2]*prev[2][2];
      if ((trace - 1) / 2 < 0.996) { // > ~5° rotation relative to vehicle
        this.x[I.A_BIAS_X] = 0;
        this.x[I.G_BIAS_Z] = 0;
        const ic = this.config.initialCovariance;
        this.S[I.A_BIAS_X][I.A_BIAS_X] = Math.sqrt(ic.accelBias!);
        this.S[I.G_BIAS_Z][I.G_BIAS_Z] = Math.sqrt(ic.gyroBias!);
        for (let i = 0; i < I.A_BIAS_X; i++) {
          this.S[I.A_BIAS_X][i] = 0;
          this.S[I.G_BIAS_Z][i] = 0;
        }
        this.S[I.G_BIAS_Z][I.A_BIAS_X] = 0;
      }
    }
    this.deviceToEnu = R;
    this.deviceToVehicle = Rdv;
  }

  predict(ax: number, ay: number, gz: number, dt: number, timestampMs: number, az?: number, gx?: number, gy?: number): void {
    if (dt <= 0) return;
    if (timestampMs <= this.lastImuTimeMs) return;
    // During coasting (no GPS), immediately decay lastGpsSpeed toward 0 so ZUPT
    // can engage if the car stops in a tunnel/basement.  Without this, a car that
    // enters a basement at >2 m/s has gpsMoving=true forever, permanently
    // blocking ZUPT — the car icon "shoots forward" because velocity persists.
    // Decay starts immediately (no delay) with 2s time constant for fast response.
    if (this.coasting && this.lastGpsSpeed > 0) {
      const decay = Math.exp(-dt / 2);
      this.lastGpsSpeed *= decay;
      if (this.lastGpsSpeed < 0.1) this.lastGpsSpeed = 0;
    }
    // Capture speed at GPS loss as the velocity prior target.  This must happen
    // BEFORE lastGpsSpeed decays to 0.  The prior pulls v toward coastSpeed
    // during continuous-motion coasting (e.g. driving underground), preventing
    // accel bias drift from integrating into velocity error.
    if (this.coasting && !this.coastSpeedReady && this.gpsInitialized) {
      this.coastSpeed = Math.abs(this.x[I.V]);
      if (this.coastSpeed < 0.5) this.coastSpeed = 0;  // was stationary — no prior needed
      this.coastSpeedReady = true;
    }
    const rawAx = ax, rawAy = ay, rawGz = gz;
    const rawAz = az ?? 0, rawGx = gx ?? 0, rawGy = gy ?? 0;
    // When orientation alignment is active, a 6-axis IMU (including az) is
    // required. Without az, gravity leaks into the horizontal ENU axes and
    // causes velocity drift. We still accept the call for backwards
    // compatibility but issue a NaN to help catch bugs in testing.
    if (this.deviceToEnu && az === undefined) {
      throw new Error('6-axis IMU required when device orientation is set — az parameter is mandatory');
    }
    if (this.deviceToEnu) {
      const R = this.deviceToEnu;
      ax = R[0][0] * rawAx + R[0][1] * rawAy + R[0][2] * rawAz;
      ay = R[1][0] * rawAx + R[1][1] * rawAy + R[1][2] * rawAz;
      gz = R[2][0] * rawGx + R[2][1] * rawGy + R[2][2] * rawGz;
      const psiBody = this.x[I.PSI];
      const forwardA = ax * Math.cos(psiBody) + ay * Math.sin(psiBody);
      const lateralA = -ax * Math.sin(psiBody) + ay * Math.cos(psiBody);
      ax = forwardA;
      ay = lateralA;
    }
    if (!isFinite(ax) || !isFinite(gz)) {
      this.lastImuTimeMs = timestampMs;
      this.saveBuf(timestampMs, rawAx, rawAy, rawGz, dt, rawAz, rawGx, rawGy);
      return;
    }
    const a = ax - this.x[I.A_BIAS_X];
    const omega = gz - this.x[I.G_BIAS_Z];
    const oldOmega = this.prevOmega;
    const psi = this.x[I.PSI], beta = this.x[I.BETA], v = this.x[I.V];
    const absV = Math.abs(v), absOmega = Math.abs(omega);

    // Angular acceleration: |dω/dt| for transient detection (corner entry/exit)
    if (this.lastPredictTimeMs > 0) {
      this.angAccel = Math.min(Math.abs(omega - oldOmega) / Math.max(dt, 1e-6), 10);
      // Sustain angAccel boost for ~0.5s after the transient ends
      const sustainAlpha = Math.exp(-dt / 0.5);
      this.smoothAngAccel = Math.max(this.angAccel, this.smoothAngAccel * sustainAlpha);
    } else {
      this.angAccel = 0;
      this.smoothAngAccel = 0;
    }

    const omegaAvg = (this.lastPredictTimeMs > 0) ? 0.5 * (oldOmega + omega) : omega;
    this.prevOmega = omega;

    this.lastOmega = omega;
    this.lastPredictTimeMs = timestampMs;
    this.lastImuTimeMs = timestampMs;
    this.updateImuWindow(ax, ay, gz, az, gx, gy, timestampMs);
    this.updateStepDetection(ax, timestampMs);
    this.computeAdaptiveQ(dt, a, omega);
    const stillness = this.getStillness();

    // Sideslip time constant: shorter during angular transients (corner entry/exit)
    // so β tracks the rapidly changing slip angle instead of lagging and corrupting ψ.
    // At max angAccel (3 rad/s²), τ drops to 40% of its base value.
    const angAccelNorm = Math.min(this.smoothAngAccel / 3.0, 1);
    const betaTauBase = absV < 0.3 ? 0.1
      : absOmega > EPS ? 1.5
      : Math.max(0.5, 1.5 - (absV - 1.5) * (1.0 / 3.5));
    const betaTau = betaTauBase * (1 - 0.6 * angAccelNorm) * (1 - 0.5 * this.stepEnergy);
    const expDt50 = Math.exp(-dt / 50);
    this.computeJacobian(a, omegaAvg, dt, betaTau, this.accelEnergy < 0.05 ? expDt50 : 1);

    for (let i = 0; i < N; i++)
      for (let j = 0; j < N; j++) {
        let s = 0;
        for (let k = 0; k < N; k++) s += this.tmpF[i][k] * this.S[k][j];
        this.tmpFS[i][j] = s;
      }

    for (let i = 0; i < N; i++)
      for (let j = 0; j < N; j++)
        this.tmpQR[i][j] = this.tmpFS[j][i];

    for (let i = 0; i < N; i++)
      for (let j = 0; j < N; j++)
        this.tmpQR[N + i][j] = this.tmpSqrtQ[i][j];

    this.qrInPlace(this.tmpQR, 2 * N, N, this.tmpHouseV);

    this.copySfromQR(this.tmpQR, 0);

    const psiBeta = psi + beta;
    const vAvg = v + 0.5 * a * dt;
    const [dx, dy] = this.ctraDelta(psiBeta, vAvg, omegaAvg, dt);
    this.x[I.X] += dx;
    this.x[I.Y] += dy;
    this.x[I.V] = Math.max(0, v + a * dt);
    const vehMoving = Math.min(absV / 0.5, 1);
    const omegaScale = vehMoving > 0 ? Math.max(vehMoving, stillness) : 1;
    this.x[I.PSI] = this.wrapAngle(psi + omegaAvg * dt * omegaScale);
    this.x[I.BETA] *= Math.exp(-dt / betaTau);
    // A_BIAS_X updated only via GPS velocity corrections (no decay needed)

    if (this.config.useLateralAccel && absOmega > 0.1 && absV > 0.2)
      this.applyLateralAccel(ay, omega);
    else if (absOmega < 0.1 && absV > 0.15)
      this.applyNonholonomic(omega);
    // ZUPT engages when the IMU acceleration is CONSTANT (low variance): a
    // stationary device reads a constant ax = sensor bias + gravity projection,
    // which must be learned as aBiasX. A genuine (varying) acceleration has high
    // ax variance and correctly excludes ZUPT. The speed gate protects genuine
    // constant-velocity motion (which IMU alone cannot distinguish from rest):
    // ZUPT is disabled above ~3 m/s, so a cruising vehicle keeps its velocity.
    // The low-variance (constant-ax) gate lets ZUPT fire during the low-speed
    // approach and LEARN the bias BEFORE v escapes past the speed gate — this is
    // what prevents a constant bias from integrating into a runaway straight-line
    // drift (the reported "car icon moves on a line" symptom). Note: the old
    // velocityStillness ramp (cutoff at |v|=0.15) was removed — it killed ZUPT
    // after a single step, before the bias could be learned, and let v diverge.
    // ZUPT — dual-threshold hysteresis (ON=0.15, OFF=0.05) prevents chatter at
    // the boundary. Chi-square innovation gate rejects velocity innovations that
    // are statistically inconsistent with zero-velocity hypothesis (3σ).
    const speedGate = Math.min(Math.max((8.0 - Math.abs(this.x[I.V])) / 4.0, 0), 1);
    const gpsMoving = this.lastGpsSpeed > 2.0;
    const zuptWeight = stillness * speedGate;
    this._zuptWeight = zuptWeight;
    this._speedGate = speedGate;
    this._accelGate = 0;

    const ZUPT_ON = 0.15, ZUPT_OFF = 0.05;
    if (!this._zuptEngaged) {
      if (!gpsMoving && zuptWeight >= ZUPT_ON && this.axWindow.length >= 5 && this.zuptChiSqGate()) {
        this.applyZupt(zuptWeight, omega);
        this._zuptEngaged = true;
        this._wasZuptEngaged = true;
      }
    } else {
      if (gpsMoving || zuptWeight < ZUPT_OFF) {
        this._zuptEngaged = false;
      } else if (this.zuptChiSqGate()) {
        this.applyZupt(zuptWeight, omega);
      }
    }
    if (this._wasZuptEngaged && !this._zuptEngaged) {
      this._wasZuptEngaged = false;
      const si = this.S;
      si[I.V][I.V] = Math.max(si[I.V][I.V], 1.5);
      si[I.X][I.X] = Math.max(si[I.X][I.X], 3);
      si[I.Y][I.Y] = Math.max(si[I.Y][I.Y], 3);
      si[I.PSI][I.PSI] = Math.max(si[I.PSI][I.PSI], 0.3);
    }
    // Velocity damping during coasting: when GPS is lost and the car appears
    // stationary (low IMU energy), inject a soft pseudo-measurement v ≈ 0 via
    // scalar QR.  This prevents the "car icon shoots forward" symptom — without
    // GPS, uncorrected accel bias integrates into velocity and position.  The
    // measurement noise tightens over time (loose initially, tight after 30s) so
    // genuine motion is preserved but stopped cars converge to v = 0.
    if (this.coasting && !this._zuptEngaged) {
      const coastTimeS = Math.max(0, (this.lastImuTimeMs - this.lastGpsTimeMs)) / 1000;
      const stationarity = this.accelEnergy + this.gyroEnergy;
      if (coastTimeS > 3 && stationarity < 0.1 && Math.abs(this.x[I.V]) > 0.05) {
        const rDamp = Math.max(0.05, 10.0 / Math.max(coastTimeS, 1));
        const innov = -this.x[I.V];
        const SV = this.S[I.V];
        let sInnov = rDamp * rDamp;
        for (let j = 0; j < N; j++) sInnov += SV[j] * SV[j];
        for (let i = 0; i < N; i++) {
          let p = 0;
          const lim = Math.min(i, I.V);
          for (let k = 0; k <= lim; k++) p += this.S[i][k] * SV[k];
          this.x[i] += p / sInnov * innov;
        }
        this.x[I.PSI] = this.wrapAngle(this.x[I.PSI]);
        // QR covariance update (scalar QR, same pattern as lateral-accel/mag)
        const A = this.tmpLatPre;
        for (let i = 0; i < MAG_PRE; i++) A[i].fill(0);
        A[0][0] = rDamp;
        for (let j = 0; j < N; j++) A[1 + j][0] = SV[j];
        for (let i = 0; i < N; i++)
          for (let j = 0; j < N; j++)
            A[1 + i][1 + j] = this.S[j][i];
        this.qrInPlace(A, MAG_PRE, MAG_PRE, this.tmpHouseV);
        this.copySfromQR(A, 1);
      }
    }
    // Velocity prior during coasting: when GPS is lost and the car is moving
    // (stationarity check fails), pull v toward coastSpeed (the speed at GPS
    // loss) to prevent accel bias from integrating into velocity drift.  The
    // prior strength decays exponentially — tight early (strong anchor), loose
    // after ~30s (car may have turned/stopped).  This addresses the multi-
    // basement scenario where continuous motion without GPS causes v to drift
    // ~0.6 m/s per 60s from uncorrected aBiasX.
    if (this.coasting && this.coastSpeedReady && this.coastSpeed > 0) {
      const coastTimeS = Math.max(0, (this.lastImuTimeMs - this.lastGpsTimeMs)) / 1000;
      const stationarity = this.accelEnergy + this.gyroEnergy;
      // Only fire when car is genuinely moving (not stopped — ZUPT handles that)
      if (stationarity > 0.05 || Math.abs(this.x[I.V]) > 0.5) {
        // Prior R grows with coast time: tight at onset (R=2), loose after 30s (R≈20)
        const rPrior = Math.min(2.0 + coastTimeS * 0.6, 20.0);
        const innov = this.x[I.V] - this.coastSpeed;
        const SV = this.S[I.V];
        let sInnov = rPrior * rPrior;
        for (let j = 0; j < N; j++) sInnov += SV[j] * SV[j];
        for (let i = 0; i < N; i++) {
          let p = 0;
          const lim = Math.min(i, I.V);
          for (let k = 0; k <= lim; k++) p += this.S[i][k] * SV[k];
          this.x[i] += p / sInnov * innov;
        }
        this.x[I.PSI] = this.wrapAngle(this.x[I.PSI]);
        // QR covariance update (scalar QR, same pattern as lateral-accel/mag)
        const A = this.tmpLatPre;
        for (let i = 0; i < MAG_PRE; i++) A[i].fill(0);
        A[0][0] = rPrior;
        for (let j = 0; j < N; j++) A[1 + j][0] = SV[j];
        for (let i = 0; i < N; i++)
          for (let j = 0; j < N; j++)
            A[1 + i][1 + j] = this.S[j][i];
        this.qrInPlace(A, MAG_PRE, MAG_PRE, this.tmpHouseV);
        this.copySfromQR(A, 1);
      }
    }

    // ZARU runs independently of ZUPT (gated on |omega| inside the method) so
    // gyro bias is corrected even when GPS-moving disables ZUPT.
    this.applyZaru(omega);
    this.safeguardState();
    this.saveBuf(timestampMs, rawAx, rawAy, rawGz, dt, rawAz, rawGx, rawGy);
  }

  updateGps(x: number, y: number, vx: number, vy: number, timestampMs: number, accuracyMeters?: number): boolean {
    if (!this.gpsInitialized) {
      this.x[I.X] = x;
      this.x[I.Y] = y;
      const spd = Math.sqrt(vx * vx + vy * vy);
      this.lastGpsSpeed = spd;
      this.smoothedSpeed = undefined;
      if (spd > 0.1) {
        this.x[I.V] = spd;
        this.x[I.PSI] = this.wrapAngle(Math.atan2(vy, vx));
      }
      this.gpsInitialized = true;
      this.gpsInitTimeMs = timestampMs;
      this.lastGpsTimeMs = timestampMs + this.config.gpsTimeOffsetMs!;
      this.lastGatePassed = true;
      this.coastSpeedReady = false;  // reset velocity prior on GPS (re)init
      // Tighten covariance to reflect GPS-derived position/velocity knowledge.
      // Without this, S[X][X]=10 (σ=10m) lets IMU drift accumulate rapidly
      // between GPS fixes, causing the outlier guard to reject valid fixes.
      // Use Math.min so we never loosen covariance tighter than what reset() set.
      const gpsAccStd = accuracyMeters !== undefined ? Math.max(accuracyMeters / 3, 1) : 3;
      this.S[I.X][I.X] = Math.min(this.S[I.X][I.X], gpsAccStd);
      this.S[I.Y][I.Y] = Math.min(this.S[I.Y][I.Y], gpsAccStd);
      if (spd > 0.1) {
        const vStd = Math.min(this.S[I.V][I.V], 0.5);
        this.S[I.V][I.V] = vStd;
        const psiStd = Math.max(0.3, spd > 2 ? 0.15 : 0.3);
        this.S[I.PSI][I.PSI] = Math.min(this.S[I.PSI][I.PSI], psiStd);
      }
      this.S[I.A_BIAS_X][I.A_BIAS_X] = Math.sqrt(this.config.initialCovariance.accelBias!);
      this.S[I.G_BIAS_Z][I.G_BIAS_Z] = Math.sqrt(this.config.initialCovariance.gyroBias!);
      this.S[I.MAG_DECL][I.MAG_DECL] = Math.sqrt(this.config.initialCovariance.magDeclination ?? 0.25);
      let tr = 0;
      for (let i = 0; i < N; i++) {
        for (let j = 0; j <= i; j++) { const v = this.S[i][j]; tr += v * v; }
      }
      this._traceCache = tr;
      return true;
    }

    // GPS timestamp monotonicity: reject stale or duplicate GPS
    const effectiveGpsTime = timestampMs + this.config.gpsTimeOffsetMs!;
    if (effectiveGpsTime <= this.lastGpsTimeMs) return true;

    // GPS latency compensation: rewind to GPS time if delayed
    let replayCount = 0;
    if (this.lastImuTimeMs - effectiveGpsTime > 0.001 && this.bufLen > 0) {
      const rc = this.rewindTo(effectiveGpsTime);
      if (rc >= 0) {
        replayCount = rc;
        // Rewind restores a predict-buffer state which may PRE-DATE the first
        // GPS initialization, reverting x,y to the pre-GPS position (e.g. 0,0)
        // while gpsInitialized remains true. Detect this via a large position
        // discrepancy and re-apply the GPS position directly.
        const dx = x - this.x[I.X], dy = y - this.x[I.Y];
        if (dx * dx + dy * dy > 100) {
          this.x[I.X] = x;
          this.x[I.Y] = y;
          if (this.coasting) {
            const spd = Math.sqrt(vx * vx + vy * vy);
            this.lastGpsSpeed = spd;
            this.x[I.V] = spd;
            this.x[I.PSI] = this.wrapAngle(Math.atan2(vy, vx));
            this.x[I.BETA] = 0;
          }
        }
      }
    }

    // Auto-detect divergence: if the last GPS update was rejected (or heavily
    // down-weighted by Huber) AND the GPS position is >100m from the current
    // state, the state has clearly diverged (e.g. after basement exit). Force
    // coasting so resetFromGps triggers on the next fix.
    if (!this.coasting && (!this.lastGatePassed || this.robustWeight < 0.05)) {
      const dxDiv = x - this.x[I.X], dyDiv = y - this.x[I.Y];
      if (dxDiv * dxDiv + dyDiv * dyDiv > 10000) {
        this.coasting = true;
      }
    }

    const origPosR = Math.max(this.config.measurementNoise.position!, 1e-6);
    const origVelR = this.config.measurementNoise.velocity!;
    let posR = origPosR;
    let velR = origVelR;
    if (accuracyMeters !== undefined) {
      const sc = Math.max(accuracyMeters / origPosR, 0.1);
      posR = origPosR * sc;
      velR = origVelR * sc;
    }
    // Speed-dependent velR inflation: ramps from 5× at v=0 to 1× at v≥10
    // Prevents velocity jumps at moderate city speeds (5-10 m/s) where
    // GPS multipath noise is high but the old ramp (v/5) was already flat.
    const speedRamp = Math.max(0, 1 - Math.abs(this.x[I.V]) / 10);
    velR *= (1 + 2 * speedRamp);

    // Low-speed GPS Doppler distrust: in urban canyons, multipath can flip the
    // GPS velocity direction ~180°. If the measured velocity opposes the IMU-
    // integrated heading (ψ+β) at low speed, strongly distrust the Doppler so
    // it cannot drag v negative or rotate heading. Position innovation still
    // corrects heading through P[PSI][X/Y] cross-covariance.
    const psiBetaLow = this.x[I.PSI] + this.x[I.BETA];
    const dotLow = vx * Math.cos(psiBetaLow) + vy * Math.sin(psiBetaLow);
    if (Math.abs(this.x[I.V]) < 2.5 && dotLow < 0) velR *= 4;

    // Direction-aware position step-change outlier guard
    const dtSinceLastGps = Math.max(effectiveGpsTime - this.lastGpsTimeMs, 0) / 1000;
    this.lastGpsSpeed = Math.sqrt(vx * vx + vy * vy);
    const maxPlausibleSpeed = Math.max(Math.abs(this.x[I.V]) * 2, 1.0) + 2;
    const dxGps = x - this.x[I.X], dyGps = y - this.x[I.Y];
    const psiBeta = this.x[I.PSI] + this.x[I.BETA];
    const cp = Math.cos(psiBeta), sp = Math.sin(psiBeta);
    const forward = dxGps * cp + dyGps * sp;    // along-track
    const cross   = -dxGps * sp + dyGps * cp;   // cross-track
    const dtBase = Math.max(dtSinceLastGps, 0.1);
    // Guard thresholds account for GPS accuracy — posR already reflects accuracy.
    // Scale min thresholds by posR so the noise floor at e.g. 16m accuracy
    // (posR ≈ 5.3) doesn't trigger inflation on every fix.
    const minForward = Math.max(posR * 2, 1.0);
    const minCross   = Math.max(posR * 0.5, 0.5);
    const maxForward = Math.max(maxPlausibleSpeed * dtBase * 2, minForward);
    const maxCross   = Math.max(maxPlausibleSpeed * dtBase * 0.5, minCross);
    const preGuardPosR = posR;
    if (Math.abs(forward) > maxForward) posR *= Math.min(Math.abs(forward) / maxForward, 5);
    if (Math.abs(cross) > maxCross) posR *= Math.min(Math.abs(cross) / maxCross, 10);

    // Guard-inflation-masking: if guard inflates R by >3× and the position
    // jump is physically plausible (at highway speeds ~50 m/s + 2m GPS
    // noise), the chiSq would pass with inflated R but the Kalman gain
    // drops to ~0.04, producing minute-long convergence. Force reset.
    // Allow snap at low speeds (>0.5 m/s) or for very large jumps (>15m).
    if (posR / preGuardPosR > 3) {
      const dxNorm = Math.sqrt(dxGps * dxGps + dyGps * dyGps);
      if (dxNorm < dtSinceLastGps * 50 + 2 && (Math.abs(this.x[I.V]) > 0.5 || dxNorm > 15)) {
        this.resetFromGps(x, y, vx, vy);
        this.lastGpsTimeMs = effectiveGpsTime;
        this.lastGatePassed = true;
        if (replayCount > 0) this.replayFromScratch(replayCount);
        return true;
      }
    }

    this.tmpZ[0] = x; this.tmpZ[1] = y; this.tmpZ[2] = vx; this.tmpZ[3] = vy;

    // Smoothly blend GPS velocity toward zero as device becomes more stationary,
    // inflating velR proportionally — no hard threshold, no on/off oscillation.
    // Hybrid speed: GPS primary, 30 % of EKF v as fallback against GPS glitches
    // at sustained speed (GPS speed can briefly read near zero during an RF
    // dropout while the vehicle is still doing 5+ m/s).
    // Smoothed with a 3 s exponential moving average so stationaryWeight doesn't
    // oscillate in stop-and-go traffic.
    if (this.smoothedSpeed === undefined) this.smoothedSpeed = this.lastGpsSpeed;
    const ekfSpeed = Math.abs(this.x[I.V]);
    const hybridSpeed = Math.max(this.lastGpsSpeed, ekfSpeed * 0.3);
    this.smoothedSpeed = 0.7 * this.smoothedSpeed + 0.3 * hybridSpeed;
    // Transition zone 0–1.0 m/s (was 0–0.5): GPS course error at 0.5 m/s is
    // still ~30° (3 m accuracy) — keeping more compass influence during the
    // critical early-acceleration phase.
    const threshold = 1.0;
    const stationaryWeight = Math.max(0, (threshold - this.smoothedSpeed) / threshold);
    this.tmpZ[2] = vx * (1 - stationaryWeight);
    this.tmpZ[3] = vy * (1 - stationaryWeight);
    // Cap inflation at 5× (was 10×) so the Kalman gains don't drop below ~16 %
    // and the first trusted fix after stationary mode leaves quickly enough.
    velR *= (1 + Math.min(4, 9 * stationaryWeight));

    // Decouple heading from position/velocity innovations when stationary:
    // GPS multipath bias at a red light creates a consistent position error
    // that pulls heading through P[PSI][X/Y/V] cross-covariance and accumulates
    // over seconds. Zero the Cholesky cross-rows to prevent this. The CTRA
    // prediction dynamics recreate them naturally when motion resumes.
    const w = stationaryWeight;
    this.S[I.PSI][I.X] *= (1 - w);
    this.S[I.PSI][I.Y] *= (1 - w);
    this.S[I.PSI][I.V] *= (1 - w);

    const ok = this.gpsUpdateSingle(posR, velR);
    let result = ok;
    if (ok) {
      const dx = x - this.x[I.X], dy = y - this.x[I.Y];
      if (dx * dx + dy * dy > 1e8) {
        this.resetFromGps(x, y, vx, vy);
        this.lastGpsTimeMs = effectiveGpsTime;
        this.lastGatePassed = true;
        result = true;
      } else {
        this.lastGatePassed = true;
        this.lastGpsTimeMs = effectiveGpsTime;
      }
    } else {
      this.lastGatePassed = false;
      if (this.coasting) {
        this.resetFromGps(x, y, vx, vy);
        this.lastGpsTimeMs = effectiveGpsTime;
        this.lastGatePassed = true;
        result = true;
      }
    }

    // Re-predict IMU steps between GPS time and current time
    if (replayCount > 0) this.replayFromScratch(replayCount);

    return result;
  }

  // ─── Barometric altitude update ──────────────────────────────
  // Uses the phone's barometer to constrain forward speed on ramps via the
  // relationship vz = v × sin(pitch).  When |sin(pitch)| > 0.03 (≥ 1.7° grade,
  // typical of parking ramps), a scalar QR pseudo-measurement injects
  // v = vz_measured / sin(pitch) with noise reflecting barometer accuracy.
  //
  // On flat ground (|sin(pitch)| ≤ 0.03) the barometer provides no speed
  // information and the update is skipped.  When device orientation is not set
  // (deviceToEnu === null), pitch is unknown so the update is also skipped.
  //
  // This addresses the multi-basement scenario where continuous ramp driving
  // without GPS causes velocity drift from uncorrected accel bias.  The
  // barometric constraint is most effective on ramps (pitch 3–9°) where the
  // altitude change rate provides a strong speed observable.
  updateBaro(altitude: number, timestampMs: number): void {
    if (!this.gpsInitialized || !isFinite(altitude)) return;
    if (!this.coasting) { this.lastBaroAlt = altitude; this.lastBaroTimeMs = timestampMs; return; }
    if (!isFinite(this.lastBaroAlt) || timestampMs <= this.lastBaroTimeMs) {
      this.lastBaroAlt = altitude;
      this.lastBaroTimeMs = timestampMs;
      return;
    }
    const dtBaro = (timestampMs - this.lastBaroTimeMs) / 1000;
    if (dtBaro < 0.05 || dtBaro > 5) { this.lastBaroAlt = altitude; this.lastBaroTimeMs = timestampMs; return; }
    // Need pitch estimate — only available when device orientation is set
    const pitch = this.curPitch;
    if (!isFinite(pitch)) { this.lastBaroAlt = altitude; this.lastBaroTimeMs = timestampMs; return; }
    const sinPitch = Math.sin(pitch);
    if (Math.abs(sinPitch) < 0.03) {
      // On flat ground — barometer gives no speed information
      this.lastBaroAlt = altitude;
      this.lastBaroTimeMs = timestampMs;
      return;
    }
    // Vertical velocity from barometer (positive = ascending)
    const vzMeasured = (altitude - this.lastBaroAlt) / dtBaro;
    this.lastBaroAlt = altitude;
    this.lastBaroTimeMs = timestampMs;
    // Observation model: vz = v × sin(pitch)  →  H[V] = sin(pitch)
    // Measurement noise: barometer σ ≈ 0.5m → for 1s windows, vel σ ≈ 0.5 m/s
    // Plus pitch uncertainty: add 10% of the signal
    const hV = sinPitch;
    const rBaro = Math.max(0.5, Math.abs(vzMeasured) * 0.1);
    const SV = this.S[I.V];
    // hs = H × S = sinPitch × S[V][:]  (scalar QR, same pattern as ZUPT/ZARU)
    const hs = this.tmpLatHS;
    for (let j = 0; j < N; j++) hs[j] = hV * SV[j];
    const innov = vzMeasured - this.x[I.V] * hV;
    let sInnov = rBaro * rBaro;
    for (let j = 0; j < N; j++) sInnov += hs[j] * hs[j];
    // Chi-square gate: reject if innovation is statistically implausible (3σ ≈ chi²₁,0.99 = 11.3)
    if (innov * innov / sInnov > 11.3) return;
    // State update
    for (let i = 0; i < N; i++) {
      let p = 0;
      const lim = Math.min(i, I.V);
      for (let k = 0; k <= lim; k++) p += this.S[i][k] * hs[k];
      this.x[i] += p / sInnov * innov;
    }
    this.x[I.PSI] = this.wrapAngle(this.x[I.PSI]);
    // QR covariance update (scalar QR, same pattern as lateral-accel/mag)
    const A = this.tmpLatPre;
    for (let i = 0; i < MAG_PRE; i++) A[i].fill(0);
    A[0][0] = rBaro;
    for (let j = 0; j < N; j++) A[1 + j][0] = hs[j];
    for (let i = 0; i < N; i++)
      for (let j = 0; j < N; j++)
        A[1 + i][1 + j] = this.S[j][i];
    this.qrInPlace(A, MAG_PRE, MAG_PRE, this.tmpHouseV);
    this.copySfromQR(A, 1);
  }

  updateMag(bearing: number, timestampMs: number): void {
    this.lastMagBearing = bearing;
    this.lastMagTimeMs = timestampMs;
    this.magUpdateSingle(bearing);
  }

  coast(timeoutMs: number, currentTimeMs: number): boolean {
    this.coasting = currentTimeMs - this.lastGpsTimeMs > timeoutMs;
    for (let i = 0; i < N; i++) {
      if (!isFinite(this.x[i])) return false;
    }
    return true;
  }

  getState(): NavigationSolution {
    return {
      x: this.x[I.X], y: this.x[I.Y], v: this.x[I.V], psi: this.x[I.PSI], beta: this.x[I.BETA],
      aBiasX: this.x[I.A_BIAS_X], gBiasZ: this.x[I.G_BIAS_Z], magDeclination: this.x[I.MAG_DECL],
      p: matLowerToFull(this.S)
    };
  }

  getStateInto(out: NavigationSolution): void {
    out.x = this.x[I.X]; out.y = this.x[I.Y]; out.v = this.x[I.V];
    out.psi = this.x[I.PSI]; out.beta = this.x[I.BETA];
    out.aBiasX = this.x[I.A_BIAS_X]; out.gBiasZ = this.x[I.G_BIAS_Z];
    out.magDeclination = this.x[I.MAG_DECL];
    matLowerToFullInto(this.S, out.p);
  }

  /** Lightweight debug snapshot of internal stillness/bias state for UI panels. */
   getDebug(): {
    stillness: number; aBiasX: number; gBiasZ: number; v: number; psiDeg: number;
    zuptWeight: number; speedGate: number; accelGate: number; n: number;
    magRejectCount: number; magTrust: number; magInnovDeg: number;
    gateThreshDeg: number; magAlpha: number;
  } {
    return {
      stillness: this.getStillness(),
      aBiasX: this.x[I.A_BIAS_X],
      gBiasZ: this.x[I.G_BIAS_Z],
      v: this.x[I.V],
      psiDeg: (this.x[I.PSI] * 180) / Math.PI,
      zuptWeight: this._zuptWeight,
      speedGate: this._speedGate,
      accelGate: this._accelGate,
      n: this.axWindow.length,
      magRejectCount: this.magRejectCount,
      magTrust: this._debugMagTrust,
      magInnovDeg: this._debugInnovDeg,
      gateThreshDeg: this._debugGateThreshDeg,
      magAlpha: this._debugMagAlpha,
    };
  }

  /** Debug: bias-corrected mean and std of the 1s IMU windows, for diagnosing
   *  rest detection on a specific device. meanAxRel = (mean ax) - aBiasX, etc. */
  getImuStats(): {
    n: number; meanAxRel: number; stdAx: number;
    meanGzRel: number; stdGz: number; lastOmega: number;
  } {
    const wmean = (b: RingBuf): number => {
      if (b.length === 0) return 0; let s = 0; for (let i = 0; i < b.length; i++) s += b.get(i); return s / b.length;
    };
    const wstd = (b: RingBuf, m: number): number => {
      if (b.length < 2) return 0; let s = 0; for (let i = 0; i < b.length; i++) { const d = b.get(i) - m; s += d * d; } return Math.sqrt(s / b.length);
    };
    const mAx = wmean(this.axWindow), mGz = wmean(this.gzWindow);
    return {
      n: this.axWindow.length,
      meanAxRel: mAx - this.x[I.A_BIAS_X],
      stdAx: wstd(this.axWindow, mAx),
      meanGzRel: mGz - this.x[I.G_BIAS_Z],
      stdGz: wstd(this.gzWindow, mGz),
      lastOmega: this.lastOmega,
    };
  }

  getDiagnostics(): EkfDiagnostics {
    this._innovCache[0] = this.tmpInnov[0];
    this._innovCache[1] = this.tmpInnov[1];
    this._innovCache[2] = this.tmpInnov[2];
    this._innovCache[3] = this.tmpInnov[3];
    return {
      trace: this._traceCache,
      gpsInnovation: this._innovCache,
      gpsChiSq: this.lastChiSq, gatePassed: this.lastGatePassed,
      coasting: this.coasting, lastGpsTimeMs: this.lastGpsTimeMs,
      lastImuTimeMs: this.lastImuTimeMs, stationary: this.getStillness() > 0.7 && Math.abs(this.x[I.V]) < 3.0,
      magDeclination: this.x[I.MAG_DECL], robustWeight: this.robustWeight, adaNoiseScale: this.adaNoiseScale
    };
  }

  // ─── private helpers ────────────────────────────────────────────

  private wrapAngle(a: number): number {
    a = a % TWO_PI;
    if (a > Math.PI) a -= TWO_PI;
    if (a <= -Math.PI) a += TWO_PI;
    return a;
  }

  private copySfromQR(Q: Float64Array[], offset: number): void {
    for (let i = 0; i < N; i++) {
      for (let j = 0; j <= i; j++) this.S[i][j] = Q[offset + j][offset + i];
      for (let j = i + 1; j < N; j++) this.S[i][j] = 0;
    }
    ensureDiag(this.S);
    let tr = 0;
    for (let i = 0; i < N; i++) {
      for (let j = 0; j <= i; j++) {
        const v = this.S[i][j];
        tr += v * v;
      }
    }
    this._traceCache = tr;
  }

  private qrInPlace(A: Float64Array[], m: number, n: number, vBuf: Float64Array): void {
    for (let k = 0; k < Math.min(m, n); k++) {
      let nrm = 0;
      for (let i = k; i < m; i++) nrm += A[i][k] * A[i][k];
      nrm = Math.sqrt(nrm);
      if (nrm < 1e-15) continue;
      const sign = A[k][k] >= 0 ? 1 : -1;
      vBuf[0] = A[k][k] + sign * nrm;
      const len = m - k;
      for (let i = 1; i < len; i++) vBuf[i] = A[k + i][k];
      let beta = 0;
      for (let i = 0; i < len; i++) beta += vBuf[i] * vBuf[i];
      beta = 2 / beta;
      for (let j = k; j < n; j++) {
        let s = 0;
        for (let i = 0; i < len; i++) s += vBuf[i] * A[k + i][j];
        s *= beta;
        for (let i = 0; i < len; i++) A[k + i][j] -= s * vBuf[i];
      }
    }
  }

  private safeguardState(): void {
    let corrupted = false;
    for (let i = 0; i < N; i++) {
      if (!isFinite(this.x[i])) { this.x[i] = 0; corrupted = true; }
      if (!isFinite(this.S[i][i]) || this.S[i][i] <= 0) { this.S[i][i] = 1e-6; corrupted = true; }
      for (let j = 0; j < i; j++)
        if (!isFinite(this.S[i][j])) { this.S[i][j] = 0; corrupted = true; }
      for (let j = i + 1; j < N; j++)
        if (this.S[i][j] !== 0) { this.S[i][j] = 0; corrupted = true; }
    }
    if (corrupted) {
      const ic = this.config.initialCovariance;
      const sigmas = [
        Math.sqrt(ic.position!), Math.sqrt(ic.position!), Math.sqrt(ic.velocity!),
        Math.sqrt(ic.heading!), Math.sqrt(ic.sideslip!),
        Math.sqrt(ic.accelBias!), Math.sqrt(ic.gyroBias!)
      ];
      for (let i = 0; i < N; i++) {
        if (!isFinite(this.S[i][i]) || this.S[i][i] <= 0) {
          this.S[i].fill(0);
          this.S[i][i] = sigmas[i];
        }
      }
      if (this.gpsInitialized) {
        this.coasting = true;
        this.lastGpsTimeMs = 0;
      }
    }
    // During coasting (no GPS), clamp Cholesky diagonal to prevent unbounded
    // covariance growth in very long tunnels.  These are physical plausibility
    // limits — position uncertainty beyond 500m or velocity uncertainty beyond
    // 50 m/s is not useful and risks numerical degradation.
    if (this.coasting) {
      const MAX_POS_SIGMA = 500;
      const MAX_VEL_SIGMA = 50;
      if (this.S[I.X][I.X] > MAX_POS_SIGMA) this.S[I.X][I.X] = MAX_POS_SIGMA;
      if (this.S[I.Y][I.Y] > MAX_POS_SIGMA) this.S[I.Y][I.Y] = MAX_POS_SIGMA;
      if (this.S[I.V][I.V] > MAX_VEL_SIGMA) this.S[I.V][I.V] = MAX_VEL_SIGMA;
      let tr = 0;
      for (let i = 0; i < N; i++) {
        for (let j = 0; j <= i; j++) {
          const v = this.S[i][j];
          tr += v * v;
        }
      }
      this._traceCache = tr;
    }
    if (!isFinite(this.adaNoiseScale) || this.adaNoiseScale < 0) this.adaNoiseScale = 1;
    if (!isFinite(this.robustWeight) || this.robustWeight < 0 || this.robustWeight > 1) this.robustWeight = 1;
  }

  private resetFromGps(gpsX: number, gpsY: number, gpsVx: number, gpsVy: number): void {
    const gpsV = Math.sqrt(gpsVx * gpsVx + gpsVy * gpsVy);
    const gpsPsi = gpsV > 1e-6 ? Math.atan2(gpsVy, gpsVx) : 0;
    // Preserve learned biases (aBiasX, gBiasZ) — they were calibrated during
    // coasting via ZUPT/ZARU and resetting them to 0 would re-introduce the
    // same systematic drift the filter had already corrected.  Only reset
    // kinematic states that are directly observable from GPS.
    const prevABiasX = this.x[I.A_BIAS_X];
    const prevGBiasZ = this.x[I.G_BIAS_Z];
    const prevDecl = this.x[I.MAG_DECL];
    this.x[I.X] = gpsX;
    this.x[I.Y] = gpsY;
    this.x[I.V] = gpsV;
    this.x[I.PSI] = this.wrapAngle(gpsPsi);
    this.x[I.BETA] = 0;
    this.x[I.A_BIAS_X] = prevABiasX;
    this.x[I.G_BIAS_Z] = prevGBiasZ;
    this.x[I.MAG_DECL] = prevDecl;
    const posSigma = Math.min(this.S[I.X][I.X], 5);
    const velSigma = this.config.measurementNoise.velocity!;
    const psiSigma = gpsV > 1 ? Math.min(velSigma / gpsV, 0.5) : 0.5;
    // Preserve bias covariances — they reflect the filter's confidence in the
    // learned bias values; resetting to initialCovariance would let them drift
    // freely again, undoing the calibration.
    const biasABSigma = this.S[I.A_BIAS_X][I.A_BIAS_X];
    const biasGBSigma = this.S[I.G_BIAS_Z][I.G_BIAS_Z];
    const declSigma = this.S[I.MAG_DECL][I.MAG_DECL];
    for (let i = 0; i < N; i++)
      for (let j = 0; j <= i; j++)
        this.S[i][j] = 0;
    this.S[I.X][I.X] = posSigma;
    this.S[I.Y][I.Y] = posSigma;
    this.S[I.V][I.V] = velSigma;
    this.S[I.PSI][I.PSI] = psiSigma;
    this.S[I.PSI][I.V] = 0.3 * psiSigma;
    this.S[I.BETA][I.BETA] = 0.1;
    this.S[I.A_BIAS_X][I.A_BIAS_X] = biasABSigma;
    this.S[I.G_BIAS_Z][I.G_BIAS_Z] = biasGBSigma;
    this.S[I.MAG_DECL][I.MAG_DECL] = declSigma;
    this.coasting = false;
    this.coastSpeedReady = false;
    this.bufTail = 0;
    this.bufLen = 0;
    this.imuTS.clear();
    this.stepTS.clear();
    this.smoothedSpeed = undefined;
    this.adaConvergeCount = 0;
    this.deviceToVehicle = null;
  }

  // ─── GPS latency compensation ──────────────────────────────────

  private saveBuf(t: number, ax: number, ay: number, gz: number, dt: number, az: number, gx: number, gy: number): void {
    const i = (this.bufTail + this.bufLen) % SrEkf.BUF_CAP;
    this.bufT[i] = t;
    const xi = i * N;
    for (let k = 0; k < N; k++) this.bufX[xi + k] = this.x[k];
    const si = i * NTRI;
    let p = 0;
    for (let row = 0; row < N; row++) {
      for (let col = 0; col <= row; col++) {
        this.bufS[si + p] = this.S[row][col];
        p++;
      }
    }
    const ini = i * 10;
    this.bufIn[ini] = ax; this.bufIn[ini + 1] = ay; this.bufIn[ini + 2] = gz;
    this.bufIn[ini + 3] = dt; this.bufIn[ini + 4] = az;
    this.bufIn[ini + 5] = gx; this.bufIn[ini + 6] = gy;
    this.bufIn[ini + 7] = this.curAzimuth;
    this.bufIn[ini + 8] = this.curPitch;
    this.bufIn[ini + 9] = this.curRoll;
    if (this.bufLen < SrEkf.BUF_CAP) {
      this.bufLen++;
    } else {
      this.bufTail = (this.bufTail + 1) % SrEkf.BUF_CAP;
    }
  }

  private restoreBuf(i: number): void {
    const xi = i * N;
    for (let k = 0; k < N; k++) this.x[k] = this.bufX[xi + k];
    const si = i * NTRI;
    let p = 0;
    for (let row = 0; row < N; row++) {
      for (let col = 0; col <= row; col++) {
        this.S[row][col] = this.bufS[si + p];
        p++;
      }
      for (let col = row + 1; col < N; col++) this.S[row][col] = 0;
    }
  }

  private rewindTo(gpsTime: number): number {
    if (this.bufLen === 0) return -1;
    let restoreOffset = -1;
    for (let off = this.bufLen - 1; off >= 0; off--) {
      const idx = (this.bufTail + off) % SrEkf.BUF_CAP;
      if (this.bufT[idx] <= gpsTime) { restoreOffset = off; break; }
    }
    if (restoreOffset < 0) return -1;
    const futureCount = this.bufLen - restoreOffset - 1;
    if (futureCount > 0) {
      for (let off = 0; off < futureCount; off++) {
        const srcIdx = (this.bufTail + restoreOffset + 1 + off) % SrEkf.BUF_CAP;
        const dstOff = off * 11;
        this.bufReplay[dstOff] = this.bufT[srcIdx];
        const ini = srcIdx * 10;
        this.bufReplay[dstOff + 1] = this.bufIn[ini];
        this.bufReplay[dstOff + 2] = this.bufIn[ini + 1];
        this.bufReplay[dstOff + 3] = this.bufIn[ini + 2];
        this.bufReplay[dstOff + 4] = this.bufIn[ini + 3];
        this.bufReplay[dstOff + 5] = this.bufIn[ini + 4];
        this.bufReplay[dstOff + 6] = this.bufIn[ini + 5];
        this.bufReplay[dstOff + 7] = this.bufIn[ini + 6];
        this.bufReplay[dstOff + 8] = this.bufIn[ini + 7];
        this.bufReplay[dstOff + 9] = this.bufIn[ini + 8];
        this.bufReplay[dstOff + 10] = this.bufIn[ini + 9];
      }
    }
    const restoreIdx = (this.bufTail + restoreOffset) % SrEkf.BUF_CAP;
    this.restoreBuf(restoreIdx);
    this.lastImuTimeMs = this.bufT[restoreIdx];
    this.imuTS.clear();
    this.axWindow.clear();
    this.ayWindow.clear();
    this.azWindow.clear();
    this.gzWindow.clear();
    this.gxWindow.clear();
    this.gyWindow.clear();
    this.stepTS.clear();
    this.stepBuffer.clear();
    this.stepFreq = 0;
    this.accelEnergy = 0;
    this.gyroEnergy = 0;
    this.varAccelEnergy = 0;
    this.varGyroEnergy = 0;
    this.angAccel = 0;
    this.smoothAngAccel = 0;
    this.bufLen = restoreOffset + 1;
    // bufTail unchanged — prior history preserved for re-rewind
    return futureCount;
  }

  private replayFromScratch(count: number): void {
    const prevRot = this.deviceToEnu;
    const prevAz = this.curAzimuth, prevPi = this.curPitch, prevRo = this.curRoll;
    for (let off = 0; off < count; off++) {
      const srcOff = off * 11;
      const az = this.bufReplay[srcOff + 8];
      const pi = this.bufReplay[srcOff + 9];
      const ro = this.bufReplay[srcOff + 10];
      if (!isNaN(az)) {
        this.setOrientation(az, pi, ro);
      } else {
        this.deviceToEnu = null;
        this.curAzimuth = NaN; this.curPitch = NaN; this.curRoll = NaN;
      }
      this.predict(
        this.bufReplay[srcOff + 1],
        this.bufReplay[srcOff + 2],
        this.bufReplay[srcOff + 3],
        this.bufReplay[srcOff + 4],
        this.bufReplay[srcOff],
        this.bufReplay[srcOff + 5],
        this.bufReplay[srcOff + 6],
        this.bufReplay[srcOff + 7]
      );
    }
    this.deviceToEnu = prevRot;
    this.curAzimuth = prevAz; this.curPitch = prevPi; this.curRoll = prevRo;
  }

  // ─── CTRA kinematics ────────────────────────────────────────────

  private ctraDelta(psi: number, v: number, omega: number, dt: number): [number, number] {
    if (Math.abs(omega) > EPS) {
      const sp = Math.sin(psi), cp = Math.cos(psi);
      const spw = Math.sin(psi + omega * dt), cpw = Math.cos(psi + omega * dt);
      return [v / omega * (spw - sp), v / omega * (-cpw + cp)];
    }
    const cp = Math.cos(psi), sp = Math.sin(psi);
    return [v * cp * dt, v * sp * dt];
  }

  private computeJacobian(a: number, omega: number, dt: number, betaTau: number, aBiasDecay: number): void {
    const psi = this.x[I.PSI], beta = this.x[I.BETA];
    const psiBeta = psi + beta, v = this.x[I.V];
    const vAvg = v + 0.5 * a * dt;
    for (let i = 0; i < N; i++) this.tmpF[i].fill(0);
    for (let i = 0; i < N; i++) this.tmpF[i][i] = 1;

    const absOmega = Math.abs(omega);
if (absOmega > EPS) {
      const sp = Math.sin(psiBeta), cp = Math.cos(psiBeta);
      const spw = Math.sin(psiBeta + omega * dt), cpw = Math.cos(psiBeta + omega * dt);
      const o2 = omega * omega;
      this.tmpF[I.X][I.V] = (spw - sp) / omega;
      this.tmpF[I.X][I.PSI] = vAvg / omega * (cpw - cp);
      this.tmpF[I.X][I.BETA] = this.tmpF[I.X][I.PSI];
      this.tmpF[I.X][I.G_BIAS_Z] = -vAvg * ((spw - sp) / o2 - dt * cpw / omega);
      this.tmpF[I.Y][I.V] = (-cpw + cp) / omega;
      this.tmpF[I.Y][I.PSI] = vAvg / omega * (spw - sp);
      this.tmpF[I.Y][I.BETA] = this.tmpF[I.Y][I.PSI];
      this.tmpF[I.Y][I.G_BIAS_Z] = -vAvg * ((-cpw + cp) / o2 - dt * spw / omega);
    } else {
      const cp = Math.cos(psiBeta), sp = Math.sin(psiBeta);
      const dt2 = dt * dt;
      this.tmpF[I.X][I.V] = cp * dt;
      this.tmpF[I.X][I.PSI] = -vAvg * sp * dt;
      this.tmpF[I.X][I.BETA] = this.tmpF[I.X][I.PSI];
      this.tmpF[I.X][I.G_BIAS_Z] = -0.5 * vAvg * sp * dt2;
      this.tmpF[I.Y][I.V] = sp * dt;
      this.tmpF[I.Y][I.PSI] = vAvg * cp * dt;
      this.tmpF[I.Y][I.BETA] = this.tmpF[I.Y][I.PSI];
      this.tmpF[I.Y][I.G_BIAS_Z] = 0.5 * vAvg * cp * dt2;
    }
    this.tmpF[I.X][I.A_BIAS_X] = -0.5 * dt * this.tmpF[I.X][I.V];
    this.tmpF[I.Y][I.A_BIAS_X] = -0.5 * dt * this.tmpF[I.Y][I.V];
    this.tmpF[I.V][I.A_BIAS_X] = -dt;
    this.tmpF[I.PSI][I.G_BIAS_Z] = -dt;
    this.tmpF[I.BETA][I.BETA] = Math.exp(-dt / betaTau);
    this.tmpF[I.A_BIAS_X][I.A_BIAS_X] = aBiasDecay;
    this.tmpF[I.MAG_DECL][I.MAG_DECL] = 1;
  }

  // ─── Adaptive process noise ──────────────────────────────────────

  private updateStepDetection(ax: number, timestampMs: number): void {
    this.stepTS.push(timestampMs);
    this.stepBuffer.push(ax);
    const cutoff = timestampMs - 1500;
    while (this.stepTS.length > 0 && this.stepTS.get(0) < cutoff) {
      this.stepTS.shift();
      this.stepBuffer.shift();
    }
    if (this.stepBuffer.length < 20) return;
    let sum = 0;
    for (let i = 0; i < this.stepBuffer.length; i++) sum += this.stepBuffer.get(i);
    const mean = sum / this.stepBuffer.length;
    let crossings = 0;
    for (let i = 1; i < this.stepBuffer.length; i++)
      if ((this.stepBuffer.get(i - 1) - mean) * (this.stepBuffer.get(i) - mean) < 0) crossings++;
    const windowDur = (this.stepTS.get(this.stepTS.length - 1) - this.stepTS.get(0)) / 1000;
    this.stepFreq = crossings > 0 ? crossings / (2 * windowDur) : 0;
  }

  private updateImuWindow(ax: number, ay: number, gz: number, az: number | undefined, gx: number | undefined, gy: number | undefined, timestampMs: number): void {
    this.imuTS.push(timestampMs);
    this.axWindow.push(ax);
    this.ayWindow.push(ay);
    this.gzWindow.push(gz);
    const cutoff = timestampMs - 1000;
    if (az !== undefined) this.azWindow.push(az);
    if (gx !== undefined) this.gxWindow.push(gx);
    if (gy !== undefined) this.gyWindow.push(gy);
    while (this.imuTS.length > 0 && this.imuTS.get(0) < cutoff) {
      this.imuTS.shift();
      this.axWindow.shift();
      this.ayWindow.shift();
      this.gzWindow.shift();
      if (this.azWindow.length > 0) this.azWindow.shift();
      if (this.gxWindow.length > 0) this.gxWindow.shift();
      if (this.gyWindow.length > 0) this.gyWindow.shift();
    }
  }

  private windowVariance(buf: RingBuf): number {
    if (buf.length < 2) return 0;
    let sum = 0, sumSq = 0;
    for (let i = 0; i < buf.length; i++) { const v = buf.get(i); sum += v; sumSq += v * v; }
    const n = buf.length, mean = sum / n;
    return Math.max(0, sumSq / n - mean * mean);
  }

  getStillness(): number {
    // Exponential stillness: still = exp(-(σ_a²/T_a² + σ_ω²/T_ω²))
    // σ_a² = component-wise acceleration variance (bias-invariant — constant
    // gravity projection has zero variance and doesn't reduce stillness)
    // σ_ω² = mean ||ω||² (gyro ENERGY catches constant-rate turns)
    // T_a = 0.3 m/s², T_w = 0.05 rad/s (tuned for phone IMU at 5 Hz)
    // Returns 0.5 (neutral) when < 5 samples in the window.
    if (this.axWindow.length < 5) return 0.5;
    const T_a = 0.3, T_w = 0.05;
    return Math.exp(-(this._accel3DVar / (T_a * T_a) + this._gyro3DEnergy / (T_w * T_w)));
  }

  private zuptChiSqGate(): boolean {
    let pvv = 0;
    for (let k = 0; k <= I.V; k++) pvv += this.S[I.V][k] * this.S[I.V][k];
    const rVel = 0.01;
    const innovVSq = this.x[I.V] * this.x[I.V];
    return innovVSq / (pvv + rVel * rVel) <= 9.0;
  }

  private computeAdaptiveQ(dt: number, a: number, omega: number): void {
    const sqrtDt = Math.sqrt(dt);
    const speedScale = Math.min(Math.sqrt(Math.max(Math.abs(this.x[I.V]), 0.05) / 5), 2);
    for (let i = 0; i < N; i++) this.tmpSqrtQ[i].fill(0);

    // Single-pass window variance for all 3 buffers
    let sumAx = 0, sumAy = 0, sumGz = 0, sqAx = 0, sqAy = 0, sqGz = 0;
    let sumAz = 0, sqAz = 0, gyro3dSumSq = 0;
    const n = this.axWindow.length;
    const hasZ = this.azWindow.length === n;
    const hasGxy = this.gxWindow.length === n && this.gyWindow.length === n;
    for (let i = 0; i < n; i++) {
      const ax = this.axWindow.get(i), ay = this.ayWindow.get(i), gz = this.gzWindow.get(i);
      sumAx += ax; sumAy += ay; sumGz += gz;
      sqAx += ax * ax; sqAy += ay * ay; sqGz += gz * gz;
      if (hasZ) { const az = this.azWindow.get(i); sumAz += az; sqAz += az * az; }
      if (hasGxy) {
        const gx = this.gxWindow.get(i), gy = this.gyWindow.get(i);
        gyro3dSumSq += gx * gx + gy * gy + gz * gz;
      } else {
        gyro3dSumSq += gz * gz;
      }
    }
    this._varAx = n >= 2 ? Math.max(0, sqAx / n - (sumAx / n) ** 2) : 0;
    this._varAy = n >= 2 ? Math.max(0, sqAy / n - (sumAy / n) ** 2) : 0;
    this._varGz = n >= 2 ? Math.max(0, sqGz / n - (sumGz / n) ** 2) : 0;
    const varAx = this._varAx, varAy = this._varAy, varGz = this._varGz;
    let accel3dVar = varAx + varAy;
    if (hasZ) accel3dVar += Math.max(0, sqAz / n - (sumAz / n) ** 2);
    this._accel3DVar = accel3dVar;
    this._gyro3DEnergy = n > 0 ? gyro3dSumSq / n : 0;

    const sqrtAccelVar = Math.sqrt(Math.max(varAx + varAy, 0)) / 5;
    // variance-based energy catches transient motion; |a| catches constant
    // braking/deceleration (zero variance but sustained acceleration)
    const rawAccel = Math.max(sqrtAccelVar, Math.abs(a) / 5);
    this.accelEnergy = 0.9 * this.accelEnergy + 0.1 * Math.min(rawAccel, 5);

    // variance-only energy (for stillness/ZUPT/ZARU): a constant bias has zero
    // variance and must NOT count as motion
    this.varAccelEnergy = 0.9 * this.varAccelEnergy + 0.1 * Math.min(sqrtAccelVar, 5);

    const absOmega = Math.abs(omega);
    const rawGyro = Math.max(Math.sqrt(Math.max(varGz, 0)), absOmega * 0.1) / 0.5;
    this.gyroEnergy = 0.9 * this.gyroEnergy + 0.1 * Math.min(rawGyro, 5);
    this.varGyroEnergy = 0.9 * this.varGyroEnergy + 0.1 * Math.min(rawGyro, 5);

    this.stepEnergy = Math.min(this.stepFreq / 3, 1);
    const stepEnergy = this.stepEnergy;
    const sc = this.config.adaptiveScaling;

    const pn = this.config.processNoise;
    const qBase = pn.position! * sqrtDt * speedScale * (1 + sc.positionAccel! * this.accelEnergy);
    const qFwd2 = qBase * qBase;
    const qCross2 = (qBase * 0.3) * (qBase * 0.3);
    const psiBeta = this.x[I.PSI] + this.x[I.BETA];
    const cp = Math.cos(psiBeta), sp = Math.sin(psiBeta);
    const qxx = cp * cp * qFwd2 + sp * sp * qCross2;
    const qyy = sp * sp * qFwd2 + cp * cp * qCross2;
    const qxy = cp * sp * (qFwd2 - qCross2);
    this.tmpSqrtQ[I.X][I.X] = Math.sqrt(qxx);
    this.tmpSqrtQ[I.Y][I.X] = qxy / Math.max(this.tmpSqrtQ[I.X][I.X], 1e-12);
    this.tmpSqrtQ[I.Y][I.Y] = Math.sqrt(Math.max(qyy - this.tmpSqrtQ[I.Y][I.X] * this.tmpSqrtQ[I.Y][I.X], 0));
    // Angular acceleration Q boost: during corner entry/exit (high |dω/dt|), inflate
    // heading/sideslip/gyroBias Q so the filter trusts its state less and allows faster
    // GPS-driven correction. Prevents heading from "getting stuck" during transients.
    // Uses smoothAngAccel to sustain the boost for ~0.5s after the rotation stops.
    const angAccelBoost = Math.min(this.smoothAngAccel / 2.0, 1);  // 0→1 for 0→2 rad/s²
    this.tmpSqrtQ[I.V][I.V] = pn.velocity! * sqrtDt * speedScale * (1 + sc.velocityAccel! * this.accelEnergy + sc.velocityStep! * stepEnergy);
    this.tmpSqrtQ[I.PSI][I.PSI] = pn.heading! * sqrtDt * (1 + sc.headingGyro! * this.gyroEnergy + sc.headingStep! * stepEnergy + 0.3 * this.accelEnergy + 1.5 * angAccelBoost + 0.3 * absOmega);
    this.tmpSqrtQ[I.BETA][I.BETA] = pn.sideslip! * sqrtDt * (1 + sc.sideslipGyro! * this.gyroEnergy + sc.sideslipStep! * stepEnergy + 2.0 * angAccelBoost + 0.5 * absOmega);
    this.tmpSqrtQ[I.A_BIAS_X][I.A_BIAS_X] = pn.accelBias! * sqrtDt;
    this.tmpSqrtQ[I.G_BIAS_Z][I.G_BIAS_Z] = pn.gyroBias! * sqrtDt * (1 + 0.3 * this.gyroEnergy + 0.5 * angAccelBoost);
    this.tmpSqrtQ[I.MAG_DECL][I.MAG_DECL] = (pn.magDeclination ?? 1e-4) * sqrtDt;

    // During coasting (no GPS corrections), reduce position/velocity Q to slow
    // covariance growth.  Without GPS, bias-driven position error is already
    // captured in P via the state cross-covariance; adding full Q on top
    // over-inflates uncertainty.  Heading/bias Q are kept at full strength —
    // ZARU/ZUPT/Mag still provide corrections in these channels.
    if (this.coasting) {
      const COAST_Q_FACTOR = 0.3;
      this.tmpSqrtQ[I.X][I.X] *= COAST_Q_FACTOR;
      this.tmpSqrtQ[I.Y][I.X] *= COAST_Q_FACTOR;
      this.tmpSqrtQ[I.Y][I.Y] *= COAST_Q_FACTOR;
      this.tmpSqrtQ[I.V][I.V] *= COAST_Q_FACTOR;
    }
  }

  // ─── GPS measurement update ──────────────────────────────────────

  private computeH(cp: number, sp: number): void {
    const v = this.x[I.V];
    for (let i = 0; i < M; i++) this.tmpH[i].fill(0);
    this.tmpH[0][I.X] = 1;
    this.tmpH[1][I.Y] = 1;
    this.tmpH[2][I.V] = cp;   this.tmpH[2][I.PSI] = -v * sp;  this.tmpH[2][I.BETA] = -v * sp;
    this.tmpH[3][I.V] = sp;   this.tmpH[3][I.PSI] = v * cp;   this.tmpH[3][I.BETA] = v * cp;
    // At low speed, GPS velocity direction is unreliable (multipath, buildings,
    // slow city driving). Ramp starts at v=0.2 m/s (down from 0.5 to catch slow
    // corner exits) and reaches full at v=3.0 m/s (~11 km/h). A smoothAngAccel boost
    // temporarily raises effective speed after a turn so heading corrects faster when
    // it matters most (recovering from corner-exit lag).
    // NOTE: the GPS-velocity heading gain MUST remain non-trivial at low speed
    // (e.g. v=1 m/s) — it is the signal that lets magnetic-declination calibrate
    // (the mag-vs-GPS-velocity heading disagreement is absorbed into magDeclination).
    let headingGain = v < 0.2 ? 0 : Math.min((v - 0.2 + this.gyroEnergy * 0.5 + this.stepEnergy * 1.5 + this.smoothAngAccel * 0.5) / 2.8, 1);

    // Adaptive initialization boost: during first 30s after GPS init, if heading
    // variance is high, temporarily boost GPS correction to accelerate convergence.
    if (this.gpsInitTimeMs > 0 && this.gpsInitialized) {
      const timeSinceGpsInit = Math.max(0, this.lastImuTimeMs - this.gpsInitTimeMs);
      if (timeSinceGpsInit < 30000) {  // First 30 seconds
        let psiVarGate = 0;
        for (let k = 0; k <= I.PSI; k++) psiVarGate += this.S[I.PSI][k] * this.S[I.PSI][k];
        const psiStd = Math.sqrt(psiVarGate);
        // Boost gain by up to 2× when psiStd is high (uncertain heading)
        const initBoost = Math.max(0, (psiStd - 0.3) / 0.7);  // Ramps from 0 at 0.3 to 1 at 1.0
        headingGain = Math.min(headingGain * (1 + 2 * initBoost), 1);
      }
    }

    if (headingGain < 1) {
      this.tmpH[2][I.PSI] *= headingGain; this.tmpH[2][I.BETA] *= headingGain;
      this.tmpH[3][I.PSI] *= headingGain; this.tmpH[3][I.BETA] *= headingGain;
    }
  }

  private computeGpsInnovation(z: Float64Array): void {
    const psi = this.x[I.PSI], beta = this.x[I.BETA], v = this.x[I.V];
    const psiBeta = psi + beta;
    this.tmpInnov[0] = z[0] - this.x[I.X];
    this.tmpInnov[1] = z[1] - this.x[I.Y];
    this.tmpInnov[2] = z[2] - v * Math.cos(psiBeta);
    this.tmpInnov[3] = z[3] - v * Math.sin(psiBeta);

    const A = this.tmpA;
    for (let i = 0; i < M; i++) {
      const Ai = A[i];
      for (let j = 0; j < N; j++) {
        let s = 0;
        for (let k = j; k < N; k++) s += this.tmpH[i][k] * this.S[k][j];
        Ai[j] = s;
      }
    }

    for (let i = 0; i < M; i++) {
      const Ai = A[i];
      for (let j = 0; j <= i; j++) {
        let s = 0;
        for (let k = 0; k < N; k++) s += Ai[k] * A[j][k];
        this.tmpWork4x4[i][j] = s;
        this.tmpWork4x4[j][i] = s;
        this.tmpHPH[i][j] = s;
        this.tmpHPH[j][i] = s;
      }
    }
  }

  private computeGpsPostFit(posR: number, velR: number, cp: number, sp: number): number {
    // Anisotropic position R: along-track σ = 0.5·posR, cross-track σ = 1.33·posR
    const ra = posR * 0.5, rc = posR * 1.33;
    const ra2 = ra * ra, rc2 = rc * rc;
    this.tmpWork4x4[0][0] += cp * cp * ra2 + sp * sp * rc2;
    this.tmpWork4x4[1][1] += sp * sp * ra2 + cp * cp * rc2;
    this.tmpWork4x4[0][1] += cp * sp * (ra2 - rc2);
    this.tmpWork4x4[1][0] += cp * sp * (ra2 - rc2);
    this.tmpWork4x4[2][2] += velR * velR;
    this.tmpWork4x4[3][3] += velR * velR;

    if (!chol4x4(this.tmpL4x4, this.tmpWork4x4)) return Infinity;

    const y = this.tmpW;
    for (let i = 0; i < M; i++) y[i] = this.tmpInnov[i];
    cholSolve4(this.tmpL4x4, y);
    let chiSq = 0;
    for (let i = 0; i < M; i++) chiSq += this.tmpInnov[i] * y[i];
    return chiSq;
  }

  private gpsUpdateSingle(posR: number, velR: number): boolean {
    const psiBeta0 = this.x[I.PSI] + this.x[I.BETA];
    const cp0 = Math.cos(psiBeta0), sp0 = Math.sin(psiBeta0);
    this.computeH(cp0, sp0);
    // Doppler heading direction is noisy/ambiguous. Smoothly block the heading
    // columns of H so velocity innovation cannot rotate ψ.  The magnetometer
    // continues to correct heading at these speeds (magHeadingTrust > 0 for
    // v<1.5), creating a clean handoff: mag owns heading at rest, GPS velocity
    // direction owns heading at speed, with smooth transition between them.
    if (this.accelEnergy + this.gyroEnergy < 0.05) {
      const restW = Math.max(0, Math.min((1.5 - Math.abs(this.x[I.V])) / 1.0, 1));
      this.tmpH[2][I.PSI] *= (1 - restW); this.tmpH[2][I.BETA] *= (1 - restW);
      this.tmpH[3][I.PSI] *= (1 - restW); this.tmpH[3][I.BETA] *= (1 - restW);
    }
    this.computeGpsInnovation(this.tmpZ);

    // Compute separate position (2-DOF) and velocity (2-DOF) χ² for gating
    // Uses the same anisotropic model as computeGpsPostFit:
    // along-track σ = 0.5·posR, cross-track σ = 1.33·posR, rotated by heading.
    // Compute inline 2×2 Cholesky solve for each sub-block of tmpHPH.
    const sgCp = cp0, sgSp = sp0;
    const sgRa = posR * 0.5, sgRc = posR * 1.33;
    const sgRa2 = sgRa * sgRa, sgRc2 = sgRc * sgRc;
    const sgRxx = sgCp * sgCp * sgRa2 + sgSp * sgSp * sgRc2;
    const sgRxy = sgCp * sgSp * (sgRa2 - sgRc2);
    const sgRyy = sgSp * sgSp * sgRa2 + sgCp * sgCp * sgRc2;
    const a00 = this.tmpHPH[0][0] + sgRxx;
    const a01 = this.tmpHPH[0][1] + sgRxy;
    const a11 = this.tmpHPH[1][1] + sgRyy;
    const l0 = Math.sqrt(a00);
    const l1 = a01 / l0;
    const l2 = Math.sqrt(Math.max(a11 - l1 * l1, 0));
    let chiSqPos = Infinity;
    if (l2 > 1e-30) {
      const r0 = this.tmpInnov[0] / l0;
      const t1 = (this.tmpInnov[1] - l1 * r0) / l2;
      chiSqPos = r0 * r0 + t1 * t1;
    }

    const b00 = this.tmpHPH[2][2] + velR * velR;
    const b01 = this.tmpHPH[2][3];
    const b11 = this.tmpHPH[3][3] + velR * velR;
    const m0 = Math.sqrt(b00);
    const m1 = b01 / m0;
    const m2 = Math.sqrt(Math.max(b11 - m1 * m1, 0));
    let chiSqVel = Infinity;
    if (m2 > 1e-30) {
      const s0 = this.tmpInnov[2] / m0;
      const u1 = (this.tmpInnov[3] - m1 * s0) / m2;
      chiSqVel = s0 * s0 + u1 * u1;
    }

    let chiSq = this.computeGpsPostFit(posR, velR, cp0, sp0);
    this.lastChiSq = chiSq;

    // 180° flip recovery (when the GPS velocity direction is clearly reversed
    // relative to the current heading). Lowered from v>2.0 to v>0.8: city
    // stop-and-go speeds are often below 2.0 m/s. The strict anti-parallel
    // dot-product check below already prevents 90° pedestrian-turn misfires.
    // The χ² gate is intentionally NOT required here: a pure velocity-direction
    // reversal (position still consistent, e.g. urban multipath on Doppler) has
    // low total χ² yet must still be corrected, otherwise the Kalman velocity
    // update resolves the v/heading sign ambiguity into the wrong (negative-v)
    // branch — the reported "U-turn → negative velocity" bug.
    if (this.x[I.V] > 0.8) {
      const psiBeta = this.x[I.PSI] + this.x[I.BETA];
      const vxP = this.x[I.V] * Math.cos(psiBeta);
      const vyP = this.x[I.V] * Math.sin(psiBeta);
      const vxM = this.tmpZ[2], vyM = this.tmpZ[3];
      const speedM2 = vxM * vxM + vyM * vyM;
      if (speedM2 > 1.0 && vxP * vxM + vyP * vyM < -0.5 * this.x[I.V] * this.x[I.V]) {
        this.x[I.PSI] = this.wrapAngle(this.x[I.PSI] + Math.PI);
        for (let j = 0; j < I.PSI; j++) this.S[I.PSI][j] = 0;
        for (let j = I.PSI + 1; j < N; j++) this.S[j][I.PSI] = 0;
        this.S[I.PSI][I.PSI] = Math.max(this.S[I.PSI][I.PSI], 1.5);
        const fp = this.x[I.PSI] + this.x[I.BETA];
        const cpF = Math.cos(fp), spF = Math.sin(fp);
        this.computeH(cpF, spF);
        this.computeGpsInnovation(this.tmpZ);
        chiSq = this.computeGpsPostFit(posR, velR, cpF, spF);
        this.lastChiSq = chiSq;
      }
    }

    let robustW = 1;
    const rw = this.config.robustWeight;
    if (rw.enabled) {
      const thr = Math.max(rw.threshold!, 1e-3);
      robustW = rw.type === 'huber'
        ? (chiSq > thr ? thr / chiSq : 1)
        : 1 / (1 + (chiSq / thr) * (chiSq / thr));
    }
    this.robustWeight = isFinite(robustW) ? robustW : 1;

    // Per-component velocity robustness (noisy city GPS heading/speed).
    // GPS heading is inferred from Doppler velocity via H[2:3][PSI,BETA]; in
    // urban canyons multipath corrupts velocity while position stays usable.
    // The joint weight above conflates the two: bad velocity would down-weight
    // the good position fix, and moderate velocity noise below the joint gate
    // passes through at full strength and jitters heading. Derive a velocity-
    // only weight from the 2-DOF velocity χ² (threshold scaled 4-DOF→2-DOF by
    // 0.63 = χ²₀.₉₅,₂ / χ²₀.₉₅,₄) and apply it to velR alone, so noisy Doppler
    // is distrusted for heading without sacrificing position. No-op on clean
    // GPS (velRobustW = 1).
    let velRobustW = 1;
    if (rw.enabled && isFinite(chiSqVel)) {
      const thrVel = Math.max(rw.threshold! * 0.63, 1e-3);
      velRobustW = rw.type === 'huber'
        ? (chiSqVel > thrVel ? thrVel / chiSqVel : 1)
        : 1 / (1 + (chiSqVel / thrVel) * (chiSqVel / thrVel));
    }

    let adaScale = 1;
    if (this.config.adaptiveNoise.enabled && isFinite(chiSq)) {
      const innovRatio = chiSq / M;
      const α = this.config.adaptiveNoise.smoothing!;
      this.adaNoiseScale = α * innovRatio + (1 - α) * this.adaNoiseScale;
      if (innovRatio < 0.5) {
        this.adaConvergeCount++;
        if (this.adaConvergeCount >= 10) { this.adaNoiseScale = 1; this.adaConvergeCount = 0; }
      } else {
        this.adaConvergeCount = 0;
      }
      adaScale = Math.max(1, Math.min(this.adaNoiseScale, this.config.adaptiveNoise.maxScale!));
    }

    const totalWeight = Math.max(robustW / adaScale, 1e-8);
    if (totalWeight < 1 || velRobustW < 1) {
      for (let i = 0; i < M; i++)
        for (let j = 0; j <= i; j++) {
          this.tmpWork4x4[i][j] = this.tmpHPH[i][j];
          this.tmpWork4x4[j][i] = this.tmpHPH[i][j];
        }
      posR /= Math.sqrt(totalWeight);
      velR /= Math.sqrt(totalWeight * velRobustW);
      chiSq = this.computeGpsPostFit(posR, velR, cp0, sp0);
    }

    // Separate position/velocity gating: if joint gate fails, position alone
    // still deserves a pass (position has higher trust than GPS velocity).
    if (chiSq > this.config.gateThreshold) {
      if (!rw.enabled || totalWeight < 0.1) {
        if (chiSqPos <= this.config.gateThreshold) {
          // Position sub-gate passes → accept despite noisy velocity
          // BUT if coasting after GPS loss, force full reset to fix
          // heading/velocity that may be corrupted after extended IMU-only prediction
          if (this.coasting) {
            this.resetFromGps(this.tmpZ[0], this.tmpZ[1], this.tmpZ[2], this.tmpZ[3]);
            return true;
          }
        } else {
          // Both joint and position fail → reject or reset
          if (this.coasting) {
            this.resetFromGps(this.tmpZ[0], this.tmpZ[1], this.tmpZ[2], this.tmpZ[3]);
            return true;
          }
          if (!rw.enabled) return false;
        }
      }
    }

    // state update: Δx = S · Aᵀ · L⁻ᵀ · L⁻¹ · innov
    const w = this.tmpW;
    const z = this.tmpBuf;
    for (let j = 0; j < N; j++) {
      let s = 0;
      for (let i = 0; i < M; i++) s += this.tmpA[i][j] * w[i];
      z[j] = s;
    }
    for (let i = 0; i < N; i++) {
      let s = 0;
      for (let j = 0; j <= i; j++) s += this.S[i][j] * z[j];
      this.x[i] += s;
    }
    this.x[I.PSI] = this.wrapAngle(this.x[I.PSI]);

    // Velocity sign-ambiguity backstop. The state (v<0, ψ) and (v>0, ψ+π)
    // represent the identical velocity vector. After a U-turn with heading
    // lag (or noisy/ambiguous GPS velocity), the filter can settle into the
    // negative-v branch, which the user reads as "velocity went negative and
    // heading snapped to the car's nose". Reparameterize into the positive-v
    // branch so that v ≥ 0 and ψ is the direction of motion. This keeps a
    // consistent covariance: flipping the sign of one state coordinate negates
    // the corresponding Cholesky row and column (the diagonal is flipped twice
    // and stays positive). Genuine sustained reverse (v already < 0.8 when the
    // reversal occurs) is handled upstream by the 180° flip-recovery, which
    // flips ψ before v can go negative — so this block only catches residual
    // artifacts and does not oscillate.
    // IMPORTANT: do NOT reparameterize at negligible speed. When the device is
    // stationary, the GPS velocity measurement is pure noise and the Kalman
    // update can drive v to a tiny negative value (e.g. −1e-4). Flipping ψ by
    // π there would make the heading shake 180° every step as the sign of v
    // jitters around zero. So only the genuine (meaningfully non-zero) branch
    // is reparameterized; at |v| < 0.5 m/s we simply clamp v to 0. The
    // higher threshold (was 0.15) prevents spurious 180° flips when GPS
    // noise briefly pushes v to −0.2..−0.4 at stationary, where the true
    // heading was already correct (aligned with compass). Only truly
    // significant negative v (>0.5, e.g. extended IMU drift) gets flipped.
    if (this.x[I.V] < 0) {
      if (-this.x[I.V] < 0.5) {
        this.x[I.V] = 0;
      } else {
        this.x[I.V] = -this.x[I.V];
        this.x[I.PSI] = this.wrapAngle(this.x[I.PSI] + Math.PI);
        this.x[I.BETA] = this.wrapAngle(this.x[I.BETA] + Math.PI);
        for (let j = 0; j <= I.V; j++) this.S[I.V][j] = -this.S[I.V][j];
        for (let i = I.V; i < N; i++) this.S[i][I.V] = -this.S[i][I.V];
      }
    }

    // QR covariance update
    for (let i = 0; i < PRE; i++) this.tmpPreA[i].fill(0);

    // Anisotropic position R Cholesky factor (2×2 upper triangular)
    const psiBetaR = this.x[I.PSI] + this.x[I.BETA];
    const cpR = Math.cos(psiBetaR), spR = Math.sin(psiBetaR);
    const raL = posR * 0.5, rcL = posR * 1.33;
    const raL2 = raL * raL, rcL2 = rcL * rcL;
    const rxx = cpR * cpR * raL2 + spR * spR * rcL2;
    const rxy = cpR * spR * (raL2 - rcL2);
    const ryy = spR * spR * raL2 + cpR * cpR * rcL2;
    const l00 = Math.sqrt(rxx);
    const l01 = rxy / l00;
    const l11 = Math.sqrt(Math.max(ryy - l01 * l01, 0));
    this.tmpPreA[0][0] = l00; this.tmpPreA[0][1] = l01;
    this.tmpPreA[1][1] = l11;
    this.tmpPreA[2][2] = velR; this.tmpPreA[3][3] = velR;

    for (let i = 0; i < M; i++)
      for (let j = 0; j < N; j++)
        this.tmpPreA[i][M + j] = this.tmpA[i][j];

    for (let i = 0; i < N; i++)
      for (let j = 0; j < N; j++)
        this.tmpPreA[M + i][M + j] = this.S[i][j];

    for (let i = 0; i < PRE; i++)
      for (let j = 0; j < PRE; j++)
        this.tmpPreAT[i][j] = this.tmpPreA[j][i];

    this.qrInPlace(this.tmpPreAT, PRE, PRE, this.tmpHouseV);
    this.copySfromQR(this.tmpPreAT, M);
    this.safeguardState();
    return true;
  }

  // ─── Magnetometer update ─────────────────────────────────────────

  private magUpdateSingle(bearing: number): void {
    // Rest-rate cross-check (drifting-compass guard): the magnetometer only makes
    // sense as a heading reference if its angular motion is confirmed by the gyro.
    // If the compass bears a rotation the bias-corrected gyro rate does NOT
    // confirm, the compass is drifting (uncalibrated mag, hard-iron, or table
    // interference) — NOT real motion. Fusing it would rotate ψ continuously on a
    // stationary table (the reported "heading rotates" bug, with GPS present: at
    // rest v≈0 so GPS velocity carries no heading info to counter the mag). Skip
    // such updates. Genuine turns pass because there the gyro rate is high and the
    // compass rotation matches it. This runs before the init snap / trust gate so
    // a brief stillness dip during a transition cannot let a drift step through.
    // The check compares the compass angular RATE between consecutive calls against
    // the gyro's rate; a stable compass has ~0 rate and passes, while any drift
    // above a small margin (well above gyro noise ≈5e-4 rad/s) is blocked.
    if (this.prevCallMagTimeMs > 0) {
      const dtCall = (this.lastMagTimeMs - this.prevCallMagTimeMs) / 1000;
      if (dtCall > 1e-3 && dtCall < 60.0) {
        let dBearing = bearing - this.prevCallMagBearing;
        while (dBearing > Math.PI) dBearing -= 2 * Math.PI;
        while (dBearing < -Math.PI) dBearing += 2 * Math.PI;
        const magRate = Math.abs(dBearing / dtCall);
        const gyroRate = Math.abs(this.lastOmega);
        // Speed-dependent margin: at v≈0 the gyro yaw rate (rr.alpha) may be
        // unavailable on some Android devices, so phone rotation is only visible
        // to the magnetometer. GPS velocity direction carries no heading info
        // at rest, making the mag the sole heading reference. Widen the margin
        // to let genuine phone rotation pass.  At speed the gyro is reliable
        // and the standard tight margin protects against mag drift.
        const v = Math.abs(this.x[I.V]);
        const margin = v < 0.1 ? 2.0 : 0.05;
        if (magRate > gyroRate + margin) {
          this.prevCallMagBearing = bearing;
          this.prevCallMagTimeMs = this.lastMagTimeMs;
          return; // drifting compass — do not apply
        }
      }
    }
    this.prevCallMagBearing = bearing;
    this.prevCallMagTimeMs = this.lastMagTimeMs;

    // Gyro-confirmed rotation guard: when the gyro measures a genuine rotation
    // (|ω| > 0.05 rad/s), the gyro is the reliable heading reference and a stale
    // or lagging compass must not yank ψ back toward the corner-entry direction —
    // the reported "heading stuck at corner entry, then snapping to the moving
    // direction" bug. At low speed the EKF velocity can read ≈0 (ZUPT drags v to
    // ~0 during a crawl corner), so the v-based protections above (trust gate,
    // restMag direct blend, widened drift margin) are all ineffective exactly when
    // the car is genuinely rotating. Gate the mag's ψ-correction on the GYRO rate
    // instead: skip the update (blend AND Kalman branches) while rotating and let
    // the gyro integrate the turn; the compass re-anchors once rotation stops.
    // This only fires when the gyro actually reports rotation — on devices where
    // the gyro yaw rate is unavailable at rest, |ω|≈0 and the mag keeps full
    // authority for phone-rotation tracking on a table.
    if (Math.abs(this.lastOmega) > 0.05) return;

    // magHeadingTrust: 1 at low speed (mag owns ψ for init and backup), → 0 at
    // speed (GPS owns live heading ψ via its velocity-direction measurement).
    // At speed the magnetometer is SKIPPED entirely (including the init snap),
    // so it cannot fight the GPS velocity-direction update for ψ. The
    // magnetometer has no β term in its observation model, so whenever sideslip
    // β ≠ 0 the two references disagree and mag would oscillate ψ. Declination
    // calibration then happens at rest/low-speed.
    let magHeadingTrust = 1;
    if (this.gpsInitialized) {
      magHeadingTrust = Math.max(0, 1 - Math.abs(this.x[I.V]) / 1.5);
    }
    this._debugMagTrust = magHeadingTrust;
    // GPS owns heading → skip the mag correction AND the init snap (keeps
    // declination at its last calibrated value; it is a slow random walk).
    if (magHeadingTrust < 0.2) return;

    // Init-snap: bootstrap heading from the compass ONLY while heading is still
    // genuinely unlearned (pre-GPS-init, heading covariance at its initial value).
    // Previously the snap also fired when mag-declination covariance was still at
    // its initial value — which is true right after ANY GPS init (GPS init resets
    // S[MAG_DECL] to its initial covariance). That let the first low-speed compass
    // reading override a GPS-learned heading with the raw bearing, veering the car
    // icon off the road. Declination is calibrated via the normal cross-covariance
    // update below, not by re-snapping an already-learned ψ.
    const initCovHeading = this.config.initialCovariance.heading!;
    const psiCov = this.S[I.PSI][I.PSI] * this.S[I.PSI][I.PSI];
    const psiMag = this.wrapAngle(bearing + this.x[I.MAG_DECL]);
    if (psiCov > initCovHeading * 0.99)
      this.x[I.PSI] = psiMag;

    // GPS owns heading → skip the mag correction (keeps declination at its last
    // calibrated value; it is a slow random walk). Only re-engage once mag is
    // meaningfully trusted again.
    if (magHeadingTrust < 0.2) { this._debugMagTrust = magHeadingTrust; return; }

    const innov = this.wrapAngle(psiMag - this.x[I.PSI]);
    this._debugInnovDeg = Math.abs(innov) * 180 / Math.PI;
    const stillness = this.getStillness();
    // Heading variance from Cholesky factor — used in gate threshold so that
    // the gate naturally widens with heading uncertainty (e.g. during phone
    // rotation at low speed) and tightens when heading is well-converged.
    let psiVarMag = 0;
    for (let k = 0; k <= I.PSI; k++) psiVarMag += this.S[I.PSI][k] * this.S[I.PSI][k];
    const psiStd = Math.sqrt(psiVarMag);
    // Reject mag readings implausibly far from current heading.
    // Gate tightens when GPS is recent (urban mag interference common
    // at red lights, near traffic infrastructure) and widens when GPS
    // is stale (no heading cross-reference available e.g. phone on a table).
    // At rest (high stillness), skip the gate entirely — compass IS the
    // heading reference and large innovations are real phone rotations.
    // At low speed GPS heading is unreliable — compass IS the heading reference.
    // Skip the innovation gate and use direct blend regardless of stillness,
    // so phone rotations at v≈0 converge immediately via compass.
    const restMag = Math.abs(this.x[I.V]) < 0.5;
    this._debugGateThreshDeg = 0;
    if (this.gpsInitialized && !restMag) {
      const speedFactor = Math.max(0, 1 - Math.abs(this.x[I.V]) / 1.5);
      const secSinceGps = this.lastGpsTimeMs > 0
        ? (this.lastMagTimeMs - this.lastGpsTimeMs) / 1000 : 10;
      const gpsStale = Math.max(0, Math.min(secSinceGps / 5, 1));
      const magGateThresh = Math.max(0.5, 3 * psiStd) + speedFactor * 1.5 * gpsStale;
      this._debugGateThreshDeg = magGateThresh * 180 / Math.PI;
      if (Math.abs(innov) > magGateThresh) {
        this.magRejectCount = Math.min(this.magRejectCount + 1, 100);
        // Lockout recovery: if mag was rejected many times, heading may have
        // drifted beyond the gate. Inflate S[PSI][PSI] so the next attempt
        // has wider effective innovation covariance and can pass the gate.
        if (this.magRejectCount >= 10) {
          this.S[I.PSI][I.PSI] = Math.max(
            this.S[I.PSI][I.PSI],
            Math.min(Math.abs(innov) / 3, 1.5),
          );
        }
        if (this.magRejectCount < 30) return;
      }
    }
    this.magRejectCount = 0;
    // Adaptive measurement noise: for large innovations, inflate r so the
    // Kalman gain is reduced and psi converges smoothly instead of jumping.
    // Baseline: 10° innovation doubles r; 30° caps at 4× baseline.
    // At rest (stillness → 1), r drops to 0.5× the baseline so psi converges
    // to the compass faster — at v=0 the compass IS the only heading reference
    // and the tighter gain prevents the ~1s tracking lag during phone handling.
    const rBase = this.config.measurementNoise.heading!;
    const innovScale = Math.min(Math.abs(innov) / 0.175, 1);
    const stillnessFactor = Math.max(0.5, 1 - stillness * 0.5);
    const r = rBase * stillnessFactor * (1 + innovScale);

    // H = [0,..,1(PSI),0,..,-1(MAG_DECL),0,..]
    // tmpMagHS = P · H^T = P[:,PSI] - P[:,MAG_DECL]
    for (let i = 0; i < N; i++) {
      let sPsi = 0, sDecl = 0;
      const limPsi = Math.min(i, I.PSI);
      for (let k = 0; k <= limPsi; k++) sPsi += this.S[i][k] * this.S[I.PSI][k];
      const limDecl = Math.min(i, I.MAG_DECL);
      for (let k = 0; k <= limDecl; k++) sDecl += this.S[i][k] * this.S[I.MAG_DECL][k];
      this.tmpMagHS[i] = sPsi - sDecl;
    }
    // Innovation variance = H·P·H^T + R = P[PSI,PSI] + P[MAG_DECL,MAG_DECL] - 2*P[PSI,MAG_DECL] + r²
    let magDeclVar = 0, psiDeclCov = 0;
    for (let k = 0; k <= I.MAG_DECL; k++) magDeclVar += this.S[I.MAG_DECL][k] * this.S[I.MAG_DECL][k];
    for (let k = 0; k <= I.PSI; k++) psiDeclCov += this.S[I.PSI][k] * this.S[I.MAG_DECL][k];
    const S = psiVarMag + magDeclVar - 2 * psiDeclCov + r * r;

    if (Math.abs(this.x[I.V]) < 0.5) {
      // At rest: Kalman gain is near zero (P[PSI]≈0 from low gyro noise). Blend directly.
      // Adaptive blend: small innovations → 50% (smooth), large innovations → 90% (fast catch-up).
      // During rapid rotation at 2Hz compass rate, fixed 50% creates ~45° steady-state lag.
      // At α=0.9 the lag drops to ~5° while convergence after rotation stops is still <1.5s.
      const alpha = 0.5 + 0.4 * Math.min(Math.abs(innov), 1);
      this._debugMagAlpha = alpha;
      this.x[I.PSI] = this.wrapAngle(this.x[I.PSI] + innov * alpha);
    } else {
      this._debugMagAlpha = 0;
      for (let i = 0; i < N; i++) this.x[i] += this.tmpMagHS[i] / S * innov * magHeadingTrust;
      this.x[I.PSI] = this.wrapAngle(this.x[I.PSI]);
    }

    for (let i = 0; i < MAG_PRE; i++) this.tmpMagAT[i].fill(0);

    this.tmpMagAT[0][0] = r;
    // H·S row = S[PSI][:] - S[MAG_DECL][:]
    for (let i = 0; i < N; i++) this.tmpMagAT[1 + i][0] = this.S[I.PSI][i] - this.S[I.MAG_DECL][i];
    for (let i = 0; i < N; i++)
      for (let j = 0; j < N; j++)
        this.tmpMagAT[1 + i][1 + j] = this.S[j][i];

    this.qrInPlace(this.tmpMagAT, MAG_PRE, MAG_PRE, this.tmpHouseV);
    this.copySfromQR(this.tmpMagAT, 1);
    this.safeguardState();
  }

  // ─── Lateral acceleration pseudo-measurement ─────────────────────

  private applyLateralAccel(ay: number, omega: number): void {
    const v = this.x[I.V];
    const r = 1.0;

    // H = [0, 0, ω, 0, 0, 0, -v]  →  hs = H·S = ω·S[V] - v·S[GZ]
    const hs = this.tmpLatHS;
    const SV = this.S[I.V], SG = this.S[I.G_BIAS_Z];
    for (let j = 0; j < N; j++)
      hs[j] = omega * SV[j] - v * SG[j];

    const innov = ay - v * omega;
    let Si = r * r;
    for (let j = 0; j < N; j++) Si += hs[j] * hs[j];

    for (let i = 0; i < N; i++) {
      let ph = 0;
      for (let k = 0; k <= i; k++) ph += this.S[i][k] * hs[k];
      this.x[i] += ph / Si * innov;
    }
    this.x[I.PSI] = this.wrapAngle(this.x[I.PSI]);

    // QR covariance update (scalar QR, same pattern as mag)
    const A = this.tmpLatPre;
    for (let i = 0; i < MAG_PRE; i++) A[i].fill(0);

    A[0][0] = r;
    for (let i = 0; i < N; i++) A[1 + i][0] = hs[i];
    for (let i = 0; i < N; i++)
      for (let j = 0; j < N; j++)
        A[1 + i][1 + j] = this.S[j][i];

    this.qrInPlace(A, MAG_PRE, MAG_PRE, this.tmpHouseV);
    this.copySfromQR(A, 1);
  }

  // ─── Nonholonomic constraint (β ≈ 0, vehicle mode) ─────────────────
  // Cars cannot move sideways: lateral velocity = v·sin(β) ≈ v·β = 0.
  // Only active during straight-line driving (|ω| < 0.1) when the
  // lateral accel constraint (active during turns) is not triggered.
  private applyNonholonomic(omega: number): void {
    const r = 0.1;
    const hs = this.tmpLatHS;
    const SB = this.S[I.BETA];
    for (let j = 0; j < N; j++) hs[j] = SB[j];
    const innov = -this.x[I.BETA];
    let sInnov = r * r;
    for (let j = 0; j < N; j++) sInnov += hs[j] * hs[j];
    for (let i = 0; i < N; i++) {
      let p = 0;
      const lim = Math.min(i, I.BETA);
      for (let k = 0; k <= lim; k++) p += this.S[i][k] * hs[k];
      this.x[i] += p / sInnov * innov;
    }
    // QR covariance update (scalar QR)
    const A = this.tmpLatPre;
    for (let i = 0; i < MAG_PRE; i++) A[i].fill(0);
    A[0][0] = r;
    for (let i = 0; i < N; i++) A[1 + i][0] = hs[i];
    for (let i = 0; i < N; i++)
      for (let j = 0; j < N; j++)
        A[1 + i][1 + j] = this.S[j][i];
    this.qrInPlace(A, MAG_PRE, MAG_PRE, this.tmpHouseV);
    this.copySfromQR(A, 1);
  }

  // ─── ZUPT ────────────────────────────────────────────────────────

  private applyZupt(weight: number, omega: number = 0): void {
    const w = Math.max(weight, 0.01);
    const rVel = 0.01 / w, rPos = 1 / w;

    // compute S_innov = H·P·Hᵀ + R from Cholesky factor directly
    const S0 = this.S;
    const ly1 = S0[1][0];
    const lv0 = S0[I.V][I.V], lv1 = S0[2][1], lv2 = S0[2][0];

    const pvv = lv2*lv2 + lv1*lv1 + lv0*lv0;
    const pxv = S0[0][0] * lv2;
    const pyv = ly1 * lv2 + S0[1][1] * lv1;
    const pxx = S0[0][0] * S0[0][0];
    const pxy = S0[0][0] * ly1;
    const pyy = ly1*ly1 + S0[1][1]*S0[1][1];

    const innovV = -this.x[I.V];
    const a00 = pvv + rVel * rVel;
    const a01 = pxv, a02 = pyv;
    const a11 = pxx + rPos * rPos, a12 = pxy, a22 = pyy + rPos * rPos;

    const det = a00 * (a11 * a22 - a12 * a12)
              - a01 * (a01 * a22 - a12 * a02)
              + a02 * (a01 * a12 - a11 * a02);

    const savedPsi = this.x[I.PSI];
    if (Math.abs(det) > 1e-30) {
      const α0 = innovV * (a11 * a22 - a12 * a12) / det;
      const α1 = -innovV * (a01 * a22 - a12 * a02) / det;
      const α2 = innovV * (a01 * a12 - a11 * a02) / det;

      // P[i][V] = Σ S[i][k]·S[V][k], P[i][X] = S[i][0]·S[0][0], P[i][Y] = S[i][0]·ly1 + S[i][1]·S[1][1]
      for (let i = 0; i < N; i++) {
        if (i === I.PSI) continue;
        let pv = 0;
        for (let k = 0; k <= Math.min(i, I.V); k++) pv += S0[i][k] * S0[I.V][k];
        const px = S0[i][0] * S0[I.X][I.X];
        const py = S0[i][0] * ly1 + (i >= I.Y ? S0[i][1] * S0[1][1] : 0);
        this.x[i] += pv * α0 + px * α1 + py * α2;
      }
    }

    // ZARU: restore heading — V/X/Y innovation must not rotate heading when stationary
    this.x[I.PSI] = savedPsi;
    // Intentionally NOT clamping v to 0 — IMU acceleration can build v even
    // during ZUPT, preventing the "stuck at red light" problem. ZUPT naturally
    // disengages when |v| > 0.15 (via velFactor in predict()).

    const Mz = 3, ps = Mz + N;
    for (let i = 0; i < ps; i++) this.tmpPreA[i].fill(0);

    for (let j = 0; j < N; j++) {
      this.tmpPreA[0][Mz + j] = S0[I.V][j];
      this.tmpPreA[1][Mz + j] = S0[I.X][j];
      this.tmpPreA[2][Mz + j] = S0[I.Y][j];
    }
    this.tmpPreA[0][0] = rVel; this.tmpPreA[1][1] = rPos;
    this.tmpPreA[2][2] = rPos;

    for (let i = 0; i < N; i++)
      for (let j = 0; j < N; j++)
        this.tmpPreA[Mz + i][Mz + j] = S0[i][j];

    for (let i = 0; i < ps; i++)
      for (let j = 0; j < ps; j++)
        this.tmpPreAT[i][j] = this.tmpPreA[j][i];

    this.qrInPlace(this.tmpPreAT, ps, ps, this.tmpHouseV);
    this.copySfromQR(this.tmpPreAT, Mz);
    // General heading covariance floor — keeps heading responsive after
    // long stationary periods.
    this.S[I.PSI][I.PSI] = Math.max(this.S[I.PSI][I.PSI], 0.05);
  }

  // ─── ZARU: Zero Angular Rate Update ──────────────────────────────
  // When the device is angularly still, gz ≈ gBiasZ (true angular rate = 0).
  // Add a scalar QR pseudo-measurement to shrink gyro bias uncertainty.
  // H = [0,0,0,0,0,0,-1] → H·S = -S[G_BIAS_Z]
  // NOTE: ZARU engages whenever the IMU is angularly still, INDEPENDENT of
  // ZUPT and of the raw-IMU stillness metric. It is gated on the FILTER'S OWN
  // bias-corrected rate |omega| (the same quantity the heading integrates):
  // when the device is truly not rotating, |omega| is small whenever gBiasZ is
  // even roughly right, and the update pulls gBiasZ the rest of the way to gz.
  // A lenient 0.1 rad/s threshold lets it bootstrap from a wrong initial bias
  // (gBiasZ=0, true bias 0.0583 → |omega|≈0.0583 < 0.1 → ZARU fires). It is a
  // standalone method called directly from predict() so it runs even when GPS
  // confirms motion and ZUPT is disabled — otherwise gyro bias stays
  // uncorrected and ψ integrates the residual rate (heading drift at rest).
  private applyZaru(omega: number): void {
    const zaruGate = 1 - Math.min(Math.abs(omega) / 0.2, 1);
    if (zaruGate <= 0.5) return;
    const rGz = 0.01 / Math.max(zaruGate, 0.01);
    const hs = this.tmpLatHS;
    for (let j = 0; j < N; j++) hs[j] = -this.S[I.G_BIAS_Z][j];
    const innov = -omega;
    let sInnov = rGz * rGz;
    for (let j = 0; j < N; j++) sInnov += hs[j] * hs[j];

    // State update
    for (let i = 0; i < N; i++) {
      let p = 0;
      const lim = Math.min(i, I.G_BIAS_Z);
      for (let k = 0; k <= lim; k++) p += this.S[i][k] * hs[k];
      this.x[i] += p / sInnov * innov;
    }
    this.x[I.PSI] = this.wrapAngle(this.x[I.PSI]);

    // QR covariance update (scalar QR, same pattern as mag)
    const A = this.tmpLatPre;
    for (let i = 0; i < MAG_PRE; i++) A[i].fill(0);
    A[0][0] = rGz;
    for (let i = 0; i < N; i++) A[1 + i][0] = hs[i];
    for (let i = 0; i < N; i++)
      for (let j = 0; j < N; j++)
        A[1 + i][1 + j] = this.S[j][i];
    this.qrInPlace(A, MAG_PRE, MAG_PRE, this.tmpHouseV);
    this.copySfromQR(A, 1);
  }

}
