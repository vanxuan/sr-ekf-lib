const N = 7;
const M = 4;
const PRE = M + N;
const MAG_PRE = 1 + N;

const enum I {
  X, Y, V, PSI, BETA, A_BIAS_X, G_BIAS_Z
}

export interface EkfConfig {
  dt?: number
  processNoise?: {
    position?: number
    velocity?: number
    heading?: number
    sideslip?: number
    accelBias?: number
    gyroBias?: number
  }
  measurementNoise?: {
    position?: number
    velocity?: number
    heading?: number
  }
  magneticDeclination?: number
  initialCovariance?: {
    position?: number
    velocity?: number
    heading?: number
    sideslip?: number
    accelBias?: number
    gyroBias?: number
  }
  gateThreshold?: number
  coastTimeoutMs?: number
  gpsTimeOffsetMs?: number
  robustWeight?: {
    enabled?: boolean
    type?: 'huber' | 'cauchy'
    threshold?: number
  }
  adaptiveNoise?: {
    enabled?: boolean
    smoothing?: number
    maxScale?: number
  }
  adaptiveScaling?: {
    positionAccel?: number
    velocityAccel?: number
    velocityStep?: number
    headingGyro?: number
    headingStep?: number
    sideslipGyro?: number
    sideslipStep?: number
  }
}

export interface NavigationSolution {
  x: number; y: number;
  v: number; psi: number;
  beta: number;
  aBiasX: number;
  gBiasZ: number;
  p: Float64Array[];
}

export interface EkfDiagnostics {
  trace: number;
  gpsInnovation: number[];
  gpsChiSq: number;
  gatePassed: boolean;
  coasting: boolean;
  lastGpsTimeMs: number;
  lastImuTimeMs: number;
  stationary: boolean;
  magDeclination: number;
  robustWeight: number;
  adaNoiseScale: number;
}

const DEFAULTS = {
  dt: 0.01,
  processNoise: { position: 1.0, velocity: 0.5, heading: 0.05, sideslip: 0.1, accelBias: 1e-4, gyroBias: 1e-5 },
  measurementNoise: { position: 3.0, velocity: 0.5, heading: 0.1 },
  magneticDeclination: 0,
  initialCovariance: { position: 100, velocity: 10, heading: Math.PI * Math.PI, sideslip: 0.25, accelBias: 0.1, gyroBias: 0.01 },
  gateThreshold: 9.488,
  coastTimeoutMs: 5000,
  gpsTimeOffsetMs: 0,
  robustWeight: { enabled: false, type: 'cauchy', threshold: 9.488 },
  adaptiveNoise: { enabled: false, smoothing: 0.1, maxScale: 10 },
  adaptiveScaling: { positionAccel: 2.0, velocityAccel: 1.0, velocityStep: 2.5, headingGyro: 0.5, headingStep: 4.8, sideslipGyro: 0.5, sideslipStep: 1.8 }
};

const EPS = 1e-6;
const TWO_PI = 2 * Math.PI;
const SIDESLIP_DECAY = 0.98;

function matCreate(rows: number, cols: number): Float64Array[] {
  const m: Float64Array[] = new Array(rows);
  for (let i = 0; i < rows; i++) m[i] = new Float64Array(cols);
  return m;
}

function matIdentity(n: number): Float64Array[] {
  const m = matCreate(n, n);
  for (let i = 0; i < n; i++) m[i][i] = 1;
  return m;
}

function matCopy(A: Float64Array[]): Float64Array[] {
  const rows = A.length;
  const cols = A[0].length;
  const m = matCreate(rows, cols);
  for (let i = 0; i < rows; i++) m[i].set(A[i]);
  return m;
}

function matMul(A: Float64Array[], B: Float64Array[]): Float64Array[] {
  const rows = A.length, cols = B[0].length, inner = B.length;
  const C = matCreate(rows, cols);
  for (let i = 0; i < rows; i++) {
    const Ai = A[i], Ci = C[i];
    for (let k = 0; k < inner; k++) {
      const aik = Ai[k];
      if (aik === 0) continue;
      const Bk = B[k];
      for (let j = 0; j < cols; j++) Ci[j] += aik * Bk[j];
    }
  }
  return C;
}

function matMulABt(A: Float64Array[], B: Float64Array[]): Float64Array[] {
  const rows = A.length, cols = B.length, inner = A[0].length;
  const C = matCreate(rows, cols);
  for (let i = 0; i < rows; i++) {
    const Ai = A[i], Ci = C[i];
    for (let k = 0; k < inner; k++) {
      const aik = Ai[k];
      if (aik === 0) continue;
      for (let j = 0; j < cols; j++) Ci[j] += aik * B[j][k];
    }
  }
  return C;
}

function matTranspose(A: Float64Array[]): Float64Array[] {
  const rows = A.length, cols = A[0].length;
  const T = matCreate(cols, rows);
  for (let i = 0; i < rows; i++)
    for (let j = 0; j < cols; j++)
      T[j][i] = A[i][j];
  return T;
}

function matAdd(A: Float64Array[], B: Float64Array[]): Float64Array[] {
  const rows = A.length, cols = A[0].length;
  const C = matCreate(rows, cols);
  for (let i = 0; i < rows; i++)
    for (let j = 0; j < cols; j++)
      C[i][j] = A[i][j] + B[i][j];
  return C;
}

function matSub(A: Float64Array[], B: Float64Array[]): Float64Array[] {
  const rows = A.length, cols = A[0].length;
  const C = matCreate(rows, cols);
  for (let i = 0; i < rows; i++)
    for (let j = 0; j < cols; j++)
      C[i][j] = A[i][j] - B[i][j];
  return C;
}

function matScale(A: Float64Array[], s: number): Float64Array[] {
  const rows = A.length, cols = A[0].length;
  const C = matCreate(rows, cols);
  for (let i = 0; i < rows; i++)
    for (let j = 0; j < cols; j++)
      C[i][j] = A[i][j] * s;
  return C;
}

function traceOfP(S: Float64Array[]): number {
  let t = 0;
  for (let i = 0; i < S.length; i++) {
    const Si = S[i];
    for (let j = 0; j <= i; j++) t += Si[j] * Si[j];
  }
  return t;
}

function matTrace(A: Float64Array[]): number {
  const n = Math.min(A.length, A[0].length);
  let t = 0;
  for (let i = 0; i < n; i++) t += A[i][i];
  return t;
}

function matLowerToFull(L: Float64Array[], out?: Float64Array[]): Float64Array[] {
  const n = L.length;
  const P = out ?? matCreate(n, n);
  for (let i = 0; i < n; i++)
    for (let j = 0; j <= i; j++) {
      let s = 0;
      const lim = Math.min(i, j);
      for (let k = 0; k <= lim; k++) s += L[i][k] * L[j][k];
      P[i][j] = s; P[j][i] = s;
    }
  return P;
}

function cholDecomposition(A: Float64Array[]): Float64Array[] {
  const n = A.length;
  const L = matCreate(n, n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j <= i; j++) {
      let s = 0;
      for (let k = 0; k < j; k++) s += L[i][k] * L[j][k];
      if (i === j) {
        const v = A[i][i] - s;
        L[i][i] = v <= 0 ? 1e-8 : Math.sqrt(v);
      } else {
        L[i][j] = (A[i][j] - s) / L[j][j];
      }
    }
  }
  return L;
}

function chol4x4(A: Float64Array[]): Float64Array[] | null {
  const L = matCreate(4, 4);
  for (let i = 0; i < 4; i++) {
    for (let j = 0; j <= i; j++) {
      let s = 0;
      for (let k = 0; k < j; k++) s += L[i][k] * L[j][k];
      if (i === j) {
        const v = A[i][i] - s;
        if (v <= 0) return null;
        L[i][i] = Math.sqrt(v);
      } else {
        L[i][j] = (A[i][j] - s) / L[j][j];
      }
    }
  }
  return L;
}

function cholSolve4(L: Float64Array[], b: Float64Array): void {
  for (let i = 0; i < 4; i++) {
    let s = b[i];
    for (let k = 0; k < i; k++) s -= L[i][k] * b[k];
    b[i] = s / L[i][i];
  }
  for (let i = 3; i >= 0; i--) {
    let s = b[i];
    for (let k = i + 1; k < 4; k++) s -= L[k][i] * b[k];
    b[i] = s / L[i][i];
  }
}

function qrDecomposition(A: Float64Array[], m: number, n: number): Float64Array[] {
  const R: Float64Array[] = new Array(m);
  for (let i = 0; i < m; i++) R[i] = new Float64Array(A[i]);
  for (let k = 0; k < Math.min(m, n); k++) {
    let nrm = 0;
    for (let i = k; i < m; i++) nrm += R[i][k] * R[i][k];
    nrm = Math.sqrt(nrm);
    if (nrm < 1e-15) continue;
    const sign = R[k][k] >= 0 ? 1 : -1;
    const u0 = R[k][k] + sign * nrm;
    const v = new Float64Array(m - k);
    v[0] = u0;
    for (let i = k + 1; i < m; i++) v[i - k] = R[i][k];
    let beta = 0;
    for (let i = 0; i < m - k; i++) beta += v[i] * v[i];
    beta = 2 / beta;
    for (let j = k; j < n; j++) {
      let s = 0;
      for (let i = 0; i < m - k; i++) s += v[i] * R[k + i][j];
      s *= beta;
      for (let i = 0; i < m - k; i++) R[k + i][j] -= s * v[i];
    }
  }
  return R;
}

function matCopyBlock(src: Float64Array[], srcRow: number, srcCol: number,
                       dst: Float64Array[], dstRow: number, dstCol: number,
                       rows: number, cols: number): void {
  for (let i = 0; i < rows; i++)
    for (let j = 0; j < cols; j++)
      dst[dstRow + i][dstCol + j] = src[srcRow + i][srcCol + j];
}

export class SrEkf {
  private x: Float64Array;
  private S: Float64Array[];
  private config: Required<EkfConfig> & { adaptiveScaling: Required<EkfConfig['adaptiveScaling']> };
  private lastGpsTimeMs: number = 0;
  private lastImuTimeMs: number = 0;
  private lastInnovation: Float64Array = new Float64Array(M);
  private lastChiSq: number = 0;
  private lastGatePassed: boolean = false;
  private coasting: boolean = false;
  private gpsInitialized: boolean = false;
  private accelEnergy: number = 0;
  private gyroEnergy: number = 0;

  private tmpF: Float64Array[];
  private tmpFS: Float64Array[];
  private tmpSqrtQ: Float64Array[];
  private tmpH: Float64Array[];
  private tmpA: Float64Array[];
  private tmpInnov: Float64Array;
  private tmpWork4x4: Float64Array[];
  private tmpSinv: Float64Array[];
  private tmpQR: Float64Array[];
  private tmpHouseV: Float64Array;
  private tmpZ: Float64Array;
  private tmpPreA: Float64Array[];
  private tmpPreAT: Float64Array[];
  private tmpMagHS: Float64Array;
  private tmpMagAT: Float64Array[];

  private stepBuffer: number[] = [];
  private stepFreq: number = 0;

  private axWindow: number[] = [];
  private ayWindow: number[] = [];
  private gzWindow: number[] = [];
  private stationary: boolean = false;
  private prevStationary: boolean = false;
  private stationaryHysteresis: number = 0;
  private minVarAx: number = Infinity;
  private minVarAy: number = Infinity;
  private minVarGz: number = Infinity;
  private magDeclination: number = 0;
  private lastMagBearing: number = 0;
  private lastMagTimeMs: number = 0;
  private deviceToEnu: Float64Array[] | null = null;
  private adaNoiseScale: number = 1;
  private robustWeight: number = 1;

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
      robustWeight: { ...DEFAULTS.robustWeight, ...config?.robustWeight } as { enabled: boolean; type: 'huber' | 'cauchy'; threshold: number },
      adaptiveNoise: { ...DEFAULTS.adaptiveNoise, ...config?.adaptiveNoise },
      adaptiveScaling: { ...DEFAULTS.adaptiveScaling, ...config?.adaptiveScaling }
    };
    this.magDeclination = this.config.magneticDeclination;

    this.x = new Float64Array(N);
    this.S = matCreate(N, N);
    this.tmpF = matCreate(N, N);
    this.tmpFS = matCreate(N, N);
    this.tmpSqrtQ = matCreate(N, N);
    this.tmpH = matCreate(M, N);
    this.tmpA = matCreate(M, N);
    this.tmpInnov = new Float64Array(M);
    this.tmpWork4x4 = matCreate(M, M);
    this.tmpSinv = matCreate(M, M);
    this.tmpQR = matCreate(2 * N, N);
    this.tmpHouseV = new Float64Array(2 * N);
    this.tmpZ = new Float64Array(M);
    this.tmpPreA = matCreate(PRE, PRE);
    this.tmpPreAT = matCreate(PRE, PRE);
    this.tmpMagHS = new Float64Array(N);
    this.tmpMagAT = matCreate(MAG_PRE, MAG_PRE);

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

    for (let i = 0; i < N; i++)
      for (let j = 0; j < N; j++)
        this.S[i][j] = 0;

    this.S[I.X][I.X] = Math.sqrt(ic.position!);
    this.S[I.Y][I.Y] = Math.sqrt(ic.position!);
    this.S[I.V][I.V] = Math.sqrt(ic.velocity!);
    this.S[I.PSI][I.PSI] = Math.sqrt(ic.heading!);
    this.S[I.BETA][I.BETA] = Math.sqrt(ic.sideslip!);
    this.S[I.A_BIAS_X][I.A_BIAS_X] = Math.sqrt(ic.accelBias!);
    this.S[I.G_BIAS_Z][I.G_BIAS_Z] = Math.sqrt(ic.gyroBias!);

    this.lastGpsTimeMs = 0;
    this.lastImuTimeMs = 0;
    this.lastInnovation.fill(0);
    this.lastChiSq = 0;
    this.lastGatePassed = false;
    this.coasting = false;
    this.gpsInitialized = false;
    this.accelEnergy = 0;
    this.gyroEnergy = 0;
    this.stepBuffer = [];
    this.stepFreq = 0;
    this.axWindow = [];
    this.ayWindow = [];
    this.gzWindow = [];
    this.stationary = false;
    this.prevStationary = false;
    this.stationaryHysteresis = 0;
    this.minVarAx = Infinity;
    this.minVarAy = Infinity;
    this.minVarGz = Infinity;
    this.magDeclination = this.config.magneticDeclination ?? 0;
    this.lastMagBearing = 0;
    this.lastMagTimeMs = 0;
    this.deviceToEnu = null;
    this.adaNoiseScale = 1;
    this.robustWeight = 1;
  }

  private wrapAngle(a: number): number {
    a = a % TWO_PI;
    if (a > Math.PI) a -= TWO_PI;
    if (a <= -Math.PI) a += TWO_PI;
    return a;
  }

  setOrientation(azimuth: number, pitch: number, roll: number): void {
    const ca = Math.cos(azimuth), sa = -Math.sin(azimuth);
    const cp = Math.cos(pitch), sp = Math.sin(pitch);
    const cr = Math.cos(roll), sr = Math.sin(roll);
    const R = matCreate(3, 3);
    R[0][0] = ca * cp;
    R[0][1] = ca * sp * sr - sa * cr;
    R[0][2] = ca * sp * cr + sa * sr;
    R[1][0] = sa * cp;
    R[1][1] = sa * sp * sr + ca * cr;
    R[1][2] = sa * sp * cr - ca * sr;
    R[2][0] = -sp;
    R[2][1] = cp * sr;
    R[2][2] = cp * cr;
    this.deviceToEnu = R;
  }

  predict(ax: number, ay: number, gz: number, dt: number, timestampMs: number, az?: number, gx?: number, gy?: number): void {
    if (this.deviceToEnu) {
      const ax0 = ax, ay0 = ay, az0 = az ?? 0;
      const gx0 = gx ?? 0, gy0 = gy ?? 0, gz0 = gz;
      const R = this.deviceToEnu;
      ax = R[0][0] * ax0 + R[0][1] * ay0 + R[0][2] * az0;
      ay = R[1][0] * ax0 + R[1][1] * ay0 + R[1][2] * az0;
      gz = R[2][0] * gx0 + R[2][1] * gy0 + R[2][2] * gz0;
    }
    const a = ax - this.x[I.A_BIAS_X];
    const omega = gz - this.x[I.G_BIAS_Z];
    const psi = this.x[I.PSI];
    const beta = this.x[I.BETA];
    const v = this.x[I.V];

    this.lastImuTimeMs = timestampMs;

    this.updateImuWindow(ax, ay, gz);
    this.detectStationary();
    this.updateStepDetection(ax, dt);
    this.computeAdaptiveQ(dt);
    this.computeJacobian(a, omega, dt);

    for (let i = 0; i < N; i++) {
      for (let j = 0; j < N; j++) {
        let s = 0;
        for (let k = 0; k < N; k++) s += this.tmpF[i][k] * this.S[k][j];
        this.tmpFS[i][j] = s;
      }
    }

    for (let i = 0; i < N; i++) {
      const row = this.tmpQR[i];
      for (let j = 0; j < N; j++) row[j] = this.tmpFS[j][i];
    }
    for (let i = 0; i < N; i++) {
      const dst = this.tmpQR[N + i];
      const src = this.tmpSqrtQ[i];
      for (let j = 0; j < N; j++) dst[j] = src[j];
    }

    this.qrInPlace(this.tmpQR, 2 * N, N, this.tmpHouseV);

    for (let i = 0; i < N; i++) {
      for (let j = 0; j <= i; j++) this.S[i][j] = this.tmpQR[j][i];
      for (let j = i + 1; j < N; j++) this.S[i][j] = 0;
    }
    for (let i = 0; i < N; i++) {
      if (this.S[i][i] < 0) {
        for (let k = i; k < N; k++) this.S[k][i] = -this.S[k][i];
      }
    }

    const psiBeta = psi + beta;
    this.x[I.X] += this.ctraDeltaX(psiBeta, v, omega, dt);
    this.x[I.Y] += this.ctraDeltaY(psiBeta, v, omega, dt);
    this.x[I.V] = Math.max(0, this.x[I.V] + a * dt);
    this.x[I.PSI] = this.wrapAngle(psi + omega * dt);
    this.x[I.BETA] *= SIDESLIP_DECAY;

    if (this.stationary) {
      this.applyZupt();
    } else if (this.prevStationary) {
      for (let j = 0; j < N; j++) {
        this.S[I.V][j] = 0;
        this.S[I.X][j] = 0;
        this.S[I.Y][j] = 0;
        this.S[I.PSI][j] *= 3;
      }
      this.S[I.V][I.V] = 5.0;
      this.S[I.X][I.X] = 10.0;
      this.S[I.Y][I.Y] = 10.0;
    }
    this.prevStationary = this.stationary;
    this.safeguardState();
  }

  private updateStepDetection(ax: number, dt: number): void {
    this.stepBuffer.push(ax);
    if (this.stepBuffer.length > Math.floor(1.5 / dt)) {
      this.stepBuffer.shift();
    }
    if (this.stepBuffer.length < 20) return;
    const mean = this.stepBuffer.reduce((s, v) => s + v, 0) / this.stepBuffer.length;
    let crossings = 0;
    for (let i = 1; i < this.stepBuffer.length; i++) {
      if ((this.stepBuffer[i - 1] - mean) * (this.stepBuffer[i] - mean) < 0) crossings++;
    }
    const windowDur = (this.stepBuffer.length - 1) * dt;
    this.stepFreq = crossings > 0 ? crossings / (2 * windowDur) : 0;
  }

  private updateImuWindow(ax: number, ay: number, gz: number): void {
    this.axWindow.push(ax);
    this.ayWindow.push(ay);
    this.gzWindow.push(gz);
    const maxLen = 100;
    if (this.axWindow.length > maxLen) {
      this.axWindow.shift();
      this.ayWindow.shift();
      this.gzWindow.shift();
    }
  }

  private windowVariance(buf: number[]): number {
    if (buf.length < 10) return Infinity;
    let sum = 0, sumSq = 0;
    for (let i = 0; i < buf.length; i++) {
      sum += buf[i];
      sumSq += buf[i] * buf[i];
    }
    const n = buf.length;
    const mean = sum / n;
    return sumSq / n - mean * mean;
  }

  private detectStationary(): void {
    const varAx = this.windowVariance(this.axWindow);
    const varAy = this.windowVariance(this.ayWindow);
    const varGz = this.windowVariance(this.gzWindow);

    const thrAx = Math.max(this.minVarAx * 10, 1e-4);
    const thrAy = Math.max(this.minVarAy * 10, 1e-4);
    const thrGz = Math.max(this.minVarGz * 10, 1e-4);

    const lowImuNoise = varAx < thrAx && varAy < thrAy && varGz < thrGz;
    const lowSpeed = this.x[I.V] < 0.3;

    if (lowImuNoise && lowSpeed) {
      this.stationaryHysteresis = Math.min(this.stationaryHysteresis + 1, 20);
    } else {
      this.stationaryHysteresis = Math.max(this.stationaryHysteresis - 1, 0);
    }

    this.stationary = this.stationaryHysteresis >= 10;

    if (this.stationary && this.axWindow.length >= 10) {
      this.minVarAx = Math.min(this.minVarAx, varAx);
      this.minVarAy = Math.min(this.minVarAy, varAy);
      this.minVarGz = Math.min(this.minVarGz, varGz);
    }
  }

  private applyZupt(): void {
    const rVel = 0.01;
    const rPos = 1.0;
    const P = matLowerToFull(this.S);

    const innovV = 0 - this.x[I.V];
    const sv = P[I.V][I.V] + rVel * rVel;
    for (let i = 0; i < N; i++) {
      this.x[i] += (P[i][I.V] / sv) * innovV;
    }
    this.x[I.PSI] = this.wrapAngle(this.x[I.PSI]);
    this.x[I.V] = 0;

    const M_z = 3;
    const preSize = M_z + N;
    const A = matCreate(preSize, preSize);
    for (let i = 0; i < preSize; i++)
      for (let j = 0; j < preSize; j++)
        A[i][j] = 0;

    for (let j = 0; j < N; j++) {
      A[0][M_z + j] = this.S[I.V][j];
      A[1][M_z + j] = this.S[I.X][j];
      A[2][M_z + j] = this.S[I.Y][j];
    }
    A[0][0] = rVel;
    A[1][1] = rPos;
    A[2][2] = rPos;

    for (let i = 0; i < N; i++)
      for (let j = 0; j < N; j++)
        A[M_z + i][M_z + j] = this.S[i][j];

    const AT = matTranspose(A);
    const R_qr = qrDecomposition(AT, preSize, preSize);
    for (let i = 0; i < N; i++) {
      for (let j = 0; j <= i; j++) this.S[i][j] = R_qr[M_z + j][M_z + i];
      for (let j = i + 1; j < N; j++) this.S[i][j] = 0;
    }
    for (let i = 0; i < N; i++) {
      if (this.S[i][i] < 0) {
        for (let k = i; k < N; k++) this.S[k][i] = -this.S[k][i];
      }
    }
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
      if (!isFinite(this.x[i])) {
        this.x[i] = 0;
        corrupted = true;
      }
      if (!isFinite(this.S[i][i]) || this.S[i][i] <= 0) {
        this.S[i][i] = 1e-6;
        corrupted = true;
      }
      for (let j = 0; j < i; j++) {
        if (!isFinite(this.S[i][j])) {
          this.S[i][j] = 0;
          corrupted = true;
        }
      }
      for (let j = i + 1; j < N; j++) {
        if (this.S[i][j] !== 0) {
          this.S[i][j] = 0;
          corrupted = true;
        }
      }
    }
    if (corrupted && this.gpsInitialized) {
      this.coasting = true;
    }
  }

  private resetFromGps(gpsX: number, gpsY: number, gpsVx: number, gpsVy: number): void {
    const gpsV = Math.sqrt(gpsVx * gpsVx + gpsVy * gpsVy);
    const gpsPsi = Math.atan2(gpsVy, gpsVx);
    this.x[I.X] = gpsX;
    this.x[I.Y] = gpsY;
    this.x[I.V] = gpsV;
    this.x[I.PSI] = this.wrapAngle(gpsPsi);
    this.x[I.BETA] = 0;
    const posSigma = Math.max(this.S[I.X][I.X], 7);
    const velSigma = this.config.measurementNoise.velocity!;
    const psiSigma = gpsV > 1 ? Math.min(velSigma / gpsV, 0.5) : 0.5;
    for (let i = 0; i < N; i++) {
      const lim = Math.min(i, I.BETA);
      for (let j = 0; j <= lim; j++) this.S[i][j] = 0;
    }
    this.S[I.X][I.X] = posSigma;
    this.S[I.Y][I.Y] = posSigma;
    this.S[I.V][I.V] = velSigma;
    this.S[I.PSI][I.PSI] = psiSigma;
    this.S[I.BETA][I.BETA] = 0.1;
    this.coasting = false;
  }

  private magUpdateSingle(bearing: number, timestampMs: number): number {
    const psiCov = this.S[I.PSI][I.PSI] * this.S[I.PSI][I.PSI];
    if (psiCov > this.config.initialCovariance.heading! * 0.99) {
      this.x[I.PSI] = this.wrapAngle(bearing + this.magDeclination);
    }

    const psiMag = this.wrapAngle(bearing + this.magDeclination);
    const innov = this.wrapAngle(psiMag - this.x[I.PSI]);
    const r = this.config.measurementNoise.heading!;

    let psiVar = 0;
    for (let k = 0; k <= I.PSI; k++) psiVar += this.S[I.PSI][k] * this.S[I.PSI][k];
    for (let i = 0; i < N; i++) {
      let s = 0;
      const lim = Math.min(i, I.PSI);
      for (let k = 0; k <= lim; k++) s += this.S[i][k] * this.S[I.PSI][k];
      this.tmpMagHS[i] = s;
    }

    const S = psiVar + r * r;
    for (let i = 0; i < N; i++) {
      this.x[i] += (this.tmpMagHS[i] / S) * innov;
    }
    this.x[I.PSI] = this.wrapAngle(this.x[I.PSI]);

    for (let i = 0; i < MAG_PRE; i++)
      for (let j = 0; j < MAG_PRE; j++)
        this.tmpMagAT[i][j] = 0;

    this.tmpMagAT[0][0] = r;
    for (let i = 0; i < N; i++) this.tmpMagAT[1 + i][0] = this.S[I.PSI][i];
    for (let i = 0; i < N; i++)
      for (let j = 0; j < N; j++)
        this.tmpMagAT[1 + i][1 + j] = this.S[j][i];

    this.qrInPlace(this.tmpMagAT, MAG_PRE, MAG_PRE, this.tmpHouseV);

    for (let i = 0; i < N; i++) {
      for (let j = 0; j <= i; j++) this.S[i][j] = this.tmpMagAT[1 + j][1 + i];
      for (let j = i + 1; j < N; j++) this.S[i][j] = 0;
    }
    for (let i = 0; i < N; i++) {
      if (this.S[i][i] < 0) {
        for (let k = i; k < N; k++) this.S[k][i] = -this.S[k][i];
      }
    }

    const normInnov = innov * innov / (psiVar + r * r);
    return normInnov;
  }

  updateMag(bearing: number, timestampMs: number): void {
    this.lastImuTimeMs = timestampMs;
    this.lastMagBearing = bearing;
    this.lastMagTimeMs = timestampMs;

    this.magUpdateSingle(bearing, timestampMs);
  }

  private computeAdaptiveQ(dt: number): void {
    const sqrtDt = Math.sqrt(dt);
    const v = this.x[I.V];
    const speedScale = Math.sqrt(Math.max(v, 0.5) / 5.0);
    for (let i = 0; i < N; i++)
      for (let j = 0; j < N; j++)
        this.tmpSqrtQ[i][j] = 0;

    const varAx = this.windowVariance(this.axWindow);
    const varAy = this.windowVariance(this.ayWindow);
    const varGz = this.windowVariance(this.gzWindow);

    const rawAccel = Math.sqrt(Math.max(varAx + varAy, 0)) / 5.0;
    this.accelEnergy = 0.9 * this.accelEnergy + 0.1 * Math.min(rawAccel, 5);

    const rawGyro = Math.sqrt(Math.max(varGz, 0)) / 0.5;
    this.gyroEnergy = 0.9 * this.gyroEnergy + 0.1 * Math.min(rawGyro, 5);

    const stepEnergy = Math.min(this.stepFreq / 3.0, 1.0);

    const sc = this.config.adaptiveScaling;
    const posScale = 1 + sc.positionAccel * this.accelEnergy;
    const velScale = 1 + sc.velocityAccel * this.accelEnergy + sc.velocityStep * stepEnergy;
    const headScale = 1 + sc.headingGyro * this.gyroEnergy + sc.headingStep * stepEnergy;
    const slipScale = 1 + sc.sideslipGyro * this.gyroEnergy + sc.sideslipStep * stepEnergy;

    const pn = this.config.processNoise;
    this.tmpSqrtQ[I.X][I.X] = pn.position! * sqrtDt * speedScale * posScale;
    this.tmpSqrtQ[I.Y][I.Y] = pn.position! * sqrtDt * speedScale * posScale;
    this.tmpSqrtQ[I.V][I.V] = pn.velocity! * sqrtDt * velScale;
    this.tmpSqrtQ[I.PSI][I.PSI] = pn.heading! * sqrtDt * headScale;
    this.tmpSqrtQ[I.BETA][I.BETA] = pn.sideslip! * sqrtDt * slipScale;
    this.tmpSqrtQ[I.A_BIAS_X][I.A_BIAS_X] = pn.accelBias! * sqrtDt;
    this.tmpSqrtQ[I.G_BIAS_Z][I.G_BIAS_Z] = pn.gyroBias! * sqrtDt;
  }

  private ctraDeltaX(psi: number, v: number, omega: number, dt: number): number {
    if (Math.abs(omega) > EPS) return (v / omega) * (Math.sin(psi + omega * dt) - Math.sin(psi));
    return v * Math.cos(psi) * dt;
  }

  private ctraDeltaY(psi: number, v: number, omega: number, dt: number): number {
    if (Math.abs(omega) > EPS) return (v / omega) * (-Math.cos(psi + omega * dt) + Math.cos(psi));
    return v * Math.sin(psi) * dt;
  }

  private computeJacobian(a: number, omega: number, dt: number): void {
    const psi = this.x[I.PSI];
    const beta = this.x[I.BETA];
    const psiBeta = psi + beta;
    const v = this.x[I.V];
    for (let i = 0; i < N; i++)
      for (let j = 0; j < N; j++)
        this.tmpF[i][j] = i === j ? 1 : 0;

    if (Math.abs(omega) > EPS) {
      const sp = Math.sin(psiBeta), cp = Math.cos(psiBeta);
      const spw = Math.sin(psiBeta + omega * dt), cpw = Math.cos(psiBeta + omega * dt);
      const o2 = omega * omega;
      this.tmpF[I.X][I.V] = (spw - sp) / omega;
      this.tmpF[I.X][I.PSI] = (v / omega) * (cpw - cp);
      this.tmpF[I.X][I.BETA] = this.tmpF[I.X][I.PSI];
      this.tmpF[I.X][I.G_BIAS_Z] = v * ((spw - sp) / o2 - dt * cpw / omega);
      this.tmpF[I.Y][I.V] = (-cpw + cp) / omega;
      this.tmpF[I.Y][I.PSI] = (v / omega) * (spw - sp);
      this.tmpF[I.Y][I.BETA] = this.tmpF[I.Y][I.PSI];
      this.tmpF[I.Y][I.G_BIAS_Z] = v * ((-cpw + cp) / o2 - dt * spw / omega);
    } else {
      const cp = Math.cos(psiBeta), sp = Math.sin(psiBeta);
      const dt2 = dt * dt;
      this.tmpF[I.X][I.V] = cp * dt;
      this.tmpF[I.X][I.PSI] = -v * sp * dt;
      this.tmpF[I.X][I.BETA] = this.tmpF[I.X][I.PSI];
      this.tmpF[I.X][I.G_BIAS_Z] = 0.5 * v * sp * dt2;
      this.tmpF[I.Y][I.V] = sp * dt;
      this.tmpF[I.Y][I.PSI] = v * cp * dt;
      this.tmpF[I.Y][I.BETA] = this.tmpF[I.Y][I.PSI];
      this.tmpF[I.Y][I.G_BIAS_Z] = -0.5 * v * cp * dt2;
    }
    this.tmpF[I.V][I.A_BIAS_X] = -dt;
    this.tmpF[I.PSI][I.G_BIAS_Z] = -dt;
    this.tmpF[I.BETA][I.BETA] = SIDESLIP_DECAY;
  }

  private computeH(): void {
    const psi = this.x[I.PSI], beta = this.x[I.BETA];
    const psiBeta = psi + beta, v = this.x[I.V];
    const cp = Math.cos(psiBeta), sp = Math.sin(psiBeta);
    for (let i = 0; i < M; i++)
      for (let j = 0; j < N; j++)
        this.tmpH[i][j] = 0;
    this.tmpH[0][I.X] = 1;
    this.tmpH[1][I.Y] = 1;
    this.tmpH[2][I.V] = cp;
    this.tmpH[2][I.PSI] = -v * sp;
    this.tmpH[2][I.BETA] = -v * sp;
    this.tmpH[3][I.V] = sp;
    this.tmpH[3][I.PSI] = v * cp;
    this.tmpH[3][I.BETA] = v * cp;
  }

  private computeGpsInnovation(z: Float64Array): void {
    const x = this.x[I.X], y = this.x[I.Y], v = this.x[I.V], psi = this.x[I.PSI], beta = this.x[I.BETA];
    const psiBeta = psi + beta;
    this.tmpInnov[0] = z[0] - x;
    this.tmpInnov[1] = z[1] - y;
    this.tmpInnov[2] = z[2] - v * Math.cos(psiBeta);
    this.tmpInnov[3] = z[3] - v * Math.sin(psiBeta);

    const H = this.tmpH;
    const S = this.S;
    const A = this.tmpA;

    for (let i = 0; i < M; i++) {
      const Ai = A[i];
      for (let j = 0; j < N; j++) {
        let s = 0;
        for (let k = j; k < N; k++) s += H[i][k] * S[k][j];
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
      }
    }
  }

  private computeGpsPostFit(posR: number, velR: number): number {
    this.tmpWork4x4[0][0] += posR * posR;
    this.tmpWork4x4[1][1] += posR * posR;
    this.tmpWork4x4[2][2] += velR * velR;
    this.tmpWork4x4[3][3] += velR * velR;

    const L = chol4x4(this.tmpWork4x4);
    if (!L) return Infinity;
    for (let i = 0; i < M; i++) this.tmpSinv[i].set(L[i]);

    const y = new Float64Array(M);
    for (let i = 0; i < M; i++) y[i] = this.tmpInnov[i];
    cholSolve4(L, y);
    let chiSq = 0;
    for (let i = 0; i < M; i++) chiSq += this.tmpInnov[i] * y[i];
    return chiSq;
  }

  private gpsUpdateSingle(posR: number, velR: number, x: number, y: number): boolean {
    this.computeH();
    this.computeGpsInnovation(this.tmpZ);

    let chiSq = this.computeGpsPostFit(posR, velR);
    this.lastChiSq = chiSq;

    const vSpeed = this.x[I.V];
    if (chiSq > this.config.gateThreshold && vSpeed > 0.5) {
      const psiBeta = this.x[I.PSI] + this.x[I.BETA];
      const vxPred = vSpeed * Math.cos(psiBeta);
      const vyPred = vSpeed * Math.sin(psiBeta);
      const vxMeas = this.tmpZ[2];
      const vyMeas = this.tmpZ[3];
      const speedMeasSq = vxMeas * vxMeas + vyMeas * vyMeas;
      if (speedMeasSq > 0.25) {
        const dot = vxPred * vxMeas + vyPred * vyMeas;
        if (dot < -0.5 * vSpeed * vSpeed) {
          this.x[I.PSI] = this.wrapAngle(this.x[I.PSI] + Math.PI);
          for (let j = 0; j < I.PSI; j++) this.S[I.PSI][j] = 0;
          for (let j = I.PSI + 1; j < N; j++) this.S[j][I.PSI] = 0;
          this.S[I.PSI][I.PSI] = Math.max(this.S[I.PSI][I.PSI], 0.5);
          this.computeH();
          this.computeGpsInnovation(this.tmpZ);
          chiSq = this.computeGpsPostFit(posR, velR);
          this.lastChiSq = chiSq;
        }
      }
    }

    let robustW = 1;
    const rw = this.config.robustWeight;
    if (rw.enabled) {
      const thr = rw.threshold!;
      if (rw.type === 'huber') {
        robustW = chiSq > thr ? thr / chiSq : 1;
      } else {
        robustW = 1 / (1 + (chiSq / thr) * (chiSq / thr));
      }
    }

    this.robustWeight = robustW;

    let adaScale = 1;
    if (this.config.adaptiveNoise.enabled) {
      const innovRatio = chiSq / M;
      const α = this.config.adaptiveNoise.smoothing!;
      this.adaNoiseScale = α * Math.max(1, innovRatio) + (1 - α) * this.adaNoiseScale;
      adaScale = Math.min(this.adaNoiseScale, this.config.adaptiveNoise.maxScale!);
    }

    const totalWeight = robustW / adaScale;
    if (totalWeight < 1) {
      this.tmpWork4x4[0][0] -= posR * posR;
      this.tmpWork4x4[1][1] -= posR * posR;
      this.tmpWork4x4[2][2] -= velR * velR;
      this.tmpWork4x4[3][3] -= velR * velR;
      posR /= Math.sqrt(totalWeight);
      velR /= Math.sqrt(totalWeight);
      chiSq = this.computeGpsPostFit(posR, velR);
    }

    if (!this.config.robustWeight.enabled && chiSq > this.config.gateThreshold) {
      if (this.coasting) {
        this.resetFromGps(this.tmpZ[0], this.tmpZ[1], this.tmpZ[2], this.tmpZ[3]);
        return true;
      }
      return false;
    }

    const A = this.tmpA;
    const H = this.tmpH;
    const L = this.tmpSinv;
    const w = new Float64Array(M);

    for (let i = 0; i < M; i++) w[i] = this.tmpInnov[i];
    cholSolve4(L, w);

    const z = new Float64Array(N);
    for (let j = 0; j < N; j++) {
      let s = 0;
      for (let i = 0; i < M; i++) s += A[i][j] * w[i];
      z[j] = s;
    }

    for (let i = 0; i < N; i++) {
      let s = 0;
      for (let j = 0; j <= i; j++) s += this.S[i][j] * z[j];
      this.x[i] += s;
    }
    this.x[I.PSI] = this.wrapAngle(this.x[I.PSI]);
    this.x[I.V] = Math.max(0, this.x[I.V]);

    for (let i = 0; i < PRE; i++)
      for (let j = 0; j < PRE; j++)
        this.tmpPreA[i][j] = 0;

    this.tmpPreA[0][0] = posR;
    this.tmpPreA[1][1] = posR;
    this.tmpPreA[2][2] = velR;
    this.tmpPreA[3][3] = velR;

    for (let i = 0; i < M; i++)
      for (let j = 0; j < N; j++) {
        let s = 0;
        for (let k = 0; k < N; k++) s += H[i][k] * this.S[k][j];
        this.tmpPreA[i][M + j] = s;
      }

    for (let i = 0; i < N; i++)
      for (let j = 0; j < N; j++)
        this.tmpPreA[M + i][M + j] = this.S[i][j];

    for (let i = 0; i < PRE; i++)
      for (let j = 0; j < PRE; j++)
        this.tmpPreAT[i][j] = this.tmpPreA[j][i];

    this.qrInPlace(this.tmpPreAT, PRE, PRE, this.tmpHouseV);

    for (let i = 0; i < N; i++) {
      for (let j = 0; j <= i; j++) this.S[i][j] = this.tmpPreAT[M + j][M + i];
      for (let j = i + 1; j < N; j++) this.S[i][j] = 0;
    }
    for (let i = 0; i < N; i++) {
      if (this.S[i][i] < 0) {
        for (let k = i; k < N; k++) this.S[k][i] = -this.S[k][i];
      }
    }

    this.safeguardState();
    return true;
  }

  updateGps(x: number, y: number, vx: number, vy: number, timestampMs: number, accuracyMeters?: number): boolean {
    if (!this.gpsInitialized) {
      this.x[I.X] = x;
      this.x[I.Y] = y;
      const v = Math.sqrt(vx * vx + vy * vy);
      if (v > 0.01) {
        this.x[I.V] = v;
        this.x[I.PSI] = this.wrapAngle(Math.atan2(vy, vx));
      }
      this.gpsInitialized = true;
      this.lastGpsTimeMs = timestampMs + this.config.gpsTimeOffsetMs!;
      this.lastGatePassed = true;
      return true;
    }

    const origPosR = this.config.measurementNoise.position!;
    const origVelR = this.config.measurementNoise.velocity!;
    let posR = origPosR;
    let velR = origVelR;
    if (accuracyMeters !== undefined) {
      const scale = Math.max(accuracyMeters / origPosR, 0.1);
      posR = origPosR * scale;
      velR = origVelR * scale;
    }

    const gpsV = Math.sqrt(vx * vx + vy * vy);
    if (gpsV > 3.0 && this.lastMagTimeMs > 0 && (accuracyMeters === undefined || accuracyMeters < 15) && timestampMs - this.lastMagTimeMs < 5000) {
      const gpsPsi = Math.atan2(vy, vx);
      const observedDec = this.wrapAngle(gpsPsi - this.lastMagBearing);
      this.magDeclination += 0.05 * this.wrapAngle(observedDec - this.magDeclination);
    }

    this.tmpZ[0] = x; this.tmpZ[1] = y; this.tmpZ[2] = vx; this.tmpZ[3] = vy;

    const result = this.gpsUpdateSingle(posR, velR, x, y);
    if (result) {
      this.lastGatePassed = true;
      this.lastGpsTimeMs = timestampMs + this.config.gpsTimeOffsetMs!;
    } else {
      this.lastGatePassed = false;
      if (this.coasting) {
        this.resetFromGps(x, y, vx, vy);
        this.lastGpsTimeMs = timestampMs + this.config.gpsTimeOffsetMs!;
        this.lastGatePassed = true;
        return true;
      }
    }
    return result;
  }

  coast(timeoutMs: number, currentTimeMs: number): boolean {
    if (currentTimeMs - this.lastGpsTimeMs > timeoutMs) {
      this.coasting = true;
      return true;
    }
    this.coasting = false;
    return true;
  }

  getState(): NavigationSolution {
    return {
      x: this.x[I.X], y: this.x[I.Y],
      v: this.x[I.V], psi: this.x[I.PSI],
      beta: this.x[I.BETA],
      aBiasX: this.x[I.A_BIAS_X],
      gBiasZ: this.x[I.G_BIAS_Z],
      p: matLowerToFull(this.S)
    };
  }

  getDiagnostics(): EkfDiagnostics {
    return {
      trace: traceOfP(this.S),
      gpsInnovation: [this.tmpInnov[0], this.tmpInnov[1], this.tmpInnov[2], this.tmpInnov[3]],
      gpsChiSq: this.lastChiSq,
      gatePassed: this.lastGatePassed,
      coasting: this.coasting,
      lastGpsTimeMs: this.lastGpsTimeMs,
      lastImuTimeMs: this.lastImuTimeMs,
      stationary: this.stationary,
      magDeclination: this.magDeclination,
      robustWeight: this.robustWeight,
      adaNoiseScale: this.adaNoiseScale
    };
  }
}
