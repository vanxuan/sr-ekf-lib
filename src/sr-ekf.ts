const N = 8;
const M = 4;

const enum I {
  X, Y, V, PSI, BETA, A_BIAS_X, A_BIAS_Y, G_BIAS_Z
}

export type Mode = 'walk' | 'drive' | 'auto';

export interface EkfConfig {
  dt?: number
  mode?: Mode
  processNoise?: {
    position?: number
    velocity?: number
    heading?: number
    accelBias?: number
    gyroBias?: number
  }
  walkingProcessNoise?: {
    position?: number
    velocity?: number
    heading?: number
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
}

export interface NavigationSolution {
  x: number; y: number;
  v: number; psi: number;
  beta: number;
  aBiasX: number; aBiasY: number;
  gBiasZ: number;
  p: Float64Array[];
  mode: Mode;
}

export interface EkfDiagnostics {
  trace: number;
  gpsInnovation: number[];
  gpsChiSq: number;
  gatePassed: boolean;
  coasting: boolean;
  lastGpsTimeMs: number;
  lastImuTimeMs: number;
  mode: Mode;
  walkLikelihood: number;
  stationary: boolean;
  magDeclination: number;
}

const DEFAULTS = {
  dt: 0.01,
  mode: 'auto' as Mode,
  processNoise: { position: 1.0, velocity: 0.5, heading: 0.05, accelBias: 1e-4, gyroBias: 1e-5 },
  walkingProcessNoise: { position: 2.0, velocity: 2.0, heading: 0.3, accelBias: 1e-3, gyroBias: 1e-4 },
  measurementNoise: { position: 3.0, velocity: 0.5, heading: 0.1 },
  magneticDeclination: 0,
  initialCovariance: { position: 100, velocity: 10, heading: Math.PI * Math.PI, sideslip: 0.25, accelBias: 0.1, gyroBias: 0.01 },
  gateThreshold: 9.488,
  coastTimeoutMs: 5000,
  gpsTimeOffsetMs: 0
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

function matTrace(A: Float64Array[]): number {
  const n = Math.min(A.length, A[0].length);
  let t = 0;
  for (let i = 0; i < n; i++) t += A[i][i];
  return t;
}

function matLowerToFull(L: Float64Array[]): Float64Array[] {
  const n = L.length;
  const P = matCreate(n, n);
  for (let i = 0; i < n; i++)
    for (let j = 0; j <= i; j++) {
      let s = 0;
      for (let k = 0; k <= i; k++) s += L[i][k] * L[j][k];
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

function matInvert4x4(A: Float64Array[]): Float64Array[] {
  const inv = matCreate(4, 4);
  for (let i = 0; i < 4; i++) inv[i][i] = 1;
  const m = matCreate(4, 4);
  for (let i = 0; i < 4; i++) m[i].set(A[i]);
  for (let col = 0; col < 4; col++) {
    let maxRow = col;
    let maxVal = Math.abs(m[col][col]);
    for (let row = col + 1; row < 4; row++) {
      const v = Math.abs(m[row][col]);
      if (v > maxVal) { maxVal = v; maxRow = row; }
    }
    if (maxVal < 1e-15) continue;
    if (maxRow !== col) {
      [m[col], m[maxRow]] = [m[maxRow], m[col]];
      [inv[col], inv[maxRow]] = [inv[maxRow], inv[col]];
    }
    const piv = m[col][col];
    for (let j = 0; j < 4; j++) { m[col][j] /= piv; inv[col][j] /= piv; }
    for (let row = 0; row < 4; row++) {
      if (row === col) continue;
      const factor = m[row][col];
      for (let j = 0; j < 4; j++) { m[row][j] -= factor * m[col][j]; inv[row][j] -= factor * inv[col][j]; }
    }
  }
  return inv;
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
  private config: Required<EkfConfig> & { walkingProcessNoise: Required<EkfConfig['processNoise']> };
  private lastGpsTimeMs: number = 0;
  private lastImuTimeMs: number = 0;
  private lastInnovation: Float64Array = new Float64Array(M);
  private lastChiSq: number = 0;
  private lastGatePassed: boolean = false;
  private coasting: boolean = false;
  private gpsInitialized: boolean = false;
  private effectiveMode: 'walk' | 'drive' = 'drive';
  private walkLikelihood: number = 0.5;

  private tmpF: Float64Array[];
  private tmpFS: Float64Array[];
  private tmpSqrtQ: Float64Array[];
  private tmpH: Float64Array[];
  private tmpP: Float64Array[];
  private tmpInnov: Float64Array;
  private tmpWork4x4: Float64Array[];

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
  private positionFrozen: boolean = false;

  constructor(config?: EkfConfig) {
    const wPN = { ...DEFAULTS.walkingProcessNoise, ...config?.walkingProcessNoise };
    this.config = {
      dt: config?.dt ?? DEFAULTS.dt,
      mode: config?.mode ?? DEFAULTS.mode,
      processNoise: { ...DEFAULTS.processNoise, ...config?.processNoise },
      walkingProcessNoise: wPN,
      measurementNoise: { ...DEFAULTS.measurementNoise, ...config?.measurementNoise },
      initialCovariance: { ...DEFAULTS.initialCovariance, ...config?.initialCovariance },
      gateThreshold: config?.gateThreshold ?? DEFAULTS.gateThreshold,
      coastTimeoutMs: config?.coastTimeoutMs ?? DEFAULTS.coastTimeoutMs,
      gpsTimeOffsetMs: config?.gpsTimeOffsetMs ?? DEFAULTS.gpsTimeOffsetMs,
      magneticDeclination: config?.magneticDeclination ?? DEFAULTS.magneticDeclination
    };
    this.magDeclination = this.config.magneticDeclination;

    this.effectiveMode = this.config.mode === 'auto' ? 'drive' : this.config.mode;

    this.x = new Float64Array(N);
    this.S = matCreate(N, N);
    this.tmpF = matCreate(N, N);
    this.tmpFS = matCreate(N, N);
    this.tmpSqrtQ = matCreate(N, N);
    this.tmpH = matCreate(M, N);
    this.tmpP = matCreate(N, N);
    this.tmpInnov = new Float64Array(M);
    this.tmpWork4x4 = matCreate(M, M);

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
    this.x[I.A_BIAS_Y] = 0;
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
    this.S[I.A_BIAS_Y][I.A_BIAS_Y] = Math.sqrt(ic.accelBias!);
    this.S[I.G_BIAS_Z][I.G_BIAS_Z] = Math.sqrt(ic.gyroBias!);

    this.lastGpsTimeMs = 0;
    this.lastImuTimeMs = 0;
    this.lastInnovation.fill(0);
    this.lastChiSq = 0;
    this.lastGatePassed = false;
    this.coasting = false;
    this.gpsInitialized = false;
    this.effectiveMode = this.config.mode === 'auto' ? 'drive' : this.config.mode;
    this.walkLikelihood = 0.5;
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
  }

  private wrapAngle(a: number): number {
    a = a % TWO_PI;
    if (a > Math.PI) a -= TWO_PI;
    if (a <= -Math.PI) a += TWO_PI;
    return a;
  }

  setOrientation(azimuth: number, pitch: number, roll: number): void {
    const ca = Math.cos(azimuth), sa = Math.sin(azimuth);
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

    if (this.config.mode === 'auto') {
      this.computeWalkLikelihood(v, omega, gz, a, dt);
      this.effectiveMode = this.walkLikelihood > 0.5 ? 'walk' : 'drive';
    } else {
      this.effectiveMode = this.config.mode;
    }

    this.computeQ(dt);

    this.computeJacobian(a, omega, dt);

    for (let i = 0; i < N; i++) {
      for (let j = 0; j < N; j++) {
        let s = 0;
        for (let k = 0; k < N; k++) s += this.tmpF[i][k] * this.S[k][j];
        this.tmpFS[i][j] = s;
      }
    }

    const mRows = 2 * N;
    const M_qr: Float64Array[] = new Array(mRows);
    for (let i = 0; i < N; i++) {
      M_qr[i] = new Float64Array(N);
      for (let j = 0; j < N; j++) M_qr[i][j] = this.tmpFS[j][i];
    }
    for (let i = 0; i < N; i++) M_qr[N + i] = new Float64Array(this.tmpSqrtQ[i]);

    const R = qrDecomposition(M_qr, mRows, N);

    for (let i = 0; i < N; i++) {
      for (let j = 0; j <= i; j++) this.S[i][j] = R[j][i];
      for (let j = i + 1; j < N; j++) this.S[i][j] = 0;
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
        this.S[I.V][j] *= 5;
        this.S[I.PSI][j] *= 5;
        this.S[I.X][j] *= 5;
        this.S[I.Y][j] *= 5;
      }
    }
    this.prevStationary = this.stationary;
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

  private computeWalkLikelihood(v: number, omega: number, gz: number, a: number, dt: number): void {
    const walkSpeedScore = this.sigmoid(2.0 - v, 0.5);
    const stepFreqScore = this.stepFreq >= 1.0 && this.stepFreq <= 3.0 ? 0.9 : 0.1;
    const accelNoiseScore = this.sigmoid(Math.abs(a) * 5, 0.2);
    const gyroScore = this.sigmoid(Math.abs(gz) * 5 - 0.5, 0.3);
    const speedHard = v > 5.0 ? 0.0 : 1.0;

    const raw = (walkSpeedScore * 0.35 + stepFreqScore * 0.35 + accelNoiseScore * 0.15 + gyroScore * 0.15) * speedHard;
    this.walkLikelihood = 0.9 * this.walkLikelihood + 0.1 * raw;
  }

  private sigmoid(x: number, k: number): number {
    return 1 / (1 + Math.exp(-x / k));
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

    if (this.axWindow.length >= 10 && this.x[I.V] < 0.3) {
      this.minVarAx = Math.min(this.minVarAx, varAx);
      this.minVarAy = Math.min(this.minVarAy, varAy);
      this.minVarGz = Math.min(this.minVarGz, varGz);
    }

    const thrAx = Math.max(this.minVarAx * 3, 1e-8);
    const thrAy = Math.max(this.minVarAy * 3, 1e-8);
    const thrGz = Math.max(this.minVarGz * 3, 1e-8);

    const lowImuNoise = varAx < thrAx && varAy < thrAy && varGz < thrGz;
    const lowSpeed = this.x[I.V] < 0.3;

    if (lowImuNoise && lowSpeed) {
      this.stationaryHysteresis = Math.min(this.stationaryHysteresis + 1, 20);
    } else {
      this.stationaryHysteresis = Math.max(this.stationaryHysteresis - 1, 0);
    }

    this.stationary = this.stationaryHysteresis >= 10;
  }

  private applyZupt(): void {
    const rVel = 0.01;
    const rPos = 0.01;
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
  }

  updateMag(bearing: number, timestampMs: number): void {
    this.lastImuTimeMs = timestampMs;
    this.lastMagBearing = bearing;
    this.lastMagTimeMs = timestampMs;

    const psiCov = this.S[I.PSI][I.PSI] * this.S[I.PSI][I.PSI];
    if (psiCov > this.config.initialCovariance.heading! * 0.99) {
      this.x[I.PSI] = this.wrapAngle(bearing + this.magDeclination);
    }

    const psiMag = this.wrapAngle(bearing + this.magDeclination);
    const innov = this.wrapAngle(psiMag - this.x[I.PSI]);

    const r = this.config.measurementNoise.heading!;
    const P = matLowerToFull(this.S);

    let S = P[I.PSI][I.PSI] + r * r;
    for (let i = 0; i < N; i++) {
      this.x[i] += (P[i][I.PSI] / S) * innov;
    }
    this.x[I.PSI] = this.wrapAngle(this.x[I.PSI]);

    const H = new Float64Array(N);
    H[I.PSI] = 1;

    const HS = new Float64Array(N);
    for (let j = 0; j < N; j++) {
      let s = 0;
      for (let k = 0; k < N; k++) s += H[k] * this.S[k][j];
      HS[j] = s;
    }

    const preSize = 1 + N;
    const A = matCreate(preSize, preSize);
    A[0][0] = r;
    for (let j = 0; j < N; j++) A[0][1 + j] = HS[j];
    for (let i = 0; i < N; i++)
      for (let j = 0; j < N; j++)
        A[1 + i][1 + j] = this.S[i][j];

    const AT = matTranspose(A);
    const R_qr = qrDecomposition(AT, preSize, preSize);
    for (let i = 0; i < N; i++) {
      for (let j = 0; j <= i; j++) this.S[i][j] = R_qr[1 + j][1 + i];
      for (let j = i + 1; j < N; j++) this.S[i][j] = 0;
    }
  }

  private computeQ(dt: number): void {
    const pn = this.effectiveMode === 'walk' ? this.config.walkingProcessNoise : this.config.processNoise;
    const sqrtDt = Math.sqrt(dt);
    const v = this.x[I.V];
    const speedScale = Math.sqrt(Math.max(v, 0.5) / 5.0);
    for (let i = 0; i < N; i++)
      for (let j = 0; j < N; j++)
        this.tmpSqrtQ[i][j] = 0;

    if (this.config.mode === 'auto' && this.effectiveMode === 'drive' && this.walkLikelihood > 0.15) {
      const pnDrive = this.config.processNoise;
      const pnWalk = this.config.walkingProcessNoise;
      const w = this.walkLikelihood;
      const blend = (d: number, wd: number) => Math.sqrt((1 - w) * d * d + w * wd * wd);
      this.tmpSqrtQ[I.X][I.X] = blend(pnDrive.position!, pnWalk.position!) * sqrtDt * speedScale;
      this.tmpSqrtQ[I.Y][I.Y] = blend(pnDrive.position!, pnWalk.position!) * sqrtDt * speedScale;
      this.tmpSqrtQ[I.V][I.V] = blend(pnDrive.velocity!, pnWalk.velocity!) * sqrtDt * speedScale;
      this.tmpSqrtQ[I.PSI][I.PSI] = blend(pnDrive.heading!, pnWalk.heading!) * sqrtDt;
      this.tmpSqrtQ[I.BETA][I.BETA] = blend(0.1, 0.3) * sqrtDt;
      this.tmpSqrtQ[I.A_BIAS_X][I.A_BIAS_X] = blend(pnDrive.accelBias!, pnWalk.accelBias!) * sqrtDt;
      this.tmpSqrtQ[I.A_BIAS_Y][I.A_BIAS_Y] = blend(pnDrive.accelBias!, pnWalk.accelBias!) * sqrtDt;
      this.tmpSqrtQ[I.G_BIAS_Z][I.G_BIAS_Z] = blend(pnDrive.gyroBias!, pnWalk.gyroBias!) * sqrtDt;
    } else {
      this.tmpSqrtQ[I.X][I.X] = pn.position! * sqrtDt * speedScale;
      this.tmpSqrtQ[I.Y][I.Y] = pn.position! * sqrtDt * speedScale;
      this.tmpSqrtQ[I.V][I.V] = pn.velocity! * sqrtDt * speedScale;
      this.tmpSqrtQ[I.PSI][I.PSI] = pn.heading! * sqrtDt;
      this.tmpSqrtQ[I.BETA][I.BETA] = 0.1 * sqrtDt;
      this.tmpSqrtQ[I.A_BIAS_X][I.A_BIAS_X] = pn.accelBias! * sqrtDt;
      this.tmpSqrtQ[I.A_BIAS_Y][I.A_BIAS_Y] = pn.accelBias! * sqrtDt;
      this.tmpSqrtQ[I.G_BIAS_Z][I.G_BIAS_Z] = pn.gyroBias! * sqrtDt;
    }
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
      this.tmpF[I.X][I.V] = cp * dt;
      this.tmpF[I.X][I.PSI] = -v * sp * dt;
      this.tmpF[I.X][I.BETA] = this.tmpF[I.X][I.PSI];
      this.tmpF[I.Y][I.V] = sp * dt;
      this.tmpF[I.Y][I.PSI] = v * cp * dt;
      this.tmpF[I.Y][I.BETA] = this.tmpF[I.Y][I.PSI];
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

  private computeInnovationAndChiSq(z: Float64Array): void {
    const x = this.x[I.X], y = this.x[I.Y], v = this.x[I.V], psi = this.x[I.PSI], beta = this.x[I.BETA];
    const psiBeta = psi + beta;
    this.tmpInnov[0] = z[0] - x;
    this.tmpInnov[1] = z[1] - y;
    this.tmpInnov[2] = z[2] - v * Math.cos(psiBeta);
    this.tmpInnov[3] = z[3] - v * Math.sin(psiBeta);

    const P = matLowerToFull(this.S);
    this.tmpP = P;
    const H = this.tmpH;

    for (let i = 0; i < M; i++)
      for (let j = 0; j < M; j++) {
        let s = 0;
        for (let k = 0; k < N; k++)
          for (let l = 0; l < N; l++)
            s += H[i][k] * P[k][l] * H[j][l];
        this.tmpWork4x4[i][j] = s;
      }

    const mn = this.config.measurementNoise;
    const mp = mn.position!;
    const mv = mn.velocity!;
    this.tmpWork4x4[0][0] += mp * mp;
    this.tmpWork4x4[1][1] += mp * mp;
    this.tmpWork4x4[2][2] += mv * mv;
    this.tmpWork4x4[3][3] += mv * mv;

    const S_inv = matInvert4x4(this.tmpWork4x4);

    let chiSq = 0;
    for (let i = 0; i < M; i++)
      for (let j = 0; j < M; j++)
        chiSq += this.tmpInnov[i] * S_inv[i][j] * this.tmpInnov[j];
    this.lastChiSq = chiSq;
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

    const origPosR = this.config.measurementNoise.position;
    const origVelR = this.config.measurementNoise.velocity;
    if (accuracyMeters !== undefined) {
      const scale = Math.max(accuracyMeters / origPosR!, 0.1);
      this.config.measurementNoise.position = origPosR! * scale;
      this.config.measurementNoise.velocity = origVelR! * scale;
    }

    if (this.config.mode === 'auto' && this.stepFreq >= 0.8 && this.stepFreq <= 3.5) {
      const v = Math.sqrt(vx * vx + vy * vy);
      if (v < 3.0) {
        this.walkLikelihood = Math.min(0.95, this.walkLikelihood + 0.1);
      }
    }

    const gpsV = Math.sqrt(vx * vx + vy * vy);
    if (gpsV > 3.0 && this.lastMagTimeMs > 0 && (accuracyMeters === undefined || accuracyMeters < 15) && timestampMs - this.lastMagTimeMs < 5000) {
      const gpsPsi = Math.atan2(vy, vx);
      const observedDec = this.wrapAngle(gpsPsi - this.lastMagBearing);
      this.magDeclination += 0.05 * this.wrapAngle(observedDec - this.magDeclination);
    }

    const z = new Float64Array([x, y, vx, vy]);
    this.computeH();
    this.computeInnovationAndChiSq(z);

    if (this.lastChiSq > this.config.gateThreshold) {
      this.lastGatePassed = false;
      if (accuracyMeters !== undefined) {
        this.config.measurementNoise.position = origPosR;
        this.config.measurementNoise.velocity = origVelR;
      }
      if (this.coasting) {
        this.x[I.X] = x;
        this.x[I.Y] = y;
        this.S[I.X][I.X] = Math.max(this.S[I.X][I.X], 7);
        this.S[I.Y][I.Y] = Math.max(this.S[I.Y][I.Y], 7);
        this.lastGpsTimeMs = timestampMs + this.config.gpsTimeOffsetMs!;
        this.positionFrozen = false;
        this.coasting = false;
        this.lastGatePassed = true;
        return true;
      }
      return false;
    }
    this.lastGatePassed = true;
    this.positionFrozen = false;
    this.lastGpsTimeMs = timestampMs + this.config.gpsTimeOffsetMs!;

    const P = this.tmpP;
    const H = this.tmpH;

    for (let i = 0; i < M; i++)
      for (let j = 0; j < M; j++) {
        let s = 0;
        for (let k = 0; k < N; k++)
          for (let l = 0; l < N; l++)
            s += H[i][k] * P[k][l] * H[j][l];
        this.tmpWork4x4[i][j] = s;
      }

    const mn = this.config.measurementNoise;
    const mp2 = mn.position!;
    const mv2 = mn.velocity!;
    this.tmpWork4x4[0][0] += mp2 * mp2;
    this.tmpWork4x4[1][1] += mp2 * mp2;
    this.tmpWork4x4[2][2] += mv2 * mv2;
    this.tmpWork4x4[3][3] += mv2 * mv2;

    const S_inv = matInvert4x4(this.tmpWork4x4);

    const K = matCreate(N, M);
    for (let i = 0; i < N; i++)
      for (let j = 0; j < M; j++) {
        let sum = 0;
        for (let k = 0; k < N; k++)
          for (let l = 0; l < M; l++)
            sum += P[i][k] * H[l][k] * S_inv[l][j];
        K[i][j] = sum;
      }

    for (let i = 0; i < N; i++)
      for (let j = 0; j < M; j++)
        this.x[i] += K[i][j] * this.tmpInnov[j];
    this.x[I.PSI] = this.wrapAngle(this.x[I.PSI]);
    this.x[I.V] = Math.max(0, this.x[I.V]);

    const mnStd = this.config.measurementNoise;
    const R_chol = matCreate(M, M);
    R_chol[0][0] = mnStd.position!;
    R_chol[1][1] = mnStd.position!;
    R_chol[2][2] = mnStd.velocity!;
    R_chol[3][3] = mnStd.velocity!;

    const HS = matCreate(M, N);
    for (let i = 0; i < M; i++)
      for (let j = 0; j < N; j++) {
        let s = 0;
        for (let k = 0; k < N; k++) s += H[i][k] * this.S[k][j];
        HS[i][j] = s;
      }

    const preSize = M + N;
    const A = matCreate(preSize, preSize);
    for (let i = 0; i < preSize; i++)
      for (let j = 0; j < preSize; j++)
        A[i][j] = 0;
    matCopyBlock(R_chol, 0, 0, A, 0, 0, M, M);
    matCopyBlock(HS, 0, 0, A, 0, M, M, N);
    matCopyBlock(this.S, 0, 0, A, M, M, N, N);

    const AT = matTranspose(A);
    const R_qr = qrDecomposition(AT, preSize, preSize);

    for (let i = 0; i < N; i++) {
      for (let j = 0; j <= i; j++) this.S[i][j] = R_qr[M + j][M + i];
      for (let j = i + 1; j < N; j++) this.S[i][j] = 0;
    }

    const gpsHeadingSpeed = Math.sqrt(vx * vx + vy * vy);
    if (gpsHeadingSpeed > 0.5) {
      const gpsHeading = Math.atan2(vy, vx);
      const rH = 0.2 / Math.max(gpsHeadingSpeed, 0.5);
      const psiBeta = this.x[I.PSI] + this.x[I.BETA];
      const innovH = this.wrapAngle(gpsHeading - psiBeta);

      const P_H = matLowerToFull(this.S);
      const SH = P_H[I.PSI][I.PSI] + 2 * P_H[I.PSI][I.BETA] + P_H[I.BETA][I.BETA] + rH * rH;
      for (let i = 0; i < N; i++) {
        this.x[i] += ((P_H[i][I.PSI] + P_H[i][I.BETA]) / SH) * innovH;
      }
      this.x[I.PSI] = this.wrapAngle(this.x[I.PSI]);

      const HH = new Float64Array(N);
      HH[I.PSI] = 1;
      HH[I.BETA] = 1;
      const HSH = new Float64Array(N);
      for (let j = 0; j < N; j++) {
        let s = 0;
        for (let k = 0; k < N; k++) s += HH[k] * this.S[k][j];
        HSH[j] = s;
      }
      const preH = 1 + N;
      const AH = matCreate(preH, preH);
      AH[0][0] = rH;
      for (let j = 0; j < N; j++) AH[0][1 + j] = HSH[j];
      for (let i = 0; i < N; i++)
        for (let j = 0; j < N; j++)
          AH[1 + i][1 + j] = this.S[i][j];
      const AHT = matTranspose(AH);
      const R_qrH = qrDecomposition(AHT, preH, preH);
      for (let i = 0; i < N; i++) {
        for (let j = 0; j <= i; j++) this.S[i][j] = R_qrH[1 + j][1 + i];
        for (let j = i + 1; j < N; j++) this.S[i][j] = 0;
      }
    }

    if (accuracyMeters !== undefined) {
      this.config.measurementNoise.position = origPosR;
      this.config.measurementNoise.velocity = origVelR;
    }

    return true;
  }

  coast(timeoutMs: number, currentTimeMs: number): boolean {
    if (currentTimeMs - this.lastGpsTimeMs > timeoutMs) {
      this.coasting = true;
      const P = matLowerToFull(this.S);
      const trace = matTrace(P);
      if (trace > 10000) return false;
      return true;
    }
    this.coasting = false;
    this.positionFrozen = false;
    return true;
  }

  getState(): NavigationSolution {
    return {
      x: this.x[I.X], y: this.x[I.Y],
      v: this.x[I.V], psi: this.x[I.PSI],
      beta: this.x[I.BETA],
      aBiasX: this.x[I.A_BIAS_X], aBiasY: this.x[I.A_BIAS_Y],
      gBiasZ: this.x[I.G_BIAS_Z],
      p: matLowerToFull(this.S),
      mode: this.effectiveMode
    };
  }

  getDiagnostics(): EkfDiagnostics {
    const P = matLowerToFull(this.S);
    return {
      trace: matTrace(P),
      gpsInnovation: [this.tmpInnov[0], this.tmpInnov[1], this.tmpInnov[2], this.tmpInnov[3]],
      gpsChiSq: this.lastChiSq,
      gatePassed: this.lastGatePassed,
      coasting: this.coasting,
      lastGpsTimeMs: this.lastGpsTimeMs,
      lastImuTimeMs: this.lastImuTimeMs,
      mode: this.effectiveMode,
      walkLikelihood: this.walkLikelihood,
      stationary: this.stationary,
      magDeclination: this.magDeclination
    };
  }
}
