import { RingBuf } from './ring-buf';
import { EkfDiagnostics } from './config';

export function wmean(b: RingBuf): number {
  if (b.length === 0) return 0;
  let s = 0;
  for (let i = 0; i < b.length; i++) s += b.get(i);
  return s / b.length;
}

export function wstd(b: RingBuf, m: number): number {
  if (b.length < 2) return 0;
  let s = 0;
  for (let i = 0; i < b.length; i++) { const d = b.get(i) - m; s += d * d; }
  return Math.sqrt(s / b.length);
}

export function buildDiagnostics(
  traceCache: number, gpsInnovation: Float64Array,
  lastChiSq: number, gatePassed: boolean, coasting: boolean,
  lastGpsTimeMs: number, lastImuTimeMs: number,
  v: number, motionStillness: number, magDeclination: number,
  robustWeight: number, adaNoiseScale: number
): EkfDiagnostics {
  return {
    trace: traceCache,
    gpsInnovation,
    gpsChiSq: lastChiSq, gatePassed,
    coasting, lastGpsTimeMs,
    lastImuTimeMs,
    stationary: motionStillness > 0.7 && Math.abs(v) < 3.0,
    motionStillness,
    magDeclination,
    robustWeight,
    adaNoiseScale
  };
}

export function buildDebug(
  stillness: number, aBiasX: number, gBiasZ: number, v: number,
  psiDeg: number, zuptWeight: number, speedGate: number, accelGate: number,
  n: number, magRejectCount: number, magTrust: number,
  magInnovDeg: number, gateThreshDeg: number, magAlpha: number
): {
  stillness: number; aBiasX: number; gBiasZ: number; v: number; psiDeg: number;
  zuptWeight: number; speedGate: number; accelGate: number; n: number;
  magRejectCount: number; magTrust: number; magInnovDeg: number;
  gateThreshDeg: number; magAlpha: number;
} {
  return {
    stillness, aBiasX, gBiasZ, v, psiDeg,
    zuptWeight, speedGate, accelGate, n,
    magRejectCount, magTrust, magInnovDeg, gateThreshDeg, magAlpha
  };
}

export function buildImuStats(
  axWindow: RingBuf, gzWindow: RingBuf, aBiasX: number, gBiasZ: number, lastOmega: number
): {
  n: number; meanAxRel: number; stdAx: number;
  meanGzRel: number; stdGz: number; lastOmega: number;
} {
  const mAx = wmean(axWindow);
  const mGz = wmean(gzWindow);
  return {
    n: axWindow.length,
    meanAxRel: mAx - aBiasX,
    stdAx: wstd(axWindow, mAx),
    meanGzRel: mGz - gBiasZ,
    stdGz: wstd(gzWindow, mGz),
    lastOmega
  };
}
