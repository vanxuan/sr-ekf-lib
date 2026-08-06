export const N = 8, M = 4, PRE = M + N, MAG_PRE = 1 + N;
export const NTRI = N * (N + 1) / 2;

export const I = { X: 0, Y: 1, V: 2, PSI: 3, BETA: 4, A_BIAS_X: 5, G_BIAS_Z: 6, MAG_DECL: 7 } as const;

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
  lastImuTimeMs: number; stationary: boolean; motionStillness: number; magDeclination: number
  robustWeight: number; adaNoiseScale: number
}

export const DEFAULTS = {
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

export const EPS = 1e-4;

export const MOTION_V_CUT = 1.0;
export const GPS_REST_NOISE = 1.0;
export const MOTION_V_CUT_STALE = 9.0;
export const DEVICE_VAR_SCALE = 8.0;
export const DEVICE_GYRO_SCALE = 0.5;
export const DEVICE_ACTIVITY_FLOOR = 0.01;
export const COAST_DAMP_STILL = 0.5;
