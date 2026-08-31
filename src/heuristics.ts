/**
 * heuristics.ts — Named thresholds and pure gate functions for SrEkf.
 *
 * Extracts magic numbers scattered across sr-ekf.ts into a single,
 * documented table so each heuristic can be unit-tested in isolation.
 * All values are the CURRENT tuned defaults; change here, not inline.
 */

// ─── EPS — split from ambiguous `EPS` for distinct uses ───────────
/** CTRA singularity threshold (|ω| ≤ EPS uses small-ω branch). */
export const EPS_CTRA = 1e-4;
/** Yaw-rate threshold for β time-constant base selection. */
export const EPS_TAU_YAW = 1e-4;

// ─── Process-noise / covariance limits ────────────────────────────
export const COAST_Q_FACTOR = 0.3;
export const MAX_POS_SIGMA = 500; // m
export const MAX_VEL_SIGMA = 50;  // m/s
export const COAST_POS_SIGMA_CAP = 5; // min(current,5) in resetFromGps

// ─── Bias time constants ──────────────────────────────────────────
export const ABIAS_REVERSION_TAU_S = 50;
export const ABIAS_RAW_ENERGY_GATE = 0.05;

// ─── CTRA integration ─────────────────────────────────────────────
export const VEL_AVG_FACTOR = 0.5; // vAvg = v + 0.5*a*dt
export const V_CLAMP_MIN = 0 as const;

// ─── Heading integration ──────────────────────────────────────────
export const VEH_MOVING_V_SCALE = 0.5; // vehMoving = min(|v|/0.5,1)

// ─── Angular acceleration ─────────────────────────────────────────
export const ANG_ACCEL_CLAMP = 10; // rad/s²
export const ANG_ACCEL_SUSTAIN_TAU_S = 0.5;
export const ANG_ACCEL_BOOST_SCALE = 2.0; // angAccelBoost = min(smooth/2,1)
export const ANG_ACCEL_NORM_SCALE = 3.0; // angAccelNorm = min(smooth/3,1)

// ─── β (sideslip) ─────────────────────────────────────────────────
export const BETA_TAU_LOW_V_THRESH = 0.3; // |v|<0.3 → τ=0.1
export const BETA_TAU_LOW_V = 0.1;
export const BETA_TAU_TURN = 1.5;
export const BETA_TAU_MIN = 0.5;
// betaTau = base * (1 - BETA_TAU_ANG_SCALE*angAccelNorm) * (1 - BETA_TAU_STEP_SCALE*stepEnergy)
export const BETA_TAU_ANG_SCALE = 0.6;
export const BETA_TAU_STEP_SCALE = 0.5;

// ─── ZUPT ─────────────────────────────────────────────────────────
export const ZUPT_ON = 0.15;
export const ZUPT_OFF = 0.05;
export const ZUPT_R_VEL = 0.01; // m/s
export const ZUPT_R_POS = 1.0;
export const ZUPT_CHI_SQ_GATE = 9.0; // 3σ (1-DOF)  — innovation²/(pvv+R²) ≤ 9
export const ZUPT_SPEED_GATE_V0 = 8.0; // speedGate = (8-|v|)/4
export const ZUPT_SPEED_GATE_SCALE = 4.0;
export const ZUPT_HEADING_FLOOR = 0.05; // rad — post-ZUPT S[PSI][PSI] clamp
export const ZUPT_V_INFLATE = 1.5;
export const ZUPT_XY_INFLATE = 3.0;
export const ZUPT_PSI_INFLATE = 0.3;

// ─── Coasting (tunnel) ────────────────────────────────────────────
export const COAST_DAMP_TIME_S = 3; // seconds before damping engages
export const COAST_DAMP_MIN_V = 0.05; // |v| must exceed
export const COAST_DAMP_R_BASE = 10.0; // rDamp = max(0.05, 10/coastTime)
export const COAST_DAMP_R_MIN = 0.05;
export const COAST_PRIOR_R_BASE = 2.0;
export const COAST_PRIOR_R_SLOPE = 0.6; // R = base + slope*coastTime
export const COAST_PRIOR_R_CAP = 20.0;
export const COAST_STALE_DECAY_TAU_S = 2.0;
export const COAST_SPEED_CAPTURE_MIN = 0.5; // |v| <0.5 → coastSpeed=0 (was stationary)

// ─── ZARU ─────────────────────────────────────────────────────────
export const ZARU_GATE_OMEGA_SCALE = 0.2; // zaruGate = 1 - min(|ω|/0.2,1)
export const ZARU_GATE_THRESH = 0.5; // gate >0.5 → |ω|<0.1
export const ZARU_R_BASE = 0.01;

// ─── Lateral accel / nonholonomic ─────────────────────────────────
export const LAT_OMEGA_THRESH = 0.1;
export const LAT_V_THRESH = 0.2;
export const LAT_R = 1.0; // m/s²
export const NONHOL_OMEGA_THRESH = 0.1;
export const NONHOL_V_THRESH = 0.15;
export const NONHOL_R = 0.1; // rad

// ─── Baro ─────────────────────────────────────────────────────────
export const BARO_DT_MIN = 0.05;
export const BARO_DT_MAX = 5.0;
export const BARO_SIN_PITCH_THRESH = 0.03; // ~1.7°
export const BARO_R_MIN = 0.5;
export const BARO_R_VZ_SCALE = 0.1;
export const BARO_CHI_SQ_GATE = 11.3; // 3σ, 1-DOF ≈ 0.99 quantile
export const BARO_NON_COAST_R_SCALE = 2.0; // inflate R when GPS fresh

// ─── GPS ──────────────────────────────────────────────────────────
export const GPS_INIT_SPD_THRESH = 0.1; // m/s
export const GPS_INIT_DIR_STD_MAX = 0.2; // rad
export const GPS_INIT_COV_RELAX = 0.99;
export const GPS_INIT_POS_STD_FLOOR = 1.0;
export const GPS_INIT_POS_STD_DEFAULT = 3.0;
export const GPS_INIT_VEL_STD_CAP = 0.5;
export const GPS_INIT_PSI_STD_CAP = 0.5;
export const GPS_INIT_PSI_STD_TURN = 0.15;
export const GPS_SPEED_RAMP_V0 = 10.0; // ramp = 1 - |v|/10
export const GPS_SPEED_RAMP_SCALE = 2.0;
export const GPS_SPEED_RAMP_CREEP_BAND = 2.5;
export const GPS_VEL_DOT_LOW_V = 2.5;
export const GPS_VEL_OPPOSE_SCALE = 4.0;
export const GPS_MAX_PLAUSIBLE_K = 2.0;
export const GPS_MAX_PLAUSIBLE_FLOOR = 1.0;
export const GPS_MAX_PLAUSIBLE_OFFSET = 2.0;
export const GPS_DT_BASE_FLOOR = 0.1;
export const GPS_GUARD_POSR_FWD_FLOOR_K = 2.0;
export const GPS_GUARD_POSR_FWD_FLOOR_MIN = 1.0;
export const GPS_GUARD_POSR_CROSS_FLOOR_K = 0.5;
export const GPS_GUARD_POSR_CROSS_FLOOR_MIN = 0.5;
export const GPS_GUARD_FWD_K = 2.0;
export const GPS_GUARD_CROSS_K = 0.5;
export const GPS_GUARD_FWD_CAP = 5.0;
export const GPS_GUARD_CROSS_CAP = 10.0;
export const GPS_GUARD_INFLATE_RATIO = 3.0;
export const GPS_PLAUSIBLE_JUMP_VMAX = 50; // m/s
export const GPS_PLAUSIBLE_JUMP_OFFSET = 2; // m
export const GPS_NEAR_STOP_V_THRESH = 1.0;
export const GPS_NEAR_STOP_DX = 15; // m
export const GPS_DIVERGENCE_DIST2 = 10000; // 100m squared
export const GPS_REWIND_GAP2 = 100; // 10m squared
export const GPS_LARGE_POS_RESET2 = 1e8; // 10km squared
export const GPS_STATIONARY_WEIGHT_CAP_SCALE = 4.0; // velR *=1+min(4,9*w)
export const GPS_STATIONARY_WEIGHT_SLOPE = 9.0;
export const GPS_REWIND_PSI_INFLATE = 1.5;

// ─── Rest-exit snap ───────────────────────────────────────────────
export const REST_EXIT_V_EKF = 0.5;
export const REST_EXIT_SPEED_GPS = 2.5;
export const REST_EXIT_OMEGA = 0.1;
export const REST_EXIT_GAP2 = 400; // 20m squared
export const REST_EXIT_PSI_DELTA = 0.26; // ~15°

// ─── Heading gain (computeH) ──────────────────────────────────────
export const HDG_GAIN_V_FLOOR = 0.3;
export const HDG_GAIN_DENOM = 3.6;
export const HDG_GAIN_GYRO_K = 0.5;
export const HDG_GAIN_STEP_K = 1.5;
export const HDG_GAIN_ANG_K = 0.5;
export const HDG_GAIN_INIT_WINDOW_S = 30000;
export const HDG_GAIN_INIT_PSI_STD_FLOOR = 0.3;
export const HDG_GAIN_INIT_SCALE = 0.7;
export const HDG_GAIN_INIT_BOOST = 2.0;

// ─── Rest heading block ───────────────────────────────────────────
export const REST_HDG_V = 1.5;
export const REST_HDG_BLOCK_V0 = 1.5;
export const REST_HDG_BLOCK_SCALE = 1.0;

// ─── Robust / adaptive ────────────────────────────────────────────
export const ROBUST_POS_SCALE = 0.63; // 2-DOF vs 4-DOF χ² 0.95 ratio
export const VEL_ROBUST_SCALE = ROBUST_POS_SCALE;
export const ADA_CONVERGE_THRESH = 0.5;
export const ADA_CONVERGE_COUNT = 10;

// ─── Chi-square positional sub-gate ───────────────────────────────
// 4-DOF gateThreshold (9.488) is χ² 0.95,4. Position sub-gate is 2-DOF →
// threshold = 5.991 (χ² 0.95,2) or gateThreshold * 0.63, whichever the user
// configured is tighter. Using the 4-DOF value here is overly permissive.
export const CHI2_95_2DOF = 5.991;
export function posGateThreshold(gateThreshold: number): number {
  return Math.min(CHI2_95_2DOF, gateThreshold * ROBUST_POS_SCALE);
}

// ─── Mag ──────────────────────────────────────────────────────────
export const MAG_HEADING_TRUST_V0 = 2.5; // trust = max(0,1 - v/2.5)
export const MAG_HEADING_TRUST_CUTOFF = 0.2;
export const MAG_GYRO_GUARD_THRESH = 0.05; // |ω|>0.05 → skip mag
export const MAG_REST_V = 0.5;
export const MAG_GPS_STATIONARY_V = 0.3;
export const MAG_DRIFT_SMOOTH_V0 = 1.5; // start of ramp
export const MAG_DRIFT_SMOOTH_V1 = 2.5; // end of ramp
export const MAG_DRIFT_MARGIN_LOW = 2.0;
export const MAG_DRIFT_MARGIN_HIGH = 0.05;
export const MAG_GATE_BASE_K = 3.0;
export const MAG_GATE_FLOOR = 0.5;
export const MAG_GATE_SPEED_K = 1.5;
export const MAG_GATE_STALE_TAU_S = 5.0;
export const MAG_LOCKOUT_WARN = 10;
export const MAG_LOCKOUT_HARD = 30;
export const MAG_R_INNOV_SCALE_DEN = 0.175; // ~10°
export const MAG_R_STILLNESS_FLOOR = 0.5;
export const MAG_REST_ALPHA_BASE = 0.5;
export const MAG_REST_ALPHA_SCALE = 0.5;
export const MAG_REST_ALPHA_CLAMP = 1.0;
export const MAG_REST_PRIOR_DECAY_S = 30;
export const MAG_REST_PRIOR_DELTA_SCALE = 0.3;

// ─── Pure functions (table-driven, unit-testable) ─────────────────

export function magDriftMargin(vAbs: number): number {
  if (vAbs < MAG_DRIFT_SMOOTH_V0) return MAG_DRIFT_MARGIN_LOW;
  if (vAbs < MAG_DRIFT_SMOOTH_V1) {
    const t = (vAbs - MAG_DRIFT_SMOOTH_V0) / (MAG_DRIFT_SMOOTH_V1 - MAG_DRIFT_SMOOTH_V0);
    return MAG_DRIFT_MARGIN_LOW + t * (MAG_DRIFT_MARGIN_HIGH - MAG_DRIFT_MARGIN_LOW);
  }
  return MAG_DRIFT_MARGIN_HIGH;
}

export function magHeadingTrust(magSpeed: number): number {
  return Math.max(0, 1 - magSpeed / MAG_HEADING_TRUST_V0);
}

export function shouldSkipMagForTrust(trust: number): boolean {
  return trust < MAG_HEADING_TRUST_CUTOFF;
}

export function headingGainRaw(
  gainSpeed: number,
  gyroEnergy: number,
  stepEnergy: number,
  smoothAngAccel: number,
): number {
  if (gainSpeed < HDG_GAIN_V_FLOOR) return 0;
  return Math.min(
    (gainSpeed - HDG_GAIN_V_FLOOR +
      gyroEnergy * HDG_GAIN_GYRO_K +
      stepEnergy * HDG_GAIN_STEP_K +
      smoothAngAccel * HDG_GAIN_ANG_K) /
      HDG_GAIN_DENOM,
    1,
  );
}

export function betaTau(
  absV: number,
  absOmega: number,
  smoothAngAccel: number,
  stepEnergy: number,
  eps: number = EPS_CTRA,
): number {
  const angAccelNorm = Math.min(smoothAngAccel / ANG_ACCEL_NORM_SCALE, 1);
  const base =
    absV < BETA_TAU_LOW_V_THRESH
      ? BETA_TAU_LOW_V
      : absOmega > eps
        ? BETA_TAU_TURN
        : Math.max(BETA_TAU_MIN, BETA_TAU_TURN - (absV - 1.5) * (1.0 / 3.5));
  return base * (1 - BETA_TAU_ANG_SCALE * angAccelNorm) * (1 - BETA_TAU_STEP_SCALE * stepEnergy);
}

export function coastDampingR(coastTimeS: number): number {
  return Math.max(COAST_DAMP_R_MIN, COAST_DAMP_R_BASE / Math.max(coastTimeS, 1));
}
export function coastPriorR(coastTimeS: number): number {
  return Math.min(COAST_PRIOR_R_BASE + coastTimeS * COAST_PRIOR_R_SLOPE, COAST_PRIOR_R_CAP);
}

export function speedRampFactor(absV: number): number {
  return Math.max(0, 1 - absV / GPS_SPEED_RAMP_V0);
}
export function isCoastDampingActive(motionStillness: number, coastStillThresh: number): boolean {
  return motionStillness >= coastStillThresh; // inclusive to close 0.5 deadband
}
export function isVelocityPriorActive(motionStillness: number, coastStillThresh: number): boolean {
  return motionStillness < coastStillThresh;
}

export function traceOfP(S: Float64Array[]): number {
  let tr = 0;
  for (let i = 0; i < S.length; i++) for (let j = 0; j <= i; j++) tr += S[i][j] * S[i][j];
  return tr;
}
