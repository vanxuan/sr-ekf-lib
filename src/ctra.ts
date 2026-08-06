import { EPS } from './config';

export function ctraDelta(psi: number, v: number, omega: number, dt: number): [number, number] {
  if (Math.abs(omega) > EPS) {
    const sp = Math.sin(psi), cp = Math.cos(psi);
    const spw = Math.sin(psi + omega * dt), cpw = Math.cos(psi + omega * dt);
    return [v / omega * (spw - sp), v / omega * (-cpw + cp)];
  }
  const cp = Math.cos(psi), sp = Math.sin(psi);
  return [v * cp * dt, v * sp * dt];
}

export function computeJacobian(
  psiBeta: number, v: number, a: number, omega: number, dt: number,
  betaTau: number, aBiasDecay: number, EPS: number,
  F: Float64Array[]
): void {
  const N = 8;
  for (let i = 0; i < N; i++) F[i].fill(0);
  for (let i = 0; i < N; i++) F[i][i] = 1;

  const vAvg = v + 0.5 * a * dt;
  const absOmega = Math.abs(omega);
  if (absOmega > EPS) {
    const sp = Math.sin(psiBeta), cp = Math.cos(psiBeta);
    const spw = Math.sin(psiBeta + omega * dt), cpw = Math.cos(psiBeta + omega * dt);
    const o2 = omega * omega;
    F[0][2] = (spw - sp) / omega;
    F[0][3] = vAvg / omega * (cpw - cp);
    F[0][4] = F[0][3];
    F[0][6] = -vAvg * ((spw - sp) / o2 - dt * cpw / omega);
    F[1][2] = (-cpw + cp) / omega;
    F[1][3] = vAvg / omega * (spw - sp);
    F[1][4] = F[1][3];
    F[1][6] = -vAvg * ((-cpw + cp) / o2 - dt * spw / omega);
  } else {
    const cp = Math.cos(psiBeta), sp = Math.sin(psiBeta);
    const dt2 = dt * dt;
    F[0][2] = cp * dt;
    F[0][3] = -vAvg * sp * dt;
    F[0][4] = F[0][3];
    F[0][6] = -0.5 * vAvg * sp * dt2;
    F[1][2] = sp * dt;
    F[1][3] = vAvg * cp * dt;
    F[1][4] = F[1][3];
    F[1][6] = 0.5 * vAvg * cp * dt2;
  }
  F[0][5] = -0.5 * dt * F[0][2];
  F[1][5] = -0.5 * dt * F[1][2];
  F[2][5] = -dt;
  F[3][6] = -dt;
  F[4][4] = Math.exp(-dt / betaTau);
  F[5][5] = aBiasDecay;
  F[7][7] = 1;
}
