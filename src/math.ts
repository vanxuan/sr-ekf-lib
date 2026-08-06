export const TWO_PI = 2 * Math.PI;

export function matCreate(rows: number, cols: number): Float64Array[] {
  const m: Float64Array[] = new Array(rows);
  for (let i = 0; i < rows; i++) m[i] = new Float64Array(cols);
  return m;
}

export function matLowerToFullInto(L: Float64Array[], out: Float64Array[]): void {
  const n = L.length;
  for (let i = 0; i < n; i++) {
    const oi = out[i];
    for (let j = 0; j <= i; j++) {
      let s = 0;
      for (let k = 0; k <= j; k++) s += L[i][k] * L[j][k];
      oi[j] = s;
      out[j][i] = s;
    }
  }
}

export function matLowerToFull(L: Float64Array[]): Float64Array[] {
  const n = L.length;
  const P = new Array<Float64Array>(n);
  for (let i = 0; i < n; i++) {
    P[i] = new Float64Array(n);
    for (let j = 0; j <= i; j++) {
      let s = 0;
      for (let k = 0; k <= j; k++) s += L[i][k] * L[j][k];
      P[i][j] = s; P[j][i] = s;
    }
  }
  return P;
}

export function chol4x4(out: Float64Array[], A: Float64Array[]): boolean {
  for (let i = 0; i < 4; i++) {
    for (let j = 0; j <= i; j++) {
      let s = 0;
      for (let k = 0; k < j; k++) s += out[i][k] * out[j][k];
      if (i === j) {
        let v = A[i][i] - s;
        if (v <= 0) {
          if (v < -1e-12) return false;
          v = 1e-15;
        }
        out[i][i] = Math.sqrt(v);
      } else {
        out[i][j] = (A[i][j] - s) / out[j][j];
      }
    }
  }
  return true;
}

export function cholSolve4(L: Float64Array[], b: Float64Array): void {
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

export function ensureDiag(S: Float64Array[]): void {
  for (let i = 0; i < S.length; i++) {
    if (S[i][i] < 0) {
      for (let k = i; k < S.length; k++) S[k][i] = -S[k][i];
    }
  }
}

export function traceOfP(S: Float64Array[]): number {
  let t = 0;
  for (let i = 0; i < S.length; i++)
    for (let j = 0; j <= i; j++)
      t += S[i][j] * S[i][j];
  return t;
}

export function qrInPlace(A: Float64Array[], m: number, n: number, vBuf: Float64Array): void {
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

export function wrapAngle(a: number): number {
  a = a % TWO_PI;
  if (a > Math.PI) a -= TWO_PI;
  if (a <= -Math.PI) a += TWO_PI;
  return a;
}

export function copySfromQR(Q: Float64Array[], offset: number, S: Float64Array[], N: number): number {
  for (let i = 0; i < N; i++) {
    for (let j = 0; j <= i; j++) S[i][j] = Q[offset + j][offset + i];
    for (let j = i + 1; j < N; j++) S[i][j] = 0;
  }
  ensureDiag(S);
  let tr = 0;
  for (let i = 0; i < N; i++) {
    for (let j = 0; j <= i; j++) {
      const v = S[i][j];
      tr += v * v;
    }
  }
  return tr;
}
