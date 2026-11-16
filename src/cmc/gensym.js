// The generalized Sym formula (cw_candidates.md, F7): a conformal immersion into R^4 = H from two
// parallel sections of the associated family at arbitrary spectral parameters.
//
// With F_x = F X(lam), F(0) = I, the parallel sections are the rows psi(z) = u^T F(z, lam), and a row
// (a, b) is the quaternion Q(a, b) = [[a, b], [-conj b, conj a]]. For two rows
//
//   f = Q(psi_1) Q(psi_2)^{-1} in H = R^4,   coordinates (Re a, Im a, Re b, Im b) of the first row.
//
// The sections close up (trivial multiplier) when u_j is a left eigenvector of zeta(0, lam_j) at a
// point of the spectral curve with trivial monodromy; here u_j is the left eigenvector of zeta(0, lam_j)
// for the eigenvalue nu_j / lam_j, nu_j^2 = lam_j a(lam_j), on the sheet chosen by sign_j.
// For lam_1, lam_2 on S^1 this is the S^3 Sym–Bobenko formula F(lam_1) F(lam_2)^{-1}, up to a Möbius
// transformation.

import {
  mulAcc, generator, commutator, projectKilling, FrameIntegrator, killingField, kappa0, shapeFlows,
  flowKilling, aPoly,
} from './cmc.js';
import { polyEval } from './poly.js';

// complex helpers on [re, im]
const cmul = ([ar, ai], [br, bi]) => [ar * br - ai * bi, ar * bi + ai * br];
const cdiv = ([ar, ai], [br, bi]) => { const d = br * br + bi * bi; return [(ar * br + ai * bi) / d, (ai * br - ar * bi) / d]; };
const csqrt = ([r, i]) => {
  const m = Math.hypot(r, i);
  const sr = Math.sqrt(Math.max(0, (m + r) / 2)), si = Math.sqrt(Math.max(0, (m - r) / 2));
  return [sr, i < 0 ? -si : si];
};

/** zeta(lam) = sum_k zeta_k lam^k (k = -1..g) as 8 doubles. */
export function zetaAt(Z, [lr, li]) {
  const n = Z.length / 8, M = new Float64Array(8);
  const d2 = lr * lr + li * li;
  let pr = lr / d2, pi = -li / d2;
  for (let i = 0; i < n; i++) {
    for (let e = 0; e < 4; e++) {
      const zr = Z[8 * i + 2 * e], zi = Z[8 * i + 2 * e + 1];
      M[2 * e] += zr * pr - zi * pi;
      M[2 * e + 1] += zr * pi + zi * pr;
    }
    const t = pr * lr - pi * li; pi = pr * li + pi * lr; pr = t;
  }
  return M;
}

/** nu = sign * sqrt(lam a(lam)) (principal root), the point (lam, nu) of the spectral curve. */
export function nuAt(alphas, lam, sign = 1) {
  const a = aPoly(alphas);
  const [s, t] = csqrt(cmul(lam, polyEval(a, lam[0], lam[1])));
  return [sign * s, sign * t];
}

/**
 * Left eigenvector u (u^T zeta = mu u^T) of zeta(lam) = [[A, B], [C, -A]] for mu = nu / lam,
 * normalized to |u| = 1. Returns [u1r, u1i, u2r, u2i].
 */
export function leftEigenvector(Z, lam, nu) {
  const M = zetaAt(Z, lam);
  const A = [M[0], M[1]], B = [M[2], M[3]], C = [M[4], M[5]];
  const mu = cdiv(nu, lam);
  // zeta^T u = mu u: candidates (C, mu - A) and (A + mu, B)
  const c1 = [C, [mu[0] - A[0], mu[1] - A[1]]];
  const c2 = [[A[0] + mu[0], A[1] + mu[1]], B];
  const nrm = (c) => Math.hypot(c[0][0], c[0][1], c[1][0], c[1][1]);
  const c = nrm(c1) >= nrm(c2) ? c1 : c2;
  const s = 1 / nrm(c);
  return [c[0][0] * s, c[0][1] * s, c[1][0] * s, c[1][1] * s];
}

/**
 * Integrator for (zeta, F(lam_1), F(lam_2)): the layout of FrameIntegrator, with F(lam_2) in place of
 * G = dF/dlam, so marching (RK4 or adaptive) is inherited.
 */
export class PairIntegrator extends FrameIntegrator {
  constructor(alphas, lams) {
    super(alphas, 0);
    this.lams = lams.map(([r, i]) => [r, i]);
    this.inv = this.lams.map(([r, i]) => { const d = r * r + i * i; return [r / d, -i / d]; });
    this.Xs = [new Float64Array(8), new Float64Array(8)];
  }

  initialState(Z) {
    const y = new Float64Array(this.len);
    y.set(Z);
    const f = 8 * this.n;
    y[f] = 1; y[f + 6] = 1; y[f + 8] = 1; y[f + 14] = 1;
    return y;
  }

  /** X(lam_j) = G_{-1}/lam_j + G_0 + G_1 lam_j for the generator Gen currently in this.Gen. */
  evalGen(j) {
    const Gen = this.Gen, X = this.Xs[j];
    const [lr, li] = this.lams[j], [ir, ii] = this.inv[j];
    for (let e = 0; e < 4; e++) {
      const ar = Gen[2 * e], ai = Gen[2 * e + 1];
      const br = Gen[8 + 2 * e], bi = Gen[8 + 2 * e + 1];
      const cr = Gen[16 + 2 * e], ci = Gen[16 + 2 * e + 1];
      X[2 * e] = ar * ir - ai * ii + br + cr * lr - ci * li;
      X[2 * e + 1] = ar * ii + ai * ir + bi + cr * li + ci * lr;
    }
    return X;
  }

  rhs(y, dy, cr, ci) {
    const n = this.n;
    generator(y, this.kap, 0, cr, ci, this.Gen);
    commutator(y, n, this.Gen, 0, dy);
    const f = 8 * n;
    for (let e = 0; e < 16; e++) dy[f + e] = 0;
    for (let j = 0; j < 2; j++) mulAcc(y, f + 8 * j, this.evalGen(j), 0, dy, f + 8 * j, 1);
  }

  /** Re-impose the Killing-field structure, and det F_j = 1 (F_j is not unitary off S^1). */
  project(y) {
    projectKilling(y.subarray(0, 8 * this.n));
    const f = 8 * this.n;
    for (let j = 0; j < 2; j++) {
      const o = f + 8 * j;
      const det = cmul([y[o], y[o + 1]], [y[o + 6], y[o + 7]]);
      const bc = cmul([y[o + 2], y[o + 3]], [y[o + 4], y[o + 5]]);
      const s = csqrt([det[0] - bc[0], det[1] - bc[1]]);
      const [sr, si] = cdiv([1, 0], s);
      for (let e = 0; e < 4; e++) {
        const r = y[o + 2 * e], i = y[o + 2 * e + 1];
        y[o + 2 * e] = r * sr - i * si; y[o + 2 * e + 1] = r * si + i * sr;
      }
    }
  }

  /** The rows r_j = u_j^T F_j, as [ar, ai, br, bi]. */
  rows(y, us) {
    const f = 8 * this.n;
    return us.map((u, j) => {
      const o = f + 8 * j;
      const row = (c) => [
        u[0] * y[o + c] - u[1] * y[o + c + 1] + u[2] * y[o + 4 + c] - u[3] * y[o + 4 + c + 1],
        u[0] * y[o + c + 1] + u[1] * y[o + c] + u[2] * y[o + 4 + c + 1] + u[3] * y[o + 4 + c],
      ];
      return [...row(0), ...row(2)];
    });
  }

  /** f = Q(r_1) Q(r_2)^{-1} as 4 reals, written to out at offset o. */
  evaluate4(y, us, out, o = 0) {
    const [r1, r2] = this.rows(y, us);
    qdiv(r1, r2, out, o);
  }

  /**
   * f and its exact first derivatives at state y: [f, f_x, f_y], each 4 reals.
   * Uses r_x = r X(lam) for the x-generator U + V and r_y = r X for i(U - V).
   */
  jet(y, us) {
    const r = this.rows(y, us);
    const d = (cr, ci) => {
      generator(y, this.kap, 0, cr, ci, this.Gen);
      return r.map((ri, j) => rowTimes(ri, this.evalGen(j)));
    };
    const rx = d(1, 0), ry = d(0, 1);
    const f = new Float64Array(4);
    qdiv(r[0], r[1], f, 0);
    // f' = (Q(r1') - f Q(r2')) Q(r2)^{-1}
    const der = (dr) => {
      const t = qmul(f, dr[1]);
      const a = [dr[0][0] - t[0], dr[0][1] - t[1], dr[0][2] - t[2], dr[0][3] - t[3]];
      const out = new Float64Array(4);
      qdiv(a, r[1], out, 0);
      return out;
    };
    return [f, der(rx), der(ry)];
  }
}

/** row (a, b) times the 2x2 matrix X (8 doubles). */
function rowTimes([ar, ai, br, bi], X) {
  return [
    ar * X[0] - ai * X[1] + br * X[4] - bi * X[5],
    ar * X[1] + ai * X[0] + br * X[5] + bi * X[4],
    ar * X[2] - ai * X[3] + br * X[6] - bi * X[7],
    ar * X[3] + ai * X[2] + br * X[7] + bi * X[6],
  ];
}

/** Quaternion product in first-row form: Q(a, b) Q(c, d) = Q(ac - b conj d, ad + b conj c). */
export function qmul([ar, ai, br, bi], [cr, ci, dr, di]) {
  return [
    ar * cr - ai * ci - (br * dr + bi * di),
    ar * ci + ai * cr - (bi * dr - br * di),
    ar * dr - ai * di + (br * cr + bi * ci),
    ar * di + ai * dr + (bi * cr - br * ci),
  ];
}

/** out[o..o+3] = Q(p) Q(q)^{-1}, with Q(c, d)^{-1} = Q(conj c, -d) / (|c|^2 + |d|^2). */
export function qdiv(p, [cr, ci, dr, di], out, o) {
  const n = cr * cr + ci * ci + dr * dr + di * di;
  const m = qmul(p, [cr, -ci, -dr, -di]);
  for (let e = 0; e < 4; e++) out[o + e] = m[e] / n;
}

/**
 * Projections R^4 -> R^3. For two Sym points on S^1, f lies on the unit S^3 (|f| = |u_1| / |u_2| = 1),
 * and stereographic projection from e_4 is the usual picture of the S^3 torus in R^3.
 */
export const projections = {
  drop4: (x, y, z, w) => [x, y, z],
  stereo4: (x, y, z, w) => { const d = 1 - w; return [x / d, y / d, z / d]; },
};

/**
 * Project points of R^4 (Float64Array, 4 per point) to R^3 (Float32Array, 3 per point).
 * kind 'auto': if the points lie on a round 3-sphere (relative residual below 1e-6), stereographic
 * projection of that sphere from its point farthest from the surface (the antipode of the direction
 * of the centroid from the centre), onto the hyperplane through the centre; otherwise orthogonal
 * projection onto the three principal axes of the point cloud. 'drop4' and 'stereo4' as in projections.
 * Returns { pos, used } with used one of 'sphere', 'pca', 'drop4', 'stereo4'.
 */
export function projectPoints(P, kind = 'auto') {
  const N = P.length / 4, pos = new Float32Array(3 * N);
  if (kind === 'drop4' || kind === 'stereo4') {
    const f = projections[kind];
    for (let v = 0; v < N; v++) pos.set(f(P[4 * v], P[4 * v + 1], P[4 * v + 2], P[4 * v + 3]), 3 * v);
    return { pos, used: kind };
  }
  const fin = [];
  for (let v = 0; v < N; v++) if ([0, 1, 2, 3].every((e) => Number.isFinite(P[4 * v + e]))) fin.push(v);
  const sph = sphereFit(P, fin);
  if (sph && sph.res < 1e-6) {
    const { c, R } = sph;
    const n = farthestPole(P, fin, c, R);
    const B = complementBasis(n);
    for (let v = 0; v < N; v++) {
      const y = [0, 1, 2, 3].map((e) => P[4 * v + e] - c[e]);
      const yn = y[0] * n[0] + y[1] * n[1] + y[2] * n[2] + y[3] * n[3];
      const k = R / (R - yn);
      for (let a = 0; a < 3; a++) pos[3 * v + a] = k * (y[0] * B[a][0] + y[1] * B[a][1] + y[2] * B[a][2] + y[3] * B[a][3]);
    }
    return { pos, used: 'sphere' };
  }
  // principal axes
  const m = [0, 0, 0, 0];
  for (const v of fin) for (let e = 0; e < 4; e++) m[e] += P[4 * v + e] / fin.length;
  const C = [0, 1, 2, 3].map(() => [0, 0, 0, 0]);
  for (const v of fin) {
    const y = [0, 1, 2, 3].map((e) => P[4 * v + e] - m[e]);
    for (let a = 0; a < 4; a++) for (let b = 0; b < 4; b++) C[a][b] += y[a] * y[b];
  }
  const { vals, vecs } = jacobi4(C);
  const order = [0, 1, 2, 3].sort((a, b) => vals[b] - vals[a]);
  const B = order.slice(0, 3).map((k) => vecs.map((row) => row[k]));
  for (let v = 0; v < N; v++) {
    const y = [0, 1, 2, 3].map((e) => P[4 * v + e] - m[e]);
    for (let a = 0; a < 3; a++) pos[3 * v + a] = y[0] * B[a][0] + y[1] * B[a][1] + y[2] * B[a][2] + y[3] * B[a][3];
  }
  return { pos, used: 'pca' };
}

/**
 * A unit direction n such that the pole c + R n of the sphere is far from the surface: among the
 * antipode of the centroid, the coordinate directions ±e_k and the directions (±1, ±1, ±1, ±1)/2, the
 * one maximizing the minimum distance to (a subsample of) the points.
 */
function farthestPole(P, idx, c, R) {
  const cand = [];
  const m = [0, 0, 0, 0];
  for (const v of idx) for (let e = 0; e < 4; e++) m[e] += P[4 * v + e] - c[e];
  const ml = Math.hypot(...m);
  if (ml > 0) cand.push(m.map((x) => -x / ml));
  for (let k = 0; k < 4; k++) for (const sg of [1, -1]) { const d = [0, 0, 0, 0]; d[k] = sg; cand.push(d); }
  for (let b = 0; b < 16; b++) cand.push([0, 1, 2, 3].map((k) => ((b >> k) & 1 ? 0.5 : -0.5)));
  const step = Math.max(1, Math.floor(idx.length / 4000));
  let best = cand[0], bd = -1;
  for (const d of cand) {
    let mn = Infinity;
    for (let q = 0; q < idx.length; q += step) {
      const v = idx[q];
      let s2 = 0;
      for (let e = 0; e < 4; e++) s2 += (P[4 * v + e] - c[e] - R * d[e]) ** 2;
      mn = Math.min(mn, s2);
    }
    if (mn > bd) { bd = mn; best = d; }
  }
  return best;
}

/** Least-squares sphere |x - c|^2 = R^2 through the points idx; res = max | |x-c|^2 - R^2 | / R^2. */
function sphereFit(P, idx) {
  if (idx.length < 6) return null;
  // 2 c.x + d = |x|^2, d = R^2 - |c|^2: normal equations in (c, d)
  const A = [0, 1, 2, 3, 4].map(() => [0, 0, 0, 0, 0, 0]);
  for (const v of idx) {
    const x = [P[4 * v], P[4 * v + 1], P[4 * v + 2], P[4 * v + 3]];
    const row = [2 * x[0], 2 * x[1], 2 * x[2], 2 * x[3], 1], rhs = x[0] ** 2 + x[1] ** 2 + x[2] ** 2 + x[3] ** 2;
    for (let a = 0; a < 5; a++) { for (let b = 0; b < 5; b++) A[a][b] += row[a] * row[b]; A[a][5] += row[a] * rhs; }
  }
  const sol = solve(A);
  if (!sol) return null;
  const c = sol.slice(0, 4), R2 = sol[4] + c[0] ** 2 + c[1] ** 2 + c[2] ** 2 + c[3] ** 2;
  if (!(R2 > 0)) return null;
  let res = 0;
  for (const v of idx) {
    const d2 = [0, 1, 2, 3].reduce((acc, e) => acc + (P[4 * v + e] - c[e]) ** 2, 0);
    res = Math.max(res, Math.abs(d2 - R2) / R2);
  }
  return { c, R: Math.sqrt(R2), res };
}

/** Gaussian elimination with partial pivoting on an augmented n x (n+1) matrix. */
function solve(M) {
  const n = M.length, A = M.map((r) => r.slice());
  for (let k = 0; k < n; k++) {
    let p = k;
    for (let i = k + 1; i < n; i++) if (Math.abs(A[i][k]) > Math.abs(A[p][k])) p = i;
    if (Math.abs(A[p][k]) < 1e-300) return null;
    [A[k], A[p]] = [A[p], A[k]];
    for (let i = k + 1; i < n; i++) {
      const f = A[i][k] / A[k][k];
      for (let j = k; j <= n; j++) A[i][j] -= f * A[k][j];
    }
  }
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let t = A[i][n];
    for (let j = i + 1; j < n; j++) t -= A[i][j] * x[j];
    x[i] = t / A[i][i];
  }
  return x;
}

/** An orthonormal basis (3 vectors) of the complement of the unit vector n in R^4. */
function complementBasis(n) {
  const B = [];
  for (let k = 0; k < 4 && B.length < 3; k++) {
    let w = [0, 0, 0, 0]; w[k] = 1;
    for (const b of [n, ...B]) {
      const d = w[0] * b[0] + w[1] * b[1] + w[2] * b[2] + w[3] * b[3];
      w = w.map((x, e) => x - d * b[e]);
    }
    const l = Math.hypot(...w);
    if (l > 1e-6) B.push(w.map((x) => x / l));
  }
  return B;
}

/** Eigen-decomposition of a symmetric 4x4 matrix by Jacobi rotations; vecs[:, k] belongs to vals[k]. */
function jacobi4(S) {
  const A = S.map((r) => r.slice()), V = [0, 1, 2, 3].map((i) => [0, 1, 2, 3].map((j) => (i === j ? 1 : 0)));
  for (let sweep = 0; sweep < 50; sweep++) {
    let off = 0;
    for (let p = 0; p < 4; p++) for (let q = p + 1; q < 4; q++) off += A[p][q] ** 2;
    if (off < 1e-30) break;
    for (let p = 0; p < 4; p++) for (let q = p + 1; q < 4; q++) {
      if (Math.abs(A[p][q]) < 1e-300) continue;
      const th = 0.5 * Math.atan2(2 * A[p][q], A[q][q] - A[p][p]);
      const c = Math.cos(th), s = Math.sin(th);
      for (let k = 0; k < 4; k++) { const a = A[k][p], b = A[k][q]; A[k][p] = c * a - s * b; A[k][q] = s * a + c * b; }
      for (let k = 0; k < 4; k++) { const a = A[p][k], b = A[q][k]; A[p][k] = c * a - s * b; A[q][k] = s * a + c * b; }
      for (let k = 0; k < 4; k++) { const a = V[k][p], b = V[k][q]; V[k][p] = c * a - s * b; V[k][q] = s * a + c * b; }
    }
  }
  return { vals: [0, 1, 2, 3].map((i) => A[i][i]), vecs: V };
}

/**
 * The immersion f into R^4 on an nx x ny grid, from two points of the spectral curve.
 *
 * params:
 *   alphas, tau, z0, width, height, nx, ny (or sCoords, tCoords), hmax, tol   as for computeSurface
 *   lams     [[re, im], [re, im]]   the two spectral parameters (in C^*)
 *   signs    [s1, s2]               sheets: nu_j = s_j sqrt(lam_j a(lam_j))
 *   phi      grid rotation (default 0)
 *   projection   'auto' (default), 'stereo4' or 'drop4' (see projectPoints)
 *   cell     optional [e1, e2]: the grid runs along e1 (s) and e2 (t) instead of an orthogonal frame;
 *            with sCoords in [0, |e1|], tCoords in [0, |e2|] and z0 = 0 it covers the parallelogram
 *
 * Returns { nx, ny, pos4 (Float32 4N), pos, nrm (3N, normals of the projected mesh), u, st, bbox, phi,
 *   us (the eigenvectors), projection (the one used), stats }.
 */
export function computeGenSurface(params) {
  const t0 = (typeof performance !== 'undefined' ? performance : Date).now();
  const { alphas, tau = [], z0 = [0, 0], width, height, hmax = 0.02, tol = 0, lams, signs = [1, 1] } = params;
  const nx = params.sCoords ? params.sCoords.length : params.nx;
  const ny = params.tCoords ? params.tCoords.length : params.ny;
  const g = alphas.length;
  const kap = kappa0(alphas);
  const phi = params.phi ?? 0;

  const Z = killingField(alphas);
  shapeFlows(g).forEach(([k, c], i) => { if (tau[i]) flowKilling(Z, kap, k, c, tau[i], 0.002); });
  // eigenvectors at z = 0, where F = I
  const us = lams.map((lam, j) => leftEigenvector(Z, lam, nuAt(alphas, lam, signs[j])));

  const I = new PairIntegrator(alphas, lams);
  const y0 = I.initialState(Z);
  const r0 = Math.hypot(z0[0], z0[1]);
  if (r0 > 0) I.march(y0, z0[0] / r0, z0[1] / r0, r0, hmax, tol);
  I.project(y0);

  const N = nx * ny;
  const P4 = new Float64Array(4 * N), nrm = new Float32Array(3 * N);
  const u = new Float32Array(N), st = new Float32Array(2 * N);
  const s = params.sCoords ? Array.from(params.sCoords)
    : Array.from({ length: nx }, (_, i) => -width / 2 + (nx > 1 ? (i * width) / (nx - 1) : 0));
  const t = params.tCoords ? Array.from(params.tCoords)
    : Array.from({ length: ny }, (_, j) => -height / 2 + (ny > 1 ? (j * height) / (ny - 1) : 0));
  let csr = Math.cos(phi), csi = Math.sin(phi);
  let ctr = -csi, cti = csr;
  let w0 = [csr * z0[0] + csi * z0[1], csr * z0[1] - csi * z0[0]];
  if (params.cell) {
    // a parallelogram domain: s along e1, t along e2 (coordinates are lengths along them)
    const [e1, e2] = params.cell, l1 = Math.hypot(...e1), l2 = Math.hypot(...e2);
    [csr, csi, ctr, cti] = [e1[0] / l1, e1[1] / l1, e2[0] / l2, e2[1] / l2];
    w0 = [0, 0];
  }
  const F4 = new Float64Array(4);

  const sweep = (ystart, coords, cr, ci, visit) => {
    for (const dir of [1, -1]) {
      const y = ystart.slice();
      let at = 0;
      const idx = coords.map((_, i) => i).filter((i) => (dir > 0 ? coords[i] >= 0 : coords[i] < 0));
      if (dir < 0) idx.reverse();
      for (const i of idx) {
        I.march(y, cr, ci, coords[i] - at, hmax, tol);
        at = coords[i];
        I.project(y);
        visit(i, y);
      }
    }
  };

  sweep(y0, s, csr, csi, (i, ys) => {
    sweep(ys, t, ctr, cti, (j, yt) => {
      const v = j * nx + i;
      I.evaluate4(yt, us, F4, 0);
      for (let e = 0; e < 4; e++) P4[4 * v + e] = F4[e];
      u[v] = Math.log((2 * Math.hypot(yt[2], yt[3])) / kap);
      st[2 * v] = s[i] + w0[0]; st[2 * v + 1] = t[j] + w0[1];
    });
  });

  const { pos, used } = projectPoints(P4, params.projection ?? 'auto');
  const pos4 = Float32Array.from(P4);
  meshNormals(pos, nx, ny, nrm);
  const bbox = [Infinity, Infinity, Infinity, -Infinity, -Infinity, -Infinity];
  for (let v = 0; v < N; v++) {
    for (let e = 0; e < 3; e++) {
      const x = pos[3 * v + e];
      if (Number.isFinite(x)) { if (x < bbox[e]) bbox[e] = x; if (x > bbox[e + 3]) bbox[e + 3] = x; }
    }
  }
  const t1 = (typeof performance !== 'undefined' ? performance : Date).now();
  return { nx, ny, pos4, pos, nrm, u, st, bbox, phi, us, projection: used, stats: { ms: t1 - t0 } };
}

/** Unit normals of a grid mesh by central differences (one-sided at the edges). */
export function meshNormals(pos, nx, ny, nrm) {
  const P = (i, j, e) => pos[3 * (j * nx + i) + e];
  for (let j = 0; j < ny; j++) {
    for (let i = 0; i < nx; i++) {
      const i0 = Math.max(0, i - 1), i1 = Math.min(nx - 1, i + 1);
      const j0 = Math.max(0, j - 1), j1 = Math.min(ny - 1, j + 1);
      const a = [0, 1, 2].map((e) => P(i1, j, e) - P(i0, j, e));
      const b = [0, 1, 2].map((e) => P(i, j1, e) - P(i, j0, e));
      const c = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
      const l = Math.hypot(c[0], c[1], c[2]) || 1;
      const v = 3 * (j * nx + i);
      nrm[v] = c[0] / l; nrm[v + 1] = c[1] / l; nrm[v + 2] = c[2] / l;
    }
  }
}
