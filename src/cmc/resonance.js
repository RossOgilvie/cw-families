// The map q and resonance pairs (cw_candidates.md, F3 and F6).
//
// q = (q^1, q^2) with q^1 = log mu_1, q^2 = log mu_i: the primitives of the real-normalized differentials
// Theta_w = p_w(lam) dlam / (lam nu) for w = 1 and w = i (periods.js), nu^2 = lam a(lam). Their periods
// are imaginary. For a period w of zeta, log mu_w(P) = Re(w) q^1(P) + Im(w) q^2(P) (mod 2 pi i), where
// mu_w is the eigenvalue of the monodromy on the eigenline of P; with F_x = F X the eigenline of the
// point (lam, nu) is the left eigenvector of zeta(0, lam) for nu / lam (gensym.js), and in genus 0
// q^1, q^2 are its eigenvalues for X = U + V and X = i(U - V).
//
// A resonance pair for gamma is a pair of points P != Q with q(P) - q(Q) = gamma; for a torus with
// lattice Lambda, gamma runs over 2 pi i Lambda^* (as vectors (k_1, k_2), q-space being C^2).
//
// Base point: a finite branch point (alpha_1), so that q o sigma = -q (sigma: nu -> -nu), modulo periods.
// In genus 0 there is none, and q_w = (w nu / lam + conj(w) nu) / kappa0 is the odd primitive in closed form.
// q is computed by Gauss–Legendre quadrature along straight segments, continuing nu by continuity.

import { aPoly, kappa0 } from './cmc.js';
import { basicPeriods, thetaPoly } from './periods.js';
import { polyEval } from './poly.js';

const cmul = ([ar, ai], [br, bi]) => [ar * br - ai * bi, ar * bi + ai * br];
const cdiv = ([ar, ai], [br, bi]) => { const d = br * br + bi * bi; return [(ar * br + ai * bi) / d, (ai * br - ar * bi) / d]; };
const cadd = (a, b) => [a[0] + b[0], a[1] + b[1]];
const csub = (a, b) => [a[0] - b[0], a[1] - b[1]];
const cabs = ([r, i]) => Math.hypot(r, i);
const csqrt = ([r, i]) => {
  const m = Math.hypot(r, i);
  const sr = Math.sqrt(Math.max(0, (m + r) / 2)), si = Math.sqrt(Math.max(0, (m - r) / 2));
  return [sr, i < 0 ? -si : si];
};
const pev = (p, [lr, li]) => {
  let r = 0, i = 0;
  for (let k = p.length - 1; k >= 0; k--) { const t = r * lr - i * li + p[k][0]; i = r * li + i * lr + p[k][1]; r = t; }
  return [r, i];
};

// 8-point Gauss–Legendre on [0, 1]
const GL_X = [-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
  0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363].map((x) => (x + 1) / 2);
const GL_W = [0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
  0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763].map((w) => w / 2);

export class QMap {
  constructor(alphas) {
    this.alphas = alphas.map((a) => a.slice());
    this.g = alphas.length;
    const a = aPoly(alphas);
    this.a = Array.from({ length: a.length / 2 }, (_, k) => [a[2 * k], a[2 * k + 1]]);
    this.kap = kappa0(alphas);
    if (this.g > 0) {
      const P = basicPeriods(alphas);
      this.p = [thetaPoly(alphas, [1, 0], P), thetaPoly(alphas, [0, 1], P)];
      this.periods = P;
    } else {
      this.p = [thetaPoly(alphas, [1, 0]), thetaPoly(alphas, [0, 1])];
      this.periods = [];
    }
    // branch points in C^* (0 and infinity are poles of Theta, not usable as base points)
    this.branch = [];
    for (const [re, im] of alphas) { const r2 = re * re + im * im; this.branch.push([re, im], [re / r2, im / r2]); }
  }

  /** nu^2 = lam a(lam). */
  nu2(lam) { return cmul(lam, pev(this.a, lam)); }

  /** The two values of nu over lam, the first being the principal square root. */
  nus(lam) { const s = csqrt(this.nu2(lam)); return [s, [-s[0], -s[1]]]; }

  /** theta(lam, nu) = (p_1(lam), p_i(lam)) / (lam nu): the derivative of q with respect to lam. */
  theta(lam, nu) {
    const d = cmul(lam, nu);
    return this.p.map((p) => cdiv(pev(p, lam), d));
  }

  /** Distance from lam to the nearest singular point (0 or a branch point other than skip). */
  _dist(lam, skip = null) {
    let d = cabs(lam);
    for (const e of this.branch) if (!skip || cabs(csub(e, skip)) > 1e-14) d = Math.min(d, cabs(csub(lam, e)));
    return d;
  }

  /**
   * Integrate along the straight segment from the point P = { lam, nu } to lam1, continuing nu.
   * Returns { lam: lam1, nu, dq: [dq1, dq2] }. If either end is a branch point, the segment is
   * parametrised by lam = e + D (3 t^2 - 2 t^3), whose derivative vanishes at both ends and cancels the
   * square-root singularities. From a branch point (nu = 0) the first node takes the principal sheet;
   * the caller compares the sheet at the end with the one it wants (q o sigma = -q).
   */
  segment(P, lam1) {
    const e = P.lam, D = csub(lam1, e), L = cabs(D);
    const dq = [[0, 0], [0, 0]];
    if (L === 0) return { lam: lam1.slice(), nu: P.nu.slice(), dq };
    const atStart = cabs(P.nu) < 1e-14;
    const atEnd = this.branch.some((b) => cabs(csub(b, lam1)) < 1e-13 * Math.max(1, cabs(b)));
    const smooth = atStart || atEnd;
    // number of subintervals: resolve the nearest singular point along the segment (other than the ends
    // when they are branch points, which the parametrisation takes care of)
    let dmin = Infinity;
    for (let s = 1; s < 16; s++) {
      const t = s / 16, z = cadd(e, [D[0] * t, D[1] * t]);
      let d = cabs(z);
      for (const b of this.branch) {
        if ((atStart && cabs(csub(b, e)) < 1e-13) || (atEnd && cabs(csub(b, lam1)) < 1e-13 * Math.max(1, cabs(b)))) continue;
        d = Math.min(d, cabs(csub(z, b)));
      }
      dmin = Math.min(dmin, d);
    }
    const m = Math.min(4000, Math.max(smooth ? 4 : 2, Math.ceil((4 * L) / Math.max(dmin, 1e-6))));
    let nuPrev = atStart ? null : P.nu;
    for (let k = 0; k < m; k++) {
      for (let q = 0; q < 8; q++) {
        const t = (k + GL_X[q]) / m, w = GL_W[q] / m;
        const ph = smooth ? t * t * (3 - 2 * t) : t, dph = smooth ? 6 * t * (1 - t) : 1;
        const lam = cadd(e, [D[0] * ph, D[1] * ph]), dl = [D[0] * dph, D[1] * dph];
        let nu = csqrt(this.nu2(lam));
        if (nuPrev && cabs(csub(nu, nuPrev)) > cabs(cadd(nu, nuPrev))) nu = [-nu[0], -nu[1]];
        nuPrev = nu;
        const th = this.theta(lam, nu);
        for (let j = 0; j < 2; j++) { const v = cmul(th[j], dl); dq[j][0] += w * v[0]; dq[j][1] += w * v[1]; }
      }
    }
    let nu = atEnd ? [0, 0] : csqrt(this.nu2(lam1));
    if (!atEnd && cabs(csub(nu, nuPrev)) > cabs(cadd(nu, nuPrev))) nu = [-nu[0], -nu[1]];
    return { lam: lam1.slice(), nu, dq };
  }

  /**
   * q at the point (lam, nu): closed form in genus 0; otherwise integrated from the base branch point
   * alpha_1 along a straight segment (q o sigma = -q decides the sign on the other sheet).
   */
  q(lam, nu) {
    if (this.g === 0) {
      return [[1, 0], [0, 1]].map((w) => {
        const a = cdiv(cmul(w, nu), lam), b = cmul([w[0], -w[1]], nu);
        return [(a[0] + b[0]) / this.kap, (a[1] + b[1]) / this.kap];
      });
    }
    const base = { lam: this.alphas[0], nu: [0, 0] };
    const r = this.segment(base, lam);
    const flip = cabs(csub(r.nu, nu)) > cabs(cadd(r.nu, nu));
    return flip ? r.dq.map(([x, y]) => [-x, -y]) : r.dq;
  }

  /**
   * Solve q(P) - q(Q) = gamma for points P, Q near the given ones, by Newton in (lam_P, lam_Q) with q
   * carried along incrementally (so the sheets and the period ambiguity stay continuous).
   * P, Q: { lam, nu }; gamma: [g1, g2] complex. Returns { P, Q, qP, qQ, res, iters, ok }.
   */
  resonance(P, Q, gamma, { tol = 1e-12, maxIter = 60 } = {}) {
    let p = { lam: P.lam.slice(), nu: P.nu.slice() }, qq = { lam: Q.lam.slice(), nu: Q.nu.slice() };
    let qP = this.q(p.lam, p.nu), qQ = this.q(qq.lam, qq.nu);
    const resid = () => [0, 1].map((j) => csub(csub(qP[j], qQ[j]), gamma[j]));
    let r = resid(), nr = Math.max(cabs(r[0]), cabs(r[1])), iters = 0;
    while (nr > tol && iters < maxIter) {
      iters++;
      const tP = this.theta(p.lam, p.nu), tQ = this.theta(qq.lam, qq.nu);
      // [tP, -tQ] (dP, dQ)^T = -r
      const a = tP[0], b = [-tQ[0][0], -tQ[0][1]], c = tP[1], d = [-tQ[1][0], -tQ[1][1]];
      const det = csub(cmul(a, d), cmul(b, c));
      if (cabs(det) < 1e-300) break;
      const r0 = [-r[0][0], -r[0][1]], r1 = [-r[1][0], -r[1][1]];
      let dP = cdiv(csub(cmul(d, r0), cmul(b, r1)), det);
      let dQ = cdiv(csub(cmul(a, r1), cmul(c, r0)), det);
      // damping: stay well clear of the singular points
      const lim = 0.3 * Math.min(this._dist(p.lam), this._dist(qq.lam));
      const big = Math.max(cabs(dP), cabs(dQ));
      if (big > lim) { const s = lim / big; dP = [dP[0] * s, dP[1] * s]; dQ = [dQ[0] * s, dQ[1] * s]; }
      const sp = this.segment(p, cadd(p.lam, dP)), sq = this.segment(qq, cadd(qq.lam, dQ));
      p = { lam: sp.lam, nu: sp.nu }; qq = { lam: sq.lam, nu: sq.nu };
      qP = [0, 1].map((j) => cadd(qP[j], sp.dq[j]));
      qQ = [0, 1].map((j) => cadd(qQ[j], sq.dq[j]));
      r = resid();
      nr = Math.max(cabs(r[0]), cabs(r[1]));
    }
    return { P: p, Q: qq, qP, qQ, res: nr, iters, ok: nr <= tol };
  }
}

/**
 * Genus 0 resonance pairs in closed form (cw_candidates.md, F3, in the code's normalization nu^2 = -lam,
 * kappa0 = 2): for k = k_1 + i k_2 (gamma = 2 pi i (k_1, k_2)), nu_P - nu_Q = 2 pi i k and
 * nu_P nu_Q = k / conj(k). Returns the two ordered solutions [{ P, Q }, { P, Q }], points as { lam, nu }.
 * Over S^1 (|nu| = 1) iff |k| <= 1/pi; at |k| = 1/pi the two solutions coincide in a sigma-node.
 */
export function genus0Resonance([k1, k2]) {
  const k = [k1, k2], D = [-2 * Math.PI * k2, 2 * Math.PI * k1]; // 2 pi i k
  const prod = cdiv(k, [k1, -k2]);
  // nu_Q^2 + D nu_Q - prod = 0
  const disc = csqrt(cadd(cmul(D, D), [4 * prod[0], 4 * prod[1]]));
  return [1, -1].map((s) => {
    const nuQ = [(-D[0] + s * disc[0]) / 2, (-D[1] + s * disc[1]) / 2];
    const nuP = cadd(nuQ, D);
    const pt = (nu) => ({ lam: [-(nu[0] * nu[0] - nu[1] * nu[1]), -2 * nu[0] * nu[1]], nu });
    return { P: pt(nuP), Q: pt(nuQ) };
  });
}

/**
 * Classify a resonance pair: 'off' (a rho-pair off S^1), 'on' (over S^1, over two values of lam) or
 * 'node' (a sigma-node over S^1: lam_P = lam_Q, nu_Q = -nu_P). dist is the signed distance to arrival
 * at S^1: log|lam_P| for a pair off S^1 (it vanishes as the pair arrives), and for a pair on S^1 the
 * (positive) angle between lam_P and lam_Q, which vanishes when it closes up into a sigma-node.
 */
export function classifyPair(P, Q, tol = 1e-8) {
  const rP = cabs(P.lam), rQ = cabs(Q.lam);
  const onCircle = Math.abs(rP - 1) < tol && Math.abs(rQ - 1) < tol;
  if (!onCircle) return { kind: 'off', dist: Math.log(rP) };
  const ang = Math.abs(Math.atan2(P.lam[0] * Q.lam[1] - P.lam[1] * Q.lam[0], P.lam[0] * Q.lam[0] + P.lam[1] * Q.lam[1]));
  return { kind: ang < Math.sqrt(tol) ? 'node' : 'on', dist: ang };
}
