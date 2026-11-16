// The tree of cw_candidates.md: edges of the tau-preserving flow, events where a marked resonance pair
// arrives at S^1, and candidates (milestones 3-5 of PLAN.md).
//
// Design: the internal tau is an input of everything here, never a constant, so that the whole tree can
// be recomputed along a path of tau (warm starts: every solver accepts previous solutions as initial
// guesses). The external tau of a candidate (of the lattice on which it closes) is derived and reported.
//
// Conventions (code normalization, see resonance.js). The torus lattice Lambda lies in the z-plane, and
// the resonance vectors are gamma = 2 pi i k for k in Lambda^*, k = (k_1, k_2) in R^2. The internal lattice
// is Lambda^* = t i <1, tau> (as complex numbers k_1 + i k_2): the factor i puts the sigma-node of the
// pair for 1 at lam = -1, where the opened genus 1 curves (alpha = -r, r -> 1) continue the genus 0 curve
// nu^2 = -lam without a rotation of the z-plane (at a double root lam_0 the nodal curve normalizes to
// nu^2 = conj(lam_0) lam, which is the genus 0 curve exactly when lam_0 = -1).
//
// Genus 0 edge: the curve is fixed and t decreases from infinity (the vacuum, W = 0); the pair of k
// arrives at S^1 when pi |k| = 1. Its sigma-node is at lam = pi^2 k^2.
// Genus 1 edge: alpha = -r, r from 1 (the node) towards 0. Theta_1 is exact (rotation fixed by alpha < 0),
// the B-period of q is 2 pi i k_1, which gives t(r) = |B-period| / (2 pi); the tau and third pairs are
// continued in r, and an event is the arrival of one of them at S^1 (a sigma-node).

import { QMap, genus0Resonance, classifyPair } from './resonance.js';
import { basicPeriods, thetaPoly } from './periods.js';
import { killingField } from './cmc.js';
import { PairIntegrator, leftEigenvector, nuAt } from './gensym.js';
import { leafData, willmoreFrame } from './leaf.js';

const cmul = ([a, b], [c, d]) => [a * c - b * d, a * d + b * c];
const dist = (p, q) => Math.hypot(p[0] - q[0], p[1] - q[1]);

/** The marked vectors 1, tau, third (the shortest representatives of the nonzero classes of Gamma / 2 Gamma). */
export function markedVectors(tau) {
  const m = [tau[0] - 1, tau[1]], p = [tau[0] + 1, tau[1]];
  return { one: [1, 0], tau: tau.slice(), third: Math.hypot(...m) <= Math.hypot(...p) ? m : p };
}

/** k = t i z for z in the internal lattice <1, tau> (complex), as a vector of R^2. */
export const kOf = (t, [x, y]) => [-t * y, t * x];

/** Positively oriented basis (w_1, w_2) of Lambda dual to (k_1, k_2) = (k(1), k(tau)): <w_i, k_j> = delta_ij. */
export function torusBasis(t, tau) {
  const k1 = kOf(t, [1, 0]), k2 = kOf(t, tau);
  const det = k1[0] * k2[1] - k1[1] * k2[0];
  let w1 = [k2[1] / det, -k2[0] / det], w2 = [-k1[1] / det, k1[0] / det];
  if (w1[0] * w2[1] - w1[1] * w2[0] < 0) w2 = [-w2[0], -w2[1]];
  return [w1, w2];
}

// ---------------------------------------------------------------- genus 0 edge

/**
 * The genus 0 edge for internal tau: arrivals and W(t) = area of Lambda = 1 / (t^2 Im tau) (the residue
 * formula on Lambda; at t = infinity the vacuum, W = 0). The edge ends at the arrival of the pair for 1,
 * t = 1/pi (the starting scale, F4). Returns { tStart, W(t), arrivals: [{ name, t, lam, W }] }.
 */
export function genus0Edge(tau) {
  const mk = markedVectors(tau);
  const W = (t) => 1 / (t * t * tau[1]);
  const arrivals = Object.entries(mk).map(([name, z]) => {
    const t = 1 / (Math.PI * Math.hypot(...z));
    const k = kOf(t, z), lam = [Math.PI ** 2 * (k[0] ** 2 - k[1] ** 2), Math.PI ** 2 * 2 * k[0] * k[1]];
    return { name, t, lam, W: W(t) };
  }).sort((a, b) => b.t - a.t);
  return { tStart: 1 / Math.PI, W, arrivals };
}

// ---------------------------------------------------------------- genus 1 edge

/** t(r) for alpha = -r: the B-period of q over 2 pi (Theta_1 is exact there). */
export function genus1Scale(r) {
  const al = [[-r, 0]], P = basicPeriods(al);
  const B = [[1, 0], [0, 1]].map((w) => {
    const p = thetaPoly(al, w, P);
    return p.reduce((a, c, k) => [a[0] + c[0] * P[0][k][0] - c[1] * P[0][k][1], a[1] + c[0] * P[0][k][1] + c[1] * P[0][k][0]], [0, 0]);
  });
  return Math.hypot(B[0][1], B[1][1]) / (2 * Math.PI);
}

const nearestSheet = (Q, p) => {
  const ns = Q.nus(p.lam);
  return { lam: p.lam.slice(), nu: dist(ns[0], p.nu) <= dist(ns[1], p.nu) ? ns[0] : ns[1] };
};

/** Solve the pair named name at r, continuing from start = { P, Q }. */
function pairAt(tau, name, r, start) {
  const al = [[-r, 0]], Q = new QMap(al), t = genus1Scale(r);
  const z = markedVectors(tau)[name], k = kOf(t, z);
  const gam = [[0, 2 * Math.PI * k[0]], [0, 2 * Math.PI * k[1]]];
  const R = Q.resonance(nearestSheet(Q, start.P), nearestSheet(Q, start.Q), gam, { lattice: [kOf(t, [1, 0]), kOf(t, tau)] });
  return { ...R, t, cls: classifyPair(R.P, R.Q) };
}

/**
 * The genus 1 edge for internal tau, from the opening of the pair for 1 at lam = -1 (r = 1) towards r = 0.
 * Tracks the tau and third pairs from their genus 0 positions at the starting scale (nu_1 = nu_0 (lam + 1)
 * near the node), samples { r, t, W } (W by the residue formula on Lambda(t)), and refines each arrival at
 * S^1 by bisection (continuing from the side where the pair is off S^1).
 * Options: rMin (stop), factor (step in r), maxEvents.
 * Returns { samples, events: [{ name, r, t, W, P, Q, lamNode }], end: { r, why } }.
 */
export function genus1Edge(tau, { rStart = 0.9999, rMin = 0.03, factor = 0.97, maxEvents = 1 } = {}) {
  const t0 = 1 / Math.PI, names = ['tau', 'third'];
  const mk = markedVectors(tau);
  const lift = (p) => ({ lam: p.lam, nu: cmul(p.nu, [p.lam[0] + 1, p.lam[1]]) });
  const cur = {};
  for (const n of names) { const s = genus0Resonance(kOf(t0, mk[n]))[0]; cur[n] = { P: lift(s.P), Q: lift(s.Q) }; }
  const samples = [], events = [];
  let r = rStart, end = { r: rMin, why: 'rMin' };
  const Wint = (rr, t) => Math.abs(willmoreFrame(leafData([[-rr, 0]]), ...torusBasis(t, tau)));
  let first = true;
  while (r > rMin) {
    let next = first ? r : Math.max(rMin, r * factor);
    first = false;
    const res = {};
    let ok = true;
    for (const n of names) {
      res[n] = pairAt(tau, n, next, cur[n]);
      if (!res[n].ok) ok = false;
    }
    if (!ok) { end = { r: next, why: 'continuation failed' }; break; }
    // an arrival between r and next: the pair was off and is no longer off
    for (const n of names) {
      if (r !== next && cur[n].cls && cur[n].cls.kind === 'off' && res[n].cls.kind !== 'off') {
        let hi = r, lo = next, sHi = cur[n];
        for (let i = 0; i < 50 && hi - lo > 1e-13; i++) {
          const m = (hi + lo) / 2, R = pairAt(tau, n, m, sHi);
          if (R.ok && R.cls.kind === 'off') { hi = m; sHi = { P: R.P, Q: R.Q, cls: R.cls }; } else lo = m;
        }
        const R = pairAt(tau, n, hi, sHi);
        const mid = [(R.P.lam[0] + R.Q.lam[0]) / 2, (R.P.lam[1] + R.Q.lam[1]) / 2], l = Math.hypot(...mid);
        events.push({ name: n, r: hi, t: R.t, W: Wint(hi, R.t), P: R.P, Q: R.Q, lamNode: [mid[0] / l, mid[1] / l], gap: dist(R.P.lam, R.Q.lam) });
      }
    }
    for (const n of names) cur[n] = { P: res[n].P, Q: res[n].Q, cls: res[n].cls };
    r = next;
    const t = res.tau.t;
    samples.push({ r, t, W: Wint(r, t), pairs: Object.fromEntries(names.map((n) => [n, res[n].cls])) });
    if (events.length >= maxEvents) { end = { r, why: 'event' }; break; }
  }
  return { samples, events, end };
}

// ---------------------------------------------------------------- candidates

/**
 * A candidate from two generalized Sym points on a genus g curve with internal lattice Lambda(t, tau):
 * syms = [{ lam, nu }, { lam, nu }] (a branch point has nu = 0). Computes the monodromy at both by
 * shooting along the basis of Lambda, the sign characters, and the lattice on which f closes: Lambda if
 * the characters agree, otherwise the index 2 sublattice on which they do (F8).
 * Returns { lattice: [e1, e2], doubled, chars, us, W (residue formula on the closing lattice), tauExt,
 *   defects: { node (|M - (±I)| at sigma-nodes), closing (|f(z + e) - f(z)| for e = e1, e2) } }.
 */
export function candidate(alphas, t, tau, syms, { h = 0.002 } = {}) {
  const Z = killingField(alphas);
  const lams = syms.map((s) => s.lam);
  const I = new PairIntegrator(alphas, lams);
  const us = syms.map((s) => leftEigenvector(Z, s.lam, s.nu));
  const [w1, w2] = torusBasis(t, tau);
  const shoot = (w) => { const y = I.initialState(Z); const L = Math.hypot(...w); I.march(y, w[0] / L, w[1] / L, L, h); return y; };
  const f = 8 * I.n;
  let nodeDefect = 0;
  // character: the eigenvalue of M_w(lam_j) on the eigenline u_j (u M = mu u), real ±1
  const chars = [w1, w2].map((w) => {
    const y = shoot(w);
    return [0, 1].map((j) => {
      const M = y.slice(f + 8 * j, f + 8 * j + 8), u = us[j];
      // (u M)_c / u_c for the larger component
      const c = Math.hypot(u[0], u[1]) >= Math.hypot(u[2], u[3]) ? 0 : 1;
      const uc = c === 0 ? [u[0], u[1]] : [u[2], u[3]];
      const um = cmul([u[0], u[1]], [M[2 * c], M[2 * c + 1]]);
      const um2 = cmul([u[2], u[3]], [M[4 + 2 * c], M[5 + 2 * c]]);
      const mu = [(um[0] + um2[0]), (um[1] + um2[1])];
      const d = uc[0] ** 2 + uc[1] ** 2, val = (mu[0] * uc[0] + mu[1] * uc[1]) / d;
      if (Math.hypot(syms[j].nu[0], syms[j].nu[1]) > 1e-12) {
        // a sigma-node: M = ±I
        const s = Math.sign(val);
        nodeDefect = Math.max(nodeDefect, Math.hypot(M[0] - s, M[1], M[2], M[3], M[4], M[5], M[6] - s, M[7]));
      }
      return Math.sign(val);
    });
  });
  const agree = chars.map(([a, b]) => a === b);
  let lattice, doubled = false;
  if (agree[0] && agree[1]) lattice = [w1, w2];
  else {
    doubled = true;
    const add = (a, b) => [a[0] + b[0], a[1] + b[1]], dbl = (a) => [2 * a[0], 2 * a[1]];
    lattice = agree[0] ? [w1, dbl(w2)] : agree[1] ? [dbl(w1), w2] : [add(w1, w2), dbl(w1)];
  }
  const [e1, e2] = lattice;
  // closing check: f(z + e) = f(z) at a sample point
  const fAt = (z) => { const y = I.initialState(Z); const L = Math.hypot(...z); if (L > 0) I.march(y, z[0] / L, z[1] / L, L, h); const out = new Float64Array(4); I.evaluate4(y, us, out, 0); return out; };
  const z0 = [0.3 * e1[0] + 0.2 * e2[0], 0.3 * e1[1] + 0.2 * e2[1]], f0 = fAt(z0);
  const closing = Math.max(...[e1, e2].map((e) => Math.max(...fAt([z0[0] + e[0], z0[1] + e[1]]).map((x, i) => Math.abs(x - f0[i])))));
  const W = Math.abs(willmoreFrame(leafData(alphas), e1, e2));
  const den = e1[0] ** 2 + e1[1] ** 2;
  const tauExt = [(e2[0] * e1[0] + e2[1] * e1[1]) / den, (e2[1] * e1[0] - e2[0] * e1[1]) / den];
  return { lattice, doubled, chars, us, W, tauExt: reduceTau(tauExt), defects: { node: nodeDefect, closing } };
}

/** Reduce tau (upper half plane) into the fundamental domain F_1 of SL(2, Z). */
export function reduceTau([x, y]) {
  if (y < 0) { x = -x; y = -y; }
  for (let i = 0; i < 100; i++) {
    x -= Math.round(x);
    const n = x * x + y * y;
    if (n >= 1 - 1e-12) break;
    x = -x / n; y = y / n;
  }
  return [x, y];
}

/** The genus 0 nu of a point is irrelevant here; for a branch point nu = 0, for a node either sheet. */
export const symPoint = (alphas, lam, nu = null) => ({ lam, nu: nu ?? nuAt(alphas, lam, 1) });

/**
 * Candidate 1 for internal tau: the genus 1 edge to the first arrival (of the tau or the third pair),
 * with the opened branch point alpha and the arriving sigma-node as generalized Sym points.
 * Returns { tau, edge, event, alphas, syms, ...candidate } or { tau, edge, event: null } if none was found.
 */
export function candidate1(tau, opts = {}) {
  const edge = genus1Edge(tau, opts);
  const event = edge.events[0] || null;
  if (!event) return { tau, edge, event: null };
  const alphas = [[-event.r, 0]];
  const syms = [{ lam: alphas[0], nu: [0, 0] }, symPoint(alphas, event.lamNode)];
  return { tau, edge, event, alphas, syms, ...candidate(alphas, event.t, tau, syms) };
}
