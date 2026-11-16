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
import { WhithamCurve, ellVector } from './whitham.js';

const cmul = ([a, b], [c, d]) => [a * c - b * d, a * d + b * c];
const dist = (p, q) => Math.hypot(p[0] - q[0], p[1] - q[1]);

/**
 * The shortest vector of each nonzero class of <1, tau> / 2 <1, tau>, with the labels of the tau pad
 * (classOrder in leaf.js): class 1 = [1] (p odd, q even), 2 = [tau] (p even, q odd), 3 = [1 + tau].
 * Returns [{ v: [x, y], len, pq: [p, q] }] indexed by class - 1.
 */
export function classReps(tau) {
  const best = [null, null, null];
  for (let p = -8; p <= 8; p++) for (let q = -8; q <= 8; q++) {
    const pp = Math.abs(p) % 2, qq = Math.abs(q) % 2;
    if (!pp && !qq) continue;
    const c = pp && !qq ? 0 : !pp && qq ? 1 : 2;
    const v = [p + q * tau[0], q * tau[1]], len = Math.hypot(...v);
    if (!best[c] || len < best[c].len - 1e-14) best[c] = { v, len, pq: [p, q] };
  }
  return best;
}

/**
 * The internal tau as a point of F_6 (a lattice with labelled classes): the pair of the shortest class
 * is opened first, so the lattice is rescaled to put its shortest vector v at 1, and the shortest vector
 * u of the second class gives the reduced tau1 = u / v (in F_1). The computation runs with tau1, whose
 * marked vectors one, tau, third are the shortest vectors of the classes classes.one, .tau, .third.
 * Returns { tau1, classes, order } (order as the pad's labels, e.g. '123').
 */
export function internalFrame(tau) {
  const reps = classReps(tau);
  const order = [0, 1, 2].sort((a, b) => reps[a].len - reps[b].len);
  const v = reps[order[0]].v, u = reps[order[1]].v, d = v[0] * v[0] + v[1] * v[1];
  let tau1 = [(u[0] * v[0] + u[1] * v[1]) / d, (u[1] * v[0] - u[0] * v[1]) / d];
  if (tau1[1] < 0) tau1 = [-tau1[0], -tau1[1]];
  return { tau1, classes: { one: order[0] + 1, tau: order[1] + 1, third: order[2] + 1 }, order: order.map((c) => c + 1).join('') };
}

/** The marked vectors 1, tau, third (the shortest representatives of the nonzero classes of Gamma / 2 Gamma). */
export function markedVectors(tau) {
  const m = [tau[0] - 1, tau[1]], p = [tau[0] + 1, tau[1]];
  return { one: [1, 0], tau: tau.slice(), third: Math.hypot(...m) <= Math.hypot(...p) ? m : p };
}

/** k = t i z for z in the internal lattice <1, tau> (complex), as a vector of R^2. */
export const kOf = (t, [x, y]) => [-t * y, t * x];

/** Positively oriented basis of the lattice dual to the lattice spanned by k1, k2 (vectors of R^2). */
export function dualBasis(k1, k2) {
  const det = k1[0] * k2[1] - k1[1] * k2[0];
  const w1 = [k2[1] / det, -k2[0] / det];
  let w2 = [-k1[1] / det, k1[0] / det];
  if (w1[0] * w2[1] - w1[1] * w2[0] < 0) w2 = [-w2[0], -w2[1]];
  return [w1, w2];
}

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
export function genus1Edge(tau, { rStart = 0.9999, rMin = 0.003, factor = 0.97, maxEvents = 1 } = {}) {
  const t0 = 1 / Math.PI, names = ['tau', 'third'];
  const mk = markedVectors(tau);
  const lift = (p) => ({ lam: p.lam, nu: cmul(p.nu, [p.lam[0] + 1, p.lam[1]]) });
  const cur = {};
  for (const n of names) { const s = genus0Resonance(kOf(t0, mk[n]))[0]; cur[n] = { P: lift(s.P), Q: lift(s.Q) }; }
  const samples = [], events = [];
  let r = rStart, end = { r: rMin, why: 'rMin' };
  const Wint = (rr, t) => Math.abs(willmoreFrame(leafData([[-rr, 0]]), ...torusBasis(t, tau)));
  let first = true, jumped = false;
  let step = 1 - factor; // relative step in r, halved on failure and grown back after successes
  while (r > rMin) {
    let next, res, ok = false;
    for (let tries = 0; tries < 8 && !ok; tries++) {
      next = first ? r : Math.max(rMin, r * (1 - step));
      res = {};
      ok = true;
      for (const n of names) {
        res[n] = pairAt(tau, n, next, cur[n]);
        // reject a step on which a point jumps (Newton found another branch)
        const jump = (a, b) => dist(a.lam, b.lam) > 0.15 * Math.max(0.05, Math.hypot(...b.lam));
        if (!res[n].ok || (!first && (jump(res[n].P, cur[n].P) || jump(res[n].Q, cur[n].Q)))) ok = false;
      }
      if (!ok) step /= 2;
    }
    first = false;
    if (!ok) { end = { r: next, why: 'continuation failed' }; break; }
    step = Math.min(1 - factor, step * 1.5);
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
        const gap = dist(R.P.lam, R.Q.lam);
        // an arrival is a sigma-node: the pair must close up (otherwise the branch jumped)
        if (gap < 1e-2) {
          // the other pairs at the same point, for continuing them beyond the event
          const others = {};
          for (const o of names) if (o !== n) { const Ro = pairAt(tau, o, hi, cur[o]); others[o] = { P: Ro.P, Q: Ro.Q, ok: Ro.ok }; }
          events.push({ name: n, r: hi, t: R.t, W: Wint(hi, R.t), P: R.P, Q: R.Q, lamNode: [mid[0] / l, mid[1] / l], gap, others });
        }
        else { end = { r: hi, why: `the ${n} pair jumped branches (gap ${gap.toFixed(3)})` }; jumped = true; }
      }
    }
    if (jumped) break;
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
 * A candidate from two generalized Sym points on a genus g curve with internal torus lattice Lambda,
 * given by a basis [w1, w2] (torusBasis or dualBasis): syms = [{ lam, nu }, { lam, nu }] (a branch point
 * has nu = 0). Computes the monodromy at both by
 * shooting along the basis of Lambda, the sign characters, and the lattice on which f closes: Lambda if
 * the characters agree, otherwise the index 2 sublattice on which they do (F8).
 * Returns { lattice: [e1, e2], doubled, chars, us, W (residue formula on the closing lattice), tauExt,
 *   defects: { node (|M - (±I)| at sigma-nodes), closing (|f(z + e) - f(z)| for e = e1, e2) } }.
 */
export function candidate(alphas, basis, syms, { h = 0.002 } = {}) {
  const Z = killingField(alphas);
  const lams = syms.map((s) => s.lam);
  const I = new PairIntegrator(alphas, lams);
  const us = syms.map((s) => leftEigenvector(Z, s.lam, s.nu));
  const [w1, w2] = basis;
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
 * Candidate 1 for internal tau (any point of the upper half plane, read as F_6): the pair of the shortest
 * class is opened (internalFrame), then the genus 1 edge is followed to the first arrival (of the second
 * or third class), with the opened branch point alpha and the arriving sigma-node as generalized Sym points.
 * Returns { tau, frame, edge, event (with cls, the class label of the arriving pair), alphas, syms,
 * ...candidate } or { tau, frame, edge, event: null } if none was found.
 */
export function candidate1(tau, opts = {}) {
  // tau may be any point of the upper half plane (F_6 with labelled classes): normalize first
  const frame = internalFrame(tau), tau1 = frame.tau1;
  const edge = genus1Edge(tau1, opts);
  const event = edge.events[0] || null;
  if (!event) return { tau, frame, edge, event: null };
  event.cls = frame.classes[event.name];
  const alphas = [[-event.r, 0]];
  const syms = [{ lam: alphas[0], nu: [0, 0] }, symPoint(alphas, event.lamNode)];
  return { tau, frame, edge, event, alphas, syms, ...candidate(alphas, torusBasis(event.t, tau1), syms) };
}

// ---------------------------------------------------------------- genus 2 edge and candidate 2

/** k_j = conj(ell_j) for the cycles around {0, alpha_j} (ell of whitham.js), signs aligned to ref (ell vectors). */
function kVectors(alphas, ref = null) {
  const ell = ellVector(alphas, ref);
  return { ell, k: ell.map(([a, b]) => [a, -b]) };
}
const cdivc = ([a, b], [c, d]) => { const n = c * c + d * d; return [(a * c + b * d) / n, (b * c - a * d) / n]; };
const csqrtc = ([r, i]) => { const m = Math.hypot(r, i); const sr = Math.sqrt(Math.max(0, (m + r) / 2)), si = Math.sqrt(Math.max(0, (m - r) / 2)); return [sr, i < 0 ? -si : si]; };

/**
 * The anchor of the genus 2 edge: the node of candidate 1 opened by eps, alpha_1 = -r_1, alpha_2 =
 * (1 - eps) e^{i theta}, with (r_1, theta) solved by Newton so that the period ratio ell_2 / ell_1 is
 * exactly that of the internal lattice, conj(s tau1 + n) (s = ±1, n from the nearest candidate at the start:
 * the B-period through the node is the arriving vector modulo k_1).
 * Returns { alphas, ratio, target, iters }.
 */
export function genus2Anchor(r, lamNode, tau1, eps = 1e-3) {
  const th0 = Math.atan2(lamNode[1], lamNode[0]);
  const al = ([r1, th]) => [[-r1, 0], [(1 - eps) * Math.cos(th), (1 - eps) * Math.sin(th)]];
  const ratio = (x) => { const e = ellVector(al(x)); return cdivc(e[1], e[0]); };
  let x = [r, th0];
  const r0 = ratio(x);
  let target = null, bd = Infinity;
  for (const s of [1, -1]) for (let n = -3; n <= 3; n++) {
    const c = [s * tau1[0] + n, -s * tau1[1]];
    const d = Math.hypot(c[0] - r0[0], c[1] - r0[1]);
    if (d < bd) { bd = d; target = c; }
  }
  let iters = 0, F = [r0[0] - target[0], r0[1] - target[1]];
  for (; iters < 20 && Math.hypot(...F) > 1e-12; iters++) {
    const h = 1e-7, J = [[0, 0], [0, 0]];
    for (let j = 0; j < 2; j++) {
      const xp = x.slice(); xp[j] += h;
      const rp = ratio(xp);
      J[0][j] = (rp[0] - target[0] - F[0]) / h; J[1][j] = (rp[1] - target[1] - F[1]) / h;
    }
    const det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    x = [x[0] - (J[1][1] * F[0] - J[0][1] * F[1]) / det, x[1] - (J[0][0] * F[1] - J[1][0] * F[0]) / det];
    const rr = ratio(x);
    F = [rr[0] - target[0], rr[1] - target[1]];
  }
  return { alphas: al(x), ratio: [F[0] + target[0], F[1] + target[1]], target, iters, residual: Math.hypot(...F) };
}

/**
 * The genus 2 edge from candidate 1 (the result of candidate1): open its node into the anchor, follow the
 * Whitham curve (the period plane fixed up to scale: the internal lattice, scaled) in the direction in
 * which alpha_2 leaves the circle, and continue the pair of the remaining class from its genus 1 position
 * at the event (lifted by nu_2 = nu_1 sqrt(-conj lam_0) (lam - lam_0)) until it arrives at S^1.
 * Returns { anchor, samples: [{ s, alphas, W }], events: [{ s, alphas, k, W, lamNode, P, Q, gap }], end }.
 */
export function genus2Edge(c1, { eps = 1e-3, ds = 0.02, sMax = 6, maxEvents = 1 } = {}) {
  const tau1 = c1.frame.tau1, ev = c1.event;
  const otherName = ev.name === 'tau' ? 'third' : 'tau';
  const other = ev.others[otherName];
  const anchor = genus2Anchor(ev.r, ev.lamNode, tau1, eps);
  // a_{g+1} = a_g (lam - alpha)(1 - conj(alpha) lam) and (lam - lam0)(1 - conj(lam0) lam) = -conj(lam0) (lam - lam0)^2
  const lam0 = ev.lamNode, c = csqrtc([-lam0[0], lam0[1]]);
  const lift = (p) => ({ lam: p.lam, nu: cmul(cmul(p.nu, c), [p.lam[0] - lam0[0], p.lam[1] - lam0[1]]) });
  let cur = { P: lift(other.P), Q: lift(other.Q), cls: null };
  const curve = new WhithamCurve(anchor.alphas);
  // direction: alpha_2 moves into the disc
  // (very short steps from the anchor can fail in WhithamCurve.at, so probe at ±0.05)
  const rad = (s) => { const p = curve.at(s); return Math.abs(p.s - s) < 1e-9 ? Math.hypot(...p.alphas[1]) : 1; };
  const dir = rad(0.05) < rad(-0.05) ? 1 : -1;
  // a point of the curve at s, trying slightly further along when the corrector fails on a short step
  const pointAt = (sv) => {
    for (const extra of [0, 0.25, 0.5, 1]) {
      const target = sv + dir * extra * ds, p = curve.at(target);
      if (Math.abs(p.s - target) < 1e-12) return p;
    }
    return null;
  };
  let ref = ellVector(anchor.alphas);
  const pairAt2 = (alphas, start, refEll) => {
    const { ell, k } = kVectors(alphas, refEll);
    const Q = new QMap(alphas);
    const R = Q.resonance(nearestSheet(Q, start.P), nearestSheet(Q, start.Q), [[0, 0], [0, 0]], { lattice: k });
    return { ...R, ell, k, cls: classifyPair(R.P, R.Q) };
  };
  const Wat = (alphas, k) => Math.abs(willmoreFrame(leafData(alphas), ...dualBasis(k[0], k[1])));
  const samples = [], events = [];
  let s = 0, step = ds, end = { s: dir * sMax, why: 'sMax' };
  // the starting point
  {
    const R = pairAt2(anchor.alphas, cur, ref);
    if (!R.ok) return { anchor, samples, events, end: { s: 0, why: 'the pair could not be continued into genus 2' } };
    cur = { P: R.P, Q: R.Q, cls: R.cls };
    ref = R.ell;
    samples.push({ s: 0, alphas: anchor.alphas, W: Wat(anchor.alphas, R.k), pair: R.cls });
  }
  while (Math.abs(s) < sMax) {
    let R, p, next, ok = false;
    for (let tries = 0; tries < 8 && !ok; tries++) {
      next = s + dir * step;
      p = pointAt(next);
      if (!p) { end = { s: next, why: `the Whitham curve ends (${curve.fail || 'corrector'})` }; break; }
      next = p.s;
      R = pairAt2(p.alphas, cur, ref);
      const jump = (a, b) => dist(a.lam, b.lam) > 0.15 * Math.max(0.05, Math.hypot(...b.lam));
      ok = R.ok && !jump(R.P, cur.P) && !jump(R.Q, cur.Q);
      if (!ok) step /= 2;
    }
    if (!ok) { if (end.why === 'sMax') end = { s: next, why: 'continuation failed' }; break; }
    step = Math.min(ds, step * 1.5);
    if (cur.cls && cur.cls.kind === 'off' && R.cls.kind !== 'off') {
      let hi = s, lo = next, sHi = cur, refHi = ref;
      for (let i = 0; i < 50 && Math.abs(hi - lo) > 1e-12; i++) {
        const m = (hi + lo) / 2, pm = curve.at(m);
        if (Math.abs(pm.s - m) > 1e-12) break;
        const Rm = pairAt2(pm.alphas, sHi, refHi);
        if (Rm.ok && Rm.cls.kind === 'off') { hi = m; sHi = { P: Rm.P, Q: Rm.Q }; refHi = Rm.ell; } else lo = m;
      }
      const pe = curve.at(hi), Re = pairAt2(pe.alphas, sHi, refHi);
      const mid = [(Re.P.lam[0] + Re.Q.lam[0]) / 2, (Re.P.lam[1] + Re.Q.lam[1]) / 2], l = Math.hypot(...mid);
      const gap = dist(Re.P.lam, Re.Q.lam);
      if (gap < 1e-2) {
        events.push({ s: hi, alphas: pe.alphas, k: Re.k, W: Wat(pe.alphas, Re.k), lamNode: [mid[0] / l, mid[1] / l], P: Re.P, Q: Re.Q, gap, name: otherName });
        samples.push({ s: hi, alphas: pe.alphas, W: Wat(pe.alphas, Re.k), pair: Re.cls });
      } else end = { s: hi, why: `the pair jumped branches (gap ${gap.toFixed(3)})` };
      if (events.length >= maxEvents || gap >= 1e-2) { if (events.length >= maxEvents) end = { s: hi, why: 'event' }; break; }
    }
    s = next;
    cur = { P: R.P, Q: R.Q, cls: R.cls };
    ref = R.ell;
    samples.push({ s, alphas: p.alphas, W: Wat(p.alphas, R.k), pair: R.cls });
  }
  return { anchor, samples, events, end };
}

/**
 * Candidate 2 for internal tau: candidate 1, then the genus 2 edge to the arrival of the remaining pair.
 * At the arrival there are three trivial-multiplier rho-pairs: the branch points alpha_1, alpha_2 (from the
 * two openings) and the new node; each choice of two gives a torus (choices: [{ pair, ...candidate }]).
 * Returns { tau, c1, edge, event, alphas, choices } (event null if none was found).
 */
export function candidate2(tau, opts = {}) {
  const c1 = candidate1(tau, opts.genus1 || {});
  if (!c1.event) return { tau, c1, edge: null, event: null };
  const edge = genus2Edge(c1, opts.genus2 || {});
  const event = edge.events[0] || null;
  if (!event) return { tau, c1, edge, event: null };
  const al = event.alphas;
  const pts = { a1: { lam: al[0], nu: [0, 0] }, a2: { lam: al[1], nu: [0, 0] }, node: symPoint(al, event.lamNode) };
  const basis = dualBasis(event.k[0], event.k[1]);
  const choices = [['a1', 'a2'], ['a1', 'node'], ['a2', 'node']].map(([x, y]) => ({ pair: [x, y], syms: [pts[x], pts[y]], ...candidate(al, basis, [pts[x], pts[y]]) }));
  return { tau, frame: c1.frame, c1, edge, event, alphas: al, choices };
}
