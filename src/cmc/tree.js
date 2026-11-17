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
import { PairIntegrator, leftEigenvector, nuAt, computeGenSurface, roundnessS4 } from './gensym.js';
import { leafData, willmoreFrame } from './leaf.js';
import { WhithamCurve, ellVector } from './whitham.js';
import { openEndTau } from './genus2.js';

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

/**
 * Refine an arrival (a passage of a pair through a sigma-node) between hi, where the pair is off S^1, and
 * lo, where it is not. Near the node the pair moves like a square root of the edge parameter on either
 * side, so g = -dist^2 (off S^1: dist = log|lam_P|) or +dist^2 (on S^1: dist = the angle between P and Q)
 * is close to linear through the arrival: regula falsi with the Illinois modification, keeping the bracket
 * (bisection if a solve fails or the step leaves the bracket). evalAt(m, s) solves at m continuing from
 * the state s of the off side and returns { ok, off, dist, s }. gHi, gLo: g at hi and lo. Stops when the
 * pair is within 1e-7 of S^1 on the off side (classifyPair counts a pair within 1e-8 of S^1 as on it, so
 * the off-side distance cannot fall much below that; the bisection before this ran on to a 1e-13
 * bracket, 50 solves, without getting closer) or the bracket is below 1e-13.
 * Returns { hi, s, evals }: the off-side end and its state.
 */
function refineArrival(hi, lo, gHi, gLo, s, evalAt) {
  let side = 0, bisect = false, evals = 0;
  for (let it = 0; it < 60; it++) {
    if (Math.abs(hi - lo) <= 1e-13 || Math.sqrt(-gHi) < 1e-7) break;
    let m = bisect || !(gLo > gHi) ? (hi + lo) / 2 : (hi * gLo - lo * gHi) / (gLo - gHi);
    const a = Math.min(hi, lo), b = Math.max(hi, lo), margin = 1e-3 * (b - a);
    if (!(m > a + margin && m < b - margin)) m = Math.min(b - margin, Math.max(a + margin, m));
    bisect = false;
    const R = evalAt(m, s);
    evals++;
    if (!R.ok) { lo = m; bisect = true; continue; }
    const g = R.off ? -R.dist * R.dist : R.dist * R.dist;
    if (R.off) { hi = m; gHi = g; s = R.s; if (side === -1) gLo /= 2; side = -1; } else { lo = m; gLo = g; if (side === 1) gHi /= 2; side = 1; }
  }
  return { hi, s, evals };
}
const gOf = (cls) => (cls.kind === 'off' ? -cls.dist * cls.dist : cls.dist * cls.dist);

/** Solve the pair named name at r, continuing from start = { P, Q }. */
let g1Cache = null;
/** The curve data at r shared by the pairs (the q map and the scale), cached for the last r. */
function genus1At(r) {
  if (!g1Cache || g1Cache.r !== r) g1Cache = { r, Q: new QMap([[-r, 0]]), t: genus1Scale(r) };
  return g1Cache;
}
/**
 * A starting guess at x for a pair continued along an edge: linear extrapolation of lam and nu in the edge
 * parameter u (log r or log t) from the last two solutions h = [{ u, P, Q }, { u, P, Q }] (or the last one).
 */
function extrapolate(h, u) {
  if (h.length < 2) return { P: h[h.length - 1].P, Q: h[h.length - 1].Q };
  const [a, b] = h.slice(-2), s = (u - b.u) / (b.u - a.u);
  const ex = (p, q) => ({ lam: [q.lam[0] + s * (q.lam[0] - p.lam[0]), q.lam[1] + s * (q.lam[1] - p.lam[1])], nu: [q.nu[0] + s * (q.nu[0] - p.nu[0]), q.nu[1] + s * (q.nu[1] - p.nu[1])] });
  return { P: ex(a.P, b.P), Q: ex(a.Q, b.Q) };
}

function pairAt(tau, name, r, start) {
  const { Q, t } = genus1At(r);
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
export function genus1Edge(tau, { rStart = 0.9999, rMin = 0.003, factor = 0.8, maxEvents = 1, predict = true } = {}) {
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
  const hist = Object.fromEntries(names.map((n) => [n, []])); // the last solutions, for the predictor
  while (r > rMin) {
    let next, res, ok = false;
    for (let tries = 0; tries < 8 && !ok; tries++) {
      next = first ? r : Math.max(rMin, r * (1 - step));
      res = {};
      ok = true;
      for (const n of names) {
        res[n] = pairAt(tau, n, next, predict && hist[n].length ? extrapolate(hist[n], Math.log(next)) : cur[n]);
        // reject a step on which a point jumps (Newton found another branch)
        const jump = (a, b) => dist(a.lam, b.lam) > 0.15 * Math.max(0.05, Math.hypot(...b.lam));
        if (!res[n].ok || (!first && (jump(res[n].P, cur[n].P) || jump(res[n].Q, cur[n].Q)))) ok = false;
      }
      if (!ok) step /= 2;
    }
    first = false;
    if (!ok) { end = { r: next, why: 'continuation failed' }; break; }
    step = Math.min(1 - factor, step * 1.5);
    // an arrival between r and next: the pair passes through a sigma-node, from off S^1 onto it (or from
    // on S^1 off it); refined by bisection, continuing from the side where the pair is off S^1
    for (const n of names) {
      if (r !== next && cur[n].cls && (cur[n].cls.kind === 'off') !== (res[n].cls.kind === 'off')) {
        const fromOff = cur[n].cls.kind === 'off';
        const offSide = fromOff ? cur[n] : res[n], other = fromOff ? res[n] : cur[n];
        const ref = refineArrival(fromOff ? r : next, fromOff ? next : r, gOf(offSide.cls), gOf(other.cls), offSide, (m, st) => {
          const Rm = pairAt(tau, n, m, st);
          return { ok: Rm.ok, off: Rm.ok && Rm.cls.kind === 'off', dist: Rm.ok ? Rm.cls.dist : 0, s: { P: Rm.P, Q: Rm.Q, cls: Rm.cls } };
        });
        const hi = ref.hi, sHi = ref.s;
        const R = pairAt(tau, n, hi, sHi);
        const mid = [(R.P.lam[0] + R.Q.lam[0]) / 2, (R.P.lam[1] + R.Q.lam[1]) / 2], l = Math.hypot(...mid);
        const gap = dist(R.P.lam, R.Q.lam);
        // an arrival is a sigma-node: the pair must close up (otherwise the branch jumped)
        if (gap < 1e-2) {
          // the other pairs at the same point, for continuing them beyond the event
          const others = {};
          for (const o of names) if (o !== n) { const Ro = pairAt(tau, o, hi, cur[o]); others[o] = { P: Ro.P, Q: Ro.Q, ok: Ro.ok }; }
          // a double arrival (on the walls of F_6 where the other two classes are equally long): another
          // pair is a sigma-node at the same point of the edge
          let double = null;
          for (const [o, p] of Object.entries(others)) {
            const c = classifyPair(p.P, p.Q, 1e-6), g = dist(p.P.lam, p.Q.lam);
            if (p.ok && g < 1e-3 && (c.kind !== 'off' || Math.abs(c.dist) < 1e-6)) {
              const mo = [(p.P.lam[0] + p.Q.lam[0]) / 2, (p.P.lam[1] + p.Q.lam[1]) / 2], lo2 = Math.hypot(...mo);
              double = { name: o, lamNode: [mo[0] / lo2, mo[1] / lo2] };
            }
          }
          events.push({ name: n, r: hi, t: R.t, W: Wint(hi, R.t), P: R.P, Q: R.Q, lamNode: [mid[0] / l, mid[1] / l], gap, others, double, refineEvals: ref.evals });
        }
        else { end = { r: hi, why: `the ${n} pair jumped branches (gap ${gap.toFixed(3)})` }; jumped = true; }
      }
    }
    if (jumped) break;
    for (const n of names) {
      cur[n] = { P: res[n].P, Q: res[n].Q, cls: res[n].cls };
      hist[n].push({ u: Math.log(next), P: res[n].P, Q: res[n].Q });
      if (hist[n].length > 2) hist[n].shift();
    }
    r = next;
    const t = res.tau.t;
    samples.push({ r, t, W: Wint(r, t), pairs: Object.fromEntries(names.map((n) => [n, res[n].cls])),
      pts: Object.fromEntries(names.map((n) => [n, { P: res[n].P, Q: res[n].Q }])) });
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
export function candidate(alphas, basis, syms, { h = Infinity } = {}) {
  const Z = killingField(alphas);
  const lams = syms.map((s) => s.lam);
  const I = new PairIntegrator(alphas, lams);
  // a quarter of the step cap: the characters and W do not depend on it, and the closing defect (RK4,
  // ~h^4) is ~1e-6 (with h = 0.002 it was ~1e-9, at ten times the cost)
  h = Math.min(h, I.stepCap(Z) / 4);
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
  // the edge is followed to the second event, for the branch that does not open the node
  const edge = genus1Edge(tau1, { maxEvents: 2, ...opts });
  const event = edge.events[0] || null;
  if (!event) return { tau, frame, edge, event: null, branch: null };
  event.cls = frame.classes[event.name];
  const alphas = [[-event.r, 0]];
  const syms = [{ lam: alphas[0], nu: [0, 0] }, symPoint(alphas, event.lamNode)];
  return { tau, frame, edge, event, alphas, syms, ...candidate(alphas, torusBasis(event.t, tau1), syms), branch: genus1Branch(edge, event, frame) };
}

/**
 * The "(or not)" branch at candidate 1: the genus 1 edge continued past the arrival without opening the
 * node. The arrived pair stays over S^1 (two points), and the branch ends at the next sigma-node event of a
 * marked pair: the other pair arriving, or the arrived pair leaving S^1 again. There the genus 1 curve
 * has two ripe rho-pairs, alpha and the node: candidate 1' (a provisional name).
 * Returns null at a double arrival or when there is no first event, else { event (or null: the edge
 * ended without one, see end), alphas, syms, ...candidate }.
 */
function genus1Branch(edge, event, frame) {
  if (!event || event.double) return null;
  const ev = edge.events[1];
  if (!ev || ev.r >= event.r) return { event: null, end: edge.end };
  ev.cls = frame.classes[ev.name];
  ev.leaving = ev.name === event.name; // the arrived pair leaves S^1 again (else the other pair arrives)
  const alphas = [[-ev.r, 0]];
  const syms = [{ lam: alphas[0], nu: [0, 0] }, symPoint(alphas, ev.lamNode)];
  return { event: ev, alphas, syms, ...candidate(alphas, torusBasis(ev.t, frame.tau1), syms) };
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
  // ellVector gives each period vector up to sign (cycle orientations): align to the first evaluation, so
  // that the sign cannot flip between the evaluations of the Newton iteration
  let ref = null;
  const ratio = (x) => { const e = ellVector(al(x), ref); if (!ref) ref = e; return cdivc(e[1], e[0]); };
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
 * The genus 2 edge from candidate 1 (the result of candidate1), parametrized by the scale of the branch
 * points. Gauge (as WhithamCurve: arg(alpha_1 alpha_2) fixed at the anchor, 2 psi):
 *   alpha_1 = t beta e^{i psi},   alpha_2 = (t / beta) e^{i psi},   t = sqrt|alpha_1 alpha_2|,
 * and the edge is the set ratio(t, beta) = the internal period ratio, solved by Newton in beta (2 real
 * equations, 2 real unknowns) and continued in sigma = log t from the anchor (the node of candidate 1,
 * opened by eps) down to tMin (default 1e-3: the open-end quantities converge with error O(t)). t decreases monotonically along the edge, and beta stays of size 1 up to
 * the open end, where the rescaled roots are check-alpha = beta, 1/beta (genus 2 paper, "The Open Ends").
 * The remaining pair is lifted at the node (nu_2 = nu_1 sqrt(-conj lam_0) (lam - lam_0)) and continued;
 * its arrival at S^1 is refined by bisection in sigma.
 * Returns { anchor, samples: [{ t, s (= -log(t / t_anchor)), alphas, W, pair }], events: [{ t, s, alphas,
 *   k, W, lamNode, P, Q, gap, name }], end, openEnd (when tMin is reached without an arrival) }.
 */
export function genus2Edge(c1, { eps = 1e-3, tMin = 1e-3, dSigma = 0.1, maxEvents = 1 } = {}) {
  const tau1 = c1.frame.tau1, ev = c1.event;
  const otherName = ev.name === 'tau' ? 'third' : 'tau';
  const other = ev.others[otherName];
  const anchor = genus2Anchor(ev.r, ev.lamNode, tau1, eps);
  // the target ratio in the sign convention of the reference used below (the period vectors are only
  // defined up to sign, so the ratio is defined up to sign)
  const ref0 = ellVector(anchor.alphas), r0 = cdivc(ref0[1], ref0[0]);
  const target = dist(r0, anchor.target) <= dist(r0, anchor.target.map((x) => -x)) ? anchor.target : anchor.target.map((x) => -x);
  // a_{g+1} = a_g (lam - alpha)(1 - conj(alpha) lam) and (lam - lam0)(1 - conj(lam0) lam) = -conj(lam0) (lam - lam0)^2
  const lam0 = ev.lamNode, c = csqrtc([-lam0[0], lam0[1]]);
  const lift = (p) => ({ lam: p.lam, nu: cmul(cmul(p.nu, c), [p.lam[0] - lam0[0], p.lam[1] - lam0[1]]) });
  const [A1, A2] = anchor.alphas, prod = cmul(A1, A2);
  const psi = Math.atan2(prod[1], prod[0]) / 2, ePsi = [Math.cos(psi), Math.sin(psi)];
  const tA = Math.sqrt(Math.hypot(...prod));
  const alphasOf = (t, beta) => [cmul(beta, ePsi).map((x) => t * x), cmul(cdivc([1, 0], beta), ePsi).map((x) => t * x)];
  let beta = cmul(A1, [ePsi[0], -ePsi[1]]).map((x) => x / tA);
  let ref = ref0;
  /** beta on the edge at scale t, by Newton from a guess. Returns { beta, ell, ok }. */
  const solveBeta = (t, guess, refEll) => {
    let b = guess.slice(), r = refEll;
    const F = (bb) => { const e = ellVector(alphasOf(t, bb), r); const q = cdivc(e[1], e[0]); return { f: [q[0] - target[0], q[1] - target[1]], e }; };
    let cur = F(b);
    for (let it = 0; it < 20; it++) {
      if (Math.hypot(...cur.f) < 1e-12) return { beta: b, ell: cur.e, ok: true };
      const h = 1e-7 * Math.max(1, Math.hypot(...b)), J = [[0, 0], [0, 0]];
      for (let j = 0; j < 2; j++) {
        const bp = b.slice(); bp[j] += h;
        const fp = F(bp).f;
        J[0][j] = (fp[0] - cur.f[0]) / h; J[1][j] = (fp[1] - cur.f[1]) / h;
      }
      const det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
      if (!(Math.abs(det) > 0)) break;
      const db = [-(J[1][1] * cur.f[0] - J[0][1] * cur.f[1]) / det, -(J[0][0] * cur.f[1] - J[1][0] * cur.f[0]) / det];
      // damping: beta moves by at most a quarter of its size per iteration
      const sc = Math.min(1, (0.25 * Math.hypot(...b)) / Math.max(1e-300, Math.hypot(...db)));
      b = [b[0] + sc * db[0], b[1] + sc * db[1]];
      cur = F(b); r = cur.e;
    }
    return { beta: b, ell: cur.e, ok: Math.hypot(...cur.f) < 1e-9 };
  };
  const pairAt2 = (alphas, start, refEll) => {
    const { ell, k } = kVectors(alphas, refEll);
    const Q = new QMap(alphas);
    const R = Q.resonance(nearestSheet(Q, start.P), nearestSheet(Q, start.Q), [[0, 0], [0, 0]], { lattice: k });
    return { ...R, ell, k, cls: classifyPair(R.P, R.Q) };
  };
  const Wat = (alphas, k) => Math.abs(willmoreFrame(leafData(alphas), ...dualBasis(k[0], k[1])));
  const jump = (a, b) => dist(a.lam, b.lam) > 0.15 * Math.max(0.05, Math.hypot(...b.lam));
  const samples = [], events = [];
  let end = { t: tMin, why: 'tMin' };
  let cur = { P: lift(other.P), Q: lift(other.Q), cls: null };
  {
    const R = pairAt2(anchor.alphas, cur, ref);
    if (!R.ok) return { anchor, samples, events, end: { t: tA, why: 'the pair could not be continued into genus 2' }, openEnd: null };
    cur = { P: R.P, Q: R.Q, cls: R.cls };
    ref = R.ell;
    samples.push({ t: tA, s: 0, alphas: anchor.alphas, W: Wat(anchor.alphas, R.k), pair: R.cls, k: R.k, pts: { P: R.P, Q: R.Q } });
  }
  let sigma = Math.log(tA), step = dSigma, betaPrev = null, sigmaPrev = null;
  const sigmaMin = Math.log(tMin);
  while (sigma > sigmaMin + 1e-12) {
    let next, sol, R, ok = false;
    for (let tries = 0; tries < 10 && !ok; tries++) {
      next = Math.max(sigmaMin, sigma - step);
      // predictor: linear in sigma
      const guess = betaPrev ? beta.map((x, i) => x + ((x - betaPrev[i]) * (next - sigma)) / (sigma - sigmaPrev)) : beta;
      sol = solveBeta(Math.exp(next), guess, ref);
      ok = sol.ok && dist(sol.beta, beta) < 0.2 * Math.hypot(...beta);
      if (ok) {
        R = pairAt2(alphasOf(Math.exp(next), sol.beta), cur, sol.ell);
        ok = R.ok && !jump(R.P, cur.P) && !jump(R.Q, cur.Q);
      }
      if (!ok) step /= 2;
    }
    if (!ok) { end = { t: Math.exp(next), why: 'continuation failed' }; break; }
    step = Math.min(dSigma, step * 1.5);
    if (cur.cls && (cur.cls.kind === 'off') !== (R.cls.kind === 'off')) {
      // arrival (a passage through a sigma-node, onto S^1 or off it): bisection in sigma, continuing from the
      // side where the pair is off S^1
      const fromOff = cur.cls.kind === 'off';
      const s0 = fromOff ? { P: cur.P, Q: cur.Q, beta, ref } : { P: R.P, Q: R.Q, beta: sol.beta, ref: R.ell };
      const rf = refineArrival(fromOff ? sigma : next, fromOff ? next : sigma, gOf(fromOff ? cur.cls : R.cls), gOf(fromOff ? R.cls : cur.cls), s0, (m, st) => {
        const sm = solveBeta(Math.exp(m), st.beta, st.ref);
        if (!sm.ok) return { ok: false };
        const Rm = pairAt2(alphasOf(Math.exp(m), sm.beta), st, sm.ell);
        return { ok: Rm.ok, off: Rm.ok && Rm.cls.kind === 'off', dist: Rm.ok ? Rm.cls.dist : 0, s: { P: Rm.P, Q: Rm.Q, beta: sm.beta, ref: sm.ell } };
      });
      const hi = rf.hi, sHi = { P: rf.s.P, Q: rf.s.Q }, bHi = rf.s.beta, refHi = rf.s.ref;
      const tE = Math.exp(hi), alE = alphasOf(tE, bHi), Re = pairAt2(alE, sHi, refHi);
      const mid = [(Re.P.lam[0] + Re.Q.lam[0]) / 2, (Re.P.lam[1] + Re.Q.lam[1]) / 2], l = Math.hypot(...mid);
      const gap = dist(Re.P.lam, Re.Q.lam);
      if (gap < 1e-2) {
        const W = Wat(alE, Re.k);
        events.push({ t: tE, s: Math.log(tA / tE), alphas: alE, k: Re.k, W, lamNode: [mid[0] / l, mid[1] / l], P: Re.P, Q: Re.Q, gap, name: otherName, refineEvals: rf.evals });
        samples.push({ t: tE, s: Math.log(tA / tE), alphas: alE, W, pair: Re.cls, k: Re.k, pts: { P: Re.P, Q: Re.Q } });
        if (events.length >= maxEvents) { end = { t: tE, why: 'event' }; break; }
      } else { end = { t: tE, why: `the pair jumped branches (gap ${gap.toFixed(3)})` }; break; }
    }
    betaPrev = beta; sigmaPrev = sigma;
    beta = sol.beta; sigma = next; ref = R.ell;
    cur = { P: R.P, Q: R.Q, cls: R.cls };
    const al = alphasOf(Math.exp(sigma), beta);
    samples.push({ t: Math.exp(sigma), s: Math.log(tA / Math.exp(sigma)), alphas: al, W: Wat(al, R.k), pair: R.cls, k: R.k, pts: { P: R.P, Q: R.Q } });
  }
  let openEnd = null;
  if (!events.length && end.why === 'tMin') {
    // the open end: the rescaled roots are beta, 1/beta (product 1 in this gauge)
    const tau = reduceTau(openEndTau(beta)), refTau = reduceTau(tau1);
    const tauError = Math.min(Math.hypot(tau[0] - refTau[0], tau[1] - refTau[1]), Math.hypot(tau[0] + refTau[0], tau[1] - refTau[1]));
    openEnd = { t: Math.exp(sigma), alpha: beta, tau, tauError, WRatio: samples[samples.length - 1].W / (4 * Math.PI) };
  } else if (!events.length) openEnd = openEndInfo(samples, tau1);
  return { anchor, samples, events, end, openEnd };
}

/**
 * Where the genus 2 edge ends without an arrival: is it the open end (both branch points to 0, genus 2
 * paper, section "The open ends")? In the rescaled coordinate lam = t check-lam, t = sqrt|alpha_1 alpha_2|,
 * the roots tend to check-alpha, 1/check-alpha, the curve to the elliptic curve check-nu^2 = check-lam
 * (check-lam - check-alpha)(check-lam - 1/check-alpha), and W to 4 pi (8 pi on a doubled lattice); the limit
 * of the candidates is a branched double cover of the sphere. The open-end tau of check-alpha (genus2.js)
 * must tend to the internal tau. Returns null unless the last sample has t < 0.15, else { t, alpha
 * (check-alpha_1), tau (open-end tau, reduced into F_1), tauError (against the reduced internal tau), WRatio
 * (W_int / 4 pi at the last sample) }.
 */
export function openEndInfo(samples, tau1) {
  if (!samples.length) return null;
  const sm = samples[samples.length - 1], [a1, a2] = sm.alphas, p = cmul(a1, a2);
  const t = Math.sqrt(Math.hypot(...p));
  if (!(t < 0.15)) return null;
  const ph = Math.atan2(p[1], p[0]) / 2, rot = [Math.cos(-ph), Math.sin(-ph)];
  const alpha = cmul(a1, rot).map((x) => x / t);
  const tau = reduceTau(openEndTau(alpha)), ref = reduceTau(tau1);
  // tau and its mirror image -conj(tau) describe the same conformal class up to orientation
  const tauError = Math.min(Math.hypot(tau[0] - ref[0], tau[1] - ref[1]), Math.hypot(tau[0] + ref[0], tau[1] - ref[1]));
  return { t, alpha, tau, tauError, WRatio: sm.W / (4 * Math.PI) };
}

/**
 * Candidate 2 for internal tau: candidate 1, then the genus 2 edge to the arrival of the remaining pair.
 * At the arrival there are three trivial-multiplier rho-pairs: the branch points alpha_1, alpha_2 (from the
 * two openings) and the new node; each choice of two gives a torus (choices: [{ pair, ...candidate }]).
 * At a double arrival in genus 1 (c1.event.double) there is no genus 2 edge: the choices are the three
 * pairs of { alpha, node a, node b } on the genus 1 curve, with double: true and edge null.
 * Returns { tau, c1, edge, event, alphas, choices } (event null if none was found). If the edge reaches the
 * open end instead, also openEnd, and the choices are candidate 2 at the scales opts.drawT (default 0.03,
 * 0.01, 0.003, 0.001; with opts.fast only 0.03), each with t and alphas; with opts.shapes also shape
 * (openEndShapes, which is otherwise left to be requested separately).
 */
export function candidate2(tau, opts = {}) {
  const c1 = candidate1(tau, opts.genus1 || {});
  if (!c1.event) return { tau, c1, edge: null, event: null };
  if (c1.event.double) {
    // a double arrival (on a wall of F_6): the genus 2 edge has length 0. All three candidates lie on the
    // genus 1 curve of candidate 1, with three ripe rho-pairs: the branch point alpha and the two nodes.
    // In the limit from either side alpha_2 -> node a, so candidate 2 (alpha_1 + alpha_2) is candidate 1
    // (alpha + node a), candidate 3 is alpha + node b, and the other choice is node a + node b.
    const al = c1.alphas, ev = c1.event;
    const pts = { a1: { lam: al[0], nu: [0, 0] }, a2: symPoint(al, ev.lamNode), node: symPoint(al, ev.double.lamNode) };
    const basis = torusBasis(ev.t, c1.frame.tau1);
    const choices = [['a1', 'a2'], ['a1', 'node'], ['a2', 'node']].map(([x, y]) => ({ pair: [x, y], syms: [pts[x], pts[y]], ...candidate(al, basis, [pts[x], pts[y]]) }));
    const event = { t: ev.t, s: 0, alphas: al, W: c1.W, lamNode: ev.double.lamNode, gap: 0, name: ev.double.name, double: true };
    return { tau, frame: c1.frame, c1, edge: null, event, alphas: al, choices, double: true };
  }
  const edge = genus2Edge(c1, opts.genus2 || {});
  const event = edge.events[0] || null;
  if (!event) {
    // no arrival: near the open end, candidate 2 (the two pairs of branch points as Sym points) at the
    // scale t approximates the limit, a branched double cover of the round sphere (roundnessS4: the
    // distance from a round sphere is O(t^(1/2)), and the S^4 area tends to 2 (x 4 pi))
    if (!edge.openEnd) return { tau, frame: c1.frame, c1, edge, event: null };
    const choices = [];
    for (const tDraw of opts.drawT || (opts.fast ? [0.03] : [0.03, 0.01, 0.003, 0.001])) {
      const sm = edge.samples.reduce((b, x) => (Math.abs(Math.log(x.t / tDraw)) < Math.abs(Math.log(b.t / tDraw)) ? x : b));
      if (Math.abs(Math.log(sm.t / tDraw)) > 0.5 || choices.some((ch) => ch.t === sm.t)) continue;
      const al = sm.alphas, k = ellVector(al).map(([x, y]) => [x, -y]);
      const syms = [{ lam: al[0], nu: [0, 0] }, { lam: al[1], nu: [0, 0] }], basis = dualBasis(k[0], k[1]);
      let ch;
      if (!choices.length) ch = { pair: ['a1', 'a2'], t: sm.t, alphas: al, syms, ...candidate(al, basis, syms) };
      else {
        // the characters of the branch points are constant along the edge: the closing lattice has the
        // parity of the first scale's, and W is the residue formula there (no frame integration, which
        // costs ~1/t near the open end); the closing defects are not computed (null)
        const c0 = choices[0], par = c0.chars.map(([u, v]) => (u === v ? 0 : 1));
        const lattice = paritySublattice(basis, par), [e1, e2] = lattice, den = e1[0] ** 2 + e1[1] ** 2;
        ch = { pair: ['a1', 'a2'], t: sm.t, alphas: al, syms, lattice, doubled: c0.doubled, chars: c0.chars,
          W: Math.abs(willmoreFrame(leafData(al), e1, e2)),
          tauExt: reduceTau([(e2[0] * e1[0] + e2[1] * e1[1]) / den, (e2[1] * e1[0] - e2[0] * e1[1]) / den]), defects: null };
      }
      choices.push(ch);
    }
    // the shape diagnostics are expensive near the open end (seconds): only on request
    if (opts.shapes) openEndShapes(choices, opts);
    return { tau, frame: c1.frame, c1, edge, event: null, openEnd: edge.openEnd, alphas: choices[0].alphas, choices };
  }
  const al = event.alphas;
  const pts = { a1: { lam: al[0], nu: [0, 0] }, a2: { lam: al[1], nu: [0, 0] }, node: symPoint(al, event.lamNode) };
  const basis = dualBasis(event.k[0], event.k[1]);
  const choices = [['a1', 'a2'], ['a1', 'node'], ['a2', 'node']].map(([x, y]) => ({ pair: [x, y], syms: [pts[x], pts[y]], ...candidate(al, basis, [pts[x], pts[y]]) }));
  return { tau, frame: c1.frame, c1, edge, event, alphas: al, choices };
}

// ---------------------------------------------------------------- the tree as one path of drawable samples

/**
 * The sublattice { a w1 + b w2 : a s1 + b s2 even } of the lattice with basis [w1, w2], for a parity
 * s in (Z/2)^2: the lattice itself for s = 0, else of index 2 (positively oriented when the basis is).
 * Two Sym points of classes g1, g2 (as (p, q) mod 2 in the basis dual to [w1, w2]) close on the
 * sublattice of parity g1 + g2 (F8); this holds for their sign characters before they are ripe too.
 */
export function paritySublattice([w1, w2], [s1, s2]) {
  const add = (a, b) => [a[0] + b[0], a[1] + b[1]], dbl = (a) => [2 * a[0], 2 * a[1]];
  if (!(s1 & 1) && !(s2 & 1)) return [w1, w2];
  if (!(s2 & 1)) return [dbl(w1), w2];
  if (!(s1 & 1)) return [w1, dbl(w2)];
  return [dbl(w1), add(w1, w2)];
}

// classes of the marked vectors of tau1, as (p, q) mod 2 in the basis (1, tau1)
const CLASS_PQ = { one: [1, 0], tau: [0, 1], third: [1, 1] };

/** The Sym point at the point p = { lam, nu } of the curve of alphas: { lam, sign } with nu = sign sqrt(lam a(lam)). */
function symOf(alphas, p) {
  const n = nuAt(alphas, p.lam, 1);
  return { lam: p.lam.slice(), sign: n[0] * p.nu[0] + n[1] * p.nu[1] >= 0 ? 1 : -1 };
}

/**
 * The tree for one internal tau as a path of samples, from near the vacuum to the last candidate, with the
 * immersion drawable at every sample (before the Sym points are ripe it is a plane, equivariant under
 * the lattice: f(z + w) = mu_1(w) f(z) mu_2(w)^{-1}, and the multipliers are real, so these are homotheties
 * of R^4 about 0; the cell drawn is the one on which the eventual candidate closes). A pair off S^1 is
 * (P, rho sigma P); the Sym point is P (Q gives the other eigenline over lam_P, and both agree once ripe).
 * C is the result of candidate2 (or of candidate1).
 * The edges (genus 0, 1, 2) occupy x in [0, 1], [1, 2], [2, 3]: genus 0 by log t from t = 10 / pi to the
 * opening at t = 1 / pi, genus 1 by log r to the first arrival, genus 2 by log t to the arrival or tMin.
 * W is the residue formula on the closing (doubled) lattice: twice W on the internal lattice.
 * Each sample: { edge, x, W, t, alphas, views: [{ key, syms: [{ lam, sign }, ...], lattice, ripe }],
 *   node (a label at the candidates and the opening) }.
 * Views on the genus 2 edge, in the order of candidate2's choices: [a1, a2] (candidate 2: ripe along the
 * whole edge), [a1, b] (towards candidate 3), [a2, b], with b the remaining pair; elsewhere one view.
 * Returns { frame, samples, nodes: [{ x, W, label, index }], end }.
 */
export function treeOf(C, { n0 = 24 } = {}) {
  const c1 = C.c1 || C, frame = c1.frame, tau1 = frame.tau1;
  const edge1 = c1.edge, ev1 = c1.event;
  const a = ev1 ? ev1.name : 'tau', b = a === 'tau' ? 'third' : 'tau';
  const samples = [], nodes = [];
  const sOne = CLASS_PQ.one, sA = CLASS_PQ[a];
  const par = [(sOne[0] + sA[0]) & 1, (sOne[1] + sA[1]) & 1];
  // genus 0: the pair for 1 and the pair a, at scale t (the curve nu^2 = -lam is fixed, the lattice scales)
  const tTop = 10 / Math.PI, tBot = 1 / Math.PI;
  const mk = markedVectors(tau1);
  for (let i = 0; i <= n0; i++) {
    const t = tTop * (tBot / tTop) ** (i / n0);
    const p1 = genus0Resonance(kOf(t, mk.one))[0].P, pa = genus0Resonance(kOf(t, mk[a]))[0].P;
    const basis = torusBasis(t, tau1);
    samples.push({ edge: 0, x: i / n0, t, W: 2 / (t * t * tau1[1]), alphas: [],
      views: [{ key: `one,${a}`, syms: [symOf([], p1), symOf([], pa)], lattice: paritySublattice(basis, par), parity: par, ripe: false }],
      node: i === n0 ? 'opening of the pair for 1' : null });
  }
  nodes.push({ index: samples.length - 1, label: 'opening' });
  // genus 1: the branch point alpha = -r and the pair a, to the first arrival
  const s1 = edge1.samples.filter((sm) => !ev1 || sm.r >= ev1.r);
  const rEnd = ev1 ? ev1.r : s1.length ? s1[s1.length - 1].r : 1;
  for (const sm of s1) {
    const al = [[-sm.r, 0]], pa = sm.pts[a].P;
    samples.push({ edge: 1, x: 1 + Math.log(sm.r) / Math.log(rEnd), t: sm.t, W: 2 * sm.W, alphas: al,
      views: [{ key: `alpha,${a}`, syms: [{ lam: al[0], sign: 1 }, symOf(al, pa)], lattice: paritySublattice(torusBasis(sm.t, tau1), par), parity: par, ripe: false }] });
  }
  if (ev1 && c1.lattice && C.double) {
    // a double arrival: candidates 1 (= 2), 3 and the other choice on the genus 1 curve (keys double:...)
    samples.push({ edge: 1, x: 2, t: ev1.t, W: c1.W, alphas: c1.alphas,
      views: C.choices.map((ch) => ({ key: `double:${ch.pair.join()}`, syms: ch.syms.map((p) => symOf(c1.alphas, p)), lattice: ch.lattice, ripe: true })),
      node: 'candidates 1, 2, 3 (double arrival)' });
    nodes.push({ index: samples.length - 1, label: 'candidates 1, 2, 3' });
  } else if (ev1 && c1.lattice) {
    samples.push({ edge: 1, x: 2, t: ev1.t, W: c1.W, alphas: c1.alphas,
      views: [{ key: 'candidate 1', syms: c1.syms.map((p) => symOf(c1.alphas, p)), lattice: c1.lattice, ripe: true }], node: 'candidate 1' });
    nodes.push({ index: samples.length - 1, label: 'candidate 1' });
  }
  // genus 2: three views; the parities from the characters of the branch points (constant along the edge)
  const edge2 = C.edge && C.c1 ? C.edge : null;
  if (edge2 && C.choices && C.choices.length) {
    const ch0 = C.choices[0], chars = ch0.chars; // chars[w][j], j = 0 (a1), 1 (a2)
    const flip = (x) => (x < 0 ? 1 : 0);
    const parOf = { 'a1,a2': chars.map(([u, v]) => (u === v ? 0 : 1)), 'a1,b': chars.map(([, v]) => flip(v)), 'a2,b': chars.map(([u]) => flip(u)) };
    const ev2 = C.event, list = edge2.samples.filter((sm) => sm.k && (!ev2 || sm.t >= ev2.t));
    const tA = list[0].t, tEnd = ev2 ? ev2.t : list[list.length - 1].t;
    for (const sm of list) {
      if (ev2 && sm.t === ev2.t) continue;
      const al = sm.alphas, basis = dualBasis(sm.k[0], sm.k[1]);
      const pts = { a1: { lam: al[0], sign: 1 }, a2: { lam: al[1], sign: 1 }, b: symOf(al, sm.pts.P) };
      samples.push({ edge: 2, x: 2 + Math.log(tA / sm.t) / Math.log(tA / tEnd), t: sm.t, W: 2 * sm.W, alphas: al,
        views: ['a1,a2', 'a1,b', 'a2,b'].map((key) => ({ key, syms: key.split(',').map((k) => pts[k]), lattice: paritySublattice(basis, parOf[key]), parity: parOf[key], ripe: key === 'a1,a2' })) });
    }
    if (ev2) {
      samples.push({ edge: 2, x: 3, t: ev2.t, W: C.choices[0].W, alphas: C.alphas,
        views: C.choices.map((ch) => ({ key: ch.pair.join(), syms: ch.syms.map((p) => symOf(C.alphas, p)), lattice: ch.lattice, ripe: true })), node: 'candidates 2 and 3' });
      nodes.push({ index: samples.length - 1, label: 'candidates 2, 3' });
    } else if (C.openEnd) nodes.push({ index: samples.length - 1, label: 'open end' });
  }
  // the branch that does not open the node at candidate 1 (edge 3, drawn over the span of genus 2):
  // views alpha + the arrived pair (over S^1: unimodular multipliers, the lattice acts by rotations) and
  // alpha + the other pair; at its end candidate 1'
  const br = c1.branch;
  if (ev1 && br) {
    const r2 = br.event ? br.event.r : Math.min(...edge1.samples.map((x) => x.r));
    const parA = par, parB = [(sOne[0] + CLASS_PQ[b][0]) & 1, (sOne[1] + CLASS_PQ[b][1]) & 1];
    for (const sm of edge1.samples) {
      if (!(sm.r < ev1.r) || sm.r < r2) continue;
      const al = [[-sm.r, 0]], basis = torusBasis(sm.t, tau1);
      samples.push({ edge: 3, x: 2 + Math.log(ev1.r / sm.r) / Math.log(ev1.r / r2), t: sm.t, r: sm.r, W: 2 * sm.W, alphas: al,
        views: [
          { key: 'branch:alpha,arrived', syms: [{ lam: al[0], sign: 1 }, symOf(al, sm.pts[a].P)], lattice: paritySublattice(basis, parA), parity: parA, ripe: false },
          { key: 'branch:alpha,other', syms: [{ lam: al[0], sign: 1 }, symOf(al, sm.pts[b].P)], lattice: paritySublattice(basis, parB), parity: parB, ripe: false },
        ] });
    }
    if (br.event && br.lattice) {
      samples.push({ edge: 3, x: 3, t: br.event.t, r: br.event.r, W: br.W, alphas: br.alphas,
        views: [{ key: 'branch:candidate', syms: br.syms.map((p) => symOf(br.alphas, p)), lattice: br.lattice, ripe: true }],
        node: `candidate 1′ (the pair of class ${br.event.cls} ${br.event.leaving ? 'leaves' : 'arrives at'} S¹)` });
      nodes.push({ index: samples.length - 1, label: 'candidate 1′' });
    }
  }
  return { frame, samples, nodes, end: edge2 ? edge2.end : edge1.end };
}

/**
 * Where a genus 2 curve near the tree of the internal tau (a point of F_6) lies: its own internal tau
 * (in F_6) and its scale t = sqrt|alpha_1 alpha_2|, so that it is the point at scale t of the genus 2 edge
 * of that tau. alphas = [alpha_1, alpha_2] in the order of the edge (alpha_1 from the genus 1 edge, of the
 * class opened first). The period ratio is ell_2 / ell_1 = conj(s tau1 + n) (genus2Anchor); s and n, and
 * the class representatives inverting internalFrame, are taken from tau, so this is continuous in the
 * curve near tau (within its chamber of F_6).
 * Returns { tau, tau1, t }.
 */
export function locateOnTree(alphas, tau) {
  const fr = internalFrame(tau), reps = classReps(tau);
  const [o1, o2] = [0, 1, 2].sort((a, b) => reps[a].len - reps[b].len);
  const [p1, q1] = reps[o1].pq, [p2, q2] = reps[o2].pq;
  const ell = ellVector(alphas), r = cdivc(ell[1], ell[0]);
  // tau1' = s (conj(r) - n), in the upper half plane and nearest to tau1
  let tau1 = null, bd = Infinity;
  for (const sg of [1, -1]) for (let n = -6; n <= 6; n++) {
    const c = [sg * (r[0] - n), -sg * r[1]];
    if (c[1] <= 0) continue;
    const d = Math.hypot(c[0] - fr.tau1[0], c[1] - fr.tau1[1]);
    if (d < bd) { bd = d; tau1 = c; }
  }
  // internalFrame: tau1 = sigma u / v with u = p2 + q2 tau, v = p1 + q1 tau and sigma = ±1 (Im tau1 > 0)
  const uv = cdivc([p2 + q2 * tau[0], q2 * tau[1]], [p1 + q1 * tau[0], q1 * tau[1]]), sigma = uv[1] < 0 ? -1 : 1;
  const w = [sigma * tau1[0], sigma * tau1[1]];
  // p2 + q2 tau = w (p1 + q1 tau)  =>  tau = (w p1 - p2) / (q2 - w q1)
  const t2 = cdivc([w[0] * p1 - p2, w[1] * p1], [q2 - w[0] * q1, -w[1] * q1]);
  const pr = cmul(alphas[0], alphas[1]);
  return { tau: t2[1] < 0 ? [-t2[0], -t2[1]] : t2, tau1, t: Math.sqrt(Math.hypot(...pr)) };
}

/**
 * The open-end shape diagnostics (F12), in place: for each choice { alphas, syms, lattice } (candidate 2
 * at a scale t), shape = roundnessS4 of the surface on a shapeGrid^2 uniform grid of one cell (first and
 * last rows coincide): the distance from a round 2-sphere and the S^4 area / 4 pi after a Möbius
 * normalization. Seconds near t = 1e-3 (the step cap is ~t). Returns the choices.
 */
export function openEndShapes(choices, { shapeGrid = 36 } = {}) {
  const n = shapeGrid;
  for (const ch of choices) {
    const [l1, l2] = ch.lattice.map((e) => Math.hypot(...e));
    const S = computeGenSurface({ alphas: ch.alphas, lams: ch.syms.map((p) => p.lam), cell: ch.lattice, projection: 'drop4',
      sCoords: Array.from({ length: n }, (_, i) => (i * l1) / (n - 1)), tCoords: Array.from({ length: n }, (_, i) => (i * l2) / (n - 1)) });
    ch.shape = roundnessS4(S.pos4, n, n);
  }
  return choices;
}
