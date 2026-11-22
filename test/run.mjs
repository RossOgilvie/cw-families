// Tests for src/cmc against the Python reference values in test/fixtures.json and the checks of §8.
// Run from the repository root:  node test/run.mjs

import { readFileSync } from 'node:fs';
import {
  aPoly, kappa0, killingField, flowKilling, isospectralBasis, divisor, FrameIntegrator,
  computeSurface, hopfArg, closingInfo,
} from '../src/cmc/cmc.js';
import { periodVector, periodLattice, commonRootsOnCircle, basicPeriods, thetaPoly } from '../src/cmc/periods.js';
import { WhithamCurve } from '../src/cmc/whitham.js';
import { RootFlow, symData } from '../src/cmc/rootflow.js';
import { ellipticKE, nodalBoundary, f0, openEndTau } from '../src/cmc/genus2.js';
import { PairIntegrator, leftEigenvector, nuAt, qmul, qdiv, computeGenSurface, doubleCoverSurface } from '../src/cmc/gensym.js';
import { QMap, genus0Resonance, classifyPair } from '../src/cmc/resonance.js';
import { genus0Edge, genus1Scale, candidate1, candidate2, genus2Edge, reduceTau, torusBasis, kOf, treeOf, paritySublattice, dualBasis, locateOnTree, familyOf, superlattices, sameClass } from '../src/cmc/tree.js';
import { generator } from '../src/cmc/cmc.js';
import {
  TauFlow, leafData, periodsOf, rationalPlanes, saturate, reduceFrame, classOrder, nearbyPeriods,
  PlaneFlow, latticeFrame, predictMoves,
} from '../src/cmc/leaf.js';

let failures = 0;
function check(name, value, tol) {
  const ok = value <= tol;
  if (!ok) failures++;
  console.log(`${ok ? '  ok  ' : '  FAIL'} ${name}: ${value.toExponential(2)} (tol ${tol.toExponential(0)})`);
}
const maxDiff = (a, b) => a.reduce((m, x, i) => Math.max(m, Math.abs(x - b[i])), 0);
const setDist = (A, B) => A.reduce((m, [ar, ai]) =>
  Math.max(m, Math.min(...B.map(([br, bi]) => Math.hypot(ar - br, ai - bi)))), 0);

// ------------------------------------------------------------------ fixtures
console.log('fixtures (JS vs the Python reference):');
const fixtures = JSON.parse(readFileSync(new URL('./fixtures.json', import.meta.url)));
for (const fx of fixtures) {
  const g = fx.alphas.length;
  const tag = `g=${g}`;
  check(`${tag} a(lam)`, maxDiff(Array.from(aPoly(fx.alphas)), fx.a), 1e-14);
  const Z = killingField(fx.alphas, fx.A);
  check(`${tag} zeta base`, maxDiff(Array.from(Z), fx.Zbase), 1e-10);
  if (fx.times) {
    const kap = kappa0(fx.alphas);
    isospectralBasis(g).forEach(([k, c], i) => flowKilling(Z, kap, k, c, fx.times[i], 0.002));
    check(`${tag} zeta after isospectral flows`, maxDiff(Array.from(Z), fx.Zinit), 1e-9);
  }
  check(`${tag} divisor`, setDist(divisor(Z).map((d) => d.mu), fx.divisor), 1e-9);
  const I = new FrameIntegrator(fx.alphas, fx.theta0);
  for (const p of fx.points) {
    const y = I.initialState(Z);
    I.march(y, 1, 0, p.z[0], 0.005);
    I.march(y, 0, 1, p.z[1], 0.005);
    const P = new Float64Array(3), N = new Float64Array(3);
    const u = I.evaluate(y, P, 0, N, 0);
    const err = Math.max(maxDiff(Array.from(P), p.f), maxDiff(Array.from(N), p.N), Math.abs(u - p.u));
    check(`${tag} f, N, u at z=${p.z}`, err, 1e-9);
  }
}

// ------------------------------------------------------------------ geometric invariants on the grid
console.log('finite-difference invariants on computeSurface (errors are O(h^2)):');
function geomCheck(label, params) {
  const S = computeSurface(params);
  const { nx, ny, pos, nrm, u } = S;
  const hs = params.width / (nx - 1), ht = params.height / (ny - 1);
  const P = (i, j, e) => pos[3 * (j * nx + i) + e];
  let conf = 0, orth = 0, metric = 0, Herr = 0, Qerr = 0, gauss = 0;
  // expected Hopf differential in the rotated coordinate w = s + i t: Q_w = e^{2 i phi} Q_z
  const argQ = hopfArg(params.alphas, params.theta0) + 2 * S.phi;
  for (let j = 2; j < ny - 2; j++) {
    for (let i = 2; i < nx - 2; i++) {
      const fs = [], ft = [], fss = [], ftt = [], fst = [], Nn = [];
      for (let e = 0; e < 3; e++) {
        fs[e] = (P(i + 1, j, e) - P(i - 1, j, e)) / (2 * hs);
        ft[e] = (P(i, j + 1, e) - P(i, j - 1, e)) / (2 * ht);
        fss[e] = (P(i + 1, j, e) - 2 * P(i, j, e) + P(i - 1, j, e)) / (hs * hs);
        ftt[e] = (P(i, j + 1, e) - 2 * P(i, j, e) + P(i, j - 1, e)) / (ht * ht);
        fst[e] = (P(i + 1, j + 1, e) - P(i + 1, j - 1, e) - P(i - 1, j + 1, e) + P(i - 1, j - 1, e)) / (4 * hs * ht);
        Nn[e] = nrm[3 * (j * nx + i) + e];
      }
      const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
      const E = dot(fs, fs), G = dot(ft, ft), F = dot(fs, ft);
      const uu = u[j * nx + i];
      conf = Math.max(conf, Math.abs(E - G) / E);
      orth = Math.max(orth, Math.abs(F) / E);
      metric = Math.max(metric, Math.abs(E / (4 * Math.exp(2 * uu)) - 1));
      const H = (dot(fss, Nn) + dot(ftt, Nn)) / (2 * E);
      Herr = Math.max(Herr, Math.abs(H - 0.5));
      const Qr = 0.25 * (dot(fss, Nn) - dot(ftt, Nn)), Qi = -0.5 * dot(fst, Nn);
      Qerr = Math.max(Qerr, Math.hypot(Qr - Math.cos(argQ), Qi - Math.sin(argQ)));
      // sinh-Gordon: Delta u + 2 sinh 2u = 0
      const lap = (u[j * nx + i + 1] - 2 * uu + u[j * nx + i - 1]) / (hs * hs)
        + (u[(j + 1) * nx + i] - 2 * uu + u[(j - 1) * nx + i]) / (ht * ht);
      gauss = Math.max(gauss, Math.abs(lap + 2 * Math.sinh(2 * uu)) / (1 + Math.abs(2 * Math.sinh(2 * uu))));
    }
  }
  const tol = 3e-3;
  check(`${label}: |E-G|/E`, conf, tol);
  check(`${label}: |F|/E`, orth, tol);
  check(`${label}: E/(4e^2u) - 1`, metric, tol);
  check(`${label}: H - 1/2`, Herr, tol);
  check(`${label}: Q error`, Qerr, tol);
  check(`${label}: sinh-Gordon residual`, gauss, 2e-2);
  check(`${label}: det zeta drift`, S.stats.detDefect, 1e-8);
}
const base = { width: 1.5, height: 1.5, nx: 61, ny: 61, hmax: 0.01 };
geomCheck('g=0', { ...base, alphas: [], theta0: 0.3 });
geomCheck('g=1', { ...base, alphas: [[0.4, 0.2]], theta0: 0.7, z0: [0.3, -0.4] });
geomCheck('g=2', { ...base, alphas: [[0.5, 0.1], [-0.3, 0.6]], theta0: 1.3, phi: 0.2 });
geomCheck('g=3', { ...base, alphas: [[0.3, 0], [0, 0.5], [-0.6, -0.2]], theta0: Math.PI / 2, tau: [0.4] });
geomCheck('g=4', { ...base, alphas: [[0.6, 0], [0, 0.4], [-0.5, -0.3], [0.2, -0.7]], theta0: 2, tau: [0.3, -0.5] });

// ------------------------------------------------------------------ Delaunay radii
console.log('Delaunay surfaces (g=1): neck and bulge radii of the parallels:');
function circumradius(a, b, c) {
  const d = (p, q) => Math.hypot(p[0] - q[0], p[1] - q[1], p[2] - q[2]);
  const A = d(b, c), B = d(a, c), C = d(a, b);
  const s = (A + B + C) / 2;
  return (A * B * C) / (4 * Math.sqrt(s * (s - A) * (s - B) * (s - C)));
}
for (const r of [0.3, 0.7]) {
  for (const sign of [-1, 1]) {
    const th = 0.9;
    const S = computeSurface({
      alphas: [[r * Math.cos(th), r * Math.sin(th)]], theta0: th + (sign < 0 ? Math.PI : 0),
      phi: 0, width: 1.4, height: 12, nx: 3, ny: 1201, hmax: 0.01,
    });
    let lo = Infinity, hi = 0;
    for (let j = 0; j < S.ny; j++) {
      const p = (i) => Array.from(S.pos.subarray(3 * (j * 3 + i), 3 * (j * 3 + i) + 3));
      const R = circumradius(p(0), p(1), p(2));
      lo = Math.min(lo, R); hi = Math.max(hi, R);
    }
    const exp = sign < 0 ? [2 * r / (1 + r), 2 / (1 + r)] : [2 * r / (1 - r), 2 / (1 - r)];
    check(`|alpha|=${r} ${sign < 0 ? 'unduloid' : 'nodoid'} radii [${lo.toFixed(4)}, ${hi.toFixed(4)}]`,
      Math.max(Math.abs(lo - exp[0]), Math.abs(hi - exp[1])), 1e-4);
  }
}
{
  // curvature-line grid: one family of parameter lines is the unit circles, the other straight lines
  const S = computeSurface({ alphas: [], theta0: 1.1, width: 2, height: 2, nx: 3, ny: 3, hmax: 0.01 });
  const p = (i, j) => Array.from(S.pos.subarray(3 * (3 * j + i), 3 * (3 * j + i) + 3));
  const Rs = circumradius(p(0, 1), p(1, 1), p(2, 1)), Rt = circumradius(p(1, 0), p(1, 1), p(1, 2));
  check('g=0 cylinder radius 1', Math.min(Math.abs(Rs - 1), Math.abs(Rt - 1)), 1e-9);
}

// ------------------------------------------------------------------ commutativity and divisor
console.log('commutativity around a closed rectangle, and base-point divisor:');
{
  const alphas = [[0.5, 0.1], [-0.3, 0.6], [0.1, -0.4]];
  const I = new FrameIntegrator(alphas, 0.8);
  const Z = killingField(alphas);
  flowKilling(Z, kappa0(alphas), 1, [0, 1], 0.7);
  const y0 = I.initialState(Z);
  const y = y0.slice();
  I.march(y, 1, 0, 1.3, 0.005); I.march(y, 0, 1, 0.9, 0.005);
  I.march(y, 1, 0, -1.3, 0.005); I.march(y, 0, 1, -0.9, 0.005);
  check('g=3 rectangle holonomy of (zeta, F, G)', maxDiff(Array.from(y), Array.from(y0)), 1e-9);
  const alphas2 = [[0.3, 0], [0, 0.5], [-0.2, -0.6]];
  check('base-point divisor = branch points', setDist(divisor(killingField(alphas2)).map((d) => d.mu), alphas2), 1e-12);
}

// ------------------------------------------------------------------ periods and closing (§9)
console.log('periods of Theta_w and closing conditions:');
{
  // g = 1: the y-period of the unduloid, against shooting (closingInfo)
  const T = Math.abs(1 / periodVector([[0.45, 0]], [0, 1])[0]);
  check('g=1 y-period from Theta_w vs shooting', Math.abs(T - closingInfo({ alphas: [[0.45, 0]], theta0: Math.PI }, 't').T), 1e-6);
  // g = 2: every lattice generator is a period of zeta (§9.4)
  const al = [[0.3, 0.4], [0.3, -0.4]];
  const { generators } = periodLattice(al);
  let worst = 0;
  for (const w of generators) {
    const I = new FrameIntegrator(al, 0), Z = killingField(al), y = I.initialState(Z);
    const L = Math.hypot(...w);
    I.march(y, w[0] / L, w[1] / L, L, 0.002);
    for (let e = 0; e < 32; e++) worst = Math.max(worst, Math.abs(y[e] - Z[e]));
  }
  check(`g=2 lattice generators (${generators.length}) are periods of zeta`, generators.length === 2 ? worst : Infinity, 1e-9);
  // Wente torus: closes around 3 w1 and 2 w2 - w1
  const W = [[0.1412634686, 0.1017768953], [0.1412634686, -0.1017768953]];
  // (generator signs are arbitrary: use T1 = |w1| on the real axis and T2/2 = |Im w2|)
  const [w1, w2] = periodLattice(W).generators;
  const T1 = Math.abs(w1[0]), T2h = Math.abs(w2[1]);
  let err = 0;
  for (const w of [[3 * T1, 0], [1.5 * T1, T2h], [0, 2 * T2h]]) {
    const I = new FrameIntegrator(W, 0), Z = killingField(W), y = I.initialState(Z);
    const L = Math.hypot(...w);
    I.march(y, w[0] / L, w[1] / L, L, 0.002);
    const P = new Float64Array(3), N = new Float64Array(3);
    I.evaluate(y, P, 0, N, 0);
    for (let e = 0; e < 32; e++) err = Math.max(err, Math.abs(y[e] - Z[e]));
    err = Math.max(err, Math.hypot(...P), Math.abs(Math.abs(y[32]) - 1));
  }
  check('Wente torus closes on the brick lattice <3 T1, 3 T1/2 + i T2/2>', err, 1e-7);
  // bubbleton: closes around the cylinder at the resonance point
  const r2 = (2 - Math.sqrt(3)) ** 2;
  const c = closingInfo({ alphas: [[r2, 0], [r2, 0]], theta0: 0, hmax: 0.002 }, 's');
  check('two-lobed bubbleton closes after one period pi', c && c.q === 1 ? Math.abs(c.T - Math.PI) : Infinity, 1e-5);
}

// ------------------------------------------------------------------ Whitham deformation
console.log('Whitham deformation keeps the period lattice up to Gamma(0)/z:');
for (const [name, al] of [
  ['Wente', [[0.1412634686, 0.1017768953], [0.1412634686, -0.1017768953]]],
  ['generic g=2', [[0.4, 0.25], [-0.3, 0.5]]],
]) {
  const C = new WhithamCurve(al);
  const L0 = periodLattice(al).generators;
  const cov = (L) => { // lattice invariants: |w1|, |w2|, |w1 x w2|, up to the choice of generators
    const [a, b] = L;
    return [Math.hypot(...a), Math.hypot(...b), Math.abs(a[0] * b[1] - a[1] * b[0])];
  };
  let worst = 0, ok = true;
  for (const s of [-0.1, 0.1, 0.25]) {
    const r = C.at(s);
    if (Math.abs(r.s - s) > 1e-3) { ok = false; continue; } // s is the arclength actually travelled
    // Gamma(s) z = Gamma(0): compare the invariants of Gamma(s) scaled by |z| (and |z|^2 for the area)
    const L = periodLattice(r.alphas).generators;
    if (L.length < 2 || L0.length < 2) { ok = false; continue; }
    const m = Math.hypot(...r.z);
    const [a, b, A] = cov(L), [a0, b0, A0] = cov(L0);
    worst = Math.max(worst, Math.abs(a * m - a0) / a0, Math.abs(b * m - b0) / b0, Math.abs(A * m * m - A0) / A0);
  }
  check(`${name}: lattice of the deformed curve = Gamma(0)/z`, ok ? worst : Infinity, 1e-6);
}
{
  // the Willmore functional of the paper is real along the flow, and stationary where B_a has a
  // common root (the Wente torus has one, at lam0 = 1)
  const C = new WhithamCurve([[0.1412634686, 0.1017768953], [0.1412634686, -0.1017768953]]);
  let im = 0;
  for (const s of [-0.1, 0, 0.3]) { const W = C.willmore(C.at(s)); im = Math.max(im, Math.abs(W[1] / W[0])); }
  check('Willmore functional is real', im, 1e-10);
  const e = 1e-3, w = (s) => C.willmore(C.at(s))[0];
  check('dW/ds = 0 at the Wente torus (relative to the change over 0.2)',
    Math.abs((w(e) - w(-e)) / (2 * e)) / Math.abs((w(0.2) - w(0)) / 0.2), 1e-3);
}
{
  // common roots of B_a on S^1: the Wente torus has one at its Sym point lam0 = 1, generic data none
  const W = commonRootsOnCircle([[0.1412634686, 0.1017768953], [0.1412634686, -0.1017768953]], 1e-5);
  check('Wente: common root of the differentials at lam = 1', W.length === 1 ? Math.hypot(W[0].lam[0] - 1, W[0].lam[1]) : Infinity, 1e-9);
  check('generic genus 2: no common root', commonRootsOnCircle([[0.4, 0.25], [-0.3, 0.5]], 1e-5).length, 0);
}
{
  // genus 3 (generally no lattice): the curve continues and keeps [ell] in CP^2 fixed
  const C = new WhithamCurve([[0.4, 0.25], [-0.3, 0.5], [0.1, -0.6]]);
  const r = C.at(0.2), q = C.at(-0.2);
  check('g=3: Whitham curve continues to s = ±0.2', Math.max(Math.abs(r.s - 0.2), Math.abs(q.s + 0.2)), 1e-3);
}

// ------------------------------------------------------------------ root-preserving flow (CKKS), genus 2
{
  const PI = Math.PI, W = [[0.1412634686, 0.1017768953], [0.1412634686, -0.1017768953]];
  const phi = symData(W, 0).phi;
  check('Wente: Sym integrals phi = (pi/3, pi/3)', Math.hypot(phi[0] - PI / 3, phi[1] - PI / 3), 1e-7);
  const run = (f, target) => { let r; for (let i = 0; i < 500; i++) { r = f.moveToward(target, { budget: 50 }); if (r.reached || r.blocked) break; } return r; };
  const f = new RootFlow(W, 0);
  const r = run(f, [PI / 4, PI / 2]);
  const cr = commonRootsOnCircle(r.alphas, 1e-5);
  check('flow reaches phi = (pi/4, pi/2)', r.reached ? Math.hypot(r.phi[0] - PI / 4, r.phi[1] - PI / 2) : Infinity, 1e-9);
  check('flow keeps the common root at lam0 = 1', cr.length ? Math.min(...cr.map((c) => Math.hypot(c.lam[0] - 1, c.lam[1]))) : Infinity, 1e-8);
  run(f, [0.4 * PI, 0.2 * PI]);
  const back = run(f, [PI / 3, PI / 3]);
  check('flow around a loop returns to the Wente data (no holonomy)', Math.hypot(...back.alphas.flat().map((v, i) => v - W.flat()[i])), 1e-7);
  const end = run(new RootFlow(W, 0), [0.6 * PI, 0.6 * PI]);
  check('flow stops just short of the hypotenuse phi_1 + phi_2 = pi', end.blocked ? Math.abs(end.phi[0] + end.phi[1] - PI) / PI : Infinity, 5e-3);
}

// ------------------------------------------------------------------ root-preserving flow, genus 4
// Data on S^4 with real coefficients (the critical point of W on the Whitham curve through a symmetric
// start), lam0 = 1. The flow keeps the common root and the period plane; phi_2 = 0 by the symmetry.
{
  const PI = Math.PI;
  const A = [[0.241595657003, 0.131253423469], [0.241595657003, -0.131253423469], [-0.226370306867, 0.288849491705], [-0.226370306867, -0.288849491705]];
  const d0 = symData(A, 0);
  check('genus 4: phi_2 = 0 for real coefficients', Math.abs(d0.phi[1]), 1e-9);
  const P = d0.P.map((Pl, l) => (d0.phi[l] < 0 ? Pl.map((v) => -v) : Pl));
  const f = new RootFlow(A, 0, P);
  const target = [f.phi[0] + 0.004 * PI, f.phi[1] + 0.004 * PI];
  let r; for (let i = 0; i < 500; i++) { r = f.moveToward(target, { budget: 50, step: 0.005 }); if (r.reached || r.blocked) break; }
  check('genus 4: flow reaches a nearby target', r.reached ? Math.hypot(r.phi[0] - target[0], r.phi[1] - target[1]) : Infinity, 1e-9);
  const cr = commonRootsOnCircle(r.alphas, 1e-5);
  check('genus 4: flow keeps the common root at lam0 = 1', cr.length ? Math.min(...cr.map((c) => Math.hypot(c.lam[0] - 1, c.lam[1]))) : Infinity, 1e-8);
  const d = symData(r.alphas, 0, P, f.ref);
  const planeErr = Math.max(...[0, 1].flatMap((l) => { const c = d.Bm.map((row) => row[l]), n = Math.hypot(...c); return f.normals.map((N) => Math.abs(c.reduce((acc, v, j) => acc + v * N[j], 0) / n)); }));
  check('genus 4: flow keeps the period plane', planeErr, 1e-9);
}

// ------------------------------------------------------------------ leaves: the tau pad, LLL, locking periods
console.log('leaves of the Whitham foliation (leaf.js):');
{
  // tau at fixed psi = sum logit|alpha_j|: the periods of the frame and psi stay fixed, and a loop in tau
  // returns to the start
  const run = (f, t) => { let r; for (let i = 0; i < 2000; i++) { r = f.moveToward(t, { budget: 50 }); if (r.reached || r.blocked) break; } return r; };
  for (const al of [[[0.4, 0.25], [-0.3, 0.5]], [[0.4, 0.25], [-0.3, 0.5], [0.1, -0.6]]]) {
    const g = al.length, f = new TauFlow(al, [[1, 0], [0, 1]]);
    const r = run(f, [0.01, 1.01]);
    const d = leafData(r.alphas, f.ref);
    const perr = Math.max(...[0, 1].flatMap((k) => periodsOf(d.ell, r.w[k]).map((v, j) => Math.abs(v - f.P[k][j]))));
    check(`g=${g} tau pad reaches tau = 0.01 + 1.01 i`, r.reached ? Math.hypot(r.tau[0] - 0.01, r.tau[1] - 1.01) : Infinity, 1e-9);
    check(`g=${g} tau pad keeps the periods of the frame`, perr, 1e-9);
    const psi = r.alphas.reduce((s, [a, b]) => { const q = Math.hypot(a, b); return s + Math.log(q / (1 - q)); }, 0);
    check(`g=${g} tau pad keeps psi`, Math.abs(psi - f.psi0), 1e-9);
    run(f, [-0.01, 0.99]);
    const back = run(f, [0, 1]);
    check(`g=${g} tau pad around a loop returns to the start (no holonomy)`, Math.hypot(...back.alphas.flat().map((v, i) => v - al.flat()[i])), 1e-8);
  }
}
{
  // LLL finds a rational plane, with a basis of its lattice: span(m1, m2) with m1 = (1, 2, 0, -1), m2 = (0, 1, 3, 1)
  const m1 = [1, 2, 0, -1], m2 = [0, 1, 3, 1];
  const v1 = m1.map((v, j) => 0.7 * v + 0.2 * m2[j]), v2 = m1.map((v, j) => -0.3 * v + 1.1 * m2[j]);
  const P = rationalPlanes(v1.map((v, j) => [v, -v2[j]]));
  const best = P.find((p) => p.dist < 1e-10);
  // the found basis spans the same lattice: its coordinates in (m1, m2) form a unimodular matrix
  const coords = best ? best.m.map((m) => { const a = m[0], b = m[1] - 2 * a; return [a, b, m.every((v, j) => v === a * m1[j] + b * m2[j])]; }) : null;
  check('LLL recognises a rational plane with a basis of its lattice',
    coords && coords.every((c) => c[2]) ? Math.abs(Math.abs(coords[0][0] * coords[1][1] - coords[0][1] * coords[1][0]) - 1) : Infinity, 0);
  const S = saturate([2, 4, 0, -2], [1, 3, 3, 0]); // 2 m1 and m1 + m2: index 2 in the lattice
  const area = (a, b) => Math.sqrt(a.reduce((s, x) => s + x * x, 0) * b.reduce((s, x) => s + x * x, 0) - a.reduce((s, x, j) => s + x * b[j], 0) ** 2);
  check('saturate: the lattice of the plane (same area as m1, m2)', Math.abs(area(...S) - area(m1, m2)), 1e-9);
}
{
  // reduction into F_6, F_3, F_1 keeps the lattice and the classes mod 2
  let worst = 0;
  for (const kind of ['F6', 'F3', 'F1']) for (const t of [[3.7, 0.05], [-0.45, 0.3], [0.2, 0.01], [-7.1, 2]]) {
    const [b1, b2] = reduceFrame([{ w: [1, 0], m: [1, 0] }, { w: t, m: [0, 1] }], kind);
    const [x, y] = [(b2.w[0] * b1.w[0] + b2.w[1] * b1.w[1]) / (b1.w[0] ** 2 + b1.w[1] ** 2), (b2.w[1] * b1.w[0] - b2.w[0] * b1.w[1]) / (b1.w[0] ** 2 + b1.w[1] ** 2)];
    const inside = kind === 'F1' ? Math.abs(x) <= 0.5 + 1e-9 && x * x + y * y >= 1 - 1e-9
      : Math.abs(x) <= 1 + 1e-9 && Math.hypot(x - 0.5, y) >= 0.5 - 1e-9 && Math.hypot(x + 0.5, y) >= 0.5 - 1e-9 && (kind === 'F6' || x * x + y * y >= 1 - 1e-9);
    const det = b1.m[0] * b2.m[1] - b1.m[1] * b2.m[0];
    const parity = kind === 'F6' ? (b1.m[0] % 2 !== 0) && (b1.m[1] % 2 === 0) && (b2.m[0] % 2 === 0) && (b2.m[1] % 2 !== 0) : true;
    if (!inside || det !== 1 || !parity) worst = Infinity;
  }
  check('reduction into F_6 (by Gamma(2)), F_3 and F_1, with determinant 1', worst, 0);
  const labels = new Set([[0, 1.6], [0.78, 1.6], [0, 0.75], [0.3, 0.62], [0.75, 0.82], [0.8, 0.45]].map(classOrder));
  check('the six copies of F_1 in F_6 have six different class orders', 6 - labels.size, 0);
}
{
  // lock two integer periods into the plane (genus 3): the data become doubly periodic, and zeta closes up
  const al = [[0.4, 0.25], [-0.3, 0.5], [0.1, -0.6]];
  const go = (f) => { let s; for (let i = 0; i < 500; i++) { s = f.advance(100); if (s.done || s.blocked) break; } return s; };
  const d0 = leafData(al), pts = nearbyPeriods(d0.ell);
  const f1 = new PlaneFlow(al, [pts[0].m], d0.ell), s1 = go(f1);
  const d1 = f1.data(), second = nearbyPeriods(d1.ell).find((p) => p.dist > 1e-8);
  const pred = predictMoves(al, d0.ell, [pts[0].m])[0], moved = Math.hypot(...s1.alphas.flat().map((v, i) => v - al.flat()[i]));
  check('predicted move of a lock vs the actual move (relative)', Math.abs(pred - moved) / moved, 0.05);
  check('lock one period: it lies in the plane', s1.done ? nearbyPeriods(d1.ell).find((p) => p.m.join() === pts[0].m.join()).dist : Infinity, 1e-9);
  const f2 = new PlaneFlow(s1.alphas, [pts[0].m, second.m], d1.ell), s2 = go(f2);
  const d2 = f2.data(), L = latticeFrame(d2.ell, ...saturate(pts[0].m, second.m));
  const intErr = Math.max(...L.w.flatMap((w) => periodsOf(d2.ell, w).map((v) => Math.abs(v - Math.round(v)))));
  check('lock two periods: the lattice basis has integer periods', s2.done ? intErr : Infinity, 1e-9);
  const Z = killingField(s2.alphas), I = new FrameIntegrator(s2.alphas, 0);
  let zerr = 0;
  for (const w of L.w) {
    const y0 = I.initialState(Z), y = y0.slice(), r = Math.hypot(...w);
    I.march(y, w[0] / r, w[1] / r, r, 0.002);
    for (let k = 0; k < 8 * 5; k++) zerr = Math.max(zerr, Math.abs(y[k] - y0[k]));
  }
  check('lock two periods: zeta(z + w) = zeta(z) for both basis vectors (shooting)', zerr, 1e-8);
}

// ------------------------------------------------------------------ genus 2: the boundary of S^2 in F_3
console.log('genus 2: the boundary of S^2 (genus2.js) against the periods:');
{
  const { K, E } = ellipticKE([0.6, 0]);
  check('K(0.6), E(0.6) by the AGM', Math.max(Math.abs(K[0] - 1.7507538029157526), Math.abs(E[0] - 1.4180833944487245)), 1e-13);
  // nodal end: open the node at -e^{i phi}, t = cot(phi/2) = t0(k), slightly; simple root at -r
  const k = 0.5, { t0, tau } = nodalBoundary(k), r = (1 - k) / (1 + k), phi = 2 * Math.atan(1 / t0);
  const al = [[-(1 - 1e-4) * Math.cos(phi), -(1 - 1e-4) * Math.sin(phi)], [-r, 0]];
  const L = latticeFrame(leafData(al).ell, [1, 0], [0, 1]);
  const err = Math.min(Math.hypot(L.tau[0] - tau[0], L.tau[1] - tau[1]), Math.hypot(L.tau[0] + tau[0], L.tau[1] - tau[1]));
  check('nodal end: tau of the paper = tau of the lattice (node opened by 1e-4)', err, 1e-4);
  // open end: on the real locus of f0, rescaled roots s alpha, s / alpha with s = 1e-3
  const rho = 0.6, g = (th) => f0([rho * Math.cos(th), rho * Math.sin(th)])[1];
  let lo = 1e-6, hi = Math.PI / 2 - 1e-6;
  const glo = g(lo);
  for (let i = 0; i < 60; i++) { const m = (lo + hi) / 2; if (Math.sign(g(m)) === Math.sign(glo)) lo = m; else hi = m; }
  const a = [rho * Math.cos(lo), rho * Math.sin(lo)], s = 1e-3, a2 = a[0] ** 2 + a[1] ** 2;
  const d = leafData([[s * a[0], s * a[1]], [(s * a[0]) / a2, (-s * a[1]) / a2]]), M = latticeFrame(d.ell, [1, 0], [0, 1]);
  const t2 = openEndTau(a);
  const err2 = Math.min(Math.hypot(M.tau[0] - t2[0], M.tau[1] - t2[1]), Math.hypot(M.tau[0] + t2[0], M.tau[1] - t2[1]));
  check('open end: tau of the formula = tau of the lattice (s = 1e-3)', err2, 2e-3);
  const lin = (w) => [w[0] * d.p1[1][0] + w[1] * d.pi[1][0], w[0] * d.p1[1][1] + w[1] * d.pi[1][1]];
  const [u, v] = [lin(M.w[0]), lin(M.w[1])];
  check('open end: on the real locus of f0 the ratio of linear coefficients is real', Math.abs(u[0] * v[1] - u[1] * v[0]) / Math.hypot(...u) / Math.hypot(...v), 1e-3);
}

// ------------------------------------------------------------------ generalized Sym formula (gensym.js)
console.log('generalized Sym formula f = Q(u1 F(lam1)) Q(u2 F(lam2))^{-1} (cw_candidates.md, F7):');
{
  let seed = 7;
  const rnd = () => { seed = (seed * 16807) % 2147483647; return seed / 2147483647 - 0.5; };
  const at = (I, Z, x, y) => {
    const st = I.initialState(Z);
    I.march(st, 1, 0, x, 0.005); I.march(st, 0, 1, y, 0.005);
    return st;
  };
  const qn = (q) => q[0] ** 2 + q[1] ** 2 + q[2] ** 2 + q[3] ** 2;
  const qdot = (p, q) => p[0] * q[0] + p[1] * q[1] + p[2] * q[2] + p[3] * q[3];
  // 1. conformal for arbitrary lambda in C^* and arbitrary rows (all sections are holomorphic sections of one bundle)
  for (const alphas of [[], [[0.5, 0]], [[0.4, 0.3], [-0.2, 0.5]]]) {
    const Z = killingField(alphas);
    for (const lams of [[[Math.cos(0.4), Math.sin(0.4)], [Math.cos(2), Math.sin(2)]], [[0.57, 0.18], [-0.13, 0.99]], [[0.5, 0.2], [1.8, -0.9]]]) {
      const I = new PairIntegrator(alphas, lams);
      const us = [[rnd(), rnd(), rnd(), rnd()], [rnd(), rnd(), rnd(), rnd()]];
      let worst = 0;
      for (const [x, y] of [[0.3, 0.7], [1.1, -0.4], [-0.8, 1.3]]) {
        const [, fx, fy] = I.jet(at(I, Z, x, y), us);
        const n = qn(fx);
        worst = Math.max(worst, Math.abs(n - qn(fy)) / n, Math.abs(qdot(fx, fy)) / n);
      }
      check(`conformal, g=${alphas.length}, lam = ${lams.map(([r, i]) => `${r.toFixed(2)}${i < 0 ? '-' : '+'}${Math.abs(i).toFixed(2)}i`).join(', ')}`, worst, 1e-12);
    }
  }
  // 2. on S^1 it is the S^3 Sym–Bobenko formula: Q(u1)^{-1} f Q(u2) = F(lam1) F(lam2)^{-1}
  {
    const alphas = [[0.5 * Math.cos(0.4), 0.5 * Math.sin(0.4)]], Z = killingField(alphas);
    const lams = [[Math.cos(0.9), Math.sin(0.9)], [Math.cos(2.3), Math.sin(2.3)]];
    const I = new PairIntegrator(alphas, lams);
    const us = [[rnd(), rnd(), rnd(), rnd()], [rnd(), rnd(), rnd(), rnd()]];
    let worst = 0;
    for (const [x, y] of [[0.5, -1.2], [-1.4, 0.3]]) {
      const st = at(I, Z, x, y), f = new Float64Array(4);
      I.evaluate4(st, us, f, 0);
      const e1 = [1, 0, 0, 0], [F1, F2] = I.rows(st, [e1, e1]), sym = new Float64Array(4);
      qdiv(F1, F2, sym, 0);
      // Q(u1)^{-1} f Q(u2), with Q(u1)^{-1} = conj(Q(u1)) / |u1|^2
      const alt = qmul([us[0][0], -us[0][1], -us[0][2], -us[0][3]], qmul(f, us[1])).map((x) => x / qn(us[0]));
      worst = Math.max(worst, Math.max(...alt.map((x, e) => Math.abs(x - sym[e]))));
    }
    check('lam_1, lam_2 on S^1: Q(u1)^{-1} f Q(u2) = F(lam1) F(lam2)^{-1}', worst, 1e-9);
  }
  // 3. genus 0 homogeneous tori: two rho-pairs over S^1, closing on (pi/cos s) x (pi/sin s),
  //    W = int |Laplace f|^2 / (4 |f_x|^2) = pi^2 (t + 1/t) with t = cot s; t = 1 is the Clifford torus.
  for (const sAng of [Math.PI / 4, 0.6, 0.4]) {
    const L1 = Math.PI / Math.cos(sAng), L2 = Math.PI / Math.sin(sAng), n = 64;
    const lams = [[Math.cos(2 * sAng), Math.sin(2 * sAng)], [Math.cos(2 * sAng), -Math.sin(2 * sAng)]];
    const Z = killingField([]), I = new PairIntegrator([], lams);
    let best = null;
    for (const s1 of [1, -1]) for (const s2 of [1, -1]) {
      const us = lams.map((lam, j) => leftEigenvector(Z, lam, nuAt([], lam, [s1, s2][j])));
      const F = (x, y) => { const f = new Float64Array(4); I.evaluate4(at(I, Z, x, y), us, f, 0); return f; };
      const f00 = F(0.3, 0.2);
      const per = Math.max(...[F(0.3 + L1, 0.2), F(0.3, 0.2 + L2)].map((q) => Math.max(...q.map((x, e) => Math.abs(x - f00[e])))));
      if (!best || per < best.per) best = { per, us };
    }
    check(`s = ${sAng.toFixed(3)}: f closes on the rectangle (pi/cos s) x (pi/sin s)`, best.per, 1e-8);
    // f on a periodic grid, then fourth-order periodic differences
    const fs = [];
    for (let j = 0; j < n; j++) {
      const st = at(I, Z, 0, (j * L2) / n);
      for (let i = 0; i < n; i++) {
        const f = new Float64Array(4); I.evaluate4(st, best.us, f, 0); fs.push(f);
        I.march(st, 1, 0, L1 / n, 0.005);
      }
    }
    const h1 = L1 / n, h2 = L2 / n, G = (i, j) => fs[((j + n) % n) * n + ((i + n) % n)];
    let W = 0;
    for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) {
      const fx = [0, 1, 2, 3].map((e) => (-G(i + 2, j)[e] + 8 * G(i + 1, j)[e] - 8 * G(i - 1, j)[e] + G(i - 2, j)[e]) / (12 * h1));
      const d2 = (A, B, C, D, E, h) => (-A + 16 * B - 30 * C + 16 * D - E) / (12 * h * h);
      const lap = [0, 1, 2, 3].map((e) => d2(G(i - 2, j)[e], G(i - 1, j)[e], G(i, j)[e], G(i + 1, j)[e], G(i + 2, j)[e], h1)
        + d2(G(i, j - 2)[e], G(i, j - 1)[e], G(i, j)[e], G(i, j + 1)[e], G(i, j + 2)[e], h2));
      W += qn(lap) / (4 * qn(fx));
    }
    W *= h1 * h2;
    const t = 1 / Math.tan(sAng);
    check(`s = ${sAng.toFixed(3)}: W = pi^2 (t + 1/t) = ${(Math.PI ** 2 * (t + 1 / t)).toFixed(6)} (relative)`, Math.abs(W / (Math.PI ** 2 * (t + 1 / t)) - 1), 1e-5);
  }
  // 4. computeGenSurface runs and agrees with the integrator at a vertex
  {
    const alphas = [[0.5, 0.1]], lams = [[0.5, 0.1], [Math.cos(2.3), Math.sin(2.3)]];
    const S = computeGenSurface({ alphas, lams, signs: [1, 1], width: 2, height: 2, nx: 5, ny: 5, hmax: 0.005 });
    const I = new PairIntegrator(alphas, lams), f = new Float64Array(4);
    I.evaluate4(at(I, killingField(alphas), 1, 1), S.us, f, 0);
    const v = 4 * 5 + 4;
    check('computeGenSurface vertex (1, 1) = direct integration (Float32)', Math.max(...[0, 1, 2, 3].map((e) => Math.abs(S.pos4[4 * v + e] - f[e]))), 1e-5);
  }
  // 5. automatic projection: the Clifford torus is recognized as lying on a 3-sphere (and its stereographic image
  //    is bounded); a branch point as a Sym point gives a surface on no 3-sphere (Möbius-balanced)
  {
    const L = Math.PI * Math.SQRT2;
    const C = computeGenSurface({ alphas: [], lams: [[0, 1], [0, -1]], signs: [1, 1], width: L, height: L, nx: 40, ny: 40, hmax: 0.01 });
    const ext = Math.max(C.bbox[3] - C.bbox[0], C.bbox[4] - C.bbox[1], C.bbox[5] - C.bbox[2]);
    check(`Clifford torus: projection = ${C.projection}, stereographic extent ${ext.toFixed(2)} bounded`, C.projection === 'sphere' && ext < 20 ? 0 : 1, 0.5);
    const B = computeGenSurface({ alphas: [[0.5, 0]], lams: [[0.5, 0], [0, 1]], signs: [1, 1], width: 4, height: 4, nx: 30, ny: 30, hmax: 0.01 });
    check(`branch point as a Sym point: projection = ${B.projection} (on no 3-sphere)`, B.projection === 'balanced' ? 0 : 1, 0.5);
    // the dropped coordinate: x_4 for drop4; for the balanced projection, (pos, w4) is the stereographic image
    // in R^4, whose norm is |X| / (1 - X.p) in S^4, so |pos|^2 + w^2 is not identically 0 and w varies
    const D4 = computeGenSurface({ alphas: [[0.5, 0]], lams: [[0.5, 0], [0, 1]], signs: [1, 1], width: 4, height: 4, nx: 9, ny: 9, hmax: 0.01, projection: 'drop4' });
    check('drop4: w4 is the fourth coordinate', Math.max(...Array.from({ length: 81 }, (_, v) => Math.abs(D4.w4[v] - D4.pos4[4 * v + 3]))), 1e-6);
    let wmin = Infinity, wmax = -Infinity;
    for (const x of B.w4) { wmin = Math.min(wmin, x); wmax = Math.max(wmax, x); }
    check(`balanced: the dropped coordinate varies over the surface (${wmin.toFixed(3)} to ${wmax.toFixed(3)}), the Clifford torus on its 3-sphere has w = 0`,
      (wmax - wmin > 1e-3 ? 0 : 1) + Math.max(...Array.from(C.w4).map(Math.abs)), 1e-9);
  }
}

// ------------------------------------------------------------------ q and resonance pairs (resonance.js)
console.log('the map q = (log mu_1, log mu_i) and resonance pairs (cw_candidates.md, F3, F6):');
{
  const ab = ([r, i]) => Math.hypot(r, i);
  // 1. genus 0: q^1, q^2 are the eigenvalues of X = U + V and X = i(U - V) on the eigenline of (lam, nu)
  //    (the left eigenvector of zeta(lam) for nu/lam, as for the generalized Sym points)
  {
    const Q = new QMap([]), Z = killingField([]), G = new Float64Array(24);
    let worst = 0;
    for (const lam of [[0.3, 0.4], [-1.2, 0.7], [0.05, -2]]) for (const nu of Q.nus(lam)) {
      const u = leftEigenvector(Z, lam, nu), q = Q.q(lam, nu);
      [[1, 0], [0, 1]].forEach(([cr, ci], j) => {
        generator(Z, 2, 0, cr, ci, G);
        const d2 = lam[0] ** 2 + lam[1] ** 2, il = [lam[0] / d2, -lam[1] / d2];
        const X = [0, 1, 2, 3].map((e) => {
          const a = [G[2 * e], G[2 * e + 1]], b = [G[8 + 2 * e], G[9 + 2 * e]], c = [G[16 + 2 * e], G[17 + 2 * e]];
          return [a[0] * il[0] - a[1] * il[1] + b[0] + c[0] * lam[0] - c[1] * lam[1], a[0] * il[1] + a[1] * il[0] + b[1] + c[0] * lam[1] + c[1] * lam[0]];
        });
        // (u^T X)_1 / u_1, u = (u1, u2) complex
        const u1 = [u[0], u[1]], u2 = [u[2], u[3]];
        const mul = (x, y) => [x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0]];
        const big = ab(u1) > ab(u2);
        const num = big ? [mul(u1, X[0])[0] + mul(u2, X[2])[0], mul(u1, X[0])[1] + mul(u2, X[2])[1]]
          : [mul(u1, X[1])[0] + mul(u2, X[3])[0], mul(u1, X[1])[1] + mul(u2, X[3])[1]];
        const den = big ? u1 : u2, dd = den[0] ** 2 + den[1] ** 2;
        const m = [(num[0] * den[0] + num[1] * den[1]) / dd, (num[1] * den[0] - num[0] * den[1]) / dd];
        worst = Math.max(worst, ab([m[0] - q[j][0], m[1] - q[j][1]]));
      });
    }
    check('genus 0: q = eigenvalues of X(lam) on the eigenline (both sheets)', worst, 1e-13);
  }
  // 2. genus 0 closed forms: two solutions per gamma, over S^1 iff |k| <= 1/pi; Newton from a perturbed start
  {
    const Q = new QMap([]);
    for (const k of [[0.2, 0.1], [0.31, 0.05], [0.5, -0.3]]) {
      const gam = [[0, 2 * Math.PI * k[0]], [0, 2 * Math.PI * k[1]]];
      const inside = Math.hypot(...k) * Math.PI < 1;
      let worst = 0, circ = 0, newton = 0;
      for (const sol of genus0Resonance(k)) {
        const q1 = Q.q(sol.P.lam, sol.P.nu), q2 = Q.q(sol.Q.lam, sol.Q.nu);
        worst = Math.max(worst, ...[0, 1].map((j) => ab([q1[j][0] - q2[j][0] - gam[j][0], q1[j][1] - q2[j][1] - gam[j][1]])));
        // over S^1 (|nu| = 1) inside the threshold, a pair with |nu_P| |nu_Q| = 1 off it
        circ = Math.max(circ, inside ? Math.abs(ab(sol.P.nu) - 1) + Math.abs(ab(sol.Q.nu) - 1)
          : Math.abs(ab(sol.P.nu) * ab(sol.Q.nu) - 1) + (Math.abs(ab(sol.P.nu) - 1) < 1e-3 ? 1 : 0));
        const pert = (p) => {
          const lam = [p.lam[0] * 1.05 + 0.02, p.lam[1] * 0.97], ns = Q.nus(lam);
          return { lam, nu: ab([ns[0][0] - p.nu[0], ns[0][1] - p.nu[1]]) < ab([ns[1][0] - p.nu[0], ns[1][1] - p.nu[1]]) ? ns[0] : ns[1] };
        };
        const R = Q.resonance(pert(sol.P), pert(sol.Q), gam);
        newton = Math.max(newton, R.ok ? ab([R.P.lam[0] - sol.P.lam[0], R.P.lam[1] - sol.P.lam[1]]) : 1);
      }
      const tag = `k = (${k}), pi|k| = ${(Math.PI * Math.hypot(...k)).toFixed(3)}`;
      check(`genus 0, ${tag}: q(P) - q(Q) = gamma for both closed-form solutions`, worst, 1e-13);
      check(`genus 0, ${tag}: ${inside ? 'over S^1' : 'a rho-pair off S^1'}`, circ, 1e-12);
      check(`genus 0, ${tag}: Newton from a perturbed start finds it`, newton, 1e-10);
    }
    // at the threshold the two solutions coincide in a sigma-node: nu_Q = -nu_P
    const kt = [Math.cos(0.7) / Math.PI, Math.sin(0.7) / Math.PI];
    const [s1] = genus0Resonance(kt);
    const kinds = [[0.2, 0.1], [0.5, -0.3]].map((k) => classifyPair(genus0Resonance(k)[0].P, genus0Resonance(k)[0].Q).kind);
    check(`genus 0: classifyPair gives on, off, node (got ${kinds.join(', ')}, ${classifyPair(s1.P, s1.Q, 1e-6).kind})`,
      kinds.join() === 'on,off' && classifyPair(s1.P, s1.Q, 1e-6).kind === 'node' ? 0 : 1, 0.5);
    check('genus 0, pi|k| = 1: a sigma-node over S^1 (nu_Q = -nu_P, |nu| = 1)', ab([s1.P.nu[0] + s1.Q.nu[0], s1.P.nu[1] + s1.Q.nu[1]]) + Math.abs(ab(s1.P.nu) - 1), 1e-7);
  }
  // 3. branch points (F6): Re q = 0, q(alpha) = q(1/conj alpha), and 2 q(alpha_2) is the period of a cycle B_1 ± B_2
  for (const al of [[[0.5, 0.2]], [[0.4, 0.3], [-0.2, 0.5]]]) {
    const Q = new QMap(al), g = al.length;
    const qs = Q.branch.map((e) => Q.q(e, [0, 0]));
    check(`g=${g}: Re q = 0 at the branch points`, Math.max(...qs.flat().map(([re]) => Math.abs(re))), 1e-11);
    check(`g=${g}: q(alpha_j) = q(1/conj alpha_j)`, Math.max(...al.map((_, j) => Math.max(...[0, 1].map((c) => ab([qs[2 * j][c][0] - qs[2 * j + 1][c][0], qs[2 * j][c][1] - qs[2 * j + 1][c][1]]))))), 1e-10);
    if (g === 2) {
      const P = basicPeriods(al), per = (cyc, w) => {
        const p = thetaPoly(al, w, P);
        return p.reduce((acc, c, k) => [acc[0] + c[0] * P[cyc][k][0] - c[1] * P[cyc][k][1], acc[1] + c[0] * P[cyc][k][1] + c[1] * P[cyc][k][0]], [0, 0]);
      };
      let best = Infinity;
      for (const s1 of [1, -1]) for (const s2 of [1, -1]) {
        const d = [[1, 0], [0, 1]].map((w, c) => {
          const b1 = per(0, w), b2 = per(2, w);
          return ab([2 * qs[2][c][0] - s1 * b1[0] - s2 * b2[0], 2 * qs[2][c][1] - s1 * b1[1] - s2 * b2[1]]);
        });
        best = Math.min(best, Math.max(...d));
      }
      check('g=2: 2 q(alpha_2) = period of B_1 ± B_2 (so q(B) is a half-period)', best, 1e-9);
    }
  }
}

// ------------------------------------------------------------------ the tree (tree.js)
console.log('the tree: genus 0 edge, opening at lam = -1, genus 1 edge, candidate 1 (internal tau = 0.2 + 1.3i):');
{
  const tau = [0.2, 1.3];
  const g0 = genus0Edge(tau);
  check('genus 0: the pair for 1 arrives first, at t = 1/pi, at the sigma-node lam = -1',
    Math.abs(g0.arrivals[0].t - 1 / Math.PI) + Math.hypot(g0.arrivals[0].lam[0] + 1, g0.arrivals[0].lam[1]) + (g0.arrivals[0].name === 'one' ? 0 : 1), 1e-12);
  const [w1, w2] = torusBasis(0.3, tau), k1 = kOf(0.3, [1, 0]), k2 = kOf(0.3, tau);
  check('torus basis dual to (k(1), k(tau))', Math.max(Math.abs(w1[0] * k1[0] + w1[1] * k1[1] - 1), Math.abs(w1[0] * k2[0] + w1[1] * k2[1]),
    Math.abs(Math.abs(w2[0] * k2[0] + w2[1] * k2[1]) - 1), Math.abs(w2[0] * k1[0] + w2[1] * k1[1])), 1e-12);
  check('genus 1: the scale t(r) tends to the starting scale 1/pi at the node (r = 1 - 1e-4)', Math.abs(genus1Scale(1 - 1e-4) - 1 / Math.PI), 1e-7);
  const C = candidate1(tau);
  check(`candidate 1: the ${C.event && C.event.name} pair arrives (r* = ${C.event && C.event.r.toFixed(6)}, expected 0.092161)`, C.event ? Math.abs(C.event.r - 0.0921606) : 1, 1e-6);
  check('candidate 1: the arriving sigma-node has monodromy ±I on the internal lattice', C.defects.node, 3e-5);
  check(`candidate 1: different characters ${JSON.stringify(C.chars)}, so a double cover`, C.doubled ? 0 : 1, 0.5);
  check('candidate 1: f closes on the doubled lattice', C.defects.closing, 3e-5);
  check(`candidate 1: W = ${C.W.toFixed(4)} = ${(C.W / Math.PI ** 2).toFixed(4)} pi^2, twice W on the internal lattice`, Math.abs(C.W - 2 * C.event.W), 1e-8);
  check(`candidate 1: 2 pi^2 < W < 8 pi, external tau = ${C.tauExt.map((x) => x.toFixed(4)).join(' + ')}i in F_1`,
    (C.W > 2 * Math.PI ** 2 && C.W < 8 * Math.PI ? 0 : 1) + (Math.abs(C.tauExt[0]) <= 0.5 && Math.hypot(...C.tauExt) >= 1 ? 0 : 1), 0.5);
  // near tau = i candidate 1 tends to the Clifford torus (W = 2 pi^2, external tau = i)
  for (const im of [1.01, 1.001]) {
    const D = candidate1([0, im]);
    check(`internal tau = ${im}i: W = ${(D.W / Math.PI ** 2).toFixed(4)} pi^2 -> 2 pi^2 (Clifford), external tau -> i`,
      Math.abs(D.W / Math.PI ** 2 - 2) + Math.hypot(D.tauExt[0], D.tauExt[1] - 1) > 30 * (im - 1) ? 1 : 0, 0.5);
  }
  // F_6: the classes keep their labels; at 0.3 + 0.8i (order 213) the class 2 pair is opened and the class 1 pair arrives
  {
    const D = candidate1([0.3, 0.8]);
    check(`internal tau = 0.3 + 0.8i in F_6: order ${D.frame.order}, arriving class ${D.event && D.event.cls}`, D.frame.order === '213' && D.event && D.event.cls === 1 ? 0 : 1, 0.5);
  }
  // candidate 2 at 0.45 + 0.95i: the genus 2 edge from the node of candidate 1, the remaining pair arrives
  {
    const D = candidate2([0.45, 0.95]);
    check('candidate 2 (0.45 + 0.95i): the genus 2 anchor has the internal period ratio', D.edge ? D.edge.anchor.residual : 1, 1e-10);
    check(`candidate 2: the remaining pair arrives (s = ${D.event && D.event.s.toFixed(4)})`, D.event ? 0 : 1, 0.5);
    if (D.event) {
      check('candidate 2: W is continuous at the opening (W_int at the anchor = W_int of candidate 1)', Math.abs(D.edge.samples[0].W - D.c1.event.W) / D.c1.event.W, 1e-6);
      check('candidate 2: all three choices of Sym points close (on doubled lattices)', Math.max(...D.choices.map((c) => c.defects.closing + c.defects.node + (c.doubled ? 0 : 1))), 3e-5);
      check(`candidate 2: the same W for every choice, ${(D.choices[0].W / Math.PI ** 2).toFixed(4)} pi^2 > candidate 1's ${(D.c1.W / Math.PI ** 2).toFixed(4)} pi^2`,
        Math.max(...D.choices.map((c) => Math.abs(c.W - D.choices[0].W))) + (D.choices[0].W > D.c1.W ? 0 : 1), 1e-9);
    }
  }
  // the open end: at 0.2 + 1.3i the genus 2 edge reaches it without an arrival; the open-end tau of the
  // rescaled roots tends to the internal tau (error O(t)), and W_int to 4 pi
  {
    const D = candidate2([0.2, 1.3], { shapes: true });
    const o = D.openEnd;
    check(`open end at 0.2 + 1.3i: reached (t = ${o && o.t.toExponential(0)}), open-end tau within O(t) of the internal tau`, o ? Math.max(0, o.tauError - o.t) : 1, 1e-9);
    check(`open end: W_int / 4 pi = ${o && o.WRatio.toFixed(6)} -> 1 (error O(t))`, o ? Math.max(0, Math.abs(o.WRatio - 1) - o.t) : 1, 1e-9);
    // the error is O(t): ten times smaller at ten times smaller t
    const E5 = genus2Edge(D.c1, { tMin: 1e-4 });
    check(`open end: the tau error scales with t (${o && o.tauError.toExponential(2)} at 1e-3, ${E5.openEnd && E5.openEnd.tauError.toExponential(2)} at 1e-4)`,
      o && E5.openEnd ? Math.abs(o.tauError / E5.openEnd.tauError - 10) : 1, 0.1);
    // the limit surface: candidate 2 at t = 0.03 ... 0.001 tends to a round sphere covered twice, read in
    // S^4 after a Möbius normalization (distance O(t^(1/2)), area -> 2 x 4 pi, W -> 8 pi)
    const sh = D.choices.map((c) => c.shape), last = D.choices[D.choices.length - 1];
    check(`open end: candidate 2 at t = ${D.choices.map((c) => c.t.toExponential(0)).join(', ')}: distance from a round sphere ${sh.map((s) => s.defect.toFixed(3)).join(', ')} (ratio sqrt(t))`,
      Math.max(...sh.slice(1).map((s, i) => Math.abs(sh[i].defect / s.defect / Math.sqrt(D.choices[i].t / D.choices[i + 1].t) - 1))), 0.1);
    check(`open end: S^4 area ${sh.map((s) => s.area.toFixed(3)).join(', ')} (x 4 pi) -> 2, below W / 4 pi`,
      Math.abs(sh[sh.length - 1].area - 2) + Math.max(0, ...D.choices.map((c) => c.shape.area - c.W / (4 * Math.PI))), 0.03);
    // Möbius-balanced projection: near the limit the surface is near the unit sphere of R^3
    const n = 24, [l1, l2] = last.lattice.map((e) => Math.hypot(...e));
    const S = computeGenSurface({ alphas: last.alphas, lams: last.syms.map((s) => s.lam), cell: last.lattice, projection: 'balanced',
      sCoords: Array.from({ length: n }, (_, i) => (i * l1) / (n - 1)), tCoords: Array.from({ length: n }, (_, i) => (i * l2) / (n - 1)) });
    let rad = 0;
    for (let v = 0; v < n * n; v++) rad = Math.max(rad, Math.abs(Math.hypot(S.pos[3 * v], S.pos[3 * v + 1], S.pos[3 * v + 2]) - 1));
    check(`open end: balanced projection at t = ${last.t.toExponential(0)} within ${rad.toFixed(3)} of the unit sphere`, rad, 0.15);
  }
  // the tree as a path (treeOf): W increases from near the vacuum to the candidates; the parity rule for the
  // closing lattice gives the candidates' own lattices; before the Sym points are ripe the immersion is a
  // plane with f(z + e) = c f(z), c real (a homothety of R^4 about 0), the same c at different z
  {
    const D = candidate2([0.45, 0.95]), T = treeOf(D);
    const S = T.samples.filter((x) => x.edge !== 3), B = T.samples.filter((x) => x.edge === 3);
    const c1s = T.samples[T.nodes.find((n) => n.label === 'candidate 1').index];
    check(`tree at 0.45 + 0.95i: ${S.length} samples, W increasing (${(S[0].W / Math.PI ** 2).toFixed(3)} to ${(S[S.length - 1].W / Math.PI ** 2).toFixed(3)} pi^2)`,
      S.slice(1).reduce((a, sm, i) => a + Math.max(0, S[i].W - sm.W), 0), 1e-9);
    check(`the branch at candidate 1: ${B.length} samples, W increasing (to ${(B[B.length - 1].W / Math.PI ** 2).toFixed(3)} pi^2)`,
      [c1s, ...B].slice(1).reduce((a, sm, i) => a + Math.max(0, [c1s, ...B][i].W - sm.W), 0), 1e-9);
    const sameLattice = (A, B) => {
      const det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
      const m = B.map((v) => [(v[0] * A[1][1] - v[1] * A[1][0]) / det, (A[0][0] * v[1] - A[0][1] * v[0]) / det]);
      const d = m[0][0] * m[1][1] - m[0][1] * m[1][0];
      return m.flat().every((x) => Math.abs(x - Math.round(x)) < 1e-6) && Math.abs(Math.abs(d) - 1) < 1e-6;
    };
    const c1 = D.c1, cls = { tau: [0, 1], third: [1, 1] }[c1.event.name];
    const basis2 = dualBasis(D.event.k[0], D.event.k[1]), ch = D.choices[0].chars, fl = (x) => (x < 0 ? 1 : 0);
    const par2 = { 'a1,a2': ch.map(([u, v]) => (u === v ? 0 : 1)), 'a1,node': ch.map(([, v]) => fl(v)), 'a2,node': ch.map(([u]) => fl(u)) };
    check('parity rule = the lattices of candidates 1, 2, 3 and the other choice',
      [sameLattice(c1.lattice, paritySublattice(torusBasis(c1.event.t, c1.frame.tau1), [(1 + cls[0]) & 1, cls[1] & 1]))]
        .concat(D.choices.map((c) => sameLattice(c.lattice, paritySublattice(basis2, par2[c.pair.join()])))).filter((x) => !x).length, 0.5);
    for (const e of [0, 1, 2]) {
      const list = S.filter((sm) => sm.edge === e && !sm.node), sm = list[Math.floor(list.length / 2)];
      const v = sm.views[e === 2 ? 1 : 0], [l1, l2] = v.lattice.map((w) => Math.hypot(...w));
      const G = computeGenSurface({ alphas: sm.alphas, lams: v.syms.map((p) => p.lam), signs: v.syms.map((p) => p.sign), cell: v.lattice, projection: 'drop4',
        sCoords: [0, 0.3 * l1, l1, 1.3 * l1], tCoords: [0, 0.4 * l2], hmax: 0.002 });
      const f = (i, j) => Array.from(G.pos4.slice(4 * (j * 4 + i), 4 * (j * 4 + i) + 4));
      const ratio = (a, b) => { const n = b.reduce((x, y) => x + y * y, 0), c = a.reduce((x, y, k) => x + y * b[k], 0) / n; return { c, r: Math.hypot(...a.map((x, k) => x - c * b[k])) / Math.hypot(...a) }; };
      const A = ratio(f(2, 0), f(0, 0)), B = ratio(f(3, 1), f(1, 1));
      check(`unripe (genus ${e}, ${v.key}): f(z + e1) = c f(z) with c = ${A.c.toFixed(4)} real, the same at z = 0.3 e1 + 0.4 e2 (${B.c.toFixed(4)})`,
        Math.max(A.r, B.r, Math.abs(A.c - B.c) / Math.abs(A.c)) + (Math.abs(A.c - 1) > 1e-3 ? 0 : 1), 1e-4);
    }
  }
  // a double arrival on the wall Re tau = 1/2 of F_6 (classes 2 and 3 equally long): the genus 2 edge has
  // length 0, and the three candidates are the three choices of two of { alpha, node a, node b } on the
  // genus 1 curve (same W); candidate 2 is candidate 1, and the candidates just off the wall tend to these
  {
    const D = candidate2([0.5, 0.9]), N = candidate2([0.499, 0.9]);
    const dW = D.double ? Math.max(...D.choices.map((c) => Math.abs(c.W - D.c1.W))) / D.c1.W : 1;
    check(`double arrival at 0.5 + 0.9i: three candidates on the genus 1 curve, W = ${(D.c1.W / Math.PI ** 2).toFixed(4)} pi^2 for all`, dW + (D.choices ? Math.max(...D.choices.map((c) => c.defects.closing)) : 1), 3e-5);
    const tE = (c) => c.tauExt;
    check('double arrival: candidate 2 = candidate 1 (external tau)', Math.hypot(...tE(D.choices[0]).map((x, i) => x - tE(D.c1)[i])), 1e-9);
    check(`just off the wall (0.499 + 0.9i): W and the external taus of the three choices within 4e-3 of the double arrival`,
      Math.max(Math.abs(N.choices[0].W - D.c1.W) / D.c1.W, ...N.choices.map((c, i) => Math.hypot(...c.tauExt.map((x, k) => x - D.choices[i].tauExt[k])))), 4e-3);
  }
  // the branch at candidate 1 that does not open the node: the genus 1 edge continued to the next
  // sigma-node event (the other pair arriving, or the arrived pair leaving S^1); candidate 1' there closes,
  // and when the other pair arrives its classes are those of candidate 3 ([1] + b), so it has the same
  // external tau (the parity rule)
  {
    const D = candidate2([0.45, 0.95]), b = D.c1.branch;
    check(`branch at 0.45 + 0.95i: the ${b && b.event && b.event.name} pair arrives at alpha = -${b && b.event && b.event.r.toFixed(6)}, candidate 1' W = ${b && (b.W / Math.PI ** 2).toFixed(4)} pi^2, closes`,
      b && b.event && !b.event.leaving ? b.defects.closing + b.defects.node : 1, 3e-5);
    check("candidate 1' and candidate 3 have the same external tau (same classes)", Math.hypot(...b.tauExt.map((x, i) => x - D.choices[1].tauExt[i])), 1e-9);
    check('candidate 1 and candidate 2 have the same external tau (same classes)', Math.hypot(...D.c1.tauExt.map((x, i) => x - D.choices[0].tauExt[i])), 1e-9);
    const E = candidate2([0.2, 1.3], { fast: true }), e = E.c1.branch;
    check(`branch at 0.2 + 1.3i: the arrived pair leaves S^1 again at alpha = -${e && e.event && e.event.r.toFixed(6)}, candidate 1' closes`,
      e && e.event && e.event.leaving ? e.defects.closing + e.defects.node : 1, 3e-5);
    check('no branch at a double arrival (0.5 + 0.9i)', candidate1([0.5, 0.9]).branch === null ? 0 : 1, 0.5);
  }
  // dragging a branch point in genus 2: locateOnTree gives back the internal tau and scale of every sample
  // of the genus 2 edge, and a perturbed curve lies on the edge of its own tau at its scale (the period
  // ratio of the edge's curve at that t equals that of the perturbed curve)
  {
    const tau = [0.3, 0.8], D = candidate2(tau, { fast: true }), T = treeOf(D), e2 = T.samples.filter((x) => x.edge === 2);
    const back = Math.max(...e2.map((x) => { const L = locateOnTree(x.alphas, tau); return Math.hypot(L.tau[0] - tau[0], L.tau[1] - tau[1]) + Math.abs(L.t - x.t) / x.t; }));
    check(`locateOnTree at 0.3 + 0.8i: the ${e2.length} genus 2 samples give back tau and t`, back, 1e-9);
    const sm = e2[Math.floor(e2.length / 2)], al = sm.alphas.map(([x, y], j) => (j ? [x * 1.03, y + 0.01] : [x, y]));
    const L = locateOnTree(al, tau), E = candidate2(L.tau, { fast: true }), f2 = treeOf(E).samples.filter((x) => x.edge === 2);
    // the edge's |alpha_1|, |alpha_2| interpolated in log t at the curve's scale (rotation is a symmetry)
    const k = f2.findIndex((x, q) => q + 1 < f2.length && x.t >= L.t && f2[q + 1].t <= L.t);
    const [A, B] = [f2[k], f2[k + 1]], u = Math.log(A.t / L.t) / Math.log(A.t / B.t);
    const lerp = (j) => Math.exp((1 - u) * Math.log(Math.hypot(...A.alphas[j])) + u * Math.log(Math.hypot(...B.alphas[j])));
    const err = Math.max(...[0, 1].map((j) => Math.abs(lerp(j) / Math.hypot(...al[j]) - 1)));
    check(`a perturbed curve (tau -> ${L.tau.map((x) => x.toFixed(4)).join(' + ')}i, t = ${L.t.toFixed(4)}) lies on the genus 2 edge of its tau: |alpha_j| to ${err.toExponential(1)}`, err, 1e-3);
  }
  check('reduceTau maps 2.3 + 0.4i into F_1', (() => { const [x, y] = reduceTau([2.3, 0.4]); return Math.abs(x) <= 0.5 && x * x + y * y >= 1 - 1e-12 ? 0 : 1; })(), 0.5);
}

// ------------------------------------------------------------------ timing
console.log('the candidate family of a conformal class (familyOf; the proposal, Section 1.6, and K-PH-S):');
{
  // a rectangular class: the member with the least W is the homogeneous torus, W = pi^2 (y + 1/y)
  const y = 1.5, F = familyOf([0, y]), ws = F.members.map((m) => (m.member ? m.member.W : Infinity));
  check(`tau^ = ${y}i: least W = ${(Math.min(...ws) / Math.PI ** 2).toFixed(6)} pi^2 is that of the homogeneous torus, pi^2 (y + 1/y)`, Math.abs(Math.min(...ws) / (Math.PI ** 2 * (y + 1 / y)) - 1), 1e-5);
  // three members of genus 1 below 8 pi, one over each superlattice, each closing on the class
  const G = familyOf([0.45, 1.2]);
  check('tau^ = 0.45 + 1.2i: three members of genus 1 with W < 8 pi', 3 - G.members.filter((m) => m.member && m.member.genus === 1 && m.member.W < 8 * Math.PI).length, 0);
  check('each member closes on <1, tau^>', Math.max(...G.members.map((m) => (m.member && sameClass(m.member.tauExt, [0.45, 1.2]) ? 0 : 1))), 0);
  // the superlattices contain <1, tau^> with index 2: their covolume is half
  const sl = superlattices([0.45, 1.2]);
  check('the three superlattices are distinct classes', sl.some((a, i) => sl.some((b, j) => i < j && sameClass(a.tau, b.tau))) ? 1 : 0, 0);
  // the branched double cover: on the unit sphere and even, f(-z) = f(z)
  const D = doubleCoverSurface([0.2, 1.3], { n: 40 }), N = D.nx * D.ny;
  let off = 0, odd = 0;
  for (let v = 0; v < N; v++) {
    off = Math.max(off, Math.abs(Math.hypot(D.pos[3 * v], D.pos[3 * v + 1], D.pos[3 * v + 2]) - 1));
    const w = N - 1 - v;
    for (let e = 0; e < 3; e++) odd = Math.max(odd, Math.abs(D.pos[3 * v + e] - D.pos[3 * w + e]));
  }
  check('double cover: on the unit sphere', off, 1e-6);
  check('double cover: f(-z) = f(z) on the cell centred at 0', odd, 1e-6);
}

console.log('timing:');
for (const [g, n] of [[2, 64], [2, 256], [4, 256]]) {
  const alphas = [[0.6, 0], [0, 0.4], [-0.5, -0.3], [0.2, -0.7]].slice(0, g);
  const S = computeSurface({ alphas, theta0: 1, tau: [0.3, 0.2], width: 20, height: 20, nx: n, ny: n, hmax: 0.05 });
  console.log(`  g=${g} ${n}x${n} over 20x20, hmax=0.05: ${S.stats.ms.toFixed(0)} ms`);
}

console.log(failures ? `\n${failures} FAILED` : '\nall passed');
process.exit(failures ? 1 : 0);
