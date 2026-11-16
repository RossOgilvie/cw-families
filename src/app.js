// Glue: state, panel controls, progressive recomputation in a worker, presets, URL state, export.

import { Viewer } from './viewer.js';
import { SpectralWidget, MAX_GENUS } from './spectral-widget.js';
import { shapeFlows, kappa0, hopfArg } from './cmc/cmc.js';
import { commonRootsOnCircle } from './cmc/periods.js';
import { TauPad } from './tau-pad.js';
import { TreeView } from './tree-view.js';

const $ = (id) => document.getElementById(id);
const polar = (r, t) => [r * Math.cos(t), r * Math.sin(t)];
const PI = Math.PI;

// ------------------------------------------------------------------ panel furniture

// borderless stroke glyphs for buttons
const GLYPH = {
  chev: '<svg width="10" height="10" viewBox="0 0 10 10" fill="none" stroke="currentColor" stroke-width="1.4" stroke-linecap="round" stroke-linejoin="round"><path d="M3.5 2l3 3-3 3"/></svg>',
  play: '<svg width="14" height="14" viewBox="0 0 14 14"><path d="M4 2.5 L11.5 7 L4 11.5 Z" fill="currentColor"/></svg>',
  stop: '<svg width="14" height="14" viewBox="0 0 14 14"><rect x="3.5" y="3.5" width="7" height="7" fill="currentColor"/></svg>',
  reset: '<svg width="15" height="15" viewBox="0 0 15 15" fill="none" stroke="currentColor" stroke-width="1.4" stroke-linecap="round" stroke-linejoin="round"><path d="M12 7.5a4.5 4.5 0 1 1-1.3-3.2"/><path d="M11.2 1.8v2.7H8.5"/></svg>',
  plus: '<svg width="12" height="12" viewBox="0 0 12 12" fill="none" stroke="currentColor" stroke-width="1.4" stroke-linecap="round"><path d="M6 2v8M2 6h8"/></svg>',
  minus: '<svg width="12" height="12" viewBox="0 0 12 12" fill="none" stroke="currentColor" stroke-width="1.4" stroke-linecap="round"><path d="M2 6h8"/></svg>',
  cross: '<svg width="11" height="11" viewBox="0 0 12 12" fill="none" stroke="currentColor" stroke-width="1.4" stroke-linecap="round"><path d="M2.5 2.5l7 7M9.5 2.5l-7 7"/></svg>',
};
const glyphs = (root) => root.querySelectorAll('[data-glyph]').forEach((el) => { el.innerHTML = GLYPH[el.dataset.glyph]; });
glyphs(document);

// every section (and the help and branch-point list) collapses; the choice is remembered per browser
{
  let open = {};
  try { open = JSON.parse(localStorage.getItem('cmc-open')) || {}; } catch { /* no storage */ }
  document.querySelectorAll('.disc[aria-controls]').forEach((b) => {
    const target = $(b.getAttribute('aria-controls'));
    const setOpen = (v) => { b.setAttribute('aria-expanded', String(v)); target.hidden = !v; };
    const id = target.id;
    setOpen(id in open ? open[id] : b.getAttribute('aria-expanded') === 'true');
    b.addEventListener('click', () => {
      const v = b.getAttribute('aria-expanded') !== 'true';
      setOpen(v);
      open[id] = v;
      try { localStorage.setItem('cmc-open', JSON.stringify(open)); } catch { /* no storage */ }
      if (v) widget.draw(); // the canvas may have been laid out while hidden
    });
  });
}

/** Fixed decimals, a true minus sign, and no "−0.000". */
function fixed(x, d) {
  const t = Math.abs(x).toFixed(d);
  return (x < 0 && +t !== 0 ? '−' : '') + t;
}
/** Number typed into a field: accepts − and a trailing π (ignored: angle fields are already in π). */
function parseNum(str) {
  const t = String(str).replace(/−/g, '-').replace(/π/g, '').trim();
  if (t === '') return null;
  const v = Number(t);
  return Number.isFinite(v) ? v : null;
}

// ------------------------------------------------------------------ state

const VIEW_DEFAULTS = {
  tau: [], z0: [0, 0], width: 12, height: 12, curv: true, phi: 0,
  gridStep: 0.25, lineW: 0.6, hmax: 0.05, adaptive: false, res: 240,
  // generalized Sym points (cw_candidates.md, F7): the immersion into R^4 from two points of the
  // spectral curve, projected to R^3
  gen: { on: false, lams: [[0, 1], [0, -1]], signs: [1, 1], proj: 'auto' },
  // a parallelogram cell [e1, e2] of the closing lattice of a candidate (drawn instead of the domain)
  cell: null,
};
// Widths of the rotational examples are exact closing periods of the parallels, pi / sqrt(det X(lam0)).
// Bubbletons: a double branch point at the k-th resonance point of the cylinder, where its monodromy is ±I.
const RES = (k) => (k - Math.sqrt(k * k - 1)) ** 2;
const PRESETS = [
  { name: 'Round cylinder (g = 0)', s: { alphas: [], theta0: 0, width: PI, height: 10, gridStep: PI / 16 } },
  { name: 'Unduloid (g = 1)', s: { alphas: [polar(0.45, 0)], theta0: PI, width: 2.90682, height: 14 } },
  { name: 'Unduloid near a sphere chain', s: { alphas: [polar(0.08, 0)], theta0: PI, width: 1.645512, height: 14, gridStep: 0.1 } },
  { name: 'Nodoid (g = 1)', s: { alphas: [polar(0.45, 0)], theta0: 0, width: 10, height: 7.663434 } },
  { name: 'Twizzler (g = 1)', s: { alphas: [polar(0.45, 0)], theta0: PI / 2, width: 8, height: 8 } },
  { name: 'Bent tube (g = 2)', s: { alphas: [polar(0.49, 1), polar(0.49, -1)], theta0: PI / 2, width: 6, height: 6 } },
  { name: 'Nearly a cylinder (g = 2)', s: { alphas: [polar(0.85, 0), polar(0.85, PI / 2)], theta0: 0.7, width: 10, height: 10 } },
  { name: 'Bubbleton, two lobes', s: { alphas: [polar(RES(2), 0), polar(RES(2), 0)], theta0: 0, width: PI, height: 3.4, gridStep: PI / 48, hmax: 0.02, res: 360 } },
  { name: 'Bubbleton, three lobes', s: { alphas: [polar(RES(3), 0), polar(RES(3), 0)], theta0: 0, width: PI, height: 1.6, gridStep: PI / 60, hmax: 0.02, res: 512 } },
  { name: 'Bubbleton, four lobes', s: { alphas: [polar(RES(4), 0), polar(RES(4), 0)], theta0: 0, width: PI, height: 1.2, gridStep: PI / 80, hmax: 0.02, res: 512 } },
  // Wente torus: McIntosh's spectral data (arXiv math/0407248, Remark 1: a = 0.1413 ± 0.1018i), refined
  // by Newton on the closing conditions. The period lattice of zeta is rhombic, <w1, w2> with w1 = 2.31698
  // real and 2 w2 - w1 = 4.47988 i; E(w1) turns by 2/3 about an axis, E(2 w2 - w1) = id, so the torus
  // lattice is <3 w1, w1 + w2>, whose fundamental domain is the brick 3 w1 x (4.47988 / 2).
  { name: 'Wente torus (three lobes)', s: { alphas: [[0.1412634686, 0.1017768953], [0.1412634686, -0.1017768953]], theta0: 0, width: 6.950948617, height: 2.239937825, gridStep: 6.950948617 / 72, res: 360 } },
  { name: 'Clifford torus (two Sym points ±i, g = 0)', s: { alphas: [], theta0: 0, width: PI * Math.SQRT2, height: PI * Math.SQRT2, gridStep: PI * Math.SQRT2 / 24, gen: { on: true, lams: [[0, 1], [0, -1]], signs: [1, 1], proj: 'auto' } } },
  { name: 'Homogeneous torus, ratio 2 (two Sym points, g = 0)', s: { alphas: [], theta0: 0, width: PI / Math.cos(Math.atan(0.5)), height: PI / Math.sin(Math.atan(0.5)), gridStep: PI / Math.cos(Math.atan(0.5)) / 16, gen: { on: true, lams: [[0.6, 0.8], [0.6, -0.8]], signs: [1, 1], proj: 'auto' } } },
  { name: 'Branch point as a Sym point (g = 1, in R⁴)', s: { alphas: [polar(0.5, 0)], theta0: 0, width: 6, height: 6, gen: { on: true, lams: [polar(0.5, 0), [0, 1]], signs: [1, 1], proj: 'auto' } } },
  { name: 'Threefold (g = 3)', s: { alphas: [polar(0.4, 0), polar(0.4, 2 * PI / 3), polar(0.4, -2 * PI / 3)], theta0: 0.3, tau: [0], width: 6, height: 6 } },
  { name: 'Genus 4', s: { alphas: [polar(0.54, 0.38), polar(0.4, PI / 2), polar(0.58, -2.6), polar(0.63, -1.25)], theta0: 1, tau: [0.8, -0.5], width: 6, height: 6 } },
];

const state = {
  alphas: [polar(0.49, 1), polar(0.49, -1)], theta0: 0,
  ...VIEW_DEFAULTS, gen: structuredClone(VIEW_DEFAULTS.gen), width: 6, height: 6,
  res: 240, hmax: 0.05, colour: 'side', grid: true, transp: 0, H: 0.5, adapt: true,
  front: '#6f8fb0', back: '#e8c9a0', bg: '#f2eee3',
};

const genus = () => state.alphas.length;
function syncTau() {
  const n = shapeFlows(genus()).length;
  state.tau = Array.from({ length: n }, (_, i) => state.tau[i] || 0);
}

// ------------------------------------------------------------------ links

let pendingCandidate = null; // internal tau from a link #k=re,im (candidate 1) or #k=re,im,2 (candidate 2), computed on start
// The state goes into the URL only through "Copy link". A page opened from such a link (or with one
// pasted into its address bar) takes the state from the hash and then clears it from the address bar.

const fmt = (x, d = 4) => String(+x.toFixed(d));
function linkURL() {
  const p = new URLSearchParams();
  // a candidate is recomputed from its internal tau
  if (state.cell && cand.tau) return `${location.origin}${location.pathname}${location.search}#k=${cand.tau.map((x) => fmt(x, 10)).join(',')}${cand.which === 2 ? ',2' : ''}`;
  p.set('a', state.alphas.map(([re, im]) => `${fmt(re, 10)},${fmt(im, 10)}`).join(';'));
  p.set('l', fmt(state.theta0, 10));
  if (state.tau.length) p.set('t', state.tau.map((x) => fmt(x)).join(','));
  if (state.gen.on && state.gen.proj !== 'auto') p.set('j', state.gen.proj);
  if (state.gen.on) p.set('y', state.gen.lams.map(([re, im], i) => `${fmt(re, 10)},${fmt(im, 10)},${state.gen.signs[i]}`).join(';'));
  p.set('c', state.z0.map((x) => fmt(x)).join(','));
  p.set('d', `${fmt(state.width, 8)},${fmt(state.height, 8)}`);
  if (!state.curv) p.set('p', fmt(state.phi));
  if (state.colour !== 'side') p.set('m', state.colour);
  if (state.transp) p.set('o', fmt(state.transp, 2));
  p.set('g', fmt(state.gridStep));
  if (state.adaptive) p.set('h', 'a');
  else if (state.hmax !== 0.05) p.set('h', fmt(state.hmax));
  if (state.res !== 240) p.set('r', state.res);
  // commas and semicolons are safe in a fragment; keep the link readable
  return `${location.origin}${location.pathname}${location.search}#${p.toString().replace(/%2C/g, ',').replace(/%3B/g, ';')}`;
}
const clearHash = () => history.replaceState(null, '', location.pathname + location.search);
// Links written before lam0 kept full precision (4 decimals): a lam0 that close to a common root of the
// differentials was on it, so put it back.
function snapLinkLam() {
  if (!state.alphas.length) return;
  let roots = [];
  try { roots = commonRootsOnCircle(state.alphas, 1e-5).map((r) => Math.atan2(r.lam[1], r.lam[0])); } catch { return; }
  for (const t of roots) if (Math.abs(wrapPi(t - state.theta0)) < 1e-4) { state.theta0 = t; return; }
}
function readHash() {
  const p = new URLSearchParams(location.hash.slice(1));
  if (p.has('k')) {
    const k = (p.get('k') || '').split(',').map(Number);
    if ((k.length === 2 || k.length === 3) && k.every(Number.isFinite) && k[1] > 0) pendingCandidate = k;
    return false;
  }
  if (!p.has('a')) return false;
  const nums = (s) => (s || '').split(',').filter((x) => x !== '').map(Number);
  try {
    const a = p.get('a');
    state.alphas = a ? a.split(';').map((s) => nums(s)).filter((v) => v.length === 2 && v.every(Number.isFinite)) : [];
    state.alphas = state.alphas.slice(0, MAX_GENUS);
    if (p.has('l')) state.theta0 = +p.get('l') || 0;
    state.tau = nums(p.get('t'));
    state.gen.on = false;
    state.cell = null;
    if (p.has('y')) {
      const ys = p.get('y').split(';').map(nums);
      if (ys.length === 2 && ys.every((v) => v.length === 3 && v.every(Number.isFinite))) {
        state.gen = { ...state.gen, proj: ['drop4', 'stereo4', 'balanced', 'pca'].includes(p.get('j')) ? p.get('j') : 'auto', on: true, lams: ys.map((v) => [v[0], v[1]]), signs: ys.map((v) => (v[2] < 0 ? -1 : 1)) };
      }
    }
    if (p.has('c')) { const c = nums(p.get('c')); if (c.length === 2) state.z0 = c; }
    if (p.has('d')) { const d = nums(p.get('d')); if (d.length === 2 && d.every((x) => x > 0)) [state.width, state.height] = d; }
    state.curv = !p.has('p');
    if (p.has('p')) state.phi = +p.get('p') || 0;
    if (p.has('m')) state.colour = p.get('m');
    state.transp = Math.min(1, Math.max(0, +p.get('o') || 0));
    if (p.has('g')) state.gridStep = +p.get('g') || state.gridStep;
    state.adaptive = p.get('h') === 'a';
    if (p.has('h')) state.hmax = state.adaptive ? 0.05 : +p.get('h') || 0.05;
    if (p.has('r')) state.res = Math.min(1024, Math.max(16, +p.get('r') || 240));
  } catch (e) {
    console.warn('bad URL state', e);
  }
  return true;
}

// ------------------------------------------------------------------ computation scheduling

const viewer = new Viewer($('view'));
let worker = null, busy = false, busyFull = false, pending = null, reqId = 0, shownId = 0;
let fullTimer = null, needFrame = true, last = null;
const jobInfo = new Map(); // job id -> what the job computed

function makeWorker() {
  worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
  worker.onmessage = onResult;
  worker.onerror = (e) => showStatus(`worker error: ${e.message}`, true);
}

function gridDims(res) {
  const m = Math.max(state.width, state.height);
  return [Math.max(2, Math.round((res * state.width) / m)), Math.max(2, Math.round((res * state.height) / m))];
}

// everything that determines the surface (not the mesh or the accuracy)
const surfaceKey = () => JSON.stringify([state.alphas, state.theta0, state.tau, state.z0,
  state.curv ? null : state.phi, state.width, state.height, state.gen.on ? state.gen : null, state.cell]);

function params(full) {
  const [nx, ny] = gridDims(full ? state.res : 72);
  const p = {
    alphas: state.alphas, theta0: state.theta0, tau: state.tau, z0: state.z0,
    phi: state.curv ? null : state.phi, width: state.width, height: state.height,
    nx, ny, hmax: full ? state.hmax : Math.max(state.hmax, 0.1),
  };
  if (state.gen.on) p.gen = { lams: state.gen.lams, signs: state.gen.signs, projection: state.gen.proj };
  if (state.gen.on && state.cell) {
    // one cell of the closing lattice of a candidate, uniformly sampled (the first and last rows and
    // columns coincide on the torus)
    const [e1, e2] = state.cell, l1 = Math.hypot(...e1), l2 = Math.hypot(...e2), m = Math.max(l1, l2);
    const res = full ? state.res : 72;
    const n1 = Math.max(3, Math.round((res * l1) / m)), n2 = Math.max(3, Math.round((res * l2) / m));
    Object.assign(p, { cell: state.cell, z0: [0, 0], phi: 0,
      sCoords: Array.from({ length: n1 }, (_, i) => (i * l1) / (n1 - 1)), tCoords: Array.from({ length: n2 }, (_, j) => (j * l2) / (n2 - 1)) });
    return p;
  }
  // error-controlled steps for the full mesh: on request, or when the preview shows the conformal factor
  // bubbling (near the corner of the spectral data where the branch points go to 0, fixed steps lose it)
  // (the latest result stands in for the preview while the data move, e.g. along a flow)
  if (full && (state.adaptive || (state.hmax >= 0.05 && last && uSpread(last) > 7))) p.tol = 1e-9;
  // full meshes: put the rows and columns where the surface is, using the preview's metric
  if (full && state.adapt && last && last.key === surfaceKey() && last.preview) {
    p.sCoords = equidistribute(last, 's', nx);
    p.tCoords = equidistribute(last, 't', ny);
  }
  return p;
}

function uSpread(res) {
  let lo = Infinity, hi = -Infinity;
  for (const v of res.u) if (Number.isFinite(v)) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
  return hi - lo;
}

/**
 * Grid coordinates (relative to z0) along one axis of a result, spaced so that each interval carries
 * about the same length on the surface: density max e^u across the other axis, plus a floor.
 */
function equidistribute(res, axis, N) {
  const { nx, ny, u } = res;
  const n = axis === 's' ? nx : ny, m = axis === 's' ? ny : nx;
  const L = axis === 's' ? state.width : state.height;
  const w = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    let mx = 0;
    for (let j = 0; j < m; j++) {
      const v = axis === 's' ? u[j * nx + i] : u[i * nx + j];
      if (Number.isFinite(v)) mx = Math.max(mx, Math.exp(v));
    }
    w[i] = mx;
  }
  const mean = w.reduce((a, b) => a + b, 0) / n || 1;
  const x = (i) => -L / 2 + (n > 1 ? (i * L) / (n - 1) : 0);
  const W = new Float64Array(n); // cumulative length (trapezoid)
  for (let i = 1; i < n; i++) W[i] = W[i - 1] + 0.5 * (w[i - 1] + w[i] + 0.5 * mean) * (x(i) - x(i - 1));
  const out = new Float64Array(N);
  let k = 0;
  for (let q = 0; q < N; q++) {
    const target = N > 1 ? (q * W[n - 1]) / (N - 1) : 0;
    while (k < n - 2 && W[k + 1] < target) k++;
    const f = W[k + 1] > W[k] ? (target - W[k]) / (W[k + 1] - W[k]) : 0;
    out[q] = x(k) + Math.min(1, Math.max(0, f)) * (x(k + 1) - x(k));
  }
  out[0] = -L / 2; out[N - 1] = L / 2;
  return out;
}

function request(full) {
  pending = { id: ++reqId, full };
  if (busy && busyFull && !full) { // a slow job is in the way of interactive feedback: drop it
    worker.terminate();
    makeWorker();
    busy = false;
  }
  pump();
}

function pump() {
  if (busy || !pending) return;
  const job = pending;
  pending = null;
  busy = true;
  busyFull = job.full;
  $('busy').classList.add('on');
  const p = params(job.full);
  jobInfo.set(job.id, { key: surfaceKey(), preview: !job.full, tol: p.tol || 0 });
  worker.postMessage({ id: job.id, params: p });
}

function onResult(e) {
  busy = false;
  const { id, result, error } = e.data;
  if (error) showStatus(error, true);
  else if (id > shownId) { shownId = id; Object.assign(result, jobInfo.get(id)); show(result); }
  jobInfo.delete(id);
  if (!pending) $('busy').classList.remove('on');
  pump();
}

/** Called after every state change. dragging: more changes are coming, so only preview. */
function changed(dragging = false) {
  updateCommonRoots();
  request(false);
  clearTimeout(fullTimer);
  fullTimer = setTimeout(() => request(true), dragging || anim ? 350 : 0);
}

// ------------------------------------------------------------------ display

function colourValues(res) {
  if (state.colour === 'side' || state.colour === 'wire') return { vals: null, range: [0, 1] };
  const N = res.u.length;
  const vals = new Float32Array(N);
  for (let v = 0; v < N; v++) {
    const u = res.u[v];
    // K = (2H)^2 (1 - e^{-4u}) / 4 after the homothety by 1/(2H)
    vals[v] = state.colour === 'K' ? state.H * state.H * (1 - Math.exp(-4 * u)) : u;
  }
  // robust range: K has a long negative tail where u is small
  const sorted = vals.filter(Number.isFinite).sort();
  const q = (f) => sorted[Math.min(sorted.length - 1, Math.floor(f * sorted.length))];
  if (!sorted.length) return { vals, range: [0, 1] };
  const H2 = state.H * state.H; // K < H^2, with equality only at umbilics
  return { vals, range: state.colour === 'K' ? [Math.min(q(0.1), -1e-3 * H2), H2] : [q(0.02), q(0.98)] };
}

function show(res) {
  last = res;
  const { vals, range } = colourValues(res);
  viewer.setSurface(res, vals);
  branchMarkers(res);
  viewer.setStyle(textureStyle(range));
  if (needFrame) { viewer.frame(); needFrame = false; }
  if (res.divisor) widget.set({ divisor: res.divisor });
  const g = genus();
  const argQ = hopfArg(state.alphas, state.theta0);
  const size = Math.max(res.bbox[3] - res.bbox[0], res.bbox[4] - res.bbox[1], res.bbox[5] - res.bbox[2]) / (2 * state.H);
  if (state.gen.on && state.cell && cand.tree && res.nx > 1 && res.ny > 1) {
    const cf = cellFactors(res), ripe = cf.every((x) => Math.abs(x.c - 1) < 1e-3);
    $('treeNote').textContent = `${cand.sampleNote}: f(z + e₁) = c₁ f(z), f(z + e₂) = c₂ f(z) on the cell with c = ${cf.map((x) => fixed(x.c, 4)).join(', ')}` +
      `${ripe ? '' : ' (homotheties of ℝ⁴ about 0)'}${Math.max(...cf.map((x) => x.res)) > 1e-3 ? ` (not proportional: ${Math.max(...cf.map((x) => x.res)).toExponential(1)})` : ''}.`;
  }
  if (state.gen.on) {
    showStatus(`Genus ${g} · two Sym points, f in ℝ⁴ (${{ sphere: 'on a 3-sphere: stereographic', balanced: 'not on a 3-sphere: Möbius-balanced', pca: 'principal axes', stereo4: 'stereographic from e₄', drop4: 'drop x₄' }[res.projection]}) · κ₀ = ${fmt(kappa0(state.alphas), 3)}` +
      ` · ${res.nx}×${res.ny} in ${res.stats.ms.toFixed(0)} ms · extent ${size.toPrecision(3)}`, false);
    return;
  }
  showStatus(
    `Genus ${g} · κ₀ = ${fmt(kappa0(state.alphas), 3)} · arg Q = ${fixed(wrapPi(argQ) / PI, 3)}π` +
    ` · ${res.nx}×${res.ny} in ${res.stats.ms.toFixed(0)} ms · extent ${size.toPrecision(3)}` +
    (res.tol ? ' · adaptive steps' : '') +
    (res.stats.detDefect > 1e-6 ? ` · det drift ${res.stats.detDefect.toExponential(1)}` : ''),
    false,
  );
}

function showStatus(text, err) {
  const s = $('status');
  s.textContent = text;
  s.classList.toggle('err', !!err);
}

const wrapPi = (x) => x - 2 * PI * Math.round(x / (2 * PI));

// ------------------------------------------------------------------ spectral widget and alpha list

const widget = new SpectralWidget($('lam'), (s, dragging) => {
  const gOld = genus();
  state.alphas = s.alphas;
  state.theta0 = s.theta0;
  if (s.gen) { state.gen.lams = s.gen.lams; state.gen.signs = s.gen.signs; refreshGen(); }
  leaveCandidate();
  if (genus() !== gOld) genusChanged();
  else refreshAlphaValues();
  changed(dragging);
});

function genusChanged() {
  syncTau();
  needFrame = true;
  buildAlphaList();
  widget.set({ alphas: state.alphas, theta0: state.theta0, divisor: [] });
}

function removeAlpha(i) {
  if (i < 0 || i >= genus()) return;
  state.alphas.splice(i, 1);
  genusChanged();
  changed(false);
}

function buildAlphaList() {
  const box = $('alphaList');
  box.innerHTML = '';
  state.alphas.forEach((_, i) => {
    const row = document.createElement('div');
    row.className = 'bprow';
    row.innerHTML = `<span class="name">α<sub>${i + 1}</sub></span>
      <span class="lab">r</span><input class="num" type="text" inputmode="decimal" autocomplete="off" aria-label="Modulus of α${i + 1}">
      <span class="lab">arg</span><input class="num" type="text" inputmode="decimal" autocomplete="off" aria-label="Argument of α${i + 1} in units of π">
      <span class="muted">π</span>
      <button class="gly" type="button" data-glyph="cross" aria-label="Remove α${i + 1}" title="Remove"></button>`;
    glyphs(row);
    const [r, t] = row.querySelectorAll('input');
    const update = () => {
      const [re, im] = state.alphas[i];
      const rv = parseNum(r.value), tv = parseNum(t.value);
      const rr = Math.min(0.97, Math.max(0.03, rv ?? Math.hypot(re, im)));
      state.alphas[i] = polar(rr, tv === null ? Math.atan2(im, re) : tv * PI);
      widget.set({ alphas: state.alphas });
      refreshAlphaValues();
      changed(false);
    };
    for (const inp of [r, t]) {
      inp.addEventListener('change', update);
      inp.addEventListener('keydown', (e) => { if (e.key === 'Enter') inp.blur(); });
    }
    row.querySelector('.gly').addEventListener('click', () => removeAlpha(i));
    box.appendChild(row);
  });
  $('genus').textContent = String(genus());
  $('remAlpha').disabled = !genus();
  $('addAlpha').disabled = genus() >= MAX_GENUS;
  refreshAlphaValues();
}

function refreshAlphaValues() {
  const rows = $('alphaList').children;
  state.alphas.forEach(([re, im], i) => {
    const [r, t] = rows[i].querySelectorAll('input');
    if (document.activeElement !== r) r.value = fixed(Math.hypot(re, im), 3);
    if (document.activeElement !== t) t.value = fixed(Math.atan2(im, re) / PI, 3);
  });
  if (document.activeElement !== $('theta0')) $('theta0').value = fixed(wrapPi(state.theta0) / PI, 3);
}

// ------------------------------------------------------------------ candidates (tree.js)

// Candidate 1 for the internal tau: computed in its own worker; on success the spectral data, the two
// generalized Sym points and one cell of the closing lattice are loaded. Editing the spectral data
// afterwards leaves the candidate (the drawing goes back to the rectangular domain).
const cand = { worker: null, id: 0, last: null, tau: null, busy: false, pending: null, path: null, which: 1, tree: null, index: -1, sampleNote: '' };
function leaveCandidate() {
  if (!state.cell) return;
  state.cell = null;
  cand.path = null;
  treeView.set({ index: -1 });
  viewer.setMarkers([]);
  $('treeNote').textContent = '';
  showFamily();
  $('candNote').textContent = cand.last ? `${cand.last} (left: the data were edited)` : '';
}
// the internal tau in F_6: dragging recomputes candidate 1 (the latest request wins)
const candPad = new TauPad($('candPad'), (t, dragging) => {
  $('tauRe').value = fixed(t[0], 4); $('tauIm').value = fixed(t[1], 4);
  requestCandidate(t);
});
candPad.set({ tau: [0.2, 1.3], mode: 'F6' });
candPad.showAll();
const PAIR_NAMES = { a1: 'α₁', a2: 'α₂', node: 'the node' };
// candidates 2 and 3 share the genus 2 spectral curve: 2 takes the two pairs of branch points as Sym
// points, 3 the branch points from genus 1 and the arriving node (cw_candidates.md)
const CHOICE_NAMES = { 'a1,a2': 'Candidate 2', 'a1,node': 'Candidate 3', 'a2,node': 'Other choice' };
function candidate2Note(R) {
  const order = `Class order ${R.frame.order}.`;
  const c = (x) => `${fixed(x[0], 4)} ${x[1] < 0 ? '−' : '+'} ${fixed(Math.abs(x[1]), 4)}i`;
  const one = R.c1 ? `Candidate 1: the pair of class ${R.c1.cls} arrives at α = −${fixed(R.c1.r, 6)}; W = ${fixed(R.c1.W / PI ** 2, 4)} π².` : 'No candidate 1.';
  if (!R.event && R.openEnd) {
    const o = R.openEnd;
    const ext = R.choices[0];
    const lines = R.choices.map((ch) => `  t = ${ch.t.toExponential(0)}: W/8π = ${fixed(ch.W / (8 * PI), 4)}, distance from a round sphere ${fixed(ch.shape.defect, 3)}, area in S⁴ ${fixed(ch.shape.area, 3)} × 4π`);
    return `${order}\n${one}\nGenus 2: no arrival before the open end (both branch points go to 0, at the scale t).\n` +
      `Followed to t = ${o.t.toExponential(0)}: open-end τ of the rescaled roots = ${c(o.tau)}, the internal τ to ${o.tauError.toExponential(1)}.\n` +
      `Candidate 2 (external τ = ${c(ext.tauExt)}${ext.doubled ? ', double cover' : ''}), read in S⁴ after a Möbius normalization:\n${lines.join('\n')}\n` +
      `The limit is the round sphere covered twice (area 2 × 4π), branched at the four half-periods: the quotient of the torus by z ↦ −z. Computed in ${R.ms.toFixed(0)} ms.`;
  }
  if (!R.event) return `${order}\n${one}\nGenus 2: no arrival (${R.end.why}).`;
  if (R.double) {
    const names = { 'a1,a2': 'Candidate 1 = candidate 2 (α + node a)', 'a1,node': 'Candidate 3 (α + node b)', 'a2,node': 'Other choice (node a + node b)' };
    const lines = R.choices.map((ch) => `  ${names[ch.pair.join()]}: external τ = ${c(ch.tauExt)}, closing ${ch.defects.closing.toExponential(0)}`);
    return `${order}\n${one}\nDouble arrival: the pairs of the other two classes, equally long here, reach S¹ together, so the genus 2 edge has length 0. ` +
      `All three candidates lie on the genus 1 curve, with W = ${fixed(R.choices[0].W / PI ** 2, 4)} π², and differ only in which two of the ripe pairs (α and the two nodes) are the Sym points:\n${lines.join('\n')}\nComputed in ${R.ms.toFixed(0)} ms.`;
  }
  const lines = R.choices.map((ch) => `  ${CHOICE_NAMES[ch.pair.join()]} (${ch.pair.map((p) => PAIR_NAMES[p]).join(' + ')}): external τ = ${c(ch.tauExt)}, closing ${ch.defects.closing.toExponential(0)}`);
  return `${order}\n${one}\nGenus 2: the remaining pair arrives at t = ${fixed(R.event.t, 5)} (the scale of the branch points); W = ${fixed(R.choices[0].W, 4)} = ${fixed(R.choices[0].W / PI ** 2, 4)} π² (doubled lattice) for every choice of Sym points:\n${lines.join('\n')}\nComputed in ${R.ms.toFixed(0)} ms.`;
}
// the tree of the current candidate (treeOf in tree.js): the dot of the tree view picks the sample whose
// immersion is drawn; at unripe samples the immersion is a plane with f(z + e) = c f(z), c real (cellFactors)
const VIEW_NAMES = {
  'a1,a2': 'Candidate 2 (α₁ + α₂)', 'a1,b': 'Towards candidate 3 (α₁ + the pair)', 'a2,b': 'α₂ + the pair',
  'a1,node': 'Candidate 3 (α₁ + the node)', 'a2,node': 'Other choice (α₂ + the node)',
  'double:a1,a2': 'Candidate 1 = candidate 2 (α + node a)', 'double:a1,node': 'Candidate 3 (α + node b)', 'double:a2,node': 'Other choice (node a + node b)',
};
const viewName = (key) => VIEW_NAMES[key] || (key === 'candidate 1' ? 'Candidate 1 (α + the node)' : key.startsWith('one,') ? 'the pair for 1 + the arriving pair' : 'α + the arriving pair');
const treeView = new TreeView($('treeView'), (i, dragging) => { showSample(i); needFrame = !dragging && needFrame; changed(dragging); });
/** Draw sample i of the tree, with the view chosen in cand2Pick (kept across samples). */
function showSample(i) {
  const T = cand.tree, sm = T.samples[i];
  cand.index = i;
  const pick = $('cand2Pick'), keep = +pick.value || 0;
  pick.innerHTML = sm.views.map((v, j) => `<option value="${j}">${viewName(v.key)}</option>`).join('');
  pick.value = String(Math.min(keep, sm.views.length - 1));
  pick.hidden = sm.views.length < 2;
  const v = sm.views[+pick.value];
  if (anim) stopAnim();
  state.alphas = sm.alphas.map((a) => a.slice());
  state.gen = { ...state.gen, on: true, lams: v.syms.map((p) => p.lam.slice()), signs: v.syms.map((p) => p.sign) };
  syncTau();
  refreshAll();
  state.cell = v.lattice;
  const edge = T.samples.filter((s) => s.edge === sm.edge && s.alphas.length);
  cand.path = edge.length > 1 ? { points: edge.map((s, j) => ({ s: j, alphas: s.alphas, W: s.W })), critical: [] } : null;
  showFamily();
  cand.sampleNote = `Genus ${sm.edge}${sm.node ? `, ${sm.node}` : ''}: ${sm.edge === 1 ? `α = −${fixed(-sm.alphas[0][0], 5)}, t = ${fixed(sm.t, 5)}` : `t = ${sm.t < 0.01 ? sm.t.toExponential(2) : fixed(sm.t, 4)}`}, W = ${fixed(sm.W / PI ** 2, 4)} π² (closing lattice).\n${viewName(v.key)}: ${v.ripe ? 'ripe, a torus' : 'not ripe, a plane'}`;
  if (towardsOpenEnd()) cand.sampleNote += ' (the spikes mark the four half-periods, where the limit, the sphere covered twice, is branched)';
  $('treeNote').textContent = cand.sampleNote;
}
/** Whether the current sample is candidate 2 on a genus 2 edge that ends at an open end. */
function towardsOpenEnd() {
  const T = cand.tree, sm = T && state.cell && T.samples[cand.index];
  return !!sm && sm.edge === 2 && T.nodes.some((n) => n.label === 'open end') && sm.views[+$('cand2Pick').value || 0].key === 'a1,a2';
}
/**
 * Towards an open end, candidate 2 tends to the branched double cover T^2 -> T^2/±1 of the round sphere,
 * branched at the four half-periods of the cell (F12; they are the minima of the conformal factor all
 * along the edge). Spikes along the outward direction mark them.
 */
function branchMarkers(res) {
  if (!towardsOpenEnd() || res.nx < 3 || res.ny < 3) { viewer.setMarkers([]); return; }
  const { nx, ny, pos } = res, N = nx * ny;
  const at = (fi, fj) => {
    // bilinear interpolation at the fractional grid index (fi, fj)
    const i = Math.min(nx - 2, Math.floor(fi)), j = Math.min(ny - 2, Math.floor(fj)), a = fi - i, b = fj - j;
    const P = (ii, jj, e) => pos[3 * (jj * nx + ii) + e];
    return [0, 1, 2].map((e) => (1 - a) * (1 - b) * P(i, j, e) + a * (1 - b) * P(i + 1, j, e) + (1 - a) * b * P(i, j + 1, e) + a * b * P(i + 1, j + 1, e));
  };
  const c = [0, 0, 0];
  for (let v = 0; v < N; v++) for (let e = 0; e < 3; e++) c[e] += pos[3 * v + e] / N;
  const R = 0.5 * Math.hypot(res.bbox[3] - res.bbox[0], res.bbox[4] - res.bbox[1], res.bbox[5] - res.bbox[2]);
  const h1 = (nx - 1) / 2, h2 = (ny - 1) / 2;
  const list = [[0, 0], [h1, 0], [0, h2], [h1, h2]].map(([fi, fj]) => {
    const p = at(fi, fj), d = p.map((x, e) => x - c[e]), l = Math.hypot(...d) || 1;
    return { p, d: d.map((x) => x / l) };
  });
  viewer.setMarkers(list, { len: 0.4 * R, rad: 0.012 * R });
}
/** f(z + e_k) = c_k f(z) on the cell (c_k real at unripe Sym points; 1 when the torus closes). */
function cellFactors(res) {
  const P = res.pos4, nx = res.nx, ny = res.ny, f = (v) => [0, 1, 2, 3].map((e) => P[4 * v + e]);
  const f0 = f(0), n0 = f0.reduce((a, x) => a + x * x, 0);
  return [nx - 1, (ny - 1) * nx].map((v) => {
    const g = f(v), c = g.reduce((a, x, e) => a + x * f0[e], 0) / n0;
    return { c, res: Math.hypot(...g.map((x, e) => x - c * f0[e])) / Math.hypot(...g) };
  });
}
$('cand2Pick').addEventListener('change', () => { if (cand.tree) { showSample(cand.index); needFrame = true; changed(false); } });
$('cand2Btn').addEventListener('click', () => {
  const re = parseNum($('tauRe').value), im = parseNum($('tauIm').value);
  if (re === null || im === null || im <= 0) { $('candNote').textContent = 'Internal τ needs Im τ > 0.'; return; }
  cand.which = 2;
  candPad.set({ tau: [re, im], target: null });
  requestCandidate([re, im]);
});
function candidateNote(R) {
  const order = `Class order ${R.frame.order}: the pair of class ${R.frame.classes.one} is opened.`;
  if (!R.event) return `${order}\nNo arrival found on the genus 1 edge (it ended at α = −${fixed(R.end.r, 4)}: ${R.end.why}).`;
  const c = (x) => `${fixed(x[0], 4)} ${x[1] < 0 ? '−' : '+'} ${fixed(Math.abs(x[1]), 4)}i`;
  return `${order}\nThe pair of class ${R.event.cls} arrives at α = −${fixed(R.event.r, 6)} (t = ${fixed(R.event.t, 5)}).\n` +
    `W = ${fixed(R.W, 4)} = ${fixed(R.W / PI ** 2, 4)} π²; external τ = ${c(R.tauExt)}` +
    `${R.doubled ? ' (double cover)' : ''}.\nClosing defect ${R.defects.closing.toExponential(1)}; computed in ${R.ms.toFixed(0)} ms.`;
}
$('cand1Btn').addEventListener('click', () => {
  const re = parseNum($('tauRe').value), im = parseNum($('tauIm').value);
  if (re === null || im === null || im <= 0) { $('candNote').textContent = 'Internal τ needs Im τ > 0.'; return; }
  cand.which = 1;
  candPad.set({ tau: [re, im], target: null });
  requestCandidate([re, im]);
});
function requestCandidate(tau) {
  cand.pending = tau;
  if (cand.busy) return;
  const t = cand.pending;
  cand.pending = null;
  if (!cand.worker) {
    cand.worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
    cand.worker.onerror = (e) => { $('cand1Btn').disabled = false; $('candNote').textContent = `worker error: ${e.message}`; };
    cand.worker.onmessage = (e) => {
      const { id, result: R, error } = e.data;
      cand.busy = false;
      $('cand1Btn').disabled = false; $('cand2Btn').disabled = false;
      if (cand.pending) requestCandidate(cand.pending);
      if (error) { $('candNote').textContent = error; return; }
      cand.last = R.which === 2 ? candidate2Note(R) : candidateNote(R);
      $('candNote').textContent = cand.last;
      candPad.set({ tau: R.tau, blocked: !R.event && !R.openEnd });
      const prev = cand.tree, prevSm = prev && prev.samples[cand.index];
      cand.tree = R.tree || null;
      treeView.set({ tree: cand.tree, index: -1 });
      if (!cand.tree) { $('cand2Pick').hidden = true; $('treeNote').textContent = ''; return; }
      // the dot: kept at the same place of the path when the tree is of the same kind (e.g. while tau is
      // dragged), otherwise at the last candidate, or near the open end at t = 0.01 (closer is slow to draw)
      const S = cand.tree.samples, nearX = (x) => S.reduce((b, s, i) => (Math.abs(s.x - x) < Math.abs(S[b].x - x) ? i : b), 0);
      const kind = (T) => T.nodes.map((n) => n.label).join();
      let idx;
      if (prevSm && kind(prev) === kind(cand.tree)) idx = nearX(prevSm.x);
      else if (R.openEnd) idx = S.reduce((b, s, i) => (s.edge === 2 && Math.abs(Math.log(s.t / 0.01)) < Math.abs(Math.log(S[b].t / 0.01)) ? i : b), S.length - 1);
      else idx = cand.tree.nodes.length ? cand.tree.nodes[cand.tree.nodes.length - 1].index : S.length - 1;
      treeView.set({ index: idx });
      showSample(idx);
      if (!cand.pending) needFrame = true;
      changed(!!cand.pending);
    };
  }
  cand.busy = true;
  $('cand1Btn').disabled = true; $('cand2Btn').disabled = true;
  cand.tau = t;
  if (!state.cell) $('candNote').textContent = cand.which === 2 ? 'Following the genus 1 and genus 2 edges…' : 'Following the genus 1 edge…';
  cand.worker.postMessage({ id: ++cand.id, kind: cand.which === 2 ? 'candidate2' : 'candidate1', params: { tau: t } });
}

// generalized Sym points: checkbox, readout, projection
function refreshGen() {
  const G = state.gen;
  $('genOn').checked = G.on;
  $('genBox').hidden = !G.on;
  const show = ([re, im], s, i) => `λ<sub>${i + 1}</sub> = ${fixed(re, 3)} ${im < 0 ? '−' : '+'} ${fixed(Math.abs(im), 3)}i, sheet ${s > 0 ? '+' : '−'}`;
  $('genInfo').innerHTML = G.lams.map((l, i) => show(l, G.signs[i], i)).join('<br>');
  $('proj').value = G.proj;
  widget.set({ gen: G.on ? { lams: G.lams, signs: G.signs } : null });
}
$('genOn').addEventListener('change', () => {
  state.gen.on = $('genOn').checked;
  if (state.gen.on && state.alphas.length && state.gen.lams.every(([re, im]) => Math.hypot(re, im) > 0.999 && Math.hypot(re, im) < 1.001)) {
    // a first look: one generalized Sym point at a branch point
    state.gen.lams = [state.alphas[0].slice(), state.gen.lams[1]];
  }
  refreshGen();
  needFrame = true;
  changed(false);
});
$('proj').addEventListener('change', () => { state.gen.proj = $('proj').value; changed(false); });

$('theta0').addEventListener('change', () => {
  const v = parseNum($('theta0').value);
  if (v !== null) state.theta0 = wrapPi(v * PI);
  widget.set({ theta0: state.theta0 });
  refreshAlphaValues();
  changed(false);
});
$('theta0').addEventListener('keydown', (e) => { if (e.key === 'Enter') $('theta0').blur(); });
$('addAlpha').addEventListener('click', () => {
  if (genus() >= MAX_GENUS) return;
  // place the new point where it is well separated from the others
  let best = null, bestD = -1;
  for (let k = 0; k < 24; k++) {
    const z = polar(0.5, (2 * PI * k) / 24 + 0.1);
    const d = Math.min(...state.alphas.map((a) => Math.hypot(a[0] - z[0], a[1] - z[1])), 9);
    if (d > bestD) { bestD = d; best = z; }
  }
  state.alphas.push(best);
  genusChanged();
  changed(false);
});
$('remAlpha').addEventListener('click', () => removeAlpha(genus() - 1));

// ------------------------------------------------------------------ sliders

/**
 * A slider row: label, track, action buttons. The value floats above the thumb (display(v), default
 * fixed digits of v / unit plus suffix); clicking it opens a field for exact entry (enter(x) if given,
 * else set(x * unit)). snap() lists values that get a tick on the track and catch the thumb within a
 * few pixels. play: an animation toggle; reset: { title, onClick } for a second button.
 */
function slider(parent, { label, min: lo, max: hi, step, get, set, digits = 3, unit = 1, suffix = '', display, entry, enter,
  snap, play, reset, spacer = false, recompute = true, title }) {
  const after = recompute ? changed : () => {};
  const id = `sl${(slider.n = (slider.n || 0) + 1)}`;
  const row = document.createElement('div');
  row.className = 'sl';
  row.innerHTML = `<label for="${id}">${label}</label><div class="trk"><button type="button" class="tag" title="Click to type a value"></button>` +
    `<input type="range" class="rng" id="${id}"></div><div class="acts"></div>`;
  const rng = row.querySelector('input'), tag = row.querySelector('.tag'), trk = row.querySelector('.trk');
  const acts = row.querySelector('.acts');
  if (title) row.title = title;
  let min = lo, max = hi;
  Object.assign(rng, { min, max, step });
  const frac = (v) => Math.min(1, Math.max(0, (v - min) / (max - min)));
  const at = (f) => `calc(6.5px + (100% - 13px) * ${f.toFixed(5)})`;
  const text = (v) => (display ? display(v) : fixed(v / unit, digits) + suffix);
  let ticks = [];
  const refresh = () => {
    const v = get();
    rng.value = v;
    const f = frac(v);
    rng.style.setProperty('--p', `${(100 * f).toFixed(3)}%`);
    tag.textContent = text(v);
    tag.style.left = at(f);
    const marks = snap ? snap().filter((c) => c >= min && c <= max) : [];
    while (ticks.length > marks.length) ticks.pop().remove();
    while (ticks.length < marks.length) {
      const t = document.createElement('span');
      t.className = 'tick';
      trk.appendChild(t);
      ticks.push(t);
    }
    marks.forEach((c, k) => { ticks[k].style.left = at(frac(c)); });
  };
  rng.addEventListener('input', () => {
    let v = +rng.value;
    if (snap) {
      const px = (rng.clientWidth - 13) / (max - min);
      for (const c of snap()) if (Math.abs(v - c) * px < 3) v = c;
    }
    set(v);
    refresh();
    after(true);
  });
  tag.addEventListener('click', () => {
    const inp = document.createElement('input');
    inp.className = 'num tagin';
    inp.type = 'text';
    inp.inputMode = 'decimal';
    inp.autocomplete = 'off';
    inp.setAttribute('aria-label', `${rng.labels[0]?.textContent || 'value'}: exact value`);
    inp.value = entry ? entry() : fixed(get() / unit, digits);
    inp.style.left = tag.style.left;
    tag.hidden = true;
    trk.appendChild(inp);
    inp.focus();
    inp.select();
    let done = false;
    const finish = (commit) => {
      if (done) return;
      done = true;
      const v = parseNum(inp.value);
      inp.remove();
      tag.hidden = false;
      if (commit && v !== null) {
        if (enter) enter(v); else set(v * unit);
        after(false);
      }
      refresh();
    };
    inp.addEventListener('keydown', (e) => {
      if (e.key === 'Enter') finish(true);
      else if (e.key === 'Escape') finish(false);
    });
    inp.addEventListener('blur', () => finish(true));
  });
  const button = (glyph, tip, onClick) => {
    const b = document.createElement('button');
    b.type = 'button';
    b.className = 'gly';
    b.innerHTML = GLYPH[glyph];
    b.title = tip;
    b.setAttribute('aria-label', tip);
    b.addEventListener('click', () => onClick(b));
    acts.appendChild(b);
    return b;
  };
  if (play) button('play', 'Animate', (b) => toggleAnim(play, b));
  // spacer: the room of a play button, so the track lines up with sliders that have one
  if (spacer) { const sp = document.createElement('span'); sp.style.width = '34px'; sp.setAttribute('aria-hidden', 'true'); acts.appendChild(sp); }
  if (reset) button('reset', reset.title, reset.onClick);
  parent.appendChild(row);
  refresh();
  const setRange = (a, b) => {
    if (a === min && b === max) return;
    [min, max] = [a, b];
    Object.assign(rng, { min, max });
    refresh();
  };
  return { refresh, row, setRange };
}

const appearanceSliders = [
  slider($('appearanceRows'), {
    label: 'Line spacing', min: 0.02, max: 2, step: 0.001, get: () => state.gridStep,
    set: (v) => { state.gridStep = Math.max(0.005, v); viewer.setStyle({ gridStep: state.gridStep }); },
    recompute: false,
  }),
  slider($('appearanceRows'), {
    label: 'Line width', min: 0.2, max: 2.5, step: 0.01, digits: 2, get: () => state.lineW,
    set: (v) => { state.lineW = v; viewer.setStyle({ lineWidth: v }); },
    recompute: false,
  }),
];
const transpSlider = slider($('appearanceRows'), {
  label: 'Transparency', title: 'How much of the surface shows through (the parameter lines stay opaque)',
  min: 0, max: 1, step: 0.01, digits: 2, get: () => state.transp,
  set: (v) => { state.transp = Math.min(1, Math.max(0, v)); viewer.setStyle({ opacity: 1 - state.transp }); },
  recompute: false,
});
appearanceSliders.push(transpSlider);

// ------------------------------------------------------------------ mean curvature (a homothety in R^3)

function recolour() {
  if (!last) return;
  const { vals, range } = colourValues(last);
  if (vals) viewer.setColourValues(vals);
  viewer.setStyle(textureStyle(range));
}
function textureStyle(range) {
  const wire = state.colour === 'wire';
  // a wireframe is its parameter lines; transparency applies to the fill
  $('grid').disabled = wire;
  transpSlider.row.classList.toggle('hide', wire);
  return { mode: { side: 0, u: 1, K: 2 }[state.colour] ?? 0, range, wire, grid: wire || state.grid, opacity: 1 - state.transp };
}
// ------------------------------------------------------------------ common roots of the differentials on S^1

// where every Theta_w vanishes: a Sym point there has no translational period in any direction
function updateCommonRoots() {
  let roots = [];
  if (genus() >= 1) {
    try { roots = commonRootsOnCircle(state.alphas, 1e-5).map((r) => r.lam); } catch { roots = []; }
  }
  commonRoots = roots;
  widget.set({ commonRoots: roots });
}
let commonRoots = [];

// ------------------------------------------------------------------ paths in the λ-plane

/** A candidate shows its path along the current edge of the tree in the λ-plane. */
function showFamily() {
  widget.set({ family: state.cell && cand.path ? cand.path : null });
}

// ------------------------------------------------------------------ appearance

$('colour').addEventListener('change', () => {
  state.colour = $('colour').value;
  if (last) recolour(); else viewer.setStyle(textureStyle());
});
$('grid').addEventListener('change', () => { state.grid = $('grid').checked; viewer.setStyle(textureStyle()); });
$('front').addEventListener('input', () => viewer.setStyle({ front: $('front').value }));
$('back').addEventListener('input', () => viewer.setStyle({ back: $('back').value }));
$('bg').addEventListener('input', () => setBackground($('bg').value));
function setBackground(c) {
  state.bg = c;
  $('stage').style.background = c;
  viewer.setBackground(c);
}
$('adapt').addEventListener('change', () => { state.adapt = $('adapt').checked; changed(false); });
$('res').addEventListener('change', () => { state.res = +$('res').value; changed(false); });
$('acc').addEventListener('change', () => {
  state.adaptive = $('acc').value === 'adaptive';
  state.hmax = state.adaptive ? 0.05 : +$('acc').value;
  changed(false);
});

// ------------------------------------------------------------------ animation

let anim = null, animBtn = null, animLast = 0;
function toggleAnim(a, btn) {
  const same = anim && anim.key === a.key && anim.i === a.i;
  stopAnim();
  if (same) return;
  anim = a;
  animBtn = btn;
  btn.classList.add('on');
  btn.innerHTML = GLYPH.stop;
  animLast = performance.now();
  requestAnimationFrame(animStep);
}
function stopAnim() {
  if (animBtn) { animBtn.classList.remove('on'); animBtn.innerHTML = GLYPH.play; }
  anim = null;
  animBtn = null;
  changed(false);
}
function animStep(now) {
  if (!anim) return;
  const dt = Math.min(0.1, (now - animLast) / 1000);
  animLast = now;
  if (anim.key === 'theta0') {
    state.theta0 = wrapPi(state.theta0 + 0.3 * dt);
    widget.set({ theta0: state.theta0 });
    refreshAlphaValues();
  }
  // only ask for a new surface when the worker has caught up
  if (!busy) changed(true);
  requestAnimationFrame(animStep);
}
$('playLam').addEventListener('click', () => toggleAnim({ key: 'theta0' }, $('playLam')));

// ------------------------------------------------------------------ presets and export

const presetSel = $('preset');
presetSel.innerHTML = '<option value="">Choose a preset…</option>' +
  PRESETS.map((p, i) => `<option value="${i}">${p.name}</option>`).join('');
presetSel.addEventListener('change', () => {
  const p = PRESETS[+presetSel.value];
  if (!p) return;
  if (anim) stopAnim();
  leaveCandidate();
  Object.assign(state, structuredClone(VIEW_DEFAULTS), structuredClone(p.s));
  syncTau();
  refreshAll();
  viewer.resetOrientation();
  needFrame = true;
  changed(false);
  presetSel.value = '';
});

$('frameBtn').addEventListener('click', () => viewer.frame());
$('viewReset').addEventListener('click', () => { viewer.resetOrientation(); viewer.frame(); });

$('pngBtn').addEventListener('click', () => download(viewer.snapshot(), 'cmc-surface.png'));
$('objBtn').addEventListener('click', () => {
  const url = URL.createObjectURL(new Blob([viewer.toOBJ()], { type: 'text/plain' }));
  download(url, 'cmc-surface.obj');
  setTimeout(() => URL.revokeObjectURL(url), 1000);
});
$('linkBtn').addEventListener('click', async () => {
  const url = linkURL();
  try {
    await navigator.clipboard.writeText(url);
    $('linkBtn').textContent = 'Copied';
    setTimeout(() => { $('linkBtn').textContent = 'Copy link'; }, 1500);
  } catch {
    window.prompt('Copy this link:', url); // no clipboard access (e.g. not https)
  }
});
function download(href, name) {
  const a = document.createElement('a');
  a.href = href;
  a.download = name;
  a.click();
}

// ------------------------------------------------------------------ start

function refreshAll() {
  buildAlphaList();
  widget.set({ alphas: state.alphas, theta0: state.theta0, divisor: [] });
  refreshGen();
  for (const s of appearanceSliders) s.refresh();
  $('res').value = String(state.res);
  if (!$('res').value) { $('res').add(new Option(String(state.res), String(state.res))); $('res').value = String(state.res); }
  $('acc').value = state.adaptive ? 'adaptive' : String(state.hmax);
  if (!$('acc').value) { $('acc').add(new Option(`h ≤ ${state.hmax}`, String(state.hmax))); $('acc').value = String(state.hmax); }
  viewer.setScale(1 / (2 * state.H));
  $('colour').value = state.colour;
  viewer.setStyle({ gridStep: state.gridStep, lineWidth: state.lineW, ...textureStyle() });
}

// a link pasted into the address bar of an open page (clearing the hash with replaceState is silent)
window.addEventListener('hashchange', () => {
  const ok = readHash();
  clearHash();
  if (pendingCandidate) { runPendingCandidate(); return; }
  if (!ok) return;
  snapLinkLam();
  if (anim) stopAnim();
  syncTau();
  refreshAll();
  needFrame = true;
  changed(false);
});

if (location.hash) { readHash(); snapLinkLam(); clearHash(); }
// without a link: a typical case, candidate 2 for internal tau = 0.45 + 0.95i (candidates 2 and 3 exist)
else pendingCandidate = [0.45, 0.95, 2];
syncTau();
makeWorker();
refreshAll();
setBackground(state.bg);
// the candidate draws itself when it arrives; otherwise draw the data from the link
if (!pendingCandidate) changed(false);
function runPendingCandidate() {
  $('tauRe').value = String(pendingCandidate[0]); $('tauIm').value = String(pendingCandidate[1]);
  const two = pendingCandidate[2] === 2;
  pendingCandidate = null;
  $(two ? 'cand2Btn' : 'cand1Btn').click();
}
if (pendingCandidate) runPendingCandidate();
window.cmcStarted = true; // see index.html
