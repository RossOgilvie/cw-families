// Precompute a map of F_6: for internal tau on a grid, the kind of tree (candidates 2 and 3 at an arrival
// on the genus 2 edge, the open end, a double arrival, or failure) and the Willmore energies; and the
// boundary between "candidates" and "open end", refined by bisection along the grid edges and joined by
// marching squares. Writes data/f6map.json (read by the app for the tau pad).
//   node tools/f6map.mjs [step = 0.05] [imMax = 2.4]
import { writeFileSync } from 'node:fs';
import { candidate2 } from '../src/cmc/tree.js';

const step = +(process.argv[2] || 0.05), imMax = +(process.argv[3] || 2.4);
const inF6 = ([x, y]) => y > 0 && Math.abs(x) <= 1 + 1e-12 && Math.hypot(x - 0.5, y) >= 0.5 - 1e-12 && Math.hypot(x + 0.5, y) >= 0.5 - 1e-12;
const kindOf = (C) => (!C.c1 || !C.c1.event ? 'none' : C.double ? 'double' : C.event ? 'arrival' : C.openEnd ? 'open' : 'fail');
function classify(tau) {
  try {
    const C = candidate2(tau, { fast: true });
    return { kind: kindOf(C), C };
  } catch (e) { return { kind: 'fail', C: null, err: String(e.message || e) }; }
}

const xs = [], ys = [];
for (let x = -1; x <= 1 + 1e-9; x += step) xs.push(+x.toFixed(6));
for (let y = step; y <= imMax + 1e-9; y += step) ys.push(+y.toFixed(6));
const grid = ys.map(() => xs.map(() => null));
const t0 = Date.now();
let n = 0;
for (let j = 0; j < ys.length; j++) for (let i = 0; i < xs.length; i++) {
  const tau = [xs[i], ys[j]];
  if (!inF6(tau)) continue;
  const { kind, C } = classify(tau);
  const pt = { kind };
  if (C && C.c1 && C.c1.event) {
    pt.order = C.frame.order;
    pt.W1 = +(C.c1.W / Math.PI ** 2).toFixed(6);
    if (C.choices && C.choices.length) pt.W2 = +(C.choices[0].W / Math.PI ** 2).toFixed(6);
    if (C.event && !C.double) pt.t = +C.event.t.toExponential(6);
    if (C.c1.branch && C.c1.branch.event) pt.W1b = +(C.c1.branch.W / Math.PI ** 2).toFixed(6);
  }
  grid[j][i] = pt;
  if (++n % 100 === 0) process.stderr.write(`${n} points, ${((Date.now() - t0) / 1000).toFixed(0)} s\n`);
}

// the boundary between 'arrival' and 'open' on each grid edge, by bisection
const isA = (p) => p && p.kind === 'arrival', isO = (p) => p && p.kind === 'open';
function crossing(a, b) {
  // a: arrival side, b: open side
  let lo = a, hi = b;
  for (let k = 0; k < 14; k++) {
    const m = [(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2], { kind } = classify(m);
    if (kind === 'arrival') lo = m; else if (kind === 'open') hi = m; else break;
  }
  return [+((lo[0] + hi[0]) / 2).toFixed(6), +((lo[1] + hi[1]) / 2).toFixed(6)];
}
const edgeKey = (i1, j1, i2, j2) => `${i1},${j1},${i2},${j2}`;
const cross = new Map();
const edgeCross = (i1, j1, i2, j2) => {
  const k = edgeKey(i1, j1, i2, j2);
  if (cross.has(k)) return cross.get(k);
  const p = grid[j1][i1], q = grid[j2][i2];
  let c = null;
  if (isA(p) && isO(q)) c = crossing([xs[i1], ys[j1]], [xs[i2], ys[j2]]);
  else if (isO(p) && isA(q)) c = crossing([xs[i2], ys[j2]], [xs[i1], ys[j1]]);
  cross.set(k, c);
  return c;
};
// marching squares: each cell whose four corners are computed contributes the segments between its crossings
const segments = [];
for (let j = 0; j + 1 < ys.length; j++) for (let i = 0; i + 1 < xs.length; i++) {
  const c = [grid[j][i], grid[j][i + 1], grid[j + 1][i + 1], grid[j + 1][i]];
  if (c.some((p) => !p)) continue;
  const pts = [edgeCross(i, j, i + 1, j), edgeCross(i + 1, j, i + 1, j + 1), edgeCross(i, j + 1, i + 1, j + 1), edgeCross(i, j, i, j + 1)].filter(Boolean);
  if (pts.length === 2) segments.push(pts);
  else if (pts.length === 4) segments.push([pts[0], pts[1]], [pts[2], pts[3]]);
}
const counts = {};
for (const row of grid) for (const p of row) if (p) counts[p.kind] = (counts[p.kind] || 0) + 1;
const out = { step, xs, ys, grid, boundary: segments, counts, made: new Date().toISOString(), seconds: (Date.now() - t0) / 1000 };
writeFileSync(new URL('../data/f6map.json', import.meta.url), JSON.stringify(out));
console.log(JSON.stringify({ counts, segments: segments.length, seconds: out.seconds }));
