// The map of F_1 for the candidate family (F17 of cw_candidates.md): on a grid of the conformal class
// tau^ (Re tau^ >= 0; the map is symmetric under tau^ -> -conj(tau^), which keeps the labels h), for each
// member h = 1/2, tau/2, (1 + tau)/2 the least W of the genus 1 and of the genus 2 curves found over its
// superlattice that close on <1, tau^> (familyOf with reduce: false, so that the labels are continuous
// across the boundary of F_1; the grid reaches a little outside F_1 for the contours).
// Usage: node tools/f1map.mjs [step = 0.025] [processes = 14]   ->  data/f1map.json (several minutes)
import { fork } from 'node:child_process';
import { writeFileSync } from 'node:fs';
import { familyOf } from '../src/cmc/tree.js';

const P2 = Math.PI ** 2;
function point([x, y]) {
  const t0 = performance.now();
  const F = familyOf([x, y], { reduce: false });
  const members = F.members.map((m) => {
    const g = (k) => { const f = m.found.filter((c) => c.genus === k); return f.length ? f[0].W / P2 : null; };
    return { W1: g(1), W2: g(2) };
  });
  return { members, ms: Math.round(performance.now() - t0) };
}

if (process.argv[2] === 'child') {
  process.on('message', (pts) => {
    for (const p of pts) {
      let r;
      try { r = point(p.tau); } catch (e) { r = { error: String(e.message || e) }; }
      process.send({ ...p, ...r });
    }
    process.exit(0);
  });
} else {
  const step = +(process.argv[2] || 0.025), nproc = +(process.argv[3] || 14);
  const xs = [], ys = [];
  for (let x = 0; x <= 0.55 + 1e-9; x += step) xs.push(+x.toFixed(6));
  for (let y = 0.8; y <= 2.3 + 1e-9; y += step) ys.push(+y.toFixed(6));
  const pts = [];
  ys.forEach((y, j) => xs.forEach((x, i) => { if (x * x + y * y > 0.78) pts.push({ i, j, tau: [x, y] }); }));
  // deal the points round-robin, so that the slow ones (large Im) are spread out
  const parts = Array.from({ length: nproc }, (_, k) => pts.filter((_, q) => q % nproc === k));
  const grid = ys.map(() => xs.map(() => null));
  let done = 0;
  const t0 = Date.now();
  await Promise.all(parts.map((part) => new Promise((res) => {
    const c = fork(new URL(import.meta.url).pathname, ['child']);
    c.on('message', (m) => {
      grid[m.j][m.i] = m.error ? { error: m.error } : { members: m.members, ms: m.ms };
      if (++done % 50 === 0) console.log(`${done}/${pts.length} in ${((Date.now() - t0) / 1000).toFixed(0)} s`);
    });
    c.on('exit', res);
    c.send(part);
  })));
  writeFileSync(new URL('../data/f1map.json', import.meta.url), JSON.stringify({ step, xs, ys, grid, made: new Date().toISOString() }));
  console.log(`done: ${pts.length} points in ${((Date.now() - t0) / 1000).toFixed(0)} s`);
}
