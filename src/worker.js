// Runs the surface computation off the main thread.
// Message in:  { id, kind: 'candidate1', params: { tau } } (candidate 1 and the tree, tree.js),
//              { id, kind: 'candidate2', params: { tau, fast } } (candidates 1-3 and the tree),
//              { id, kind: 'shapes', params: { choices } } (the open-end shape diagnostics),
//              { id, params }   (params as for computeSurface; with params.gen = { lams, signs, projection }
//              the generalized Sym surface of gensym.js), { id, kind: 'close', params },
// Message out: { id, result } or { id, error }

import { computeSurface, closingInfo } from './cmc/cmc.js';
import { computeGenSurface } from './cmc/gensym.js';
import { candidate1, candidate2, treeOf, openEndShapes } from './cmc/tree.js';

/** The branch at candidate 1 that does not open the node (candidate 1'), for the note. */
const branchInfo = (b) => (b && b.event ? { name: b.event.name, cls: b.event.cls, leaving: b.event.leaving, r: b.event.r, W: b.W, tauExt: b.tauExt, closing: b.defects.closing } : null);

self.onmessage = (e) => {
  const { id, params, kind } = e.data;
  if (kind === 'candidate1') {
    try {
      const t0 = performance.now();
      const C = candidate1(params.tau);
      const { edge, event } = C;
      const out = { tau: C.tau, frame: C.frame, samples: edge.samples.map(({ r, t, W }) => ({ r, t, W })), end: edge.end, event: null, ms: 0, tree: treeOf(C), branch: branchInfo(C.branch) };
      if (event) {
        Object.assign(out, {
          event: { name: event.name, cls: event.cls, r: event.r, t: event.t, W: event.W, lamNode: event.lamNode },
          alphas: C.alphas, syms: C.syms, lattice: C.lattice, doubled: C.doubled, chars: C.chars,
          W: C.W, tauExt: C.tauExt, defects: C.defects,
        });
      }
      out.ms = performance.now() - t0;
      self.postMessage({ id, result: out });
    } catch (err) {
      self.postMessage({ id, error: String(err && err.message ? err.message : err) });
    }
    return;
  }
  if (kind === 'candidate2') {
    try {
      const t0 = performance.now();
      const C = candidate2(params.tau, { fast: !!params.fast });
      const c1 = C.c1;
      const out = { tau: C.tau, frame: c1.frame, which: 2, event: null, ms: 0,
        c1: c1.event ? { r: c1.event.r, cls: c1.event.cls, W: c1.W, tauExt: c1.tauExt, branch: branchInfo(c1.branch) } : null,
        samples: C.edge ? C.edge.samples.map(({ s, alphas, W }) => ({ s, alphas, W })) : [], end: C.edge ? C.edge.end : { why: 'no candidate 1' } };
      if (C.openEnd) out.openEnd = C.openEnd;
      if (C.double) out.double = true;
      if (c1.event) out.tree = treeOf(C);
      if (C.event || C.openEnd) {
        Object.assign(out, {
          event: C.event ? { s: C.event.s, t: C.event.t, W: C.event.W, lamNode: C.event.lamNode, gap: C.event.gap } : null,
          alphas: C.alphas,
          choices: C.choices.map((ch) => ({ pair: ch.pair, syms: ch.syms, lattice: ch.lattice, doubled: ch.doubled, chars: ch.chars, W: ch.W, tauExt: ch.tauExt, defects: ch.defects, t: ch.t, alphas: ch.alphas, shape: ch.shape })),
        });
      }
      out.ms = performance.now() - t0;
      self.postMessage({ id, result: out });
    } catch (err) {
      self.postMessage({ id, error: String(err && err.message ? err.message : err) });
    }
    return;
  }
  if (kind === 'shapes') {
    // the open-end shape diagnostics, on request (openEndShapes)
    try {
      const t0 = performance.now();
      const ch = openEndShapes(params.choices);
      self.postMessage({ id, result: { shapes: ch.map((c) => c.shape), ms: performance.now() - t0 } });
    } catch (err) {
      self.postMessage({ id, error: String(err && err.message ? err.message : err) });
    }
    return;
  }
  if (kind === 'close') {
    try {
      self.postMessage({ id, result: { s: closingInfo(params, 's'), t: closingInfo(params, 't') } });
    } catch (err) {
      self.postMessage({ id, error: String(err && err.message ? err.message : err) });
    }
    return;
  }
  try {
    if (params.gen) {
      const r = computeGenSurface({ ...params, ...params.gen });
      self.postMessage({ id, result: r }, [r.pos.buffer, r.nrm.buffer, r.u.buffer, r.st.buffer, r.pos4.buffer, r.w4.buffer]);
      return;
    }
    const r = computeSurface(params);
    self.postMessage({ id, result: r }, [r.pos.buffer, r.nrm.buffer, r.u.buffer, r.st.buffer]);
  } catch (err) {
    self.postMessage({ id, error: String(err && err.message ? err.message : err) });
  }
};
