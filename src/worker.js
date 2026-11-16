// Runs the surface computation off the main thread.
// Message in:  { id, kind: 'candidate1', params: { tau } } (candidate 1 of the tree, tree.js),
//              { id, params }   (params as for computeSurface; with params.gen = { lams, signs, projection }
//              the generalized Sym surface of gensym.js), { id, kind: 'close', params },
//              or { id, kind: 'family', params: { alphas, smax, maxSteps } } (the Whitham curve)
// Message out: { id, result } or { id, error }

import { computeSurface, closingInfo } from './cmc/cmc.js';
import { computeGenSurface } from './cmc/gensym.js';
import { candidate1 } from './cmc/tree.js';
import { WhithamCurve } from './cmc/whitham.js';

self.onmessage = (e) => {
  const { id, params, kind } = e.data;
  if (kind === 'family') {
    try {
      const t0 = performance.now();
      const tr = new WhithamCurve(params.alphas).trace(params.smax, params.maxSteps);
      self.postMessage({ id, result: { ...tr, ms: performance.now() - t0 } });
    } catch (err) {
      self.postMessage({ id, error: String(err && err.message ? err.message : err) });
    }
    return;
  }
  if (kind === 'candidate1') {
    try {
      const t0 = performance.now();
      const C = candidate1(params.tau);
      const { edge, event } = C;
      const out = { tau: C.tau, frame: C.frame, samples: edge.samples.map(({ r, t, W }) => ({ r, t, W })), end: edge.end, event: null, ms: 0 };
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
      self.postMessage({ id, result: r }, [r.pos.buffer, r.nrm.buffer, r.u.buffer, r.st.buffer, r.pos4.buffer]);
      return;
    }
    const r = computeSurface(params);
    self.postMessage({ id, result: r }, [r.pos.buffer, r.nrm.buffer, r.u.buffer, r.st.buffer]);
  } catch (err) {
    self.postMessage({ id, error: String(err && err.message ? err.message : err) });
  }
};
