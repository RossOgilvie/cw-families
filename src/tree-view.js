// The tree of one internal tau (tree.js, treeOf) as a graph of W against the path through the edges:
// the opening of the pair for 1 on the left axis, the genus 1 edge on the first third of the width (to the
// first arrival, candidate 1), the genus 2 edge on the other two thirds (to the arrival of candidates 2
// and 3, or the open end). The genus 0 edge is not shown. Path coordinates x (treeOf): genus 1 on [1, 2],
// genus 2 on [2, 3]. The branch that does not open the node at candidate 1 (edge 3, the genus 1 edge
// continued to candidate 1') is drawn dashed and blue over the same span as genus 2, from candidate 1; the
// up and down arrows switch between it and the genus 2 edge.
// W is in units of pi^2, on the closing (doubled) lattice. A dot marks the sample whose immersion is shown;
// it can be dragged along the path (snapping to the nodes) or stepped with the arrow keys.

/** Fraction of the width at path coordinate x: genus 1 [1, 2] -> [0, 1/3], genus 2 [2, 3] -> [1/3, 1]. */
const frac = (x) => (x <= 2 ? Math.max(0, (x - 1) / 3) : 1 / 3 + ((x - 2) * 2) / 3);
/** Shown: the genus 1 and 2 edges, and the opening (the last sample of genus 0). */
const shown = (s) => s.edge > 0 || s.x >= 1;

export class TreeView {
  /** onPick(index, dragging) when the user moves the dot */
  constructor(canvas, onPick) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.onPick = onPick;
    this.tree = null;   // { samples: [{ x, W, edge, node }], nodes: [{ index, label }] }
    this.index = -1;
    this.drag = false;
    this.full = false;  // W from 0 (double-click), or the band of the genus 1 and 2 edges
    this._bind();
    new ResizeObserver(() => this.draw()).observe(canvas);
    if (document.fonts) document.fonts.load("14px 'EB Garamond'", 'πW').then(() => this.draw(), () => {});
  }

  set(opts) {
    Object.assign(this, opts);
    this.draw();
  }

  _geom() {
    const w = this.canvas.clientWidth, h = this.canvas.clientHeight;
    const pad = { l: 36, r: 10, t: 16, b: 22 };
    // in units of pi^2
    let lo = 0, hi = 2.7;
    if (this.tree) {
      const Ws = this.tree.samples.map((s) => s.W / Math.PI ** 2), up = this.tree.samples.filter((s) => s.edge > 0).map((s) => s.W / Math.PI ** 2);
      const top = Math.max(8 / Math.PI, ...Ws), bot = this.full || !up.length ? 0 : Math.min(...up);
      const span = top - bot;
      hi = top + 0.12 * span;
      lo = this.full ? 0 : Math.max(0, bot - 0.25 * span);
    }
    return { w, h, pad, lo, hi, X: (x) => pad.l + frac(x) * (w - pad.l - pad.r), Y: (W) => h - pad.b - ((W / Math.PI ** 2 - lo) / (hi - lo)) * (h - pad.t - pad.b) };
  }

  /**
   * The sample nearest to the pixel (px, py), snapping to a node within 9 px. Horizontal distance counts
   * more than vertical, so that dragging along an edge follows it; the vertical part chooses between the
   * genus 2 edge and the branch where both run over the same span.
   */
  _pick(px, py) {
    const S = this.tree.samples, { X, Y } = this._geom();
    let best = 0, bd = Infinity;
    S.forEach((s, i) => { if (!shown(s)) return; const d = Math.hypot(X(s.x) - px, 0.35 * (Y(s.W) - py)); if (d < bd) { bd = d; best = i; } });
    for (const n of this.tree.nodes) if (Math.hypot(X(S[n.index].x) - px, Y(S[n.index].W) - py) < 9) return n.index;
    return best;
  }

  /** The samples along the branch of sample i, in order: the main path, or candidate 1 and edge 3. */
  _path(i) {
    const S = this.tree.samples;
    const idx = S.map((_, k) => k);
    if (S[i] && S[i].edge === 3) {
      const c1 = this.tree.nodes.find((n) => n.label === 'candidate 1');
      return [...(c1 ? [c1.index] : []), ...idx.filter((k) => S[k].edge === 3)];
    }
    return idx.filter((k) => shown(S[k]) && S[k].edge !== 3);
  }

  _bind() {
    const c = this.canvas;
    const px = (e) => e.clientX - c.getBoundingClientRect().left;
    const py = (e) => e.clientY - c.getBoundingClientRect().top;
    const emit = (i, dragging) => {
      if (i === this.index && dragging) return;
      this.index = i;
      this.draw();
      this.onPick(i, dragging);
    };
    c.addEventListener('pointerdown', (e) => {
      if (!this.tree) return;
      this.drag = true;
      c.setPointerCapture(e.pointerId);
      c.focus();
      e.preventDefault();
      emit(this._pick(px(e), py(e)), true);
    });
    c.addEventListener('pointermove', (e) => {
      c.style.cursor = this.tree ? (this.drag ? 'grabbing' : 'pointer') : 'default';
      if (this.drag) emit(this._pick(px(e), py(e)), true);
    });
    const up = () => {
      if (!this.drag) return;
      this.drag = false;
      this.onPick(this.index, false);
    };
    c.addEventListener('pointerup', up);
    c.addEventListener('pointercancel', up);
    c.addEventListener('dblclick', (e) => { e.preventDefault(); this.full = !this.full; this.draw(); });
    c.addEventListener('keydown', (e) => {
      if (!this.tree || !['ArrowLeft', 'ArrowRight', 'ArrowUp', 'ArrowDown'].includes(e.key)) return;
      e.preventDefault();
      if (e.key === 'ArrowUp' || e.key === 'ArrowDown') {
        // switch between the main path and the branch, at the nearest path coordinate x
        const S = this.tree.samples, cur = S[this.index];
        if (!cur) return;
        const onBranch = cur.edge === 3, other = this.tree.samples.map((_, k) => k).filter((k) => (onBranch ? S[k].edge === 2 : S[k].edge === 3));
        if (!other.length || (!onBranch && cur.x < 2)) return;
        emit(other.reduce((b, k) => (Math.abs(S[k].x - cur.x) < Math.abs(S[b].x - cur.x) ? k : b), other[0]), false);
        return;
      }
      // along the branch the dot is on (at candidate 1 the right arrow follows the genus 2 edge)
      const path = this._path(this.index), d = e.key === 'ArrowRight' ? 1 : -1;
      const at = Math.max(0, path.indexOf(this.index));
      emit(path[Math.max(0, Math.min(path.length - 1, at + d * (e.shiftKey ? 10 : 1)))], false);
    });
  }

  draw() {
    const c = this.canvas, ctx = this.ctx;
    const dpr = Math.min(window.devicePixelRatio || 1, 3);
    const Wd = Math.round(c.clientWidth * dpr), Hd = Math.round(c.clientHeight * dpr);
    if (!Wd || !Hd) return;
    if (c.width !== Wd || c.height !== Hd) { c.width = Wd; c.height = Hd; }
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    const css = getComputedStyle(c);
    const col = (n, f) => css.getPropertyValue(n).trim() || f;
    const C = { ink: col('--ink', '#27231E'), muted: col('--muted', '#6A6257'), rule: col('--rule', '#CFC5B0'), accent: col('--accent', '#9E3F24'), shade: col('--shade', '#EAE4D6'), blue: col('--blue', '#3D5A80') };
    const FONT = "'EB Garamond', Georgia, serif";
    const { w, h, pad, lo, hi, X, Y } = this._geom();
    ctx.clearRect(0, 0, w, h);
    ctx.globalAlpha = this.tree ? 1 : 0.45;

    // frame, genus columns, W ticks, the 8 pi line
    ctx.strokeStyle = C.rule; ctx.lineWidth = 1;
    ctx.strokeRect(pad.l + 0.5, pad.t + 0.5, w - pad.l - pad.r - 1, h - pad.t - pad.b - 1);
    ctx.font = `italic 11px ${FONT}`; ctx.fillStyle = C.muted; ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    ctx.beginPath(); ctx.moveTo(X(2) + 0.5, pad.t); ctx.lineTo(X(2) + 0.5, h - pad.b); ctx.stroke();
    ctx.fillText('genus 1', X(1.5), h - pad.b + 5);
    ctx.fillText('genus 2', X(2.5), h - pad.b + 5);
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    // ticks at a step of 1, 0.5, 0.2 or 0.1 (pi^2), about four of them
    const step = [1, 0.5, 0.2, 0.1, 0.05].find((d) => (hi - lo) / d >= 3) || 0.05;
    for (let k = Math.ceil(lo / step); k * step <= hi; k++) {
      const v = k * step;
      ctx.fillText(v ? `${+v.toFixed(2)}π²` : '0', pad.l - 4, Y(v * Math.PI ** 2));
    }
    const y8 = Y(8 * Math.PI);
    ctx.setLineDash([3, 3]); ctx.beginPath(); ctx.moveTo(pad.l, y8); ctx.lineTo(w - pad.r, y8); ctx.stroke(); ctx.setLineDash([]);
    ctx.textAlign = 'left'; ctx.textBaseline = 'top';
    ctx.fillText('8π', pad.l + 3, y8 + 1);
    if (!this.tree) return;
    ctx.save();
    ctx.beginPath(); ctx.rect(pad.l - 7, pad.t - 12, w - pad.l - pad.r + 7, h - pad.t - pad.b + 12); ctx.clip();

    // the edges
    const S = this.tree.samples;
    ctx.strokeStyle = C.ink; ctx.lineWidth = 1.5;
    for (let g = 1; g < 3; g++) {
      // genus 1 starts at the opening
      const e = S.filter((s) => s.edge === g || (g === 1 && s.edge === 0 && shown(s)));
      if (!e.length) continue;
      ctx.beginPath();
      e.forEach((s, i) => (i ? ctx.lineTo(X(s.x), Y(s.W)) : ctx.moveTo(X(s.x), Y(s.W))));
      ctx.stroke();
    }
    // the branch that does not open the node (dashed, from candidate 1)
    const br = S.filter((s) => s.edge === 3), c1n = this.tree.nodes.find((n) => n.label === 'candidate 1');
    if (br.length && c1n) {
      ctx.save();
      ctx.strokeStyle = C.blue; ctx.setLineDash([4, 3]);
      ctx.beginPath();
      [S[c1n.index], ...br].forEach((s, i) => (i ? ctx.lineTo(X(s.x), Y(s.W)) : ctx.moveTo(X(s.x), Y(s.W))));
      ctx.stroke();
      ctx.restore();
    }
    // nodes
    ctx.font = `italic 11px ${FONT}`; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom';
    for (const n of this.tree.nodes) {
      const s = S[n.index], x = X(s.x), y = Y(s.W);
      ctx.fillStyle = '#fff'; ctx.strokeStyle = C.ink; ctx.lineWidth = 1.2;
      ctx.beginPath(); ctx.arc(x, y, 3.5, 0, 2 * Math.PI); ctx.fill(); ctx.stroke();
      ctx.fillStyle = S[n.index].edge === 3 ? C.blue : C.muted;
      ctx.textAlign = x > w - 50 ? 'right' : x < pad.l + 25 ? 'left' : 'center';
      // the end of the branch is labelled below its node (it can be close to candidates 2 and 3)
      const below = S[n.index].edge === 3;
      ctx.textBaseline = below ? 'top' : 'bottom';
      ctx.fillText(n.label, x > w - 50 ? x - 4 : x < pad.l + 25 ? x + 3 : x, below ? y + 5 : y - 5);
    }
    // the dot
    if (this.index >= 0 && this.index < S.length) {
      const s = S[this.index];
      ctx.fillStyle = C.accent;
      ctx.beginPath(); ctx.arc(X(s.x), Y(s.W), 5, 0, 2 * Math.PI); ctx.fill();
    }
    ctx.restore();
  }
}
