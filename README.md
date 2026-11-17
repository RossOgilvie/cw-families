# Willmore candidates: walking the tree

An interactive page that walks a "tree" of sinh-Gordon spectral data at a fixed internal τ and draws
the candidates for constrained Willmore tori that it meets. The edges of the tree are deformations of
the spectral curve within one genus that keep τ fixed. At a branching, a marked resonance pair arrives
over the unit circle as a double point, which is inserted and opened into the next genus. The tori are
computed in ℝ⁴ by the generalized Sym formula and projected to ℝ³.

`about.html` ("Learn more" in the panel) explains the construction for readers. The mathematics,
with Ross's sketch, the findings F1–F15 and the open questions, is in `cw_candidates.md`; the code,
how the tree is computed and the next steps are in `PLAN.md`.

The repository started as a copy of cmc-immersion (`~/programming/web/cmc-immersion`), whose
`cmc_spectral_to_immersion.md` is the specification of the inherited CMC code in `src/cmc/cmc.js`.

## Running

No build step. Serve the folder statically and open `index.html`:

    python3 serve.py 8765            # then visit http://localhost:8765

`serve.py` is `python3 -m http.server` plus `Cache-Control: no-cache`, so the browser never mixes an
edited module with a stale cached one. Opening the file directly doesn't work, because browsers block
module workers on `file://`.

A link `#k=Re,Im,2` computes candidates 1, 2 and 3 for that internal τ, and `#k=Re,Im` candidate 1.
*Copy link* produces such a link for the current candidate.

## Controls

- **Viewer:** drag rotates, shift+drag rolls, wheel or `+`/`−` zooms, `WASD` moves in the screen
  plane, `Q`/`E` toward/away, `F` frames the surface, `R` resets. Touch: drag, pinch, twist.
- **Candidates:** the internal τ, as numbers or as a point of the F₆ pad (drag it; the latest
  request wins). *Find candidates 2, 3* and *Find candidate 1* compute the tree. The picker chooses
  the pair of Sym points at the current point of the tree.
- **Tree view:** W (on the closing lattice) along the path: the opening on the left axis, the genus 1
  edge on the first third, the genus 2 edge on the rest. Drag the dot (it snaps to the nodes) or use
  the arrow keys (shift: ten steps); double-click shows W from 0. The note below gives the factors
  c₁, c₂ with f(z + eₖ) = cₖ f(z) on the drawn cell (1 at the candidates).
- **λ-plane:** the branch points α; Sym points that are branch points are drawn blue, other Sym
  points are blue squares. Dragging a branch point moves along the tree (in genus 2 to the curve's own
  τ and scale). Shift+drag turns everything; double-click adds or removes a branch point (this leaves
  the tree).
- **Appearance:** the projection ℝ⁴ → ℝ³, the texture, parameter lines, transparency, colours,
  resolution and accuracy.
- **Presets:** the Clifford torus, and candidates at three internal τ.

## Layout

    src/cmc/tree.js        the tree: internal frame, genus 1 and 2 edges, candidates, treeOf
    src/cmc/resonance.js   the map q and resonance pairs
    src/cmc/gensym.js      the generalized Sym immersion, projections R^4 -> R^3, shape diagnostics
    src/cmc/cmc.js         Killing fields and frame integration (inherited)
    src/cmc/periods.js     periods of the real-normalized differentials
    src/cmc/whitham.js     genus 2 period vectors (ellVector); the Whitham curve (tests only)
    src/cmc/leaf.js        W by the residue formula (willmoreFrame)
    src/cmc/genus2.js      the open-end tau of genus 2 curves
    src/cmc/poly.js        complex polynomials and roots
    src/cmc/rootflow.js    the root-preserving flow (tests only)
    src/app.js             state, links, the panel, the candidate worker
    src/worker.js          surfaces and candidates off the main thread
    src/tree-view.js       the tree widget
    src/spectral-widget.js the λ-plane
    src/tau-pad.js         the F₆ pad
    src/viewer.js          three.js viewer
    test/run.mjs           tests: node test/run.mjs
    about.html             the "Learn more" page (maths typeset by Temml as MathML)
    vendor/                three.js; temml/: Temml (MIT); fonts/: EB Garamond and Garamond-Math (SIL OFL)
