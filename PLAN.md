# cw-families: plan and status

State at the end of 17 November 2026. The mathematics (Ross's sketch, findings F1–F15, open questions) is in `cw_candidates.md`; this file is about the code: what it does, how, and what comes next.

## 1. Goal

1. Walk the "tree" of `cw_candidates.md` through sinh-Gordon spectral data at a fixed internal τ.
   1. Edges: τ-preserving flows within one genus.
   2. Branchings: a marked resonance pair arrives at S^1 as a σ-node, which is inserted and opened into the next genus.
2. At each candidate, draw the torus.
   1. It is computed in \bbR^4 by the generalized Sym formula (F7), and projected to \bbR^3.
   2. Before the Sym points are ripe, the same formula gives an immersed plane (F13).
3. Later: vary the internal τ (and possibly the external τ), to see how the candidates change and where they exist.

## 2. Running and testing

1. The repository is a copy of cmc-immersion (at its commit 11eb1d9).
   1. Plain ES modules, no build step, three.js vendored in `vendor/`.
   2. The CMC derivation that `cmc.js` cites by section number is `cmc_spectral_to_immersion.md` in `~/programming/web/cmc-immersion/`.
2. To serve: `python3 serve.py 8765`, then open http://127.0.0.1:8765/.
   1. A link `#k=Re,Im,2` computes candidates 2 and 3 for that internal τ.
   2. `#k=Re,Im` computes candidate 1 only.
3. Tests: `node test/run.mjs`, 184 checks in about 10 s. They cover F7 (generalized Sym), q and resonance (F3, F6), the tree (candidates 1–3, the open end, double arrivals, unripe planes, `locateOnTree`), and the inherited CMC code.
4. Headless checks of the page:
   1. `chromium --virtual-time-budget` never finishes module workers that import other modules.
   2. So use a real-time DevTools-protocol script instead. Node has `WebSocket` and `fetch`, so no packages are needed.
   3. Drags need `Input.dispatchMouseEvent`, because synthetic pointer events cannot take pointer capture. `window.cmcWidget` exposes the λ-plane widget for this.

## 3. Code map

1. `src/cmc/gensym.js`: the generalized Sym immersion, and projections.
   1. `PairIntegrator`: F at two λ ∈ \bbC^*, with `stepCap` = 0.25/|X|. Near the open end |X| ≈ 0.71/t.
   2. `leftEigenvector`, `computeGenSurface`.
   3. Projections \bbR^4 → \bbR^3 (`projectPoints`):
      1. automatic: stereographic if the surface lies on a round 3-sphere, otherwise Möbius-balanced;
      2. Möbius-balanced: conformal barycenter in S^4, stereographic from the axis of least spread;
      3. principal axes;
      4. stereographic from e_4;
      5. drop x_4.
   4. `roundnessS4`: the distance from a round 2-sphere and the S^4 area after a Möbius normalization.
2. `src/cmc/resonance.js`: the map q and resonance pairs.
   1. `QMap`: q = (log μ_1, log μ_i) by graded Gauss–Legendre quadrature from α_1; `resonance()` is Newton with lattice reduction.
   2. `genus0Resonance`: the closed forms of F3.
   3. `classifyPair`: off, on or node, with a signed distance. A pair within 1e-8 of S^1 counts as on it.
3. `src/cmc/tree.js`: the tree.
   1. `internalFrame`: the internal τ as a point of F_6.
   2. `genus1Edge`, `candidate1`, `genus2Anchor`, `genus2Edge`, `candidate2`.
   3. `candidate`: the characters, the closing lattice, W and the external τ.
   4. `refineArrival`: regula falsi.
   5. `treeOf`: the drawable path of samples.
   6. `paritySublattice`.
   7. `locateOnTree`: the τ and t of a genus 2 curve.
   8. `openEndShapes`.
4. Spectral-curve modules:
   1. `src/cmc/periods.js`: periods of the real-normalized differentials Θ_w.
   2. `src/cmc/whitham.js`: `ellVector`, the genus 2 period vectors. `WhithamCurve` is only used by tests and the worker's unused `family` message.
   3. `src/cmc/leaf.js`: `leafData`, `willmoreFrame` (W by the residue formula).
   4. `src/cmc/genus2.js`: `openEndTau`.
5. Inherited CMC modules:
   1. `src/cmc/cmc.js`: Killing fields and `FrameIntegrator`.
   2. `src/cmc/poly.js`.
   3. `src/cmc/rootflow.js`: tests only.
6. The page:
   1. `src/app.js`: state, links, the candidate worker, `showSample`, `treeDrag`, `treeJump`.
   2. `src/worker.js`: the messages `candidate1`, `candidate2` (with `fast`), `shapes`, and surfaces.
   3. `src/tree-view.js`: the tree widget.
   4. `src/spectral-widget.js`: the λ-plane.
   5. `src/tau-pad.js`: the F_6 pad.
   6. `src/viewer.js`: three.js, including `setMarkers` for the spikes.
   7. `src/period-pad.js`: unused since its panel was removed.

## 4. How the tree is computed (`candidate2(τ)`)

1. **Internal frame.**
   1. The classes 1 = [1], 2 = [τ], 3 = [1 + τ] are ordered by the lengths of their shortest vectors. The shortest is opened first.
   2. The computation runs with the reduced τ1 of the rescaled lattice.
   3. Convention: Λ^* = t i ⟨1, τ1⟩ and Γ = 2πi Λ^*. The factor i puts the σ-node of the pair for 1 at λ = −1.
2. **Genus 0** (closed form).
   1. The scale t runs from the vacuum to t = 1/π, where the pair for 1 is a σ-node at λ = −1 (F4).
   2. On the closing lattice, W = 2/(t² Im τ1).
3. **Genus 1 edge** (`genus1Edge`).
   1. α = −r, starting at r = 1 and stepping down by 20%, with a linear predictor of the pairs in log r.
   2. t(r) comes from a B-period.
   3. The two remaining marked pairs are re-solved at each step.
   4. An arrival (a passage through a σ-node, onto S^1 or off it) is refined by `refineArrival`, in 10–15 solves.
   5. A double arrival, where the other pair is a σ-node at the same point (F14), is flagged.
4. **Candidate 1** (`candidate`).
   1. The monodromy characters at both Sym points are found by shooting along the lattice, at steps of stepCap/4.
   2. The closing lattice is doubled when the characters differ (F8).
   3. Then W by the residue formula, the external τ, and the closing and node defects (about 1e-6).
5. **Genus 2 edge** (`genus2Edge`).
   1. Anchor: the node is opened by ε = 1e-3, with the period ratio ell_2/ell_1 = conj(s τ1 + n) solved exactly.
   2. Gauge: α_1 = t β e^{iψ}, α_2 = (t/β) e^{iψ}, where t = sqrt|α_1 α_2|.
   3. Steps of 0.1 in log t; at each step, Newton in β for the ratio.
   4. The remaining pair is lifted at the node (ν_2 = ν_1 sqrt(−conj λ_0)(λ − λ_0)) and continued.
   5. It ends at the arrival (candidates 2 and 3) or at tMin = 1e-3 (the open end).
6. **Candidates 2 and 3.**
   1. At the arrival: the three choices of two from α_1, α_2 and the node. Candidate 2 is α_1 + α_2; candidate 3 is α_1 + the node.
   2. At a double arrival: the three choices from α, node a and node b, on the genus 1 curve.
   3. At an open end:
      1. candidate 2 at t = 0.03, 0.01, 0.003, 0.001;
      2. the frame is integrated only at 0.03, and the others reuse its parity;
      3. the shape diagnostics only on request.
7. **`treeOf`**: one path of samples.
   1. It runs from t = 10/π in genus 0, through the opening and candidate 1, to the end of the genus 2 edge.
   2. Each sample has views (Sym points, the cell of the eventual candidate, and ripe or not).
8. **Timings** (Node, warm):
   1. about 60–100 ms with an arrival;
   2. about 0.2–0.3 s at an open end;
   3. the drawn surface in the app's worker: about 0.3–1 s at resolution 240.

## 5. The app

1. **Candidates section.**
   1. Internal τ: fields and the F_6 pad. Dragging recomputes, and the latest request wins.
   2. Buttons: Find candidates 2, 3, and Find candidate 1.
   3. The view picker: the Sym-point pairs at the current sample.
   4. The tree view.
   5. The sample note, with the factors c_k: f(z + e_k) = c_k f(z) on the cell.
   6. The candidate note, and the Check the limit shape button at open ends.
2. **Tree view.**
   1. W on the closing lattice, in units of π², against the path. The opening is on the left axis, genus 1 takes the first third of the width and genus 2 the rest.
   2. The nodes: opening, candidate 1, candidates 2 and 3 (or "candidates 1, 2, 3" at a double arrival), open end.
   3. Interaction: the dot is dragged and snaps to the nodes, or stepped with the arrow keys; shift steps by 10. Double-click shows W from 0.
3. **The dot drives the drawing.**
   1. The 3D view shows the sample.
   2. The λ-plane shows the path along the current edge.
   3. Near an open end, spikes mark the four half-periods (the branch points of the limit, F12).
4. **λ-plane.**
   1. The Sym points are not draggable.
      1. A branch point that is a Sym point is drawn blue.
      2. Any other Sym point is a blue square λ_j: on S^1 a σ-node, elsewhere a pair not yet on S^1.
   2. Dragging a branch point moves along the tree.
      1. In genus 2, `locateOnTree` gives the curve's τ and t. The tree is recomputed live, and candidate 2 of the dragged curve is drawn.
      2. In genus 1, τ is not determined by the curve, so the drag moves along the edge.
5. **Changes between nearby τ** (`treeJump`), shown at the top of the note:
   1. a wall of F_6 crossed (the class order changes, and the candidates continue);
   2. a change of the kind of tree;
   3. a jump of W at a node.
6. **Appearance and presets.**
   1. The projection \bbR^4 → \bbR^3 is in Appearance.
   2. Presets:
      1. the Clifford torus;
      2. candidates at 0.45 + 0.95i;
      3. a double arrival at 0.5 + 0.9i;
      4. an open end at 0.2 + 1.3i.
   3. Without a link, the page starts at 0.45 + 0.95i.

## 6. Done (condensed log; details in the commit messages)

1. Generalized Sym immersion and projections (F7, F9, F11).
2. q and resonance pairs (F3, F6).
3. The tree to candidate 1, and the internal τ in F_6.
4. Genus 2:
   1. The genus 2 edge and candidates 2 and 3.
   2. The edge parametrized by t and followed to the open end. The open-end τ of the rescaled roots tends to the internal τ with error exactly O(t).
5. The open-end limit: a branched double cover of the round sphere (F12), with spikes at the branch points.
6. Unripe immersions along the whole path: planes with homotheties (F13), and the tree view.
7. Double arrivals on the walls of F_6 (F14). Arrivals now include passages from on S^1 to off it.
8. Interface:
   1. The CMC panels and presets were removed.
   2. The generalized formula is always on.
   3. Fixed Sym-point symbols, and dragging along the tree.
9. Speed: candidate frames at stepCap/4, regula falsi for the arrivals, 20% steps in genus 1, and the open-end diagnostics on request.
   1. A tree with an arrival: 344 → about 60 ms.
   2. An open end: 2–5 s → about 0.25 s.
10. Tried and dropped: continuing the tree in τ (starting guesses from the previous tree).
    1. The cost is the two full q integrals per solve, not the Newton iterations.
    2. Measured gain 8–15%.

## 7. Next steps

1. **The "(or not)" branches.**
   1. Continue the genus 1 edge past candidate 1 without opening, to the arrival of the third pair.
   2. Decide how branches spread horizontally in the tree view.
2. **Curves in F_6** (next direction 5). Draw on the τ pad where candidates 2 and 3 exist and where the edge ends at an open end.
   1. The transition is continuous: the arrival goes to t = 0 (F15).
   2. It can be computed directly as the curve where the arrival is at t = 0, or from the real locus of the genus 2 paper, and checked against a sparse grid.
   3. Also mark the walls of F_6 where double arrivals happen.
3. **Layout** (next direction 6): one quarter for the controls, and three quarters for candidates 1, 2 and 3. These are one canvas with three viewports.
4. **4D visualization** (next direction 3), for example:
   1. colour by the coordinate dropped by the balanced projection;
   2. a "4D rotation" drag that moves the projection pole on S^4.
5. **Smoother dragging:**
   1. a lower resolution while dragging;
   2. cheaper q integrals (fewer nodes on the graded mesh); these are now the main cost of an edge.
6. **Questions for the mathematics** (see `cw_candidates.md`):
   1. The ℝ^3 hypothesis: the choice "node a + node b" at a double arrival is a CMC torus in S^3 (F14).
   2. A sweep of W over F_6.
   3. The external τ.
7. **Housekeeping:**
   1. Rewrite the README, the page title and subtitle ("CMC surfaces"), and `about.html`; they still describe cmc-immersion.
   2. Remove `src/period-pad.js` and the worker's `family` message, both unused.
   3. Decide about the deleted `cmc_spectral_to_immersion.md`: it is deleted in the working tree but not committed, and `cmc.js` cites it.

## 8. Decisions taken

1. The surface is computed in \bbR^4, because generically the target is S^4 (F9, F11).
2. The internal τ (of Γ) is fixed and given as a point of F_6. The external τ (of the torus) is derived and reported.
3. Take the double cover whenever necessary; it is always necessary when the two Sym points are different marked pairs (F8).
4. Openings place the branch points near the node, as a ρ-pair inside the disc (F5).
5. Genus 1 has its own edge code; genus 2 uses the t-parametrization, not `WhithamCurve`.
6. Only the generalized Sym formula is used.
   1. The Sym points follow from the tree and are not editable.
   2. An unripe Sym point uses the point P of its pair, not Q = ρσP (F13).

## 9. Open questions (implementation)

1. The genus 1 edge has no global parameter other than r (or W, which is monotone along every path tried so far).
2. The tree view's horizontal axis is a path parameter: genus 1 by log r, genus 2 by log t. A more meaningful axis would help once branches are drawn.
3. Pairs tracked along S^1 after an arrival (the "on" pairs) are followed only to detect their passage through a σ-node.
