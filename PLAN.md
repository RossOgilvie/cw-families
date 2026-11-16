# Implementation plan: walking the tree of constrained Willmore candidates

Goal: an interactive page that walks the "tree" of `cw_candidates.md` through the space of sinh-Gordon spectral data at fixed conformal class τ.
The edges are τ-preserving Whitham flows within one genus, and the branchings insert and open nodes where a marked resonance pair arrives at S^1.
At each candidate the page draws the immersed torus, computed in \bbR^4 by the generalized Sym formula (F7 of `cw_candidates.md`) and projected to \bbR^3.
The mathematics is in `cw_candidates.md` (sketch, findings F1–F10 and open questions) and, for the inherited CMC pipeline, in `cmc_spectral_to_immersion.md`.

The code base is a copy of cmc-immersion at commit 11eb1d9 (its plan is `PLAN.md` in the first commit of this repository).
Architecture, conventions and tooling are inherited: plain ES modules, no build step, three.js vendored, computation in a worker, `node test/run.mjs` for tests.

---

## 1. What is reused

1. `src/cmc/cmc.js`: Killing fields, isospectral flows, grid integration of the frame.
   The frame integration is generalized from (F, ∂_λF) at one λ_0 ∈ S^1 to F at two arbitrary λ_1, λ_2 ∈ \bbC^*.
2. `src/cmc/periods.js`: periods of the differentials by quadrature; the basis of the evaluation of q along paths.
3. `src/cmc/whitham.js`, `src/cmc/leaf.js`: the τ-flow (W varies, τ fixed) for genus ≥ 2.
   In genus 2 the length of the integral period differential in the b(0) coordinates is a global parameter of the τ-flow; W is one only away from \calS.
4. `src/cmc/genus2.js`: opening a node on S^1 to a ρ-pair of branch points.
5. `src/viewer.js`, `src/worker.js`: the 3D viewer and the worker, with positions now coming from \bbR^4 through a projection.
6. The λ-plane widget, for showing branch points, marked pairs and Sym points.

The CMC-specific parts (Sym point slider, mean curvature, curvature-line grid, bubbletons, presets) stay available while the new mode is built, and are removed or hidden once they are no longer needed.

## 2. New ingredients

1. **Generalized Sym immersion** (F7).
   Input: two points of the spectral curve (λ_j with a sheet), the eigenvectors u_j of ζ(0, λ_j), and the grid.
   Output: f = Q(ψ_1) Q(ψ_2)^{-1} ∈ \bbR^4 per vertex, with ψ_j = u_j^T F(z, λ_j).
   The derivative ∂_λF is no longer needed.
   Normals in \bbR^3 are computed from the projected mesh, since the projection is not conformal in general.
2. **Projections \bbR^4 → \bbR^3.**
   First: drop the last coordinate.
   Later: stereographic projection from a chosen point of S^4, orthogonal projection along a chosen direction, and a rotation of \bbR^4 driven from the UI.
3. **The map q and resonance pairs.**
   Evaluate q = (q^1, q^2) at points of the spectral curve along paths from a branch point, modulo Γ, with the quadrature of `periods.js`.
   Solve q(P) − q(Q) = γ by Newton, starting from the genus 0 closed forms (F3) and continued along each edge.
   Classify each marked pair as off S^1 (ρ-pair), on S^1 over two values of λ, or a σ-node (the event).
4. **Genus 1 τ-flow.**
   `whitham.js` requires g ≥ 2 because it works with the period plane in Gr(2, \bbR^g); in genus 1 one differential is exact instead.
   A small separate module: one ρ-pair α, 1/conj(α), τ fixed, parametrized by W (to be checked to be global in genus 1).
5. **Events and openings.**
   Detect the arrival of a marked pair at S^1 along an edge (a sign change of a distance function, then bisection).
   Insert the node there and open it (F5) into genus g + 1, then continue with the τ-flow of genus g + 1.
6. **Candidates.**
   At an event, choose the two trivial-multiplier ρ-pairs (F7, F8), find the closing lattice (original torus or the index 2 sublattice), and draw one fundamental domain.
   Report W both from the residue formula and from the numerical integral ∫ |Δf|^2 / (4 |f_x|^2) dx dy.
7. **The tree in the UI**, drawn as on the whiteboard.
   The edges are curves in a plane whose vertical axis is the Willmore energy: W at a torus (a candidate), and the normalized Willmore energy elsewhere (the residue formula divided by the signed area of the parallelogram of B-periods of b_1, b_2, as in the Lemma "willmore" of the genus g paper).
   The horizontal axis is a parameter along the edges, laid out so that branches separate (to be designed once the edges can be computed).
   A dot moves along the tree by dragging and snaps to the bifurcation points and candidates; the 3D view and the λ-plane follow the dot.

## 3. Milestones

Each milestone ends in a check that can be verified.

1. **Generalized Sym immersion in the existing app.**
   A mode in which the user picks two points λ_1, λ_2 in the λ-plane (with sheets), drawn by dropping the last coordinate of \bbR^4.
   Tests in `test/run.mjs`, ported from the Python checks of F7: conformality for arbitrary λ, the homogeneous tori W = π^2 (t + 1/t), and agreement with the S^3 Sym formula for λ_j ∈ S^1.
2. **q and resonance pairs.**
   Tests: the genus 0 closed forms of F3 (two solutions per γ, the threshold |γ| = 4, the σ-node at the threshold), and q(B) ∈ ½Γ at branch points in genus 1 and 2 (F6).
3. **Genus 0 edge and the first opening.**
   The genus 0 scaling from the vacuum to the starting scale (F4), and the opening of the pair for 1 into genus 1.
   Check: the class of the new branch points is [1], and the genus 1 curve right after opening is close to the nodal one.
4. **Genus 1 τ-flow with event detection.**
   Checks: τ stays fixed to quadrature accuracy, W is monotone along the edge, and the formulas of `genus_one.tex` hold in the overlapping normalization (to be matched; the notes normalize the Fermi lattice, not the period lattice).
   Output: candidate 1, as an event on the edge.
5. **Candidate 1 drawn.**
   Generalized Sym points from the two trivial-multiplier pairs, the closing lattice, the torus in \bbR^4 and its projection.
   Checks: closing by shooting (f(z + w) = f(z) for the lattice generators), and W from the residue formula against the numerical integral.
6. **Second opening and candidate 2.**
   The opening of genus 1 into genus 2 at candidate 1, the τ-flow in genus 2 (`leaf.js`), and candidate 2.
   Check: the τ-flow passes through the nodal curve continuously (the paper shows the leaf passes through U^g).
7. **Tree view and projections.**
   The tree drawn with the Willmore energy on the vertical axis and a snapping dot (§2.7), more projections, and per-candidate views.
8. **Third candidate** (the doubling construction, F8), and tracking of pairs along S^1 (the "(or not)" branches), for further experiments.

## Status

1. Milestone 1 done (November 2026), apart from the README.
   `src/cmc/gensym.js` has the pair integrator (F at two λ ∈ \bbC^*, with the layout and marching of `FrameIntegrator`), the eigenvectors of ζ(0, λ) on a chosen sheet, the quaternion formula and `computeGenSurface`.
   `test/run.mjs` checks conformality for arbitrary λ (to 1e-15), agreement with the S^3 Sym formula, closing and W = π^2 (t + 1/t) for three homogeneous tori (to 4e-6, the error of the finite differences), and the grid against direct integration.
   In the app, the checkbox *Two Sym points* in Spectral data switches the surface to the generalized formula.
   λ_1, λ_2 are blue squares in the λ-plane: draggable anywhere, snapping onto the branch points and their reflections (and following a branch point they sit on), with a double-click flipping the sheet.
   Projections: automatic (the default: stereographic projection of the 3-sphere when the surface lies on one, from its point farthest from the surface, and otherwise Möbius-balanced), Möbius-balanced (f read in S^4 and moved to its conformal barycenter, then stereographic projection from the principal axis of least spread, keeping the three axes of most spread; a round 2-sphere becomes the unit sphere), principal axes of f in R^4, stereographic from e_4, and dropping x_4.
   The state is in the link (`y`, `j`), and there are three presets (the Clifford torus, a homogeneous torus of ratio 2, and a branch point as one Sym point).
   The panels Mean curvature and Root-preserving flow are removed (with `src/sym-pad.js`); `src/cmc/rootflow.js` stays, as its tests still run.

2. Milestone 2, first part done (November 2026): `src/cmc/resonance.js`.
   q = (q^1, q^2) = (log μ_1, log μ_i) is the primitive of the real-normalized Θ_1, Θ_i of `periods.js`, based at the branch point α_1 (q∘σ = −q), in closed form in genus 0 (q_w = (w ν/λ + conj(w) ν)/κ_0, in the code's normalization ν^2 = −λ, κ_0 = 2).
   It is computed by Gauss–Legendre quadrature on segments, with ν continued and the square-root singularities at branch-point ends removed.
   Resonance pairs: Newton in (λ_P, λ_Q) with q carried along incrementally; `genus0Resonance` (closed form: ν_P − ν_Q = 2πi k, ν_P ν_Q = k/conj(k), over S^1 if and only if |k| ≤ 1/π); `classifyPair` (off, on, node, with a signed distance to arrival).
   In this normalization the torus lattice Λ is in the z-plane and γ = 2πi (k_1, k_2) with k ∈ Λ^*; the sketch's Γ is 2πi Λ^*, so its scale is that of Λ.
   Tests: in genus 0, q equals the eigenvalues of X(λ) on the eigenline used by `gensym.js` (both sheets), the closed forms, Newton, the threshold and classification; in genus 1 and 2, Re q = 0 and q(α) = q(1/conj α) at branch points, and 2q(α_2) is a period of B_1 ± B_2.
   Continuing marked pairs along a deformation of the curve is done with the genus 1 edge (item 3), using the lattice reduction of `resonance()`.
3. Milestones 3–5 for candidate 1 done (November 2026): `src/cmc/tree.js`.
   Design: the internal τ is an input of every function (and the solvers take warm starts), so that the tree can later be recomputed along a path of τ; the external τ of a candidate is derived and reported.
   Conventions: Λ^* = t i ⟨1, τ⟩; the factor i puts the σ-node of the pair for 1 at λ = −1, where the opened genus 1 curves α = −r continue ν^2 = −λ without rotating the z-plane (at a double root λ_0 the nodal curve normalizes to ν^2 = conj(λ_0) λ).
   Genus 0 edge in closed form (W = area of Λ = 1/(t^2 Im τ), arrivals at π|k| = 1, the σ-node at λ = π^2 k^2).
   Genus 1 edge: t(r) from the B-period, the τ and third pairs continued from their genus 0 positions (ν_1 = ν_0 (λ + 1) at the node), arrivals refined by bisection.
   Candidates: monodromy characters at both Sym points by shooting, the closing lattice (doubled when they differ), W by the residue formula on it, the external τ, and a closing check.
   In the app: the *Candidates* section (internal τ, *Find candidate 1*), computed in a worker, loads the curve, the two Sym points and one cell of the closing lattice; a link `#k=Re,Im` recomputes it.
   Result for internal τ = 0.2 + 1.3i: the τ pair arrives at α = −0.092161; the node has monodromy ±I; the characters differ, so the torus closes on the doubled lattice (to 2e-10); W = 23.5632 = 2.3874 π^2; external τ = −0.3133 + 1.1159i.
   A finite-difference W of the drawn torus converges towards the residue value (24.84, 24.05, 23.82, 23.67 at 48, 72, 96, 144 points per side).
   Open: the genus 1 continuation fails below r ≈ 0.05 (the α → 0 corner); after an arrival the coarse steps can jump to another branch (harmless for candidate 1, which stops at the first arrival).

4. Candidate 2 and the F_6 pad (November 2026).
   The internal τ is a point of F_6 (`internalFrame`): the classes 1 = [1], 2 = [τ], 3 = [1 + τ] keep their labels, the shortest class is opened first, and the computation runs with the reduced τ' of the rescaled lattice.
   The Candidates section has an F_6 pad: dragging recomputes the last requested candidate (latest request wins).
   Genus 2 edge (`genus2Edge`): the node of candidate 1 opened by ε = 1e-3 with the period ratio solved exactly (`genus2Anchor`: ell_2/ell_1 = conj(s τ' + n), the B-period through the node being the arriving vector modulo k_1), then `WhithamCurve` in the direction where α_2 leaves S^1, the remaining pair lifted by ν_2 = ν_1 sqrt(−conj λ_0)(λ − λ_0) and continued with the lattice of the current B-periods.
   At its arrival there are three trivial-multiplier ρ-pairs (α_1, α_2, the node); every choice of two closes on a doubled lattice, with the same W and different external τ.
   The paths of the branch points along the genus 1 or genus 2 edge are drawn in the λ-plane, coloured by W.
   Results: at 0.45 + 0.95i candidate 1 has W = 2.3131 π^2 and candidate 2 W = 2.4725 π^2, with external τ −0.0871 + 1.5768i (α_1 + α_2), 0.1855 + 1.7195i (α_1 + node), −0.1 + 1.9i (α_2 + node).
   At 0.2 + 1.3i there is no candidate 2 on the genus 2 edge: every tracked pair (vectors m + nτ, |m| ≤ 2, n ≤ 2) moves away from S^1 until the edge reaches the α → 0 corner.
   Open: the arrivals in genus 2 lie close to the α → 0 corner (|α| ~ 0.05–0.3), where WhithamCurve stops (newton); at 0.4 + i Node finds the arrival at s = 0.935 and the browser does not (the curve ends first), so the corner needs better numerics (the rescaled coordinate λ/t of the genus 2 paper). The genus 2 anchor Newton diverges for some τ (0.4 + 0.97i, near i), and for τ' on the imaginary axis (0 + 1.15i, 1.6i) the remaining pair is not continued into genus 2. The third pair at candidate 1 can already be on S^1 (0.5 + 0.9i), where an arrival from off S^1 is not the right event.

5. The open ends (November 2026).
   When the genus 2 edge reaches the open end (both branch points to 0; genus 2 paper, section "The Open Ends") before an arrival, `openEndInfo` rescales by t = sqrt|α_1 α_2|: the open-end τ of the rescaled root (genus2.js) tends to the internal τ with error O(t), and W_int to 4π (8π on the doubled lattice).
   `candidate2` then returns candidate 2 at the last point of the edge; the app says that candidates 2 and 3 become a branched double cover of the sphere there (Ross), and draws it.
   Observed at 0.2 + 1.3i (t = 0.038): the drawn torus is two nearly round spheres joined by a small neck, W = 0.988 · 8π, and its fourth principal variance in R^4 tends to 0 (0.15, 0.04, 0.03 along the end of the edge). Whether this is the branched double cover after a normalization (a Möbius transformation, or a rescaling about the neck) is not yet checked.
6. The genus 2 edge parametrized by t (November 2026), replacing WhithamCurve there.
   In the gauge of WhithamCurve (arg(α_1 α_2) fixed), α_1 = t β e^{iψ}, α_2 = (t/β) e^{iψ}, and the edge is ratio(t, β) = the internal ratio, solved by Newton in β and continued in log t; t decreases monotonically along the edge and β stays of size 1, the rescaled root of the open end (check-α = β).
   Numerics made this work to t = 1e-5 in under 0.1 s: a graded quadrature mesh in `QMap.segment`; the B-cycle differentials of periods.js near the open ends (the route through all cycles divides by a_0 ~ t^2) but not near a node on S^1 (where the B-cycle route has a nonzero real A-period); one damping factor in the resonance Newton respecting each point's own distance to the singular points (P runs out like 1/t while Q stays at scale t); the period vectors' signs aligned throughout.
   Results: the open-end τ of β agrees with the internal τ to 8.3e-3, 8.3e-4, 8.3e-5, 8.3e-6 at t = 1e-2 ... 1e-5 (exactly O(t)), and W_int/4π − 1 = −3e-6 at 1e-5: the edges end at the open end of the genus 2 paper, with its τ. All seven test values of τ now end either at an arrival (0.45 + 0.95i, 0.4 + i, 0.3 + 0.8i) or at the open end (0.2 + 1.3i, 1.6i, 0.75 + 0.82i, 2i; 0.5 + 0.9i was listed here wrongly, see item 9); before, 1.6i and 0.3 + 0.8i failed and the arrival at 0.4 + i was lost in the browser.
   Near the open end the drawn surface was taken at t = 0.03 (superseded by item 7).
   Open: an arrival of a pair that is already on S^1 at the start of the edge (0.5 + 0.9i) is not detected (only arrivals from off S^1). Resolved in item 9.
7. The limit surface at the open end (November 2026).
   The "stiffness" was the step size: the generator has |X(λ_j)| ≈ 0.71/t at the Sym points, and the torus lattice shrinks like t^{1/2}. Steps are now capped at 0.25/|X| (`PairIntegrator.stepCap`, in `candidate` and `computeGenSurface`), and candidate 2 is computed at t = 0.03, 0.01, 0.003, 0.001 (closing defects ≤ 4e-6). The picker chooses the scale; the default is 0.01.
   The shape is read Möbius-invariantly (`roundnessS4`): f in S^4, moved to its conformal barycenter (uniform weights on the torus), and fitted by a round 2-sphere (two linear relations among (X, 1)).
   At 0.2 + 1.3i, 1.6i and 2i: the distance from a round sphere is 0.34, 0.20, 0.11, 0.064 at t = 0.03 ... 0.001 (∝ t^{1/2}), and the S^4 area tends to 2 × 4π (always below W/4π), while no part of the torus carries a growing share of the area. So the limit is a round sphere covered twice, not two spheres; the two spheres with a neck of item 5 were the R^4 view of a surface still far from the limit.
   The branch points (minima of the conformal factor, which tend to 0 like t) are exactly the four half-periods of the doubled lattice, and f(−z) − f(z) → 0 like t^{1/2} (0.2 + 1.3i). So the limit is the quotient T^2 → T^2/±1 ≅ S^2 of the external torus (℘ of the external lattice up to a Möbius transformation); W = 8π.
   The external τ is constant along the edge (−0.3133 + 1.1159i at 0.2 + 1.3i).
   The Möbius-balanced projection is now the automatic one for surfaces on no 3-sphere. It also changes the look of ordinary candidates: candidate 2 at 0.45 + 0.95i appears as a torus without the two-coloured halves. The candidates are still on no 3-sphere in the balanced frame (largest distance from the best hyperplane 0.22–0.52 on the unit S^4, against 4e-8 for the Clifford torus), so F11 stands.
8. Unripe data and the tree view (November 2026; next directions 1 and 2).
   `treeOf` (tree.js) turns a candidate computation into one path of samples, from near the vacuum (t = 10/π, W = 0.02π²) through the opening, candidate 1, and the genus 2 edge to candidates 2 and 3 or the open end, with the immersion drawable at every sample.
   Sym points along the path: genus 0, the pair for 1 and the arriving pair a; genus 1, the branch point α and a; genus 2, three views (α_1 + α_2, which is candidate 2 and ripe along the whole edge; α_1 + the remaining pair b, towards candidate 3; α_2 + b).
   The cell drawn is the one on which the eventual candidate closes: by F8 it depends only on the classes of the two Sym points (`paritySublattice`), and this reproduces the lattices that candidates 1, 2, 3 and the other choice find from their characters.
   Before ripeness the immersion is a plane with f(z + w) = μ_1(w) f(z) μ_2(w)^{-1} (u_j is an eigenvector of the monodromy). The multipliers turn out to be real, so the lattice acts by homotheties of R^4 about 0, f(z + w) = c(w) f(z) with c real, and ripeness is c = 1 (checked to 1e-8 on all three edges, F13).
   W increases monotonically along the whole path (all test values of τ).
   The tree view (`src/tree-view.js`, milestone 7, first version) plots W on the closing lattice against the path: genus 0, 1, 2 on [0, 1], [1, 2], [2, 3]. The nodes are the opening, candidate 1, and candidates 2 and 3 or the open end. A dot is dragged along the path (snapping to the nodes) or stepped with the arrow keys, and the 3D view, the λ-plane and a note (with c_1, c_2) follow it. The y-axis shows the band of the genus 1 and 2 edges, and double-click shows it from 0. The view picker chooses among the Sym-point pairs on the genus 2 edge.
   Not yet in the tree: the genus 1 edge continued past candidate 1 without opening (to the arrival of the third pair).
9. Double arrivals (November 2026).
   On the walls of F_6 where the two classes after the opened one are equally long (Re τ = 1/2 in the chamber 123, and its images), both pairs arrive at S^1 together on the genus 1 edge. Then the genus 2 edge has length 0: approaching the wall, the arrival on the genus 2 edge moves to its start (t_arrival / t_start = 0.40, 0.64, 0.76, 0.83, 0.93 at Re τ = 0.45, 0.48, 0.49, 0.495, 0.499 with Im τ = 0.9), and W of candidates 2 and 3 tends to W of candidate 1.
   Before, the code opened one node, missed the other (it started just on S^1 and left it through its σ-node, which was not flagged), and followed the edge on to the open end; so 0.5 + 0.9i was wrongly reported as an open end.
   Now the genus 1 edge flags a double arrival (the other pair is a σ-node at the same point), and `candidate2` returns the three choices of two of {α, node a, node b} on the genus 1 curve: candidate 2 (α_2 → node a) is candidate 1, candidate 3 is α + node b, the other choice node a + node b. All three candidates have the same spectral curve and W and differ only in the Sym points (Ross).
   At 0.5 + 0.9i: W = 2.3505π² for all three; the external τ of candidate 3 is the mirror image of that of candidate 1 (∓0.0566 + 1.6981i), and the other choice is rectangular (1.8i = 2 Im τ · i). Just off the wall (0.499 + 0.9i) the three candidates agree with these to 2e-3. The same at −0.5 + 1.2i and 0.5 + 1.4i.
   Both edges now also treat a passage of a pair from on S^1 to off S^1 (through a σ-node) as an arrival.
   The third choice there (node a + node b) lies on a 3-sphere: a CMC torus in S^3 with rectangular external τ (F14).
10. Branch points near the open end (November 2026).
   Along a genus 2 edge that ends at an open end, the drawing of candidate 2 has spikes (`Viewer.setMarkers`) at the four half-periods of the cell, along the outward direction from the centroid: they are the minima of the conformal factor all along the edge (checked at t = 0.3, 0.1, 0.03, 0.003) and the branch points of the limit, so t = 0.1 or 0.01 already shows the branched double cover (Ross).

## TODO

1. Rewrite the README for cw-families: it still describes cmc-immersion (title, controls, layout, the removed panels).
2. The label of a λ_j overlaps the label of the branch point it sits on.
3. The page title and subtitle ("CMC surfaces") and `about.html` still describe cmc-immersion.
4. Editing the domain sliders while a candidate is shown has no effect (the candidate's cell is drawn); editing the spectral data leaves the candidate.
5. Headless checks: `chromium --virtual-time-budget` does not finish loading module workers that import other modules, so screenshots of worker results need a real-time wait (a small DevTools-protocol script).

Next directions (Ross, 16 November 2026):

1. Unripe data: the immersion formula before the resonance pair reaches S^1. Expected: an immersed plane, equivariant under the lattice (see below), as for Sym–Bobenko without closing conditions.
2. The immersion all along the tree, from the vacuum to the candidate.
3. Better visualizations of 4D.
4. Only the generalized Sym formula is wanted: clean up the λ-plane interface accordingly.
5. Precompute curves in F_6 showing where candidates exist and where the edges end at open ends.
6. A new layout: the screen in four quarters, one with the controls and three with the candidates.
7. Store the whole tree, for smoother and faster updates under changes of τ (reuse the existing tree, e.g. at open ends), and to see whether candidates jump as τ changes.

## 4. Decisions taken

1. The surface is computed in \bbR^4 (generically the target is S^4, F9); the first projection drops the last coordinate.
2. Openings place the branch points near the node, as a ρ-pair inside the disc (F5).
3. Genus 1 gets its own τ-flow module instead of extending `whitham.js`.
4. Tracking resonance pairs along S^1 after arrival is a later feature; it is not needed for the candidates.

## 5. Open questions

1. Which conformal class is fixed: that of Γ, or that of the doubled torus (open question 2 of `cw_candidates.md`)?
   Class matching itself is settled: take the double cover whenever necessary.
2. Normalization of q: the constant fixing the trivial-multiplier class (the spin structure) and the identification of q with the differentials p_w of `periods.js` and of Γ with the code's period lattice.
   This is an implementation matter to settle in milestone 2, against F3.
3. A global parameter for the genus 1 τ-flow (W, or the analogue of the genus 2 integral period length).
4. The layout of the tree view: the horizontal axis, and how the branches at a bifurcation are spread out.
