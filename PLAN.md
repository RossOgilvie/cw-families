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
   Projections: automatic (the default: stereographic projection of the 3-sphere when the surface lies on one, from its point farthest from the surface, and otherwise orthogonal projection onto the three principal axes), stereographic from e_4, and dropping x_4.
   The state is in the link (`y`, `j`), and there are three presets (the Clifford torus, a homogeneous torus of ratio 2, and a branch point as one Sym point).
   The panels Mean curvature and Root-preserving flow are removed (with `src/sym-pad.js`); `src/cmc/rootflow.js` stays, as its tests still run.

2. Milestone 2, first part done (November 2026): `src/cmc/resonance.js`.
   q = (q^1, q^2) = (log μ_1, log μ_i) is the primitive of the real-normalized Θ_1, Θ_i of `periods.js`, based at the branch point α_1 (q∘σ = −q), in closed form in genus 0 (q_w = (w ν/λ + conj(w) ν)/κ_0, in the code's normalization ν^2 = −λ, κ_0 = 2).
   It is computed by Gauss–Legendre quadrature on segments, with ν continued and the square-root singularities at branch-point ends removed.
   Resonance pairs: Newton in (λ_P, λ_Q) with q carried along incrementally; `genus0Resonance` (closed form: ν_P − ν_Q = 2πi k, ν_P ν_Q = k/conj(k), over S^1 if and only if |k| ≤ 1/π); `classifyPair` (off, on, node, with a signed distance to arrival).
   In this normalization the torus lattice Λ is in the z-plane and γ = 2πi (k_1, k_2) with k ∈ Λ^*; the sketch's Γ is 2πi Λ^*, so its scale is that of Λ.
   Tests: in genus 0, q equals the eigenvalues of X(λ) on the eigenline used by `gensym.js` (both sheets), the closed forms, Newton, the threshold and classification; in genus 1 and 2, Re q = 0 and q(α) = q(1/conj α) at branch points, and 2q(α_2) is a period of B_1 ± B_2.
   Still to do: continuing marked pairs along a deformation of the curve (with the edges, milestones 3 and 4).

## TODO

1. Rewrite the README for cw-families: it still describes cmc-immersion (title, controls, layout, the removed panels).
2. The label of a λ_j overlaps the label of the branch point it sits on.
3. The page title and subtitle ("CMC surfaces") and `about.html` still describe cmc-immersion.

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
