The goal of this numerical simulator is to "walk" through a "tree" in the space of sinh-Gordon spectral data using the τ-preserving flow, both for smooth spectral data, but also inserting and opening nodes on \bbS^1.

Let us describe the "starting point". We have a genus 0 spectral curve. The spectral curve has the form ν^2 = λ ; there are no parameters. We look at the ρ-real differentials with double poles and no residues over λ=0,∞. They are exact, so we write them as dq = d(linear in λ/ν). They are completely determined by the coefficient of λ^0. Let q^1 correspond to coefficient 1, and q^2 correspond to coefficient i. In the construction of the immersion, this plane of differentials is the domain of the immersion.

Additionally we have a param τ in the upper half plane (actually the fundamental domain \calF_1), which we fix (for phase 1 of the project).
We may consider q = (q^1,q^2) as a map to \bbC^2 / Γ for Γ = \bbZ< (i,0), (i \Real τ, i \Imag τ) > .
This is a bit trivial in the genus 0 case, where q^k are well-defined functions, but in general they have purely imaginary periods, so the real parts are well defined, but the imag parts are only quasiperiodic on the spectral curve. 
The image of the spectral curve in \bbC^2 / Γ has self-intersections, which we call resonance pairs. 
That is, a resonance pair is a pair of points P ≠ Q of the spectral curve with q(P) − q(Q) ∈ Γ; in particular Re q(P) = Re q(Q).
At least in genus 0, resonance points should be in correspondence with points of the lattice Γ. That is, we can think of the image in \bbC^2 and the image shifted by an element of Γ, which should have a unique intersection point. The corresponding points on the spectral curve are the two points of the resonance pair associated to this lattice element. We mark (and will follow) the resonance pairs associated to the three shortest lattice vectors, 1, τ, and the third.

The τ, the genus zero spectral curve, and these marked resonance pairs are our starting data. 

It should be the case that the resonance pair associated to 1 is over \bbS^1. We insert a node here an open into sinh-Gordon genus 1 spectral curves. This part I don't remember the details so exactly. By symmetry the resonance pair associated to -1 should also be on \bbS^1. My recollection is that we "open" the -1 resonance pair in some way that destroys it, while the 1 resonance pair becomes the branch points. I can dig up some old notes that might help here, if you don't already understand it.

Now we do τ-preserving flow (which should be a 1 dimensional flow), deforming the spectral curve ν^2 = λ a(λ), while following the marked resonance pairs. 
We continue until another resonance pair reaches \bbS^1.
This is our first "candidate".
Then we open, like from genus 0 into genus 1, and continue until another resonance pair reaches \bbS^1.
This is our second candidate. 

There's a third candidate, but we can implement that later.

In general, the "tree" consists edges of τ-preserving flow within spectral curves on a given genus, and branching by inserting and opening nodes when resonance points reach \bbS^1 (or not). 

## Findings so far (November 2026; last updated 22 November)

Each item says how it was checked. The code and its status are in `PLAN.md` (in the code repository, `~/programming/web/cw-families`).
Sources: the notes `background/genus_zero.tex` and `background/genus_one.tex` in the papers folder (`2026-post-qft`), Bohle 2010 (Constrained Willmore tori in the 4-sphere), Bohle–Leschke–Pedit–Pinkall 2012 (BLPP, Conformal maps from a 2-torus to the 4-sphere), Knopf–Peña Hoepner–Schmidt 2025 (K-PH-S, Solutions of the sinh-Gordon equation of spectral genus two and constrained Willmore tori I), and Schmidt's DFG proposal as rewritten in November 2026 (`antrag_rewrite_11_26.pdf`).

1. F1: Sinh-Gordon curves are constrained Willmore spectral curves.
   For CMC tori (with a suitable Lagrange multiplier) the constrained Willmore spectral curve is a double cover of the spectral plane, branched over 0 and ∞, and equals the harmonic map curve ν^2 = λ a(λ) (Bohle 2010, Theorem 5.1 Case II and Proposition 6.8).
   The 4-fold cover of Case I is the generic constrained Willmore case.
   The marked points ∞^± lie over λ = 0, ∞, and the real structure η is ρ.
   Checked: read in Bohle.
2. F2: Why nodes are inserted over S^1.
   A node that stays in the sinh-Gordon class is a σ-pair {P, σP} (σ the hyperelliptic involution), and it is η-invariant exactly over S^1, because η acts as σ on the fibres over S^1.
   So a resonance pair can be opened, within sinh-Gordon, only at the moment it arrives at S^1.
   In the general Willmore notes the (0,1) node can be opened for every W below π^2 vol (genus_one.tex, line 1327); the sinh-Gordon opening is that endpoint.
   Checked: by hand.
3. F3: Genus 0 in closed form.
   Normalize ρ(ν) = −1/conj(ν) on ν^2 = λ, so that q^1 = 1/ν − ν and q^2 = i(1/ν + ν) are imaginary over S^1.
   For each γ ∈ Γ the resonance equation q(P) − q(Q) = γ has exactly two ordered solutions, (P, Q) and (σQ, σP).
   For |γ| < 4 both pairs lie over S^1, at two distinct values λ = e^{2i(φ∓ψ)} with sin ψ = |γ|/4.
   For |γ| > 4 they are ρ-pairs off S^1.
   At |γ| = 4 the pair collides with its σ-partner into a σ-node over a single λ ∈ S^1.
   This is the threshold |λ| = 4a of genus_zero.tex, lines 1117–1119.
   Checked: by hand and numerically.
4. F4: Starting scale.
   The genus 0 curve has a scale (the vacuum picture), and the starting scale is the one at which the pair for 1 has just arrived at S^1 (|1| = 4a).
   Longer vectors are then still off S^1.
   It is unique for τ in the interior of \calF_1; on the boundary |τ| = 1 there is a tie.
5. F5: Opening a node.
   Replace the double root (λ − λ_0)^2 by a ρ-pair of branch points α, 1/conj(α) near λ_0, with |α| slightly below 1, as in `src/cmc/genus2.js`.
   At the moment of arrival the pairs for 1 and −1 are the same σ-node, so opening uses up both.
6. F6: q at branch points.
   With a branch point as base point, q∘σ = −q.
   So at every branch point B, 2q(B) is a period, and q(B) ∈ ½Γ on a leaf with periods in Γ.
   The class of q(B) in ½Γ/Γ ≅ Γ/2Γ is constant along an edge.
   A σ-node of γ has q ≡ γ/2, of class [γ], and branch points born from opening it keep that class.
   The classes are defined up to one common shift; only differences are intrinsic.
   The three marked vectors 1, τ and the third are the shortest representatives of the three nonzero classes of Γ/2Γ.
   Checked: by hand.
7. F7: Generalized Sym points.
   For W < 8π, the holomorphic sections of V/L with trivial monodromy form a space of quaternionic dimension 2 (Plücker), and f is recovered from them by the Kodaira correspondence (BLPP).
   The parallel sections of the associated family at every λ ∈ \bbC^* are holomorphic sections of this one bundle (Bohle 2010, Lemmas 2.4 and 6.7).
   In our conventions (F_x = F X, F(z + w) = M_w F(z)) this gives
   f = Q(ψ_1) Q(ψ_2)^{-1},  ψ_j(z) = u_j^T F(z, λ_j),  Q(a, b) = [[a, b], [−conj b, conj a]],
   where u_j is a left eigenvector of ζ(0, λ_j) (equivalently of the monodromy), at one point of each of two ρ-pairs with trivial multiplier.
   Each ρ-pair needs only one eigenvector, so a pair of branch points works as well as a σ-node.
   Checked numerically:
   1. f is conformal (to 1e-16) for any λ_1, λ_2 ∈ \bbC^* and any rows, in genus 0, 1 and 2; controls with wrong generators fail.
   2. Genus 0 homogeneous tori: W = π^2 (t + 1/t) to 10 digits, for aspect ratio t from 1 to 3.9; t = 1 is the Clifford torus, with Sym angles ±π/2.
   3. With λ_1, λ_2 ∈ S^1 it is the S^3 Sym–Bobenko formula F(λ_1) F(λ_2)^{-1}, up to a Möbius transformation.
8. F8: Closing and the doubling trick.
   A Sym point of class [γ] has monodromy with sign character w ↦ (−1)^{⟨w, γ⟩}.
   Two Sym points of the same class close on the original torus.
   Two of different classes close on the index 2 sublattice of the torus lattice on which the character of the third class is trivial.
   This is "choose 2 of 3" with doubling.
   Note the direction: the torus lattice gets an index 2 sublattice, while the lattice in q-space gets refined.
   Checked: by hand; the genus 0 test closes with sign-matched multipliers −1, −1.
9. F9: Generically the target is S^4.
   With a branch point as one generalized Sym point, f lies in no round 3-sphere (local patch, fit residual 7e-2, against 1e-11 for two σ-nodes).
   This is consistent with Bohle 2010, Corollary 6.12: a Case II torus in a 3-sphere is CMC, which needs σ-nodes.
   So we compute in \bbR^4 and project to \bbR^3.
   The known rectangular cases presumably lie in \bbR^3 (unchecked); F14 gives one candidate that lies on a 3-sphere.
10. F10: Parameter count.
   A σ-node arriving over S^1 is a codimension 1 event on an edge (one parameter on S^1 against q ∈ ½Γ in i\bbR^2).
   The branch pair from the opened node is a trivial-multiplier ρ-pair along the whole edge, provided its class matches.
   So the candidates are isolated points on the edges.
11. F11: The candidates lie in S^4, not in a 3-sphere.
   For candidate 1 at internal τ = 0.2 + 1.3i and 0.45 + 0.95i, and for all three choices of Sym points on the genus 2 curve at 0.45 + 0.95i, the closed torus lies on no round 3-sphere and in no hyperplane: the smallest singular value of [|f|^2, f, 1] over the torus is 0.15–0.33 of the largest (for the Clifford torus it is 1e-11).
   This agrees with Bohle 2010: our tori have a harmonic left normal (im Q° constant, Lemma 6.3), and a Case II torus in a 3-sphere is CMC there (Corollary 6.12), whose constant Darboux transforms are σ-nodes over S^1 (Lemma 4.7 c); a branch point as a Sym point rules this out.
   So the two-coloured look of the drawings (one half seen from outside, one from inside) comes from the projection: a generic projection of a surface in R^4 to R^3 is an immersion with self-intersection curves ending in isolated cross-caps, and across such a curve one sees the other side.
   Checked: numerically (the fit) and by the argument above.
   Read in S^4 after a Möbius normalization (conformal barycenter), the candidates at 0.45 + 0.95i, 0.4 + i and 0.3 + 0.8i are still at distance 0.22–0.52 from the nearest 3-sphere (Clifford torus: 4e-8); in that frame the drawings lose the two-coloured look.
12. F12: The open end is a branched double cover of the round sphere.
   When the genus 2 edge reaches the open end without an arrival, candidate 2 (the two pairs of branch points as Sym points, on the doubled lattice) at scale t = sqrt|α_1 α_2| tends to the quotient map T^2 → T^2/±1 ≅ S^2 of the external torus onto a round 2-sphere, so W → 8π.
   Read in S^4 after a Möbius normalization, the distance from a round 2-sphere is O(t^{1/2}) and the area tends to 2 × 4π, with no concentration of area (no bubbling); the branch points are the four half-periods of the doubled lattice, and f(−z) − f(z) → 0.
   Checked: numerically at internal τ = 0.2 + 1.3i, 1.6i and 2i, for t = 0.03 down to 0.001; not proved.
13. F13: Before ripeness the immersion is a plane, equivariant under homotheties.
   With u_j the eigenvector of ζ(0, λ_j), which is also an eigenvector of the monodromy, ψ_j(z + w) = μ_j(w) ψ_j(z), and Q(μa, μb) = diag(μ, conj μ) Q(a, b) gives f(z + w) = μ_1(w) f(z) μ_2(w)^{-1} for w in the internal lattice.
   At the points of a marked pair off S^1 the multipliers are real: such a pair is (P, ρσP), and the resonance condition μ(P) = μ(ρσP) = conj μ(P) forces μ(P) to be real. Branch points have μ = ±1.
   So the lattice acts by homotheties of R^4 about 0, f(z + w) = c(w) f(z) with c(w) = μ_1(w)/μ_2(w) real, and the Sym points are ripe (the torus closes on the lattice) exactly when c ≡ 1. This is the analogue of Sym–Bobenko without the closing conditions, where the lattice acts by rigid motions.
   P and Q = ρσP give different planes (the two eigenlines over λ_P), which agree when the pair arrives as a σ-node.
   Checked: numerically (f(z + e) = c f(z) with the same real c at two points z, to 1e-8, on the genus 0, 1 and 2 edges at 0.45 + 0.95i), and the argument above.
14. F14: Double arrivals on the walls of F_6.
   Where the two classes after the opened one are equally long (Re τ = 1/2 in the chamber 123, and its images), both of their pairs arrive at S^1 at the same point of the genus 1 edge, and the genus 2 edge has length 0.
   Then all three candidates lie on the genus 1 curve, which has three ripe ρ-pairs (α and the two σ-nodes); they are the three choices of two of them. Candidate 2 is candidate 1, candidate 3 is α + the other node, and the third choice is the two nodes.
   They have the same W, and the candidates just off the wall tend to them. At 0.5 + 0.9i the external τ of candidate 3 is the mirror image of that of candidate 1, and the third choice is rectangular (external τ = 2 Im τ · i).
   Checked: numerically at 0.5 + 0.9i, −0.5 + 1.2i, 0.5 + 1.4i, and by continuity from 0.499 + 0.9i.
   The third choice has two σ-nodes as Sym points, so it lies on a round 3-sphere (F11): it is a CMC torus in S^3, with a rectangular external τ, and so a torus in R^3 after stereographic projection (the app's automatic projection recognizes the 3-sphere). This bears on open questions 4 and 5.
15. F15: How the tree changes with the internal τ.
   1. W increases monotonically along the whole path from the vacuum to the last candidate, at every internal τ tried.
   2. Candidates 2 and 3 merge continuously into the open end. Along Im τ = 0.95 the arrival on the genus 2 edge goes to t = 0 as Re τ decreases to about 0.33 (t = 0.27, 0.14, 0.05, 0.014 at Re τ = 0.45, 0.40, 0.36, 0.34), and W of candidates 2 and 3 tends to 8π.
   3. Across a wall of F_6 the order of the classes changes but the candidates continue: across Re τ = 1/2 candidate 1 is unchanged and the roles of classes 2 and 3 are exchanged (mirror symmetry); on the wall itself there is a double arrival (F14).
   4. No discontinuity of W has been seen in the sweeps so far (the app reports one when it happens between nearby τ).
   5. The map of F_6 (`data/f6map.json`, a grid of step 0.05 with Im τ ≤ 2.4): candidates 2 and 3 exist only in narrow regions about the walls Re τ = ±1/2 (where the double arrivals are), widest near Im τ ≈ 0.9 where the chambers 312, 123 and 213 meet, and narrowing towards i∞. Each region is bounded by two curves going up to i∞. Everywhere else the genus 2 edge reaches the open end first. Of 1679 grid points: 167 with candidates 2 and 3, 54 double arrivals, 1424 open ends, and 37 at the cusps and corners without candidate 1.
   6. Over the map, W of candidate 1 runs from 2.005π² (near the Clifford torus, τ = i) to 8π, and W of candidates 2 and 3, where they exist, from 2.3505π² (on the walls) to 8π (at the boundary).
   Checked: numerically along a few lines in F_6, and on the grid of the map (the boundary bisected along the grid edges to about 3e-6).
16. F16: The branch that does not open the node, and candidates in pairs of one conformal class.
   1. Continuing the genus 1 edge past candidate 1 without opening the node, the arrived pair stays over S^1 as two distinct points. At such a point the multipliers are unimodular, so with it as a Sym point the lattice acts by rotations of R^4 about 0 (compare F13).
   2. The branch ends at the next σ-node event: the other pair arriving (at 0.45 + 0.95i at α = −0.28555, W = 2.4857π²; at 0.3 + 0.8i at α = −0.13621, W = 2.5755π²), or the arrived pair leaving S^1 again (at 0.2 + 1.3i at α = −0.05684, W = 2.6700π²). There α and the node are ripe, and give a torus: candidate 1′ (a provisional name). There is no branch at a double arrival.
   3. The external τ depends only on the internal lattice and on the classes of the two Sym points (F8). So candidates 1 and 2 (classes [1] and a) have the same external τ, and so do candidate 3 and candidate 1′ when the other pair arrives (classes [1] and b). At 0.45 + 0.95i: candidates 1 and 2 have external τ −0.0871 + 1.5768i, with W = 2.3131π² and 2.4725π²; candidates 3 and 1′ have external τ 0.1855 + 1.7195i, with W = 2.4725π² and 2.4857π².
   4. W increases along the branch too. Unlike the other candidates, candidate 1′ mostly has W > 8π: at 1503 of the 1578 points of the map where it exists, up to 2.93π². So the argument of F7 for W < 8π does not cover it (the formula is still a conformal immersion, F7.1).
   Checked: numerically (tests at 0.45 + 0.95i, 0.2 + 1.3i and 0.5 + 0.9i), and point 3 by F8.
17. F17: The candidate family is indexed by the conformal class of the torus, with one member for each half-period (22 November; from Schmidt's rewritten proposal and K-PH-S, Sections 4–5).
   1. In K-PH-S the genus 1 curves with one σ-node over S^1 (their M^2_2) map onto the lattice of the potential (our internal τ) by a threefold covering (Remark 4.2). The Sym points are the branch points and the node (Lemma 5.1), and the torus closes on the index 2 sublattice on which their characters agree (Lemma 5.3, our F8). These are exactly the curves at the events of our genus 1 edges. At internal τ = 0.2 + 1.3i the three are candidate 1, candidate 1′, and a curve on the edge that opens class 2 first, with a node of class 1, which the tree does not reach.
   2. The three curves over one internal τ close on different sublattices. Conversely, a conformal class τ̂ has three index 2 superlattices Γ̂ + h, one for each half-period h of C/Γ̂, and over each of them one of the nodal curves found closes on Γ̂. So the genus 1 family is also a threefold covering of the conformal classes, as the proposal says, and its members are labelled by h.
   3. The least W of the three reproduces Fig. 1 of the proposal: 1.601 × 4π at τ̂ = 0.02 + 1.02i (π/2 at the Clifford torus); π^2 (y + 1/y) at the rectangular class iy, to 5e-9 for y = 1.5 (the homogeneous torus, whose internal τ lies on |τ| = 1); and 1.989 × 4π at 0.45 + 1.42i, which is on the curve W = 8π of Fig. 1.
   4. The label of a pair does not determine the character of its node. At internal 0.2 + 1.3i the pair of class 2 arrives (candidate 1, τ̂ = −0.3133 + 1.1159i) and later leaves S^1 (candidate 1′, τ̂ = −0.2312 + 1.5029i). The branch point is the same, but the closing lattices differ. The curve that opens class 2 with a node of class 1 closes on 0.4 + 2.6i. F16.3 holds for the case it states (the other pair arriving); the reason for the difference is unchecked (perhaps q at a node is only determined modulo ½Γ).
   5. The proposal's family (O-S-2) differs from that of K-PH-S where curves of genus 2 have less W. The page uses this rule: over each superlattice, among the genus 1 nodal curves and the genus 2 candidates of its tree that close on Γ̂, the member is the one with the least W if that is at most 8π, and otherwise the branched double cover of the round sphere.
      1. At τ̂ = 0.1855 + 1.7195i, over h = (1 + τ̂)/2 (internal 0.45 + 0.95i), candidate 3 (genus 2, 2.4725π^2) beats candidate 1′ (genus 1, 2.4857π^2), as in F16.3. Over h = τ̂/2 the only curve found has 2.7726π^2 > 8π, so the member is the double cover.
      2. At τ̂ = 0.45 + 1.2i all three members have genus 1 and W < 8π. Over h = τ̂/2, genus 1 (2.5222π^2) beats candidate 2 (2.5464π^2).
      3. Near the Clifford torus only one member is below 8π; the curves over the other two superlattices (internal τ near 2i) have W just above 8π.
   6. Limits of the search:
      1. Over a superlattice with Im τ̃ above about 3 no nodal curve is found. The edges end at r = 1e-4, and r = 1e-8 does not help. Example: τ̂ = 0.1855 + 1.7195i, h = ½, internal 0.371 + 3.439i.
      2. Genus 2 curves are only searched on the tree that opens the shortest class first.
   7. The map of F_1 (`data/f1map.json`, grid step 0.025, Im τ̂ ≤ 2.3; the labels h are those of the basis (1, τ̂)):
      1. The member h = (1 + τ̂)/2 has the least W almost everywhere. Its superlattice is the one that is square at τ̂ = i.
      2. Two exceptions: small patches near Re τ̂ = ±½ at Im τ̂ ≈ 1.55–1.65, where h = ½ is least, and the top corners, where all three members are double covers. These corners are bounded by the curve where the member (1 + τ̂)/2 reaches 8π, from about (±½, 1.45) to (±0.3, 2.3). This is the "relatively small set of conformal classes" of the proposal where the minimal energy is 8π.
      3. The member (1 + τ̂)/2 changes from genus 1 to genus 2 (genus 2 with less W) above a curve from about (±½, 1.3) to (0, 1.83).
      4. The least W of genus 1 reaches 8π on a curve from (0, 2.12) down to about (±0.4, 1.6), meeting Re τ̂ = ±½ at about 1.65. This has the shape of the white curve of the proposal's Fig. 1, but there the curve meets Re τ̂ = ±½ lower, at about 1.4 (read off the figure; unchecked why).
      5. The member h = ½ is not found above Im τ̂ ≈ 1.5 (F17.6).
   Checked: numerically, at the classes named above (tests in `test/run.mjs` for 1.5i and 0.45 + 1.2i), and on the grid of the map. The selection rule in point 5 is a reading of the proposal, not a theorem.

## Open questions

Settled (answers from Ross):

1. Class matching: take the double cover whenever necessary. Since the three marked pairs have three different classes, this is always the case when the two Sym pairs are different marked pairs (F8).
2. Which conformal class is fixed: the internal τ (of Γ) is fixed; the external τ (of the torus, on the index 2 sublattice of F8) is derived and generally differs. There is no good way yet to fix the external τ ab initio; to be revisited once the candidates are understood better.
   Update (F17): the proposal indexes the family by the external τ̂, and the internal τ of each member is one of the three superlattices of ⟨1, τ̂⟩. The page now starts from τ̂.
3. Candidates 2 and 3: they have the same (genus 2) spectral curve; candidate 2 takes the two pairs of branch points as Sym points, candidate 3 the pair of branch points from genus 1 and the resonance pair arriving on S^1. All candidates on one curve differ only in which resonances are the Sym points (also at a double arrival, F14). Where the edge reaches the open end first, the limit is a branched double cover of the sphere (confirmed numerically, F12).

Open:

4. The R^3 hypothesis.
   The candidates with a branch point as a Sym point lie in S^4, not in a 3-sphere (F11), so they are not tori in R^3.
   How does this fit the hypothesis that one of the three candidates minimizes W among tori in R^3 of the conformal class?
   At a double arrival the third choice (two σ-nodes) is a CMC torus in S^3 (F14).
5. Do the three candidates for a rectangular τ lie in \bbR^3, as in the known cases?
6. Where in F_6 do candidates 2 and 3 exist? Is the boundary (where the arrival reaches t = 0, F15) the curve of the genus 2 paper where the closure of \calS^2 meets the open ends? (That curve is stated for the paper's τ in \calF_3; how it corresponds to the internal τ here is unchecked.) The map (F15.5) shows two curves going up to i∞ about each wall Re τ = ±1/2, which is at least reminiscent of a curve "tending to i∞ at both ends".
7. The "(or not)" branches: what are the candidates on the genus 1 edge continued past candidate 1 without opening, and how should they enter the tree?
   Implemented as candidate 1′ (F16), drawn as a dashed branch in the tree view. Its meaning, and whether it should count among the candidates (it shares the conformal class of candidate 3), is open.
8. Can the regions of F_6 choose the candidate, instead of explicit buttons (Find candidate 1, Find candidates 2, 3, the view picker)? That is, does the position of the internal τ in F_6 (the chamber of the class order, and which side of the boundary of F15 it lies on) determine which candidate to show?
   For now the page can show all three side by side (Three candidates side by side), and the map of F_6 on the τ pad shows where candidates 2 and 3 exist.
   Partly answered by F17: the member over each superlattice is chosen by least W. The page now shows the three members of the family of τ̂, and the tree of one of them.
9. Is "least W over each superlattice, else the double cover" (F17.5) the family of Conjecture 1? In particular, is it continuous where the choice switches between genus 1, genus 2 and the double cover? And where is the third nodal curve over superlattices with large Im τ̃ (F17.6)?
