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

## Findings so far (November 2026)

Each item says how it was checked.
Sources: the notes `background/genus_zero.tex` and `background/genus_one.tex` in the papers folder (`2026-post-qft`), Bohle 2010 (Constrained Willmore tori in the 4-sphere), and Bohle–Leschke–Pedit–Pinkall 2012 (BLPP, Conformal maps from a 2-torus to the 4-sphere).

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
   The known rectangular cases presumably lie in \bbR^3 (unchecked).
10. F10: Parameter count.
   A σ-node arriving over S^1 is a codimension 1 event on an edge (one parameter on S^1 against q ∈ ½Γ in i\bbR^2).
   The branch pair from the opened node is a trivial-multiplier ρ-pair along the whole edge, provided its class matches.
   So the candidates are isolated points on the edges.

## Open questions

1. Class matching at candidates 1 and 2.
   The branch pair from opening 1 has class [1], while the pairs of τ and of the third vector have other classes.
   By F8, a candidate built from 1 and τ closes only on a double cover.
   Are candidates 1 and 2 meant to be on the double cover too, or is the second Sym pair at those candidates of class [1]?
   Answer (Ross): take the double cover whenever necessary.
   Since the three marked pairs have three different classes, this is always the case when the two Sym pairs are different marked pairs.
2. Which conformal class is fixed?
   The torus lattice of the candidate is the index 2 sublattice of F8, whose conformal class generally differs from τ (for a square lattice, one index 2 sublattice is again square and the other two are rectangles with ratio 2).
   Is τ in this sketch the class of Γ, or the class of the final (doubled) torus?
3. Do the three candidates for a rectangular τ lie in \bbR^3, as in the known cases?
