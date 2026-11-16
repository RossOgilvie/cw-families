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