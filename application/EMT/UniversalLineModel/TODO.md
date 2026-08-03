Highest suspect #1 — the sweep failures: per-mode eigenvector tracking in Gamma is structurally singular at eigenvalue crossings
Gamma.cpp advances each mode separately through the differentiated eigenproblem: for mode m it enforces (Γ − λ_m)·ṫv_m + (Γ̇ − λ̇_m)·tv_m = 0 (and the adjoint twin for ti_m), closed by the pairing ti_m^H tv_m = 1. At a simple eigenvalue this linear system has exactly the right rank: (Γ − λ_m I) is rank K−1, and the one lost rank is absorbed by the extra unknown λ̇_m. At a crossing (λ_m = λ_n), (Γ − λ_m I) drops to rank K−2, but the per-mode block still has only the single λ̇_m to absorb it — the IDA Newton matrix becomes exactly singular, and you get IDA_CONV_FAIL with h → hmin. That is precisely what we observed, with a clean experimental matrix:

Every mirror-symmetric tower (138kv-delta, 345kv, 500kv, 765kv) fails at ω ≈ 1.8 rad/s, right where the low-frequency mode clusters sit. Mirror symmetry guarantees exact crossings.
Both asymmetric geometries (69kv wood-pole, the deliberately "unsymmetric-flat" Overhead line) sweep fine, and staggering the 345kv tower rescued it.
The quad-bundle 765kv fails even staggered, because of the second aggravator: there is no bundle merging or shield-wire Kron reduction anywhere in the parameter chain (OverheadImpl.hpp:241-273 goes conductor→tower→…→gamma with no reduction element; the 345kv sweep's CSV confirms it: 8×8 Yc, 8 modes). Every subconductor and every shield is a modal DOF, so a quad bundle contributes near-identical circulating modes that are degenerate to ~machine precision — untrackable regardless of tower stagger.
Two follow-on consequences come from the same formulation even when it doesn't fail: near-crossings (veering, ~11 Hz on the Overhead line) make the tracking ill-conditioned — that's the source of the Gin/Gout fit plateau and the 50%-pointwise composition error spike at ~477 Hz — and the pairing constraint doesn't pin the per-mode phase (tv→αtv, ti→ti/ᾱ leaves every residual invariant), so the eigenvector phase drifts along the sweep at the discretization's whim — the 2.4 rad drift I had to gauge away in the ULM app.

Proposed fix (layered):

Parameter stage (small, standard, high value): merge bundles into equivalent conductors and Kron-reduce grounded shield wires before the modal stage. K drops to the phase count, and the artificial bundle/shield mode clusters — the worst degeneracies — disappear. This alone likely unblocks staggered/untransposed bundled lines end to end.
Gamma stage (the real fix): track invariant subspaces per delay group instead of individual eigenpairs — states become a basis Q_g (K×c) and a full c×c block Λ_g per group, with residuals Γ·Q_g − Q_g·Λ_g = 0 plus a biorthonormal pairing W_g^H Q_g = I. Crossings inside a group are harmless because only the subspace, which stays smooth, is tracked; the block eigenvalues can be extracted pointwise for τ/H outputs. Downstream, the ULM app then fits delay-grouped dyads P_g = Q_g·(…)·W_g^H — which are gauge-invariant by construction. That one architectural change simultaneously fixes: symmetric-tower sweeps, the veering error spike, the 4e-2 factor-fit plateau, and removes the need for any phase-gauge post-processing. It is also exactly Morched's original ULM formulation.
Highest suspect #2 — the "suspect results"/asymmetry: Yc and Zc compute the wrong matrix
Yc.cpp initializes Yc = sqrt(Z⁻¹Y) and its residual enforces Z·Yc² = Y. The characteristic admittance is Yc = Z⁻¹·sqrt(ZY), which satisfies the sandwich equation Yc·Z·Yc = Y. The two agree only when Z and Y commute — i.e., ideally transposed lines. For untransposed lines they are different matrices. Verified against the monitored sweep data:

ω [rad/s]	monitored vs sqrt(Z⁻¹Y)	monitored vs Z⁻¹sqrt(ZY)	asym of code formula	asym of true formula
8.9	5.7e-11	3.8e-3	7.5e-3	6.9e-16
2742.7	3.1e-12	4.7e-3	9.1e-3	3.2e-16
1.58e6	1.6e-15	8.3e-4	1.6e-3	4.8e-16
The monitor reproduces the code's formula to 11+ digits; the true formula is symmetric to machine precision at every frequency; the difference is a 0.1–0.6% bias in every Yc entry plus the ~1% asymmetry we measured. Zc.cpp has the mirror-image error (sqrt(Y⁻¹Z)).

Proposed fix: change the residual to the sandwich form — equation: −Y + Yc·Z·Yc = 0 (complex product (Gc+jBc)(R+jXL)(Gc+jBc)), with initialization Z⁻¹·sqrtm(Z·Y); correspondingly Zc·Y·Zc = Z with Y⁻¹·sqrtm(Y·Z). Same state layout, same sparsity class, essentially the same cost — and symmetry becomes exact by construction. Once this lands, the Yc symmetrization I added to the ULM app and fit.py becomes redundant (harmless belt-and-braces, but removable), and the Yc fits should tighten further since the fitter will no longer average down a 1% asymmetric perturbation.

One connecting observation: issue #2 is not caused by issue #1 — Yc/Zc never touch the eigenvectors — but both were exposed by the same probe, and together they account for every anomaly seen tonight: the four sweep failures (crossings), the 765kv-staggered failure (bundle clusters), the ~1% Yc/Zc asymmetry (formula), the Gin/Gout plateau and 477 Hz composition spike (veering ill-conditioning), and the eigenvector phase drift (unpinned gauge).

Where things stand meanwhile: the Overhead reorg is done and regenerated (scripts in plotting/, everything in output/, response/ gone, solver specs committed at the example root, gamma α/β plot restored), model JSONs now carry nested rows×cols matrices, and the gallery is built with the staggered 345kv working end to end — I paused it mid-run when you interrupted. Say the word and I'll finish the gallery run and stats table as-is (documenting the symmetric-line failures in the table), or hold everything until you decide on the two fixes above.



I found three major correctness problems and several system-level risks. The first two affect the generated line model even when every existing test passes.

## Major findings

1. **Critical — `Yc` and `Zc` use scalar identities that are invalid for coupled matrices.**

   [Yc.cpp](/home/lukel/GridKit/GridKit/Model/EMT/Parameters/Response/Yc/Yc.cpp:62) initializes
   \[
   Y_c=\sqrt{Z^{-1}Y}
   \]
   and its residual enforces \(Z Y_c^2=Y\). [Zc.cpp](/home/lukel/GridKit/GridKit/Model/EMT/Parameters/Response/Zc/Zc.cpp:62) does the dual operation.

   For a coupled multiconductor line with \(\Gamma^2=ZY\), the characteristic admittance is
   \[
   Y_c=Z^{-1}\Gamma,
   \]
   which satisfies
   \[
   Y_c Z Y_c=Y,
   \]
   not \(Z Y_c^2=Y\). The current equation is valid only when the matrices commute, including the scalar case.

   This explains why [UniversalLineModel.cpp](/home/lukel/GridKit/application/EMT/UniversalLineModel/UniversalLineModel.cpp:539) sees percent-level nonsymmetry and manually averages \(Y_c\). That comment attributes the asymmetry to modal reconstruction, but `Yc` is computed directly from \(Z^{-1}Y\). Symmetrizing afterward does not restore the correct Riccati equation.

   A read-only check against the tracked overhead response produced a physical Riccati residual between roughly **0.6% and 1.4%**. The tests miss this because the multiconductor test only checks the implementation’s own residual, while the independent physical reference is scalar-only in [runOverheadTests.cpp](/home/lukel/GridKit/tests/UnitTests/EMT/Parameters/runOverheadTests.cpp:570).

2. **Critical — phase, circuit, bundle, and grounded-conductor semantics disappear before model export.**

   The parser preserves `phase` and `circuit` in [OverheadDataJSONParser.hpp](/home/lukel/GridKit/GridKit/Model/EMT/Parameters/OverheadDataJSONParser.hpp:535), and the schema says conductors sharing phase and circuit have the same potential while `g` conductors are grounded in [line_schema.py](/home/lukel/GridKit/docs/_ext/gridkit/line_schema.py:272).

   However, production code never uses those fields after parsing. [UniversalLineModel.cpp](/home/lukel/GridKit/application/EMT/UniversalLineModel/UniversalLineModel.cpp:260) infers `K` from the full-conductor `Yc` matrix and exports that physical-conductor count directly.

   Consequently:

   - The twin-bundle 345 kV line is exported as `K=8`, not three electrical phase terminals.
   - The twin-bundle double-circuit line is exported as `K=14`, not six phase/circuit terminals.
   - The quad-bundle 765 kV line is exported as `K=14`, not three phase terminals.

   See [345kv-horizontal.line.json](/home/lukel/GridKit/examples/EMT/Lines/345kv-horizontal.line.json:12), [500kv-double-circuit.line.json](/home/lukel/GridKit/examples/EMT/Lines/500kv-double-circuit.line.json:12), and [765kv-horizontal.line.json](/home/lukel/GridKit/examples/EMT/Lines/765kv-horizontal.line.json:12).

   Full-conductor matrices are reasonable intermediates, but the exported system model must either perform bundle/ground reduction or carry explicit constraints and mapping. It currently does neither.

3. **High — `inner_radius` is accepted and validated but ignored by the conductor model.**

   The schema explicitly models annular conduction regions, and the parser loads the inner radius. But [SkinEffect.cpp](/home/lukel/GridKit/GridKit/Model/EMT/Parameters/Effects/SkinEffect/SkinEffect.cpp:27) uses only the outer radius:
   \[
   R_k\propto\frac{1}{\sigma r_\text{outer}^2}.
   \]

   No production call to `innerRadius()` exists. The supplied Drake and Cardinal catalog entries have nonzero inner radii in [north-american.catalog.json](/home/lukel/GridKit/examples/EMT/Lines/north-american.catalog.json:8). Even at DC, their annular resistance should be approximately **15.4% and 10.3% higher**, respectively, than the current solid-area result. Their frequency-dependent response is also different.

   All skin-effect tests force inner radius to zero in [runSkinEffectTests.cpp](/home/lukel/GridKit/tests/UnitTests/EMT/Parameters/runSkinEffectTests.cpp:24), so this path is untested.

4. **High — modal tracking assumes simple eigenvalues, but supported line geometries naturally produce degenerate modes.**

   [Gamma.cpp](/home/lukel/GridKit/GridKit/Model/EMT/Parameters/Response/Gamma/Gamma.cpp:75) independently sorts eigenvalues and directly inverts the eigenvector matrix without condition checks. The derivative tracking equations then track individual vectors. The model’s own [README](/home/lukel/GridKit/GridKit/Model/EMT/Parameters/Response/Gamma/README.md:11) acknowledges that every tracked branch must be simple.

   Symmetric or transposed three-phase lines naturally contain repeated or near-repeated modal subspaces. Individual eigenvectors are then non-unique, making the tracking Jacobian singular or ill-conditioned and allowing branch rotation/swapping. The accepted input contract does not reject this case, detect it, or provide subspace tracking.

## Additional system-level problems

5. **Absolute tolerances are too loose for small parameter states.**

   [Overhead.hpp](/home/lukel/GridKit/GridKit/Model/EMT/Parameters/Overhead.hpp:169) sets
   \[
   \mathrm{atol}_i=\mathrm{rtol}(1+|y_i|).
   \]
   At the default \(10^{-7}\), every state below one receives an absolute tolerance near \(10^{-7}\), including propagation quantities around \(10^{-8}\). This contradicts the comment claiming scaling to each variable’s initialization magnitude. IDA really uses this vector through [Ida.cpp](/home/lukel/GridKit/GridKit/Solver/Dynamic/Ida.cpp:1298).

6. **Nonpassive fitted models are still emitted as successful artifacts.**

   Passivity checking is explicitly grid-limited in [Passivity.cpp](/home/lukel/GridKit/GridKit/Solver/Optimization/Rational/Passivity.cpp:71). More importantly, [UniversalLineModel.cpp](/home/lukel/GridKit/application/EMT/UniversalLineModel/UniversalLineModel.cpp:592) only prints a warning when the fit is nonpassive, then writes the model and returns success based solely on fit error. A nonpassive EMT line can inject energy or destabilize the simulation; it should not silently become an apparently successful production artifact.

7. **Geometry validation permits physically invalid configurations, followed by release-unsafe inversion.**

   `Tower` verifies distinct conductor centers but not conductor overlap or surface clearance from ground. Yet inductance and potential use logarithms such as \(\log(2h/r)\) in [GeometricInductance.hpp](/home/lukel/GridKit/GridKit/Model/EMT/Parameters/Effects/GeometricInductance/GeometricInductance.hpp:119). The capacitance inversion in [ShuntPotential.cpp](/home/lukel/GridKit/GridKit/Model/EMT/Parameters/Effects/ShuntPotential/ShuntPotential.cpp:116) protects singular pivots only with `assert`; release builds can divide by zero and propagate NaNs.

8. **The aggregate DAE/Jacobian path will scale poorly.**

   The current layout contains approximately \(26K^2+18K+3\) states—about **5,351 states for a 14-conductor line**. Each Jacobian evaluation allocates new jet vectors in [Element.hpp](/home/lukel/GridKit/GridKit/Model/EMT/Element.hpp:260), while every sparse derivative is stored in a dynamic vector and merged through linear searches in [SparseJet.hpp](/home/lukel/GridKit/GridKit/LinearAlgebra/SparseJet.hpp:59). This is a substantial avoidable cost inside an IDA continuation sweep, especially around the dense Gamma tracking equations.

My recommended fix order is: correct the matrix `Yc/Zc` equations, define the physical-conductor-to-electrical-terminal reduction, implement or remove annular conductor support, and then address modal degeneracy. The other issues are important, but those first four determine whether the produced model represents the requested line.

No files were changed, and I did not build or run tests because of the read-only boundary. The pre-existing untracked LineGallery files remain untouched.
