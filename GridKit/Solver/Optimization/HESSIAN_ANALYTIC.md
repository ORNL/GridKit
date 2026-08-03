The previous proposal is mathematically correct, but it is not the ideal design for your requirement.
Its scalar-Lagrangian formulation, component-local assembly, lower-triangular storage, and forward-over-reverse idea are right. Its recommendation—production SecondOrderJet/manual formulas with Enzyme only as validation—directly contradicts “fully Enzyme,” preserves the manual Jacobian, and should be rejected.
The correct architecture is:
One primal equation implementation per component → Enzyme generates every gradient/Jacobian/Hessian value → topology supplies immutable sparse coordinates → System deduplicates them into CSR → Ipopt receives the exact lower-triangular Lagrangian Hessian.

No global dense Hessian, finite-difference derivative, hand-coded derivative formula, or L-BFGS default remains.
What should be retained
The OPF branch already has a useful composition layer:
Global variable and constraint ownership/binding.
Component-local contributions into shared bus-balance rows.
One-time coordinate sorting, deduplication, CSR construction, and component slot binding in System.cpp.
Fixed tap, phase, and open operating inputs. They must remain outside the optimization vector.
What should be deleted is the duplicated calculus:
The hand-written branch \(4\times4\) Jacobian and thermal chain rule in Branch.cpp.
Branch::addJacobian() at lines 353–427.
Generator objective-gradient and balance stamps at Generator.cpp.
Shunt derivatives at Shunt.cpp.
nnz_h_lag = 0, the empty eval_h(), and unconditional limited-memory mode in IpoptSolver.cpp.
PhasorDynamics style to copy
The right exemplar is Tgov1, not the PhasorDynamics branch model:
Equations exist once in an always_inline kernel: [Tgov1Impl.hpp (line 301)](/home/lukel/GridKit/GridKit/Model/PhasorDynamics/Governor/Tgov1/Tgov1Impl.hpp:301).
The Enzyme translation unit is thin: [Tgov1Enzyme.cpp (line 22)](/home/lukel/GridKit/GridKit/Model/PhasorDynamics/Governor/Tgov1/Tgov1Enzyme.cpp:22).
Generic wrappers handle activity annotations and seeding.
System assembly deduplicates local contributions once and subsequently scatter-adds values: [SystemModelImpl.hpp (line 836)](/home/lukel/GridKit/GridKit/Model/PhasorDynamics/SystemModelImpl.hpp:836).
For OPF, the layout should become:
OPF/
  Component.hpp
  OptimizationDerivatives.hpp
  System.hpp
  System.cpp

  Branch/
    Branch.hpp
    BranchImpl.hpp
    BranchEnzyme.cpp

  Generator/
    Generator.hpp
    GeneratorImpl.hpp
    GeneratorEnzyme.cpp

  Shunt/
    Shunt.hpp
    ShuntImpl.hpp
    ShuntEnzyme.cpp
Bus and Load have no derivative implementation because their NLP contributions are constant or empty.
1. Pure primal kernels
The branch should contain only power-flow equations—no derivative arrays:
enum class BranchLocalVariable : std::size_t
{
  VMF,
  VAF,
  VMT,
  VAT,
  MAXIMUM
};

enum class BranchLocalConstraint : std::size_t
{
  PF,
  QF,
  PT,
  QT,
  SF2,
  ST2,
  MAXIMUM
};

template <class ScalarT>
__attribute__((always_inline)) inline void
Branch::evaluateLocalConstraints(const ScalarT* x, ScalarT* g) const
{
  for (std::size_t i = 0; i < localConstraintCount(); ++i)
    g[i] = ScalarT{0};

  if (open_)
    return;

  const ScalarT vmf   = x[std::to_underlying(BranchLocalVariable::VMF)];
  const ScalarT vaf   = x[std::to_underlying(BranchLocalVariable::VAF)];
  const ScalarT vmt   = x[std::to_underlying(BranchLocalVariable::VMT)];
  const ScalarT vat   = x[std::to_underlying(BranchLocalVariable::VAT)];
  const ScalarT delta = vaf - vat;

  const ScalarT cos_delta = std::cos(delta);
  const ScalarT sin_delta = std::sin(delta);

  const ScalarT A = gft_ * cos_delta + bft_ * sin_delta;
  const ScalarT C = gft_ * sin_delta - bft_ * cos_delta;
  const ScalarT D = gtf_ * cos_delta - btf_ * sin_delta;
  const ScalarT E = -gtf_ * sin_delta - btf_ * cos_delta;

  const ScalarT pf = gff_ * vmf * vmf + vmf * vmt * A;
  const ScalarT qf = -bff_ * vmf * vmf + vmf * vmt * C;
  const ScalarT pt = gtt_ * vmt * vmt + vmf * vmt * D;
  const ScalarT qt = -btt_ * vmt * vmt + vmf * vmt * E;

  g[std::to_underlying(BranchLocalConstraint::PF)] = pf;
  g[std::to_underlying(BranchLocalConstraint::QF)] = qf;
  g[std::to_underlying(BranchLocalConstraint::PT)] = pt;
  g[std::to_underlying(BranchLocalConstraint::QT)] = qt;

  if (has_thermal_limit_)
  {
    g[std::to_underlying(BranchLocalConstraint::SF2)] = pf * pf + qf * qf;
    g[std::to_underlying(BranchLocalConstraint::ST2)] = pt * pt + qt * qt;
  }
}
The admittance coefficients are validated and computed once from fixed branch data. Validation and std::isfinite checks stay outside the differentiated kernel.
Generator and shunt follow the same rule:
template <class ScalarT>
__attribute__((always_inline)) inline ScalarT
Generator::evaluateLocalObjective(const ScalarT* x) const
{
  if (!online_)
    return ScalarT{0};

  const ScalarT pg = x[std::to_underlying(GeneratorLocalVariable::PG)];
  return c0_ + c1_ * pg + c2_ * pg * pg;
}

template <class ScalarT>
__attribute__((always_inline)) inline void
Generator::evaluateLocalConstraints(const ScalarT* x, ScalarT* g) const
{
  g[0] = x[std::to_underlying(GeneratorLocalVariable::PG)];
  g[1] = x[std::to_underlying(GeneratorLocalVariable::QG)];
}
These are the only equations. There is no separately maintained derivative algebra.
2. First derivatives are also Enzyme-generated
A generic OPF wrapper applies forward mode to the local constraint kernel:
for (std::size_t column = 0; column < variable_count; ++column)
{
  std::ranges::fill(seed, RealT{0});
  std::ranges::fill(dconstraints, RealT{0});
  seed[column] = RealT{1};

  __enzyme_fwddiff<void>(
      reinterpret_cast<void*>(LocalNlpWrapper<ModelT>::constraints),
      enzyme_const, model,
      enzyme_dup, x, seed.data(),
      enzyme_dupnoneed, constraints.data(), dconstraints.data());

  scatterJacobianColumn(column, dconstraints, jacobian_values);
}
The objective gradient can use reverse mode:
std::ranges::fill(local_gradient, RealT{0});

__enzyme_autodiff<void>(
    reinterpret_cast<void*>(LocalNlpWrapper<ModelT>::objective),
    enzyme_const, model,
    enzyme_dup, x, local_gradient.data());
Thus BranchEvaluation::jacobian, generator’s c1 + 2*c2*pg, shunt’s ±2*coefficient*vm, and all manual thermal derivatives disappear.
3. Exact Enzyme Lagrangian Hessian
Ipopt requests
\[
\nabla^2_{xx}L(x,\lambda,\sigma)
=
\sigma\nabla^2f(x)
+
\sum_i\lambda_i\nabla^2g_i(x).
\]It supplies the structure separately, requires that structure to remain constant, and expects only one symmetric triangle. Ipopt’s interface documentation is explicit about all three requirements.
Each component forms its own scalar contribution:
template <class ModelT>
__attribute__((always_inline)) inline
typename ModelT::RealT localLagrangian(
    ModelT*                             model,
    const typename ModelT::RealT*       x,
    typename ModelT::RealT              objective_factor,
    const typename ModelT::RealT*       lambda)
{
  using RealT = typename ModelT::RealT;

  std::array<RealT, ModelT::LocalConstraintCount> g{};
  model->evaluateLocalConstraints(x, g.data());

  RealT value = objective_factor * model->evaluateLocalObjective(x);
  for (std::size_t row = 0; row < model->localConstraintCount(); ++row)
    value += lambda[row] * g[row];

  return value;
}
Global multipliers are gathered through the existing constraint bindings. For a branch:
lambda_local = {
    lambda[DIVPF],
    lambda[DIVQF],
    lambda[DIVPT],
    lambda[DIVQT],
    has_smax ? lambda[SF2] : 0.0,
    has_smax ? lambda[ST2] : 0.0
};
This is why component-local differentiation remains mathematically exact even when many components contribute to the same bus-balance row: the global Lagrangian is the sum of the component-local Lagrangians.
The Enzyme forward-over-reverse implementation is:
template <class ModelT>
__attribute__((always_inline)) inline void lagrangianGradient(
    ModelT*                       model,
    typename ModelT::RealT*       x,
    typename ModelT::RealT        objective_factor,
    const typename ModelT::RealT* lambda,
    typename ModelT::RealT*       gradient)
{
  __enzyme_autodiff<void>(
      reinterpret_cast<void*>(localLagrangian<ModelT>),
      enzyme_const, model,
      enzyme_dup, x, gradient,
      enzyme_const, objective_factor,
      enzyme_const, lambda);
}

template <class ModelT>
void evaluateLagrangianHessian(
    ModelT*                       model,
    typename ModelT::RealT*       x,
    typename ModelT::RealT        objective_factor,
    const typename ModelT::RealT* lambda,
    typename ModelT::RealT*       global_values)
{
  using RealT = typename ModelT::RealT;

  std::array<RealT, ModelT::LocalVariableCount> seed{};
  std::array<RealT, ModelT::LocalVariableCount> gradient{};
  std::array<RealT, ModelT::LocalVariableCount> hessian_column{};

  for (std::size_t column = 0;
       column < ModelT::LocalVariableCount;
       ++column)
  {
    seed.fill(RealT{0});
    gradient.fill(RealT{0});
    hessian_column.fill(RealT{0});
    seed[column] = RealT{1};

    __enzyme_fwddiff<void>(
        reinterpret_cast<void*>(lagrangianGradient<ModelT>),
        enzyme_const, model,
        enzyme_dup, x, seed.data(),
        enzyme_const, objective_factor,
        enzyme_const, lambda,
        enzyme_dupnoneed, gradient.data(), hessian_column.data());

    model->scatterLowerHessianColumn(
        column, hessian_column.data(), global_values);
  }
}
These are design sketches, not code I compiled this turn. However, this nesting is not speculative: the exact installed Enzyme source snapshot, commit 53a46ca7, contains the same reverse-gradient-inside-forward-Hessian construction in its Sparse/ringspring.cpp integration test. Enzyme also documents forward-over-reverse as its Hessian construction. Enzyme higher-order documentation, C++ calling conventions.
4. Sparse structure comes from topology—not numerical values
This is the most important correction to a literal copy of PhasorDynamics.
[LowerSparseStorage.hpp (line 209)](/home/lukel/GridKit/GridKit/AutomaticDifferentiation/Enzyme/LowerSparseStorage.hpp:209) currently contains:
if (val == 0.0)
  return;
That cannot be used to discover Ipopt’s Hessian pattern. At a flat start, zero flow, zero multiplier, zero objective factor, or special angle, a structurally valid entry can currently evaluate to zero and disappear.
Instead, every component reports immutable support:
Component	Local Hessian support
Bus	0
Load	0
Generator	\((P_g,P_g)\): 1
Shunt	\((V_m,V_m)\): 1
Branch	Lower triangle over \((V_{mf},\theta_f,V_{mt},\theta_t)\): 10


Before deduplication,
\[
\operatorname{nnz}(H)
\le N_G + N_S + 10N_B,
\]so storage remains linear in network size.
For each local pair:
const IdxT global_i = local_variable_indices[i];
const IdxT global_j = local_variable_indices[j];

entries.emplace_back(
    std::max(global_i, global_j),
    std::min(global_i, global_j));
System then sorts, deduplicates, constructs lower-triangular CSR, and binds each raw contribution to its global value slot exactly as it already does for the Jacobian.
Slots remain present even when:
c2 == 0;
objective_factor == 0;
a multiplier is zero;
a branch is open;
a component is offline;
an operating point happens to make the value zero.
Enzyme refreshes the values; it never determines whether a coordinate exists.
The small local arrays above are stack scratch for at most four branch variables. They are not an \(n\times n\) global matrix and are not a dense approximation.
5. Component and model contracts
The component interface should distinguish Jacobian and Hessian contribution counts:
virtual IdxT jacobianContributionCount() const;
virtual void addJacobianSparsity(std::vector<SparseEntry>&) const;
virtual int  setJacobianSlots(std::span<const IdxT>);

virtual IdxT hessianContributionCount() const;
virtual void addHessianSparsity(std::vector<SparseEntry>&) const;
virtual int  setHessianSlots(std::span<const IdxT>);

virtual int addFirstDerivatives(
    const ScalarT* variables,
    ScalarT* objective_gradient,
    RealT*  jacobian_values) const;

virtual int addLagrangianHessian(
    const ScalarT*          variables,
    RealT                   objective_factor,
    std::span<const RealT>  multipliers,
    RealT*                  hessian_values) const;
Bus and Load inherit zero-contribution defaults.
For cleanliness, Hessian operations belong on an optimization-specific evaluator rather than adding more Ipopt-shaped behavior to the general dynamic Model::Evaluator:
virtual IdxT nnzLagrangianHessian() const = 0;

virtual int evaluateLagrangianHessian(
    RealT objective_factor,
    std::span<const RealT> multipliers) = 0;

virtual CsrMatrixT* getCsrLagrangianHessian() const = 0;
System::evaluateLagrangianHessian() clears only the fixed CSR values and asks every component to += its Enzyme-generated contribution.
6. Ipopt adapter
The adapter becomes a straightforward translation:
bool get_nlp_info(..., Ipopt::Index& nnz_h_lag, ...)
{
  return toIpoptIndex(model_->size(), n)
      && toIpoptIndex(model_->sizeResidual(), m)
      && toIpoptIndex(model_->nnz(), nnz_jac_g)
      && toIpoptIndex(model_->nnzLagrangianHessian(), nnz_h_lag);
}
bool eval_h(
    Ipopt::Index         n,
    const Ipopt::Number* x,
    bool                 new_x,
    Ipopt::Number        objective_factor,
    Ipopt::Index         m,
    const Ipopt::Number* lambda,
    bool,
    Ipopt::Index         nele_hess,
    Ipopt::Index*        iRow,
    Ipopt::Index*        jCol,
    Ipopt::Number*       values) override
{
  auto* hessian = model_->getCsrLagrangianHessian();
  if (!validLowerHessian(hessian, n, nele_hess))
    return false;

  if (values == nullptr)
    return copyLowerCsrStructure(*hessian, iRow, jCol);

  if (!syncVariables(x, new_x))
    return false;

  const std::span<const Ipopt::Number> multipliers{
      lambda, static_cast<std::size_t>(m)};

  if (model_->evaluateLagrangianHessian(
          objective_factor, multipliers) != 0)
    return false;

  std::copy_n(hessian->getValues(), nele_hess, values);
  return true;
}
Initially, I would recompute on each value request. If caching is later worthwhile, the key must include \(x\), lambda, and objective_factor—not merely new_x.
Remove the forced:
hessian_approximation = limited-memory
and explicitly select exact, or rely on Ipopt’s exact default. Ipopt’s own documentation recommends exact second derivatives when they are economical and provides derivative_test=second-order for validation. Ipopt options.
7. Build contract
Each nonlinear component should be compiled as an Enzyme-visible target, following current PhasorDynamics CMake style:
gridkit_add_library(
  model_opf_branch
  SOURCES BranchEnzyme.cpp
  HEADERS Branch.hpp BranchData.hpp
  LINK_LIBRARIES
    PUBLIC GridKit::model_opf_core
    PRIVATE ClangEnzymeFlags
  COMPILE_OPTIONS
    PRIVATE
    -fno-math-errno)
Because this design uses fixed structural slots rather than __enzyme_todense to discover sparsity, -enzyme-auto-sparsity=1 is not inherently required for the OPF Hessian path.
The exact-Hessian OPF solver must require Enzyme. It must never silently substitute manual derivatives, finite differences, or limited-memory approximation when Enzyme is disabled.
Verification required before calling it complete
The implementation should proceed in these reviewable slices:
Prove the generic forward-over-reverse wrapper against the pinned Enzyme/LLVM 16 environment using a two-variable scalar function.
Split OPF primal equations into Impl.hpp without changing behavior.
Replace all manual objective gradients and constraint Jacobians with Enzyme.
Add fixed component Hessian support, global deduplication, and slot binding.
Implement eval_h and exact Ipopt mode.
Delete the old derivative formulas only after parity tests pass.
The required tests are:
Enzyme objective gradient and constraint Jacobian versus independent numerical or dependency-tracking oracles.
Hessian versus finite differences of\[
\sigma\nabla f(x)+J(x)^T\lambda.
\]
Mixed nonzero multipliers, not merely objective curvature.
objective_factor == 0.
Pattern invariance across flat/nonflat states and zero/nonzero multipliers.
Duplicate accumulation from multiple branches and shunts.
Lower-triangular, sorted, unique global coordinates.
Open/offline devices.
Ipopt derivative_test=second-order.
Exact-Hessian integration solve.
Large-case assertion that Hessian storage scales with network topology, never \(n^2\).
So the answer is: use component-local Enzyme forward-over-reverse as the production Hessian, retain only explicit structural sparsity, and remove every manual derivative and the hard-coded L-BFGS path. That is the clean design consistent with your requirement.
I remained on lukel/line-params-dev@3f46c82b; origin/lukel/state-init-2-dev remains at 1b468507. The worktree is clean. I made no edits, branch switches, builds, or git mutations.