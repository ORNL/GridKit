#pragma once

#include <cstddef>
#include <vector>

#include <GridKit/Model/OPF/Component.hpp>
#include <GridKit/Model/OPF/ComponentConstraints.hpp>
#include <GridKit/Model/OPF/ComponentVariables.hpp>
#include <GridKit/Model/OPF/Generator/GeneratorData.hpp>
#include <GridKit/Model/StateData.hpp>

namespace GridKit
{
  namespace OPF
  {
    enum class GeneratorInternalVariables : std::size_t
    {
      PG, ///< Active-power generation
      QG, ///< Reactive-power generation
      MAXIMUM
    };

    enum class GeneratorExternalVariables : std::size_t
    {
      VM, ///< Bus voltage magnitude
      VA, ///< Bus voltage angle
      MAXIMUM
    };

    enum class GeneratorExternalConstraints : std::size_t
    {
      DIVP, ///< Bus active-power balance
      DIVQ, ///< Bus reactive-power balance
      MAXIMUM
    };

    /// Dispatchable generation and quadratic cost contribution.
    template <class scalar_type, typename index_type>
    class Generator final : public Component<scalar_type, index_type>
    {
    public:
      using ScalarT      = scalar_type;
      using IdxT         = index_type;
      using ComponentT   = Component<ScalarT, IdxT>;
      using RealT        = typename ComponentT::RealT;
      using DataT        = GeneratorData<RealT, IdxT>;
      using StateT       = Model::DeviceState;
      using VariablesT   = ComponentVariables<ScalarT,
                                              IdxT,
                                              GeneratorInternalVariables,
                                              GeneratorExternalVariables>;
      using ConstraintsT = ComponentConstraints<IdxT,
                                                NoConstraints,
                                                GeneratorExternalConstraints>;

      Generator(const DataT& data, const StateT& state);
      ~Generator() override = default;

      VariablesT&         variables();
      const VariablesT&   variables() const;
      ConstraintsT&       constraints();
      const ConstraintsT& constraints() const;

      IdxT sizeInternalVariables() const override;
      IdxT sizeInternalConstraints() const override;
      IdxT nnz() const override;
      void setVariableOffset(IdxT offset) override;
      void addJacobianSparsity(std::vector<typename ComponentT::JacobianEntry>& entries) const override;
      int  setJacobianSlots(const std::vector<IdxT>& slots) override;
      int  initialize(ScalarT* values) const override;
      int  setVariableBounds(RealT* lower, RealT* upper) const override;
      int  addObjective(const ScalarT* values, RealT& objective) const override;
      int  addObjectiveGradient(const ScalarT* values, ScalarT* gradient) const override;
      int  addConstraints(const ScalarT* values, ScalarT* constraints) const override;
      int  addJacobian(const ScalarT* values, RealT* jacobian_values) const override;
      int  updateSolutionState(const ScalarT*    values,
                               Model::StateData& state) const override;

    private:
      DataT             data_;
      StateT            state_;
      VariablesT        variables_;
      ConstraintsT      constraints_;
      std::vector<IdxT> jacobian_slots_;
    };

  } // namespace OPF
} // namespace GridKit
