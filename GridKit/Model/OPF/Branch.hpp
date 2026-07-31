#pragma once

#include <cstddef>

#include <GridKit/Model/OPF/Component.hpp>
#include <GridKit/Model/OPF/ComponentConstraints.hpp>
#include <GridKit/Model/OPF/ComponentVariables.hpp>
#include <GridKit/Model/OPF/SystemData.hpp>
#include <GridKit/Model/StateData.hpp>

namespace GridKit
{
  namespace OPF
  {
    enum class BranchExternalVariables : std::size_t
    {
      FROM_VOLTAGE_MAGNITUDE,
      FROM_VOLTAGE_ANGLE,
      TO_VOLTAGE_MAGNITUDE,
      TO_VOLTAGE_ANGLE,
      MAXIMUM
    };

    enum class BranchInternalConstraints : std::size_t
    {
      FROM_APPARENT_POWER,
      TO_APPARENT_POWER,
      MAXIMUM
    };

    // These rows are present only when smax is defined. Otherwise
    // sizeInternalConstraints() is zero and their bindings remain invalid.

    enum class BranchExternalConstraints : std::size_t
    {
      FROM_ACTIVE_POWER_BALANCE,
      FROM_REACTIVE_POWER_BALANCE,
      TO_ACTIVE_POWER_BALANCE,
      TO_REACTIVE_POWER_BALANCE,
      MAXIMUM
    };

    template <class scalar_type, typename index_type>
    class Branch : public Component<scalar_type, index_type>
    {
    public:
      using ScalarT      = scalar_type;
      using IdxT         = index_type;
      using RealT        = typename Component<ScalarT, IdxT>::RealT;
      using DataT        = BranchData<RealT, IdxT>;
      using StateT       = Model::DeviceState;
      using VariablesT   = ComponentVariables<ScalarT,
                                              IdxT,
                                              NoVariables,
                                              BranchExternalVariables>;
      using ConstraintsT = ComponentConstraints<IdxT,
                                                BranchInternalConstraints,
                                                BranchExternalConstraints>;

      Branch(const DataT& data, const StateT& state);
      ~Branch() override = default;

      VariablesT&         variables();
      const VariablesT&   variables() const;
      ConstraintsT&       constraints();
      const ConstraintsT& constraints() const;

      IdxT sizeInternalVariables() const override;
      IdxT sizeInternalConstraints() const override;
      IdxT nnz() const override;
      void setConstraintOffset(IdxT offset) override;
      void addJacobianSparsity(
          std::vector<typename Component<ScalarT, IdxT>::JacobianEntry>& entries) const override;
      int setJacobianSlots(const std::vector<IdxT>& slots) override;
      int setConstraintBounds(RealT* lower, RealT* upper) const override;
      int addConstraints(const ScalarT* values, ScalarT* constraints) const override;
      int addJacobian(const ScalarT* values, RealT* jacobian_values) const override;

    private:
      DataT             data_;
      StateT            state_;
      VariablesT        variables_;
      ConstraintsT      constraints_;
      std::vector<IdxT> jacobian_slots_;
    };

  } // namespace OPF
} // namespace GridKit
