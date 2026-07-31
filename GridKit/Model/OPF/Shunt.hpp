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
    enum class ShuntExternalVariables : std::size_t
    {
      VOLTAGE_MAGNITUDE,
      VOLTAGE_ANGLE,
      MAXIMUM
    };

    enum class ShuntExternalConstraints : std::size_t
    {
      ACTIVE_POWER_BALANCE,
      REACTIVE_POWER_BALANCE,
      MAXIMUM
    };

    template <class scalar_type, typename index_type>
    class Shunt : public Component<scalar_type, index_type>
    {
    public:
      using ScalarT      = scalar_type;
      using IdxT         = index_type;
      using RealT        = typename Component<ScalarT, IdxT>::RealT;
      using DataT        = ShuntData<RealT, IdxT>;
      using StateT       = Model::DeviceState;
      using VariablesT   = ComponentVariables<ScalarT,
                                              IdxT,
                                              NoVariables,
                                              ShuntExternalVariables>;
      using ConstraintsT = ComponentConstraints<IdxT,
                                                NoConstraints,
                                                ShuntExternalConstraints>;

      Shunt(const DataT& data, const StateT& state);
      ~Shunt() override = default;

      VariablesT&         variables();
      const VariablesT&   variables() const;
      ConstraintsT&       constraints();
      const ConstraintsT& constraints() const;

      IdxT sizeInternalVariables() const override;
      IdxT sizeInternalConstraints() const override;
      IdxT nnz() const override;
      void addJacobianSparsity(
          std::vector<typename Component<ScalarT, IdxT>::JacobianEntry>& entries) const override;
      int setJacobianSlots(const std::vector<IdxT>& slots) override;
      int addConstraints(const ScalarT* values, ScalarT* constraints) const override;
      int addJacobian(const ScalarT* values, RealT* jacobian_values) const override;
      int updateSolutionState(const ScalarT*    values,
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
