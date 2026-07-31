#pragma once

#include <cstddef>
#include <vector>

#include <GridKit/Model/OPF/Component.hpp>
#include <GridKit/Model/OPF/ComponentConstraints.hpp>
#include <GridKit/Model/OPF/ComponentVariables.hpp>
#include <GridKit/Model/OPF/Shunt/ShuntData.hpp>
#include <GridKit/Model/StateData.hpp>

namespace GridKit
{
  namespace OPF
  {
    enum class ShuntExternalVariables : std::size_t
    {
      VM, ///< Bus voltage magnitude
      VA, ///< Bus voltage angle
      MAXIMUM
    };

    enum class ShuntExternalConstraints : std::size_t
    {
      DIVP, ///< Bus active-power balance
      DIVQ, ///< Bus reactive-power balance
      MAXIMUM
    };

    /// Fixed shunt-admittance contribution to a bus power balance.
    template <class scalar_type, typename index_type>
    class Shunt final : public Component<scalar_type, index_type>
    {
    public:
      using ScalarT      = scalar_type;
      using IdxT         = index_type;
      using ComponentT   = Component<ScalarT, IdxT>;
      using RealT        = typename ComponentT::RealT;
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
      void addJacobianSparsity(std::vector<typename ComponentT::JacobianEntry>& entries) const override;
      int  setJacobianSlots(const std::vector<IdxT>& slots) override;
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
