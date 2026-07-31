#pragma once

#include <cstddef>

#include <GridKit/Model/OPF/Bus/BusData.hpp>
#include <GridKit/Model/OPF/Component.hpp>
#include <GridKit/Model/OPF/ComponentConstraints.hpp>
#include <GridKit/Model/OPF/ComponentVariables.hpp>
#include <GridKit/Model/StateData.hpp>

namespace GridKit
{
  namespace OPF
  {
    enum class BusInternalVariables : std::size_t
    {
      VM, ///< Voltage magnitude
      VA, ///< Voltage angle
      MAXIMUM
    };

    enum class BusInternalConstraints : std::size_t
    {
      DIVP, ///< Active-power balance
      DIVQ, ///< Reactive-power balance
      MAXIMUM
    };

    /// OPF bus owning voltage variables and power-balance constraints.
    template <class scalar_type, typename index_type>
    class Bus final : public Component<scalar_type, index_type>
    {
    public:
      using ScalarT      = scalar_type;
      using IdxT         = index_type;
      using ComponentT   = Component<ScalarT, IdxT>;
      using RealT        = typename ComponentT::RealT;
      using DataT        = BusData<RealT, IdxT>;
      using StateT       = Model::BusState;
      using VariablesT   = ComponentVariables<ScalarT,
                                              IdxT,
                                              BusInternalVariables,
                                              NoVariables>;
      using ConstraintsT = ComponentConstraints<IdxT,
                                                BusInternalConstraints,
                                                NoConstraints>;

      Bus(const DataT& data, const StateT& state);
      ~Bus() override = default;

      VariablesT&         variables();
      const VariablesT&   variables() const;
      ConstraintsT&       constraints();
      const ConstraintsT& constraints() const;

      IdxT sizeInternalVariables() const override;
      IdxT sizeInternalConstraints() const override;
      IdxT nnz() const override;
      void setVariableOffset(IdxT offset) override;
      void setConstraintOffset(IdxT offset) override;
      int  initialize(ScalarT* values) const override;
      int  setVariableBounds(RealT* lower, RealT* upper) const override;
      int  updateSolutionState(const ScalarT*    values,
                               Model::StateData& state) const override;

    private:
      DataT        data_;
      StateT       state_;
      VariablesT   variables_;
      ConstraintsT constraints_;
    };

  } // namespace OPF
} // namespace GridKit
