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
    enum class LoadExternalVariables : std::size_t
    {
      VOLTAGE_MAGNITUDE,
      VOLTAGE_ANGLE,
      MAXIMUM
    };

    enum class LoadExternalConstraints : std::size_t
    {
      ACTIVE_POWER_BALANCE,
      REACTIVE_POWER_BALANCE,
      MAXIMUM
    };

    /**
     * @brief Fixed demand read from the companion state data.
     *
     * The base OPF formulation does not dispatch load. Its p/q state is held
     * fixed while contributing to the connected bus balances. Load limits in
     * the case data are retained for validation and future controllable-load
     * formulations.
     */
    template <class scalar_type, typename index_type>
    class Load : public Component<scalar_type, index_type>
    {
    public:
      using ScalarT      = scalar_type;
      using IdxT         = index_type;
      using RealT        = typename Component<ScalarT, IdxT>::RealT;
      using DataT        = LoadData<RealT, IdxT>;
      using StateT       = Model::DeviceState;
      using VariablesT   = ComponentVariables<ScalarT,
                                              IdxT,
                                              NoVariables,
                                              LoadExternalVariables>;
      using ConstraintsT = ComponentConstraints<IdxT,
                                                NoConstraints,
                                                LoadExternalConstraints>;

      Load(const DataT& data, const StateT& state);
      ~Load() override = default;

      VariablesT&         variables();
      const VariablesT&   variables() const;
      ConstraintsT&       constraints();
      const ConstraintsT& constraints() const;

      IdxT sizeInternalVariables() const override;
      IdxT sizeInternalConstraints() const override;
      IdxT nnz() const override;
      int  addConstraints(const ScalarT* values, ScalarT* constraints) const override;
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
