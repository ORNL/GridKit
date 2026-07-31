#pragma once

#include <cstddef>

#include <GridKit/Model/OPF/Component.hpp>
#include <GridKit/Model/OPF/ComponentConstraints.hpp>
#include <GridKit/Model/OPF/ComponentVariables.hpp>
#include <GridKit/Model/OPF/Load/LoadData.hpp>
#include <GridKit/Model/StateData.hpp>

namespace GridKit
{
  namespace OPF
  {
    enum class LoadExternalVariables : std::size_t
    {
      VM, ///< Bus voltage magnitude
      VA, ///< Bus voltage angle
      MAXIMUM
    };

    enum class LoadExternalConstraints : std::size_t
    {
      DIVP, ///< Bus active-power balance
      DIVQ, ///< Bus reactive-power balance
      MAXIMUM
    };

    /**
     * @brief Fixed demand read from the companion state data.
     *
     * The base OPF formulation does not dispatch load. Its p/q state is held
     * fixed while contributing to the connected bus balances. Load limits in
     * the system data validate that operating state.
     */
    template <class scalar_type, typename index_type>
    class Load final : public Component<scalar_type, index_type>
    {
    public:
      using ScalarT      = scalar_type;
      using IdxT         = index_type;
      using ComponentT   = Component<ScalarT, IdxT>;
      using RealT        = typename ComponentT::RealT;
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
