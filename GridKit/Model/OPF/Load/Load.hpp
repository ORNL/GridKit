#pragma once

#include <array>
#include <cstddef>
#include <span>

#include <GridKit/Model/OPF/Component.hpp>
#include <GridKit/Model/OPF/Load/LoadData.hpp>

namespace GridKit
{
  namespace OPF
  {
    /// Fixed signed complex-power injection read from the input state.
    template <class scalar_type, typename index_type>
    class Load final : public Component<scalar_type, index_type>
    {
    public:
      using ScalarT        = scalar_type;
      using IdxT           = index_type;
      using ComponentT     = Component<ScalarT, IdxT>;
      using RealT          = typename ComponentT::RealT;
      using DataT          = LoadData<RealT, IdxT>;
      using StateT         = Model::DeviceState;
      using JacobianEntryT = typename ComponentT::JacobianEntryT;
      using HessianEntryT  = typename ComponentT::HessianEntryT;

      static constexpr std::size_t VARIABLE_COUNT   = 0;
      static constexpr std::size_t CONSTRAINT_COUNT = 2;

      Load(const DataT&                       data,
           const StateT&                      state,
           std::array<IdxT, CONSTRAINT_COUNT> constraint_indices,
           std::array<IdxT, 2>                bus_voltage_indices);
      ~Load() override = default;

      std::span<const IdxT>           variableIndices() const override;
      std::span<const IdxT>           constraintIndices() const override;
      std::span<const JacobianEntryT> jacobianEntries() const override;
      std::span<const HessianEntryT>  hessianEntries() const override;

      int gatherVariables(const ScalarT* global_variables) override;
      int evaluateObjective() override;
      int evaluateObjectiveGradient() override;
      int evaluateConstraints() override;
      int evaluateJacobian() override;
      int evaluateHessian(
          RealT                  objective_factor,
          std::span<const RealT> local_multipliers) override;

      RealT                    objective() const override;
      std::span<const ScalarT> objectiveGradientValues() const override;
      std::span<const ScalarT> constraintValues() const override;
      std::span<const RealT>   jacobianValues() const override;
      std::span<const RealT>   hessianValues() const override;

      bool hasJacobian() const override;
      bool hasHessian() const override;

      int updateSolutionState(const ScalarT*    global_variables,
                              Model::StateData& state) const override;

      ScalarT evaluateObjective(const ScalarT* local_variables) const;
      void    evaluateConstraints(const ScalarT* local_variables,
                                  ScalarT*       local_values) const;

    private:
      DataT                                 data_;
      RealT                                 in_service_{1};
      RealT                                 p_{0};
      RealT                                 q_{0};
      std::array<IdxT, CONSTRAINT_COUNT>    constraint_indices_;
      std::array<IdxT, 2>                   bus_voltage_indices_;
      std::array<ScalarT, CONSTRAINT_COUNT> constraint_values_{};
      RealT                                 objective_{0};
    };
  } // namespace OPF
} // namespace GridKit
