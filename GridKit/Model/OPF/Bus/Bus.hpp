#pragma once

#include <array>
#include <cstddef>
#include <span>

#include <GridKit/Model/OPF/Bus/BusData.hpp>
#include <GridKit/Model/OPF/Component.hpp>

namespace GridKit
{
  namespace OPF
  {
    /// Voltage variables and power-balance row ownership for one OPF bus.
    template <class scalar_type, typename index_type>
    class Bus final : public Component<scalar_type, index_type>
    {
    public:
      using ScalarT        = scalar_type;
      using IdxT           = index_type;
      using ComponentT     = Component<ScalarT, IdxT>;
      using RealT          = typename ComponentT::RealT;
      using DataT          = BusData<RealT, IdxT>;
      using StateT         = Model::BusState;
      using JacobianEntryT = typename ComponentT::JacobianEntryT;
      using HessianEntryT  = typename ComponentT::HessianEntryT;

      static constexpr std::size_t VARIABLE_COUNT   = 2;
      static constexpr std::size_t CONSTRAINT_COUNT = 0;

      Bus(const DataT&                     data,
          const StateT&                    state,
          std::array<IdxT, VARIABLE_COUNT> variable_indices,
          std::array<IdxT, 2>              balance_indices);
      ~Bus() override = default;

      std::span<const IdxT> variableIndices() const override;
      std::span<const IdxT> constraintIndices() const override;
      std::span<const IdxT> balanceIndices() const;

      IdxT voltageMagnitudeIndex() const;
      IdxT voltageAngleIndex() const;
      IdxT activeBalanceIndex() const;
      IdxT reactiveBalanceIndex() const;

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

      int initialize(ScalarT* global_variables) const override;
      int setVariableBounds(RealT* global_lower_bounds,
                            RealT* global_upper_bounds) const override;
      int updateSolutionState(const ScalarT*    global_variables,
                              Model::StateData& state) const override;

      ScalarT evaluateObjective(const ScalarT* local_variables) const;
      void    evaluateConstraints(const ScalarT* local_variables,
                                  ScalarT*       local_values) const;

    private:
      DataT                               data_;
      StateT                              state_;
      std::array<IdxT, VARIABLE_COUNT>    variable_indices_;
      std::array<IdxT, 2>                 balance_indices_;
      std::array<ScalarT, VARIABLE_COUNT> variables_{};
      std::array<ScalarT, VARIABLE_COUNT> objective_gradient_values_{};
      RealT                               objective_{0};
    };
  } // namespace OPF
} // namespace GridKit
