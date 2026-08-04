#pragma once

#include <array>
#include <cstddef>
#include <span>

#include <GridKit/Model/OPF/Component.hpp>
#include <GridKit/Model/OPF/Generator/GeneratorData.hpp>

namespace GridKit
{
  namespace OPF
  {
    /// Dispatchable active/reactive generation and quadratic operating cost.
    template <class scalar_type, typename index_type>
    class Generator final : public Component<scalar_type, index_type>
    {
    public:
      using ScalarT        = scalar_type;
      using IdxT           = index_type;
      using ComponentT     = Component<ScalarT, IdxT>;
      using RealT          = typename ComponentT::RealT;
      using DataT          = GeneratorData<RealT, IdxT>;
      using StateT         = Model::DeviceState;
      using JacobianEntryT = typename ComponentT::JacobianEntryT;
      using HessianEntryT  = typename ComponentT::HessianEntryT;

      static constexpr std::size_t VARIABLE_COUNT   = 2;
      static constexpr std::size_t CONSTRAINT_COUNT = 2;

      Generator(const DataT&                       data,
                const StateT&                      state,
                std::array<IdxT, VARIABLE_COUNT>   variable_indices,
                std::array<IdxT, CONSTRAINT_COUNT> constraint_indices,
                std::array<IdxT, 2>                bus_voltage_indices);
      ~Generator() override = default;

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

      int initialize(ScalarT* global_variables) const override;
      int setVariableBounds(RealT* global_lower_bounds,
                            RealT* global_upper_bounds) const override;
      int updateSolutionState(const ScalarT*    global_variables,
                              Model::StateData& state) const override;

      ScalarT evaluateObjective(const ScalarT* local_variables) const;
      void    evaluateConstraints(const ScalarT* local_variables,
                                  ScalarT*       local_values) const;

    private:
      inline static constexpr std::array<JacobianEntryT, 2>
          JACOBIAN_ENTRIES{{{IdxT{0}, IdxT{0}}, {IdxT{1}, IdxT{1}}}};
      inline static constexpr std::array<HessianEntryT, 1>
          HESSIAN_ENTRIES{{{IdxT{0}, IdxT{0}}}};

      DataT                              data_;
      StateT                             state_;
      RealT                              in_service_{1};
      std::array<IdxT, VARIABLE_COUNT>   variable_indices_;
      std::array<IdxT, CONSTRAINT_COUNT> constraint_indices_;
      std::array<IdxT, 2>                bus_voltage_indices_;

      std::array<ScalarT, VARIABLE_COUNT>        variables_{};
      std::array<ScalarT, VARIABLE_COUNT>        objective_gradient_values_{};
      std::array<ScalarT, CONSTRAINT_COUNT>      constraint_values_{};
      std::array<RealT, JACOBIAN_ENTRIES.size()> jacobian_values_{};
      std::array<RealT, HESSIAN_ENTRIES.size()>  hessian_values_{};
      std::array<RealT, CONSTRAINT_COUNT>        local_multipliers_{};
      RealT                                      objective_{0};
    };
  } // namespace OPF
} // namespace GridKit
