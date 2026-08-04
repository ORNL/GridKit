#pragma once

#include <array>
#include <cstddef>
#include <span>
#include <vector>

#include <GridKit/Model/OPF/Branch/BranchData.hpp>
#include <GridKit/Model/OPF/Component.hpp>

namespace GridKit
{
  namespace OPF
  {
    template <class scalar_type>
    struct TerminalPower
    {
      scalar_type p_from;
      scalar_type q_from;
      scalar_type p_to;
      scalar_type q_to;
    };

    /// Fixed-tap pi-model branch with optional terminal apparent-power limits.
    template <class scalar_type, typename index_type>
    class Branch final : public Component<scalar_type, index_type>
    {
    public:
      using ScalarT        = scalar_type;
      using IdxT           = index_type;
      using ComponentT     = Component<ScalarT, IdxT>;
      using RealT          = typename ComponentT::RealT;
      using DataT          = BranchData<RealT, IdxT>;
      using StateT         = Model::DeviceState;
      using JacobianEntryT = typename ComponentT::JacobianEntryT;
      using HessianEntryT  = typename ComponentT::HessianEntryT;

      static constexpr std::size_t VARIABLE_COUNT   = 4;
      static constexpr std::size_t CONSTRAINT_COUNT = 6;

      Branch(const DataT&                     data,
             const StateT&                    state,
             std::array<IdxT, VARIABLE_COUNT> variable_indices,
             std::span<const IdxT>            active_constraint_indices);
      ~Branch() override = default;

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

      ScalarT                evaluateObjective(const ScalarT* local_variables) const;
      TerminalPower<ScalarT> terminalPower(const ScalarT* local_variables) const;
      void                   evaluateConstraints(const ScalarT* local_variables,
                                                 ScalarT*       local_values) const;

    private:
      inline static constexpr std::array<JacobianEntryT, 16>
          UNRATED_JACOBIAN_ENTRIES{{{IdxT{0}, IdxT{0}},
                                    {IdxT{1}, IdxT{0}},
                                    {IdxT{2}, IdxT{0}},
                                    {IdxT{3}, IdxT{0}},
                                    {IdxT{0}, IdxT{1}},
                                    {IdxT{1}, IdxT{1}},
                                    {IdxT{2}, IdxT{1}},
                                    {IdxT{3}, IdxT{1}},
                                    {IdxT{0}, IdxT{2}},
                                    {IdxT{1}, IdxT{2}},
                                    {IdxT{2}, IdxT{2}},
                                    {IdxT{3}, IdxT{2}},
                                    {IdxT{0}, IdxT{3}},
                                    {IdxT{1}, IdxT{3}},
                                    {IdxT{2}, IdxT{3}},
                                    {IdxT{3}, IdxT{3}}}};

      inline static constexpr std::array<JacobianEntryT, 24>
          RATED_JACOBIAN_ENTRIES{{{IdxT{0}, IdxT{0}},
                                  {IdxT{1}, IdxT{0}},
                                  {IdxT{2}, IdxT{0}},
                                  {IdxT{3}, IdxT{0}},
                                  {IdxT{4}, IdxT{0}},
                                  {IdxT{5}, IdxT{0}},
                                  {IdxT{0}, IdxT{1}},
                                  {IdxT{1}, IdxT{1}},
                                  {IdxT{2}, IdxT{1}},
                                  {IdxT{3}, IdxT{1}},
                                  {IdxT{4}, IdxT{1}},
                                  {IdxT{5}, IdxT{1}},
                                  {IdxT{0}, IdxT{2}},
                                  {IdxT{1}, IdxT{2}},
                                  {IdxT{2}, IdxT{2}},
                                  {IdxT{3}, IdxT{2}},
                                  {IdxT{4}, IdxT{2}},
                                  {IdxT{5}, IdxT{2}},
                                  {IdxT{0}, IdxT{3}},
                                  {IdxT{1}, IdxT{3}},
                                  {IdxT{2}, IdxT{3}},
                                  {IdxT{3}, IdxT{3}},
                                  {IdxT{4}, IdxT{3}},
                                  {IdxT{5}, IdxT{3}}}};

      inline static constexpr std::array<HessianEntryT, 10>
          HESSIAN_ENTRIES{{{IdxT{0}, IdxT{0}},
                           {IdxT{1}, IdxT{0}},
                           {IdxT{2}, IdxT{0}},
                           {IdxT{3}, IdxT{0}},
                           {IdxT{1}, IdxT{1}},
                           {IdxT{2}, IdxT{1}},
                           {IdxT{3}, IdxT{1}},
                           {IdxT{2}, IdxT{2}},
                           {IdxT{3}, IdxT{2}},
                           {IdxT{3}, IdxT{3}}}};

      DataT                            data_;
      RealT                            in_service_{1};
      RealT                            gff_{0};
      RealT                            bff_{0};
      RealT                            gft_{0};
      RealT                            bft_{0};
      RealT                            gtf_{0};
      RealT                            btf_{0};
      RealT                            gtt_{0};
      RealT                            btt_{0};
      std::array<IdxT, VARIABLE_COUNT> variable_indices_;
      std::vector<IdxT>                constraint_indices_;

      std::array<ScalarT, VARIABLE_COUNT>              variables_{};
      std::array<ScalarT, VARIABLE_COUNT>              objective_gradient_values_{};
      std::array<ScalarT, CONSTRAINT_COUNT>            constraint_values_{};
      std::array<RealT, RATED_JACOBIAN_ENTRIES.size()> jacobian_values_{};
      std::array<RealT, HESSIAN_ENTRIES.size()>        hessian_values_{};
      std::array<RealT, CONSTRAINT_COUNT>              local_multipliers_{};
      RealT                                            objective_{0};
    };
  } // namespace OPF
} // namespace GridKit
