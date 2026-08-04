#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "Branch.hpp"

namespace GridKit
{
  namespace OPF
  {
    template <class ScalarT, typename IdxT>
    Branch<ScalarT, IdxT>::Branch(
        const DataT&                     data,
        const StateT&                    state,
        std::array<IdxT, VARIABLE_COUNT> variable_indices,
        std::span<const IdxT>            active_constraint_indices)
      : data_(data),
        in_service_(state.open.value_or(false) ? RealT{0} : RealT{1}),
        variable_indices_(variable_indices),
        constraint_indices_(active_constraint_indices.begin(),
                            active_constraint_indices.end())
    {
      const std::size_t expected_constraint_count =
          data_.smax.has_value() ? CONSTRAINT_COUNT : std::size_t{4};
      if (constraint_indices_.size() != expected_constraint_count)
      {
        throw std::invalid_argument(
            "OPF branch requires four balance rows and two rows when rated");
      }

      const RealT tap   = static_cast<RealT>(state.tap.value_or(1.0));
      const RealT phase = static_cast<RealT>(state.phase.value_or(0.0));
      if (!std::isfinite(data_.R) || !std::isfinite(data_.X)
          || !std::isfinite(data_.G) || !std::isfinite(data_.B)
          || (data_.R == RealT{0} && data_.X == RealT{0})
          || !std::isfinite(tap) || !(tap > RealT{0})
          || !std::isfinite(phase))
      {
        throw std::invalid_argument("OPF branch parameters are invalid");
      }

      const RealT impedance_scale = std::max(std::abs(data_.R),
                                             std::abs(data_.X));
      const RealT scaled_r        = data_.R / impedance_scale;
      const RealT scaled_x        = data_.X / impedance_scale;
      const RealT scaled_norm     = scaled_r * scaled_r + scaled_x * scaled_x;
      const RealT g               = (scaled_r / scaled_norm) / impedance_scale;
      const RealT b               = -(scaled_x / scaled_norm) / impedance_scale;
      const RealT inverse_tap     = RealT{1} / tap;
      const RealT c               = std::cos(phase);
      const RealT s               = std::sin(phase);

      gff_ = (g + RealT{0.5} * data_.G) * inverse_tap * inverse_tap;
      bff_ = (b + RealT{0.5} * data_.B) * inverse_tap * inverse_tap;
      gft_ = -(g * c - b * s) * inverse_tap;
      bft_ = -(g * s + b * c) * inverse_tap;
      gtf_ = -(g * c + b * s) * inverse_tap;
      btf_ = (g * s - b * c) * inverse_tap;
      gtt_ = g + RealT{0.5} * data_.G;
      btt_ = b + RealT{0.5} * data_.B;

      const std::array<RealT, 8> coefficients{{gff_, bff_, gft_, bft_, gtf_, btf_, gtt_, btt_}};
      if (!std::all_of(coefficients.begin(),
                       coefficients.end(),
                       [](RealT value)
                       { return std::isfinite(value); }))
      {
        throw std::invalid_argument(
            "OPF branch operating inputs produce non-finite admittance");
      }
    }

    template <class ScalarT, typename IdxT>
    std::span<const IdxT> Branch<ScalarT, IdxT>::variableIndices() const
    {
      return variable_indices_;
    }

    template <class ScalarT, typename IdxT>
    std::span<const IdxT> Branch<ScalarT, IdxT>::constraintIndices() const
    {
      return constraint_indices_;
    }

    template <class ScalarT, typename IdxT>
    std::span<const typename Branch<ScalarT, IdxT>::JacobianEntryT>
    Branch<ScalarT, IdxT>::jacobianEntries() const
    {
      if (data_.smax.has_value())
      {
        return RATED_JACOBIAN_ENTRIES;
      }
      return UNRATED_JACOBIAN_ENTRIES;
    }

    template <class ScalarT, typename IdxT>
    std::span<const typename Branch<ScalarT, IdxT>::HessianEntryT>
    Branch<ScalarT, IdxT>::hessianEntries() const
    {
      return HESSIAN_ENTRIES;
    }

    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::gatherVariables(const ScalarT* global_variables)
    {
      if (global_variables == nullptr)
      {
        return 1;
      }
      for (std::size_t local = 0; local < VARIABLE_COUNT; ++local)
      {
        variables_[local] = global_variables[variable_indices_[local]];
      }
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::evaluateObjective()
    {
      objective_ = RealT{0};
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::evaluateObjectiveGradient()
    {
      objective_gradient_values_.fill(ScalarT{0});
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::evaluateConstraints()
    {
      evaluateConstraints(variables_.data(), constraint_values_.data());
      for (const ScalarT value : constraint_values_)
      {
        if (!std::isfinite(static_cast<RealT>(value)))
        {
          return 1;
        }
      }
      return 0;
    }

    template <class ScalarT, typename IdxT>
    typename Branch<ScalarT, IdxT>::RealT Branch<ScalarT, IdxT>::objective() const
    {
      return objective_;
    }

    template <class ScalarT, typename IdxT>
    std::span<const ScalarT> Branch<ScalarT, IdxT>::objectiveGradientValues() const
    {
      return objective_gradient_values_;
    }

    template <class ScalarT, typename IdxT>
    std::span<const ScalarT> Branch<ScalarT, IdxT>::constraintValues() const
    {
      return {constraint_values_.data(), constraint_indices_.size()};
    }

    template <class ScalarT, typename IdxT>
    std::span<const typename Branch<ScalarT, IdxT>::RealT>
    Branch<ScalarT, IdxT>::jacobianValues() const
    {
      return {jacobian_values_.data(), jacobianEntries().size()};
    }

    template <class ScalarT, typename IdxT>
    std::span<const typename Branch<ScalarT, IdxT>::RealT>
    Branch<ScalarT, IdxT>::hessianValues() const
    {
      return hessian_values_;
    }

    template <class ScalarT, typename IdxT>
    bool Branch<ScalarT, IdxT>::hasJacobian() const
    {
      return true;
    }

    template <class ScalarT, typename IdxT>
    bool Branch<ScalarT, IdxT>::hasHessian() const
    {
      return true;
    }

    template <class ScalarT, typename IdxT>
    ScalarT Branch<ScalarT, IdxT>::evaluateObjective(
        [[maybe_unused]] const ScalarT* local_variables) const
    {
      return ScalarT{0};
    }

    template <class ScalarT, typename IdxT>
    TerminalPower<ScalarT> Branch<ScalarT, IdxT>::terminalPower(
        const ScalarT* local_variables) const
    {
      const ScalarT vf = local_variables[0];
      const ScalarT af = local_variables[1];
      const ScalarT vt = local_variables[2];
      const ScalarT at = local_variables[3];

      const ScalarT delta      = af - at;
      const ScalarT c          = std::cos(delta);
      const ScalarT s          = std::sin(delta);
      const ScalarT in_service = static_cast<ScalarT>(in_service_);

      const ScalarT pf =
          in_service
          * (static_cast<ScalarT>(gff_) * vf * vf
             + vf * vt
                   * (static_cast<ScalarT>(gft_) * c
                      + static_cast<ScalarT>(bft_) * s));
      const ScalarT qf =
          in_service
          * (-static_cast<ScalarT>(bff_) * vf * vf
             + vf * vt
                   * (static_cast<ScalarT>(gft_) * s
                      - static_cast<ScalarT>(bft_) * c));
      const ScalarT pt =
          in_service
          * (static_cast<ScalarT>(gtt_) * vt * vt
             + vf * vt
                   * (static_cast<ScalarT>(gtf_) * c
                      - static_cast<ScalarT>(btf_) * s));
      const ScalarT qt =
          in_service
          * (-static_cast<ScalarT>(btt_) * vt * vt
             + vf * vt
                   * (-static_cast<ScalarT>(gtf_) * s
                      - static_cast<ScalarT>(btf_) * c));
      return {pf, qf, pt, qt};
    }

    template <class ScalarT, typename IdxT>
    void Branch<ScalarT, IdxT>::evaluateConstraints(
        const ScalarT* local_variables,
        ScalarT*       local_values) const
    {
      const auto flow = terminalPower(local_variables);

      local_values[0] = -flow.p_from;
      local_values[1] = -flow.q_from;
      local_values[2] = -flow.p_to;
      local_values[3] = -flow.q_to;
      local_values[4] = flow.p_from * flow.p_from
                        + flow.q_from * flow.q_from;
      local_values[5] = flow.p_to * flow.p_to + flow.q_to * flow.q_to;
    }
  } // namespace OPF
} // namespace GridKit
