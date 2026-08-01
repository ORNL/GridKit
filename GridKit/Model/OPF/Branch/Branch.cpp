#include "Branch.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace GridKit
{
  namespace OPF
  {
    namespace
    {
      template <class ScalarT>
      struct BranchEvaluation
      {
        std::array<ScalarT, 4>                powers{};
        std::array<std::array<ScalarT, 4>, 4> jacobian{};
      };

      template <typename RealT>
      bool validBranchData(RealT R,
                           RealT X,
                           RealT G,
                           RealT B,
                           RealT tap,
                           RealT phase)
      {
        const RealT denominator = R * R + X * X;
        if (!std::isfinite(R) || !std::isfinite(X) || !std::isfinite(G)
            || !std::isfinite(B) || !std::isfinite(denominator)
            || !(denominator > RealT{0}) || !std::isfinite(tap)
            || !(tap > RealT{0}) || !std::isfinite(phase))
        {
          return false;
        }

        const RealT                g                   = R / denominator;
        const RealT                b                   = -X / denominator;
        const RealT                inverse_tap         = RealT{1} / tap;
        const RealT                inverse_tap_squared = inverse_tap * inverse_tap;
        const RealT                cos_phase           = std::cos(phase);
        const RealT                sin_phase           = std::sin(phase);
        const std::array<RealT, 8> admittance{
            -(g + RealT{0.5} * G) * inverse_tap_squared,
            -(b + RealT{0.5} * B) * inverse_tap_squared,
            (g * cos_phase - b * sin_phase) * inverse_tap,
            (b * cos_phase + g * sin_phase) * inverse_tap,
            (g * cos_phase + b * sin_phase) * inverse_tap,
            (b * cos_phase - g * sin_phase) * inverse_tap,
            -(g + RealT{0.5} * G),
            -(b + RealT{0.5} * B)};
        return std::isfinite(inverse_tap)
               && std::isfinite(inverse_tap_squared)
               && std::all_of(admittance.begin(),
                              admittance.end(),
                              [](RealT value)
                              { return std::isfinite(value); });
      }

      template <class ScalarT>
      bool validEvaluation(const BranchEvaluation<ScalarT>& evaluation)
      {
        for (const ScalarT value : evaluation.powers)
        {
          if (!std::isfinite(value))
          {
            return false;
          }
        }
        for (const auto& row : evaluation.jacobian)
        {
          for (const ScalarT value : row)
          {
            if (!std::isfinite(value))
            {
              return false;
            }
          }
        }
        return true;
      }

      template <class ScalarT, typename RealT>
      BranchEvaluation<ScalarT> evaluateBranch(RealT   R,
                                               RealT   X,
                                               RealT   G,
                                               RealT   B,
                                               RealT   tap,
                                               RealT   phase,
                                               bool    open,
                                               ScalarT vmf,
                                               ScalarT vaf,
                                               ScalarT vmt,
                                               ScalarT vat)
      {
        BranchEvaluation<ScalarT> result;
        if (open)
        {
          return result;
        }

        const RealT denominator = R * R + X * X;
        const RealT g           = R / denominator;
        const RealT b           = -X / denominator;
        const RealT inv_tap     = RealT{1} / tap;
        const RealT cos_phase   = std::cos(phase);
        const RealT sin_phase   = std::sin(phase);

        const RealT gff = -(g + RealT{0.5} * G) * inv_tap * inv_tap;
        const RealT bff = -(b + RealT{0.5} * B) * inv_tap * inv_tap;
        const RealT gft = (g * cos_phase - b * sin_phase) * inv_tap;
        const RealT bft = (b * cos_phase + g * sin_phase) * inv_tap;
        const RealT gtf = (g * cos_phase + b * sin_phase) * inv_tap;
        const RealT btf = (b * cos_phase - g * sin_phase) * inv_tap;
        const RealT gtt = -(g + RealT{0.5} * G);
        const RealT btt = -(b + RealT{0.5} * B);

        const ScalarT delta     = vaf - vat;
        const ScalarT cos_delta = std::cos(delta);
        const ScalarT sin_delta = std::sin(delta);
        const ScalarT A         = static_cast<ScalarT>(gft) * cos_delta
                          + static_cast<ScalarT>(bft) * sin_delta;
        const ScalarT C = static_cast<ScalarT>(gft) * sin_delta
                          - static_cast<ScalarT>(bft) * cos_delta;
        const ScalarT D = static_cast<ScalarT>(gtf) * cos_delta
                          - static_cast<ScalarT>(btf) * sin_delta;
        const ScalarT E = -static_cast<ScalarT>(gtf) * sin_delta
                          - static_cast<ScalarT>(btf) * cos_delta;

        const ScalarT pf = static_cast<ScalarT>(gff) * vmf * vmf + vmf * vmt * A;
        const ScalarT qf = -static_cast<ScalarT>(bff) * vmf * vmf + vmf * vmt * C;
        const ScalarT pt = static_cast<ScalarT>(gtt) * vmt * vmt + vmf * vmt * D;
        const ScalarT qt = -static_cast<ScalarT>(btt) * vmt * vmt + vmf * vmt * E;
        result.powers    = {pf, qf, pt, qt};

        result.jacobian[0] = {
            ScalarT{2} * static_cast<ScalarT>(gff) * vmf + vmt * A,
            -vmf * vmt * C,
            vmf * A,
            vmf * vmt * C};
        result.jacobian[1] = {
            -ScalarT{2} * static_cast<ScalarT>(bff) * vmf + vmt * C,
            vmf * vmt * A,
            vmf * C,
            -vmf * vmt * A};
        result.jacobian[2] = {
            vmt * D,
            vmf * vmt * E,
            ScalarT{2} * static_cast<ScalarT>(gtt) * vmt + vmf * D,
            -vmf * vmt * E};
        result.jacobian[3] = {
            vmt * E,
            -vmf * vmt * D,
            -ScalarT{2} * static_cast<ScalarT>(btt) * vmt + vmf * E,
            vmf * vmt * D};
        return result;
      }
    } // namespace

    template <class ScalarT, typename IdxT>
    Branch<ScalarT, IdxT>::Branch(const DataT& data, const StateT& state)
      : data_(data),
        state_(state)
    {
    }

    template <class ScalarT, typename IdxT>
    typename Branch<ScalarT, IdxT>::VariablesT& Branch<ScalarT, IdxT>::variables()
    {
      return variables_;
    }

    template <class ScalarT, typename IdxT>
    const typename Branch<ScalarT, IdxT>::VariablesT& Branch<ScalarT, IdxT>::variables() const
    {
      return variables_;
    }

    template <class ScalarT, typename IdxT>
    typename Branch<ScalarT, IdxT>::ConstraintsT& Branch<ScalarT, IdxT>::constraints()
    {
      return constraints_;
    }

    template <class ScalarT, typename IdxT>
    const typename Branch<ScalarT, IdxT>::ConstraintsT& Branch<ScalarT, IdxT>::constraints() const
    {
      return constraints_;
    }

    template <class ScalarT, typename IdxT>
    IdxT Branch<ScalarT, IdxT>::sizeInternalVariables() const
    {
      return static_cast<IdxT>(VariablesT::sizeInternal());
    }

    template <class ScalarT, typename IdxT>
    IdxT Branch<ScalarT, IdxT>::sizeInternalConstraints() const
    {
      return data_.smax.has_value()
                 ? static_cast<IdxT>(ConstraintsT::sizeInternal())
                 : IdxT{0};
    }

    template <class ScalarT, typename IdxT>
    IdxT Branch<ScalarT, IdxT>::nnz() const
    {
      return data_.smax.has_value() ? IdxT{24} : IdxT{16};
    }

    template <class ScalarT, typename IdxT>
    void Branch<ScalarT, IdxT>::setConstraintOffset(IdxT offset)
    {
      if (data_.smax.has_value())
      {
        constraints_.setInternalOffset(offset);
      }
    }

    template <class ScalarT, typename IdxT>
    void Branch<ScalarT, IdxT>::addJacobianSparsity(
        std::vector<typename ComponentT::JacobianEntry>& entries) const
    {
      const std::array<IdxT, 4> columns{
          variables_.template externalIndex<BranchExternalVariables::VMF>(),
          variables_.template externalIndex<BranchExternalVariables::VAF>(),
          variables_.template externalIndex<BranchExternalVariables::VMT>(),
          variables_.template externalIndex<BranchExternalVariables::VAT>()};
      const std::array<IdxT, 4> external_rows{
          constraints_.template externalIndex<BranchExternalConstraints::DIVPF>(),
          constraints_.template externalIndex<BranchExternalConstraints::DIVQF>(),
          constraints_.template externalIndex<BranchExternalConstraints::DIVPT>(),
          constraints_.template externalIndex<BranchExternalConstraints::DIVQT>()};

      for (const IdxT row : external_rows)
      {
        for (const IdxT column : columns)
        {
          entries.emplace_back(row, column);
        }
      }

      if (data_.smax.has_value())
      {
        const std::array<IdxT, 2> internal_rows{
            constraints_.template internalIndex<BranchInternalConstraints::SF2>(),
            constraints_.template internalIndex<BranchInternalConstraints::ST2>()};
        for (const IdxT row : internal_rows)
        {
          for (const IdxT column : columns)
          {
            entries.emplace_back(row, column);
          }
        }
      }
    }

    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::setJacobianSlots(const std::vector<IdxT>& slots)
    {
      if (slots.size() != static_cast<std::size_t>(nnz()))
      {
        return 1;
      }
      jacobian_slots_ = slots;
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::setConstraintBounds(RealT* lower,
                                                   RealT* upper) const
    {
      if (!data_.smax.has_value())
      {
        return 0;
      }
      if (!std::isfinite(*data_.smax) || *data_.smax < RealT{0})
      {
        return 1;
      }

      const RealT upper_bound = *data_.smax * *data_.smax;
      if (!std::isfinite(upper_bound))
      {
        return 1;
      }
      const IdxT sf2 = constraints_.template internalIndex<BranchInternalConstraints::SF2>();
      const IdxT st2 = constraints_.template internalIndex<BranchInternalConstraints::ST2>();
      lower[sf2]     = RealT{0};
      upper[sf2]     = upper_bound;
      lower[st2]     = RealT{0};
      upper[st2]     = upper_bound;
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::addConstraints(const ScalarT* values,
                                              ScalarT*       constraints) const
    {
      const RealT tap   = static_cast<RealT>(state_.tap.value_or(1.0));
      const RealT phase = static_cast<RealT>(state_.phase.value_or(0.0));
      if (!validBranchData(data_.R, data_.X, data_.G, data_.B, tap, phase))
      {
        return 1;
      }

      const auto evaluation = evaluateBranch<ScalarT>(
          data_.R,
          data_.X,
          data_.G,
          data_.B,
          tap,
          phase,
          state_.open.value_or(false),
          variables_.template external<BranchExternalVariables::VMF>(values),
          variables_.template external<BranchExternalVariables::VAF>(values),
          variables_.template external<BranchExternalVariables::VMT>(values),
          variables_.template external<BranchExternalVariables::VAT>(values));
      if (!validEvaluation(evaluation))
      {
        return 1;
      }

      ScalarT sf2{};
      ScalarT st2{};
      if (data_.smax.has_value())
      {
        sf2 = evaluation.powers[0] * evaluation.powers[0]
              + evaluation.powers[1] * evaluation.powers[1];
        st2 = evaluation.powers[2] * evaluation.powers[2]
              + evaluation.powers[3] * evaluation.powers[3];
        if (!std::isfinite(static_cast<RealT>(sf2))
            || !std::isfinite(static_cast<RealT>(st2)))
        {
          return 1;
        }
      }

      constraints[constraints_.template externalIndex<BranchExternalConstraints::DIVPF>()] += evaluation.powers[0];
      constraints[constraints_.template externalIndex<BranchExternalConstraints::DIVQF>()] += evaluation.powers[1];
      constraints[constraints_.template externalIndex<BranchExternalConstraints::DIVPT>()] += evaluation.powers[2];
      constraints[constraints_.template externalIndex<BranchExternalConstraints::DIVQT>()] += evaluation.powers[3];

      if (data_.smax.has_value())
      {
        constraints[constraints_.template internalIndex<BranchInternalConstraints::SF2>()] += sf2;
        constraints[constraints_.template internalIndex<BranchInternalConstraints::ST2>()] += st2;
      }
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::addJacobian(const ScalarT* values,
                                           RealT*         jacobian_values) const
    {
      if (jacobian_slots_.size() != static_cast<std::size_t>(nnz()))
      {
        return 1;
      }

      const RealT tap   = static_cast<RealT>(state_.tap.value_or(1.0));
      const RealT phase = static_cast<RealT>(state_.phase.value_or(0.0));
      if (!validBranchData(data_.R, data_.X, data_.G, data_.B, tap, phase))
      {
        return 1;
      }

      const auto evaluation = evaluateBranch<ScalarT>(
          data_.R,
          data_.X,
          data_.G,
          data_.B,
          tap,
          phase,
          state_.open.value_or(false),
          variables_.template external<BranchExternalVariables::VMF>(values),
          variables_.template external<BranchExternalVariables::VAF>(values),
          variables_.template external<BranchExternalVariables::VMT>(values),
          variables_.template external<BranchExternalVariables::VAT>(values));
      if (!validEvaluation(evaluation))
      {
        return 1;
      }

      std::array<RealT, 24> contributions{};
      std::size_t           contribution = 0;
      for (const auto& row : evaluation.jacobian)
      {
        for (const ScalarT value : row)
        {
          contributions[contribution] = static_cast<RealT>(value);
          ++contribution;
        }
      }

      if (data_.smax.has_value())
      {
        for (std::size_t column = 0; column < 4; ++column)
        {
          const ScalarT value = ScalarT{2}
                                * (evaluation.powers[0] * evaluation.jacobian[0][column]
                                   + evaluation.powers[1] * evaluation.jacobian[1][column]);
          if (!std::isfinite(static_cast<RealT>(value)))
          {
            return 1;
          }
          contributions[contribution] = static_cast<RealT>(value);
          ++contribution;
        }
        for (std::size_t column = 0; column < 4; ++column)
        {
          const ScalarT value = ScalarT{2}
                                * (evaluation.powers[2] * evaluation.jacobian[2][column]
                                   + evaluation.powers[3] * evaluation.jacobian[3][column]);
          if (!std::isfinite(static_cast<RealT>(value)))
          {
            return 1;
          }
          contributions[contribution] = static_cast<RealT>(value);
          ++contribution;
        }
      }

      for (std::size_t i = 0; i < contribution; ++i)
      {
        jacobian_values[jacobian_slots_[i]] += contributions[i];
      }
      return 0;
    }

    template class Branch<double, int>;
    template class Branch<double, long int>;
    template class Branch<double, std::size_t>;

  } // namespace OPF
} // namespace GridKit
