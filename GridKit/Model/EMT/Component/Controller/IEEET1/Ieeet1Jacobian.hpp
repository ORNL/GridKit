#pragma once

#include "Ieeet1Impl.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::evaluateJacobian()
      {
        gatherExternalVariables();
        if (J_rows_buffer_ == nullptr)
        {
          const size_t capacity = 81 * this->externalJacobianExpansion();
          J_rows_buffer_        = new IdxT[capacity];
          J_cols_buffer_        = new IdxT[capacity];
          J_vals_buffer_        = new RealT[capacity];
        }
        nnz_       = 0;
        auto entry = [&](IdxT row, IdxT column, RealT value)
        {
          J_rows_buffer_[nnz_]   = this->getResidualIndex(row);
          J_cols_buffer_[nnz_]   = column;
          J_vals_buffer_[nnz_++] = value;
        };
        auto internal = [&](IdxT row, IdxT column, RealT value)
        {
          entry(row, this->getVariableIndex(column), value);
        };
        auto external = [&](IdxT row, size_t slot, RealT value)
        {
          auto* signal = this->externalVariableSignals()[slot];
          if (signal != nullptr)
          {
            typename SignalT::GradientT gradient;
            signal->appendGradient(gradient, value);
            for (const auto& [column, derivative] : gradient)
              entry(row, column, derivative);
          }
        };
        const auto* y        = y_.getData();
        const RealT vr       = static_cast<RealT>(y[1]);
        const RealT efdp     = static_cast<RealT>(y[2]);
        const RealT rate     = (-vr + Ka_ * static_cast<RealT>(y[4])) / Ta_;
        const RealT lower    = Math::sigmoid(vr - Vrmin_);
        const RealT upper    = Math::sigmoid(Vrmax_ - vr);
        const RealT positive = Math::sigmoid(rate);
        const RealT negative = Math::sigmoid(-rate);
        const RealT mu       = Math::MU<RealT>;
        const RealT lower_x  = mu * lower * (ONE<RealT> - lower);
        const RealT upper_x  = -mu * upper * (ONE<RealT> - upper);
        const RealT gate     = lower * upper + (ONE<RealT> - upper) * negative
                           + (ONE<RealT> - lower) * positive;
        const RealT gate_x    = lower_x * (upper - positive) + upper_x * (lower - negative);
        const RealT gate_rate = -(ONE<RealT> - upper) * mu * negative * (ONE<RealT> - negative)
                                + (ONE<RealT> - lower) * mu * positive * (ONE<RealT> - positive);
        const RealT limited_rate_derivative = gate + rate * gate_rate;
        internal(0, 0, -ONE<RealT> / Tr_ - alpha_);
        internal(1, 1, rate * gate_x - limited_rate_derivative / Ta_ - alpha_);
        internal(1, 4, limited_rate_derivative * Ka_ / Ta_);
        internal(2, 1, ONE<RealT> / Te_);
        internal(2, 2, -Ke_eff_ / Te_ - alpha_);
        internal(2, 6, -ONE<RealT> / Te_);
        internal(3, 3, -alpha_);
        internal(3, 5, ONE<RealT> / Tf_);
        internal(4, 0, -ONE<RealT>);
        internal(4, 4, -ONE<RealT>);
        internal(4, 5, -ONE<RealT>);
        internal(5, 2, Kf_);
        internal(5, 3, -Tf_);
        internal(5, 5, -Tf_);
        internal(6, 6, -ONE<RealT>);
        internal(6, 8, ONE<RealT>);
        internal(7, 2, ONE<RealT> + (static_cast<RealT>(y_ext_[3]) - ONE<RealT>) *Ispdlim_);
        internal(7, 7, -ONE<RealT>);
        const RealT x     = efdp - SA_;
        const RealT sigma = Math::sigmoid(x);
        internal(8, 2, SB_ * (static_cast<RealT>(2) * x * sigma + x * x * mu * sigma * (ONE<RealT> - sigma)));
        internal(8, 8, -ONE<RealT>);
        const RealT va               = static_cast<RealT>(y_ext_[0]);
        const RealT vb               = static_cast<RealT>(y_ext_[1]);
        const RealT vc               = static_cast<RealT>(y_ext_[2]);
        const RealT norm             = std::sqrt(va * va + vb * vb + vc * vc);
        const RealT voltage_gradient = norm == ZERO<RealT> ? ZERO<RealT> : ONE<RealT> / (V_ * Tr_ * norm);
        external(0, 0, va * voltage_gradient);
        external(0, 1, vb * voltage_gradient);
        external(0, 2, vc * voltage_gradient);
        external(7, 3, efdp * Ispdlim_);
        external(4, 4, ONE<RealT>);
        external(4, 5, ONE<RealT>);
        external(4, 6, uel_on_);
        external(4, 7, oel_on_);
        this->constructCoo();
        return 0;
      }
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
