#pragma once

#include "SexsPtiImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::evaluateJacobian()
      {
        gatherExternalVariables();
        if (J_rows_buffer_ == nullptr)
        {
          const size_t capacity = 44 * this->externalJacobianExpansion();
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
        const auto* y         = y_.getData();
        const RealT efd       = static_cast<RealT>(y[1]);
        const RealT rate      = (-efd + K_ / Tb_ * (-static_cast<RealT>(y[0]) + Ta_ * static_cast<RealT>(y[2]))) / Te_;
        const RealT lower     = Math::sigmoid(efd - Efdmin_);
        const RealT upper     = Math::sigmoid(Efdmax_ - efd);
        const RealT positive  = Math::sigmoid(rate);
        const RealT negative  = Math::sigmoid(-rate);
        const RealT mu        = Math::MU<RealT>;
        const RealT lower_x   = mu * lower * (ONE<RealT> - lower);
        const RealT upper_x   = -mu * upper * (ONE<RealT> - upper);
        const RealT gate      = lower * upper + (ONE<RealT> - upper) * negative + (ONE<RealT> - lower) * positive;
        const RealT gate_x    = lower_x * (upper - positive) + upper_x * (lower - negative);
        const RealT gate_rate = -(ONE<RealT> - upper) * mu * negative * (ONE<RealT> - negative)
                                + (ONE<RealT> - lower) * mu * positive * (ONE<RealT> - positive);
        const RealT limited_rate_derivative = gate + rate * gate_rate;
        internal(0, 0, -ONE<RealT> / Tb_ - alpha_);
        internal(0, 2, Ta_ / Tb_ - ONE<RealT>);
        internal(1, 0, -limited_rate_derivative * K_ / (Tb_ * Te_));
        internal(1, 1, rate * gate_x - limited_rate_derivative / Te_ - alpha_);
        internal(1, 2, limited_rate_derivative * K_ * Ta_ / (Tb_ * Te_));
        internal(2, 2, -ONE<RealT>);
        internal(2, 3, -ONE<RealT>);
        internal(3, 3, -ONE<RealT> - Tr_ * alpha_);
        external(2, 0, ONE<RealT>);
        external(2, 1, ONE<RealT>);
        external(2, 2, uel_on_);
        external(2, 3, oel_on_);
        const RealT va       = static_cast<RealT>(y_ext_[4]);
        const RealT vb       = static_cast<RealT>(y_ext_[5]);
        const RealT vc       = static_cast<RealT>(y_ext_[6]);
        const RealT norm     = std::sqrt(va * va + vb * vb + vc * vc);
        const RealT gradient = norm == ZERO<RealT> ? ZERO<RealT> : ONE<RealT> / (V_ * norm);
        external(3, 4, va * gradient);
        external(3, 5, vb * gradient);
        external(3, 6, vc * gradient);
        this->constructCoo();
        return 0;
      }
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
