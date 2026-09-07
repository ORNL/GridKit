#pragma once

#include "GastPtiImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::assembleJacobian(RealT y_scale, RealT yp_scale)
      {
        gatherExternalVariables();
        if (J_rows_buffer_ == nullptr)
        {
          const size_t capacity = 63 * this->externalJacobianExpansion();
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
        auto internal = [&](IdxT row, IdxT column, RealT value, RealT derivative = ZERO<RealT>)
        {
          if (y_scale != ZERO<RealT> || (yp_scale != ZERO<RealT> && derivative != ZERO<RealT>) )
            entry(row, this->getVariableIndex(column), y_scale * value + yp_scale * derivative);
        };
        auto external = [&](IdxT row, size_t slot, RealT value)
        {
          if (y_scale == ZERO<RealT>)
            return;
          auto* signal = this->externalVariableSignals()[slot];
          if (signal != nullptr)
          {
            typename SignalT::GradientT gradient;
            signal->appendGradient(gradient, y_scale * value);
            for (const auto& [column, derivative] : gradient)
              entry(row, column, derivative);
          }
        };
        const auto* y         = y_.getData();
        const RealT valve     = static_cast<RealT>(y[0]);
        const RealT rate      = static_cast<RealT>(y[5]) - valve;
        const RealT lower     = Math::sigmoid(valve - Vmin_response_);
        const RealT upper     = Math::sigmoid(Vmax_response_ - valve);
        const RealT positive  = Math::sigmoid(rate);
        const RealT negative  = Math::sigmoid(-rate);
        const RealT mu        = Math::MU<RealT>;
        const RealT lower_x   = mu * lower * (ONE<RealT> - lower);
        const RealT upper_x   = -mu * upper * (ONE<RealT> - upper);
        const RealT gate      = lower * upper + (ONE<RealT> - upper) * negative + (ONE<RealT> - lower) * positive;
        const RealT gate_x    = lower_x * (upper - positive) + upper_x * (lower - negative);
        const RealT gate_rate = -(ONE<RealT> - upper) * mu * negative * (ONE<RealT> - negative)
                                + (ONE<RealT> - lower) * mu * positive * (ONE<RealT> - positive);
        const RealT rate_derivative = gate + rate * gate_rate;
        const RealT ratio           = va_machine_base_ / va_component_base_;
        internal(0, 0, s_valve_ * (rate * gate_x - rate_derivative) / T1_, -ONE<RealT>);
        internal(0, 5, s_valve_ * rate_derivative / T1_);
        internal(1, 0, ONE<RealT> / T2_);
        internal(1, 1, -ONE<RealT> / T2_, -ONE<RealT>);
        internal(2, 1, ONE<RealT> / T3_);
        internal(2, 2, -ONE<RealT> / T3_, -ONE<RealT>);
        internal(3, 3, -R_);
        internal(4, 2, -Kt_);
        internal(4, 4, -ONE<RealT>);
        const RealT selector = Math::sigmoid(static_cast<RealT>(y[3]) - static_cast<RealT>(y[4]));
        internal(5, 3, ONE<RealT> - selector);
        internal(5, 4, selector);
        internal(5, 5, -ONE<RealT>);
        internal(6, 1, ONE<RealT>);
        internal(6, 6, -ratio);
        external(3, 0, -ONE<RealT>);
        external(3, 1, R_ * ratio);
        external(6, 0, -Dturb_);
        this->constructCoo();
        return 0;
      }
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
