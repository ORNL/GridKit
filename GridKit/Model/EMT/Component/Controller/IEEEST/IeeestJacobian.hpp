#pragma once

#include "IeeestImpl.hpp"

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      template <typename scalar_type, typename index_type>
      int Ieeest<scalar_type, index_type>::evaluateJacobian()
      {
        gatherExternalVariables();
        if (J_rows_buffer_ == nullptr)
        {
          const size_t capacity = 180 * this->externalJacobianExpansion();
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
        const RealT b2   = use_2nd_order_ * safe_inv_a2_;
        const RealT b3   = use_3rd_order_ * safe_inv_a3_;
        const RealT b4   = use_4th_order_ * safe_inv_a4_;
        const RealT high = use_4th_order_ + use_3rd_order_;
        internal(0, 0, -alpha_);
        internal(0, 1, use_notch_);
        internal(1, 0, -b2);
        internal(1, 1, -a1_ * b2 - alpha_);
        internal(1, 2, high);
        internal(2, 0, -b3);
        internal(2, 1, -a1_ * b3);
        internal(2, 2, -a2_ * b3 - alpha_);
        internal(2, 3, use_4th_order_);
        internal(3, 0, -b4);
        internal(3, 1, -a1_ * b4);
        internal(3, 2, -a2_ * b4);
        internal(3, 3, -a3_ * b4 - alpha_);
        internal(4, 4, -ONE<RealT> - T2_ * alpha_);
        internal(4, 7, ONE<RealT>);
        internal(5, 5, -ONE<RealT> - T4_ * alpha_);
        internal(5, 8, ONE<RealT>);
        internal(6, 6, -ONE<RealT> - T6_ * alpha_);
        internal(6, 9, ONE<RealT>);
        internal(7, 0, use_notch_ * (ONE<RealT> - A6_ * b2));
        internal(7, 1, use_notch_ * (A5_ - A6_ * a1_ * b2));
        internal(7, 2, use_notch_ * A6_ * high);
        internal(7, 7, -ONE<RealT>);
        internal(8, 4, use_T2_block_ * (T2_ - T1_));
        internal(8, 7, use_T2_block_ * T1_ + bypass_T2_block_);
        internal(8, 8, -use_T2_block_ * T2_ - bypass_T2_block_);
        internal(9, 5, use_T4_block_ * (T4_ - T3_));
        internal(9, 8, use_T4_block_ * T3_ + bypass_T4_block_);
        internal(9, 9, -use_T4_block_ * T4_ - bypass_T4_block_);
        internal(10, 6, -use_T6_block_ * Ks_ * T5_);
        internal(10, 9, use_T6_block_ * Ks_ * T5_ + bypass_T6_block_ * Ks_);
        internal(10, 10, -use_T6_block_ * T6_ - bypass_T6_block_);
        const RealT voltage  = static_cast<RealT>(y_ext_[2]);
        const RealT lower    = Math::sigmoid(voltage - Vcl_);
        const RealT upper    = Math::sigmoid(Vcu_ - voltage);
        const RealT low_gate = ONE<RealT> - lower_cutout_ + lower_cutout_ * lower;
        const RealT up_gate  = ONE<RealT> - upper_cutout_ + upper_cutout_ * upper;
        const RealT gate     = low_gate * up_gate;
        const RealT gate_v   = Math::MU<RealT> * (lower_cutout_ * lower * (ONE<RealT> - lower) * up_gate - upper_cutout_ * upper * (ONE<RealT> - upper) * low_gate);
        const RealT v7       = static_cast<RealT>(y_.getData()[10]);
        internal(11, 10, gate * (Math::sigmoid(v7 - Lsmin_) - Math::sigmoid(v7 - Lsmax_)));
        internal(11, 11, -ONE<RealT>);
        for (size_t slot = 0; slot < 2; ++slot)
        {
          const RealT scale = slot == 0 ? ONE<RealT> - use_speed_ : use_speed_;
          external(1, slot, scale * b2);
          external(2, slot, scale * b3);
          external(3, slot, scale * b4);
          external(7, slot, scale * (bypass_notch_ + use_notch_ * A6_ * b2));
        }
        external(11, 2, gate_v * Math::clamp(v7, Lsmin_, Lsmax_));
        this->constructCoo();
        return 0;
      }
    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
