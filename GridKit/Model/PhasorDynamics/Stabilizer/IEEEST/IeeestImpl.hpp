#pragma once

/**
 * @file IeeestImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the IEEEST Power System Stabilizer.
 */

#include <cmath>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/ConnectedElementImpl.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/Ieeest.hpp>
#include <GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/IeeestData.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Stabilizer
    {
      using Log = ::GridKit::Utilities::Logger;

      template <typename ScalarP, typename IdxP>
      Ieeest<ScalarP, IdxP>::Ieeest()
      {
        size_ = 13;
      }

      template <typename ScalarP, typename IdxP>
      Ieeest<ScalarP, IdxP>::Ieeest(const ModelDataT& data)
        : ConnectedElement<Ieeest>(data)
      {
        initializeParameters(data);
        initializeMonitor();
        size_ = 13;
      }

      template <typename ScalarP, typename IdxP>
      Ieeest<ScalarP, IdxP>::~Ieeest()
      {
      }

      template <typename ScalarP, typename IdxP>
      void Ieeest<ScalarP, IdxP>::initializeParameters(const ModelDataT& data)
      {
        if (data.parameters.contains(ModelDataT::Parameters::A1))
        {
          A1_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::A1));
        }
        if (data.parameters.contains(ModelDataT::Parameters::A2))
        {
          A2_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::A2));
        }
        if (data.parameters.contains(ModelDataT::Parameters::A3))
        {
          A3_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::A3));
        }
        if (data.parameters.contains(ModelDataT::Parameters::A4))
        {
          A4_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::A4));
        }
        if (data.parameters.contains(ModelDataT::Parameters::A5))
        {
          A5_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::A5));
        }
        if (data.parameters.contains(ModelDataT::Parameters::A6))
        {
          A6_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::A6));
        }
        if (data.parameters.contains(ModelDataT::Parameters::T1))
        {
          T1_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::T1));
        }
        if (data.parameters.contains(ModelDataT::Parameters::T2))
        {
          T2_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::T2));
        }
        if (data.parameters.contains(ModelDataT::Parameters::T3))
        {
          T3_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::T3));
        }
        if (data.parameters.contains(ModelDataT::Parameters::T4))
        {
          T4_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::T4));
        }
        if (data.parameters.contains(ModelDataT::Parameters::T5))
        {
          T5_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::T5));
        }
        if (data.parameters.contains(ModelDataT::Parameters::T6))
        {
          T6_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::T6));
        }
        if (data.parameters.contains(ModelDataT::Parameters::Ks))
        {
          Ks_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::Ks));
        }
        if (data.parameters.contains(ModelDataT::Parameters::Lsmin))
        {
          Lsmin_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::Lsmin));
        }
        if (data.parameters.contains(ModelDataT::Parameters::Lsmax))
        {
          Lsmax_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::Lsmax));
        }
        if (data.parameters.contains(ModelDataT::Parameters::Vcl))
        {
          Vcl_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::Vcl));
        }
        if (data.parameters.contains(ModelDataT::Parameters::Vcu))
        {
          Vcu_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::Vcu));
        }
        if (data.parameters.contains(ModelDataT::Parameters::Tdelay))
        {
          Tdelay_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::Tdelay));
        }

        a0_ = 1;
        a1_ = A1_ + A3_;
        a2_ = A2_ + A4_ + A1_ * A3_;
        a3_ = A1_ * A4_ + A2_ * A3_;
        a4_ = A2_ * A4_;

        // Precompute masks and safe inverse coefficients so the residual stays branch-free.
        use_notch_    = static_cast<RealT>(a2_ != 0.0 || a3_ != 0.0 || a4_ != 0.0);
        bypass_notch_ = 1.0 - use_notch_;

        use_4th_order_ = static_cast<RealT>(a4_ != 0.0);
        use_3rd_order_ = static_cast<RealT>(a4_ == 0.0 && a3_ != 0.0);
        use_2nd_order_ = static_cast<RealT>(a4_ == 0.0 && a3_ == 0.0 && a2_ != 0.0);
        safe_inv_a4_   = use_4th_order_ / (a4_ + (1.0 - use_4th_order_));
        safe_inv_a3_   = use_3rd_order_ / (a3_ + (1.0 - use_3rd_order_));
        safe_inv_a2_   = use_2nd_order_ / (a2_ + (1.0 - use_2nd_order_));

        use_T2_block_    = static_cast<RealT>(T2_ != 0.0);
        bypass_T2_block_ = 1.0 - use_T2_block_;
        safe_inv_T2_     = use_T2_block_ / (T2_ + bypass_T2_block_);

        use_T4_block_    = static_cast<RealT>(T4_ != 0.0);
        bypass_T4_block_ = 1.0 - use_T4_block_;
        safe_inv_T4_     = use_T4_block_ / (T4_ + bypass_T4_block_);

        use_T6_block_    = static_cast<RealT>(T6_ != 0.0);
        bypass_T6_block_ = 1.0 - use_T6_block_;
        safe_inv_T6_     = use_T6_block_ / (T6_ + bypass_T6_block_);

        // Vcl = Vcu = 0 means "cutout disabled" (passthrough: vs = vss).
        use_cutout_    = static_cast<RealT>(Vcl_ != 0.0 || Vcu_ != 0.0);
        bypass_cutout_ = 1.0 - use_cutout_;
      }

      template <typename ScalarP, typename IdxP>
      int Ieeest<ScalarP, IdxP>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      template <typename ScalarP, typename IdxP>
      int Ieeest<ScalarP, IdxP>::allocate()
      {
        auto size = static_cast<size_t>(size_);
        f_.resize(size);
        y_.resize(size);
        yp_.resize(size);
        tag_.resize(size);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        ws_.resize(2);
        ws_indices_.resize(2);
        ws_[0]         = 0.0;
        ws_indices_[0] = INVALID_INDEX<IdxT>;
        ws_[1]         = 0.0;
        ws_indices_[1] = INVALID_INDEX<IdxT>;

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        if (signals_.template isAssigned<IeeestInternalVariables::VS>())
        {
          signals_.template getSignalNode<IeeestInternalVariables::VS>()->set(
              &y_[12], &(this->getVariableIndex(12)));
        }

        return 0;
      }

      template <typename ScalarP, typename IdxP>
      int Ieeest<ScalarP, IdxP>::verify() const
      {
        int ret = 0;

        if (signals_.template isAttached<IeeestExternalVariables::U>())
        {
          if (!signals_.template isLinked<IeeestExternalVariables::U>())
          {
            Log::error() << "Ieeest: input signal U attached with no linked source\n";
            ret += 1;
          }
        }
        else
        {
          Log::error() << "Ieeest: required input signal U is not attached\n";
          ret += 1;
        }

        if (signals_.template isAttached<IeeestExternalVariables::VCT>())
        {
          if (!signals_.template isLinked<IeeestExternalVariables::VCT>())
          {
            Log::error() << "Ieeest: cutout signal VCT attached with no linked source\n";
            ret += 1;
          }
        }

        if (a4_ == 0 && a3_ == 0 && a2_ == 0 && a1_ != 0)
        {
          Log::error() << "Ieeest: a2, a3, and a4 are all zero - no valid notch filter\n";
          ret += 1;
        }

        return ret;
      }

      template <typename ScalarP, typename IdxP>
      int Ieeest<ScalarP, IdxP>::initialize()
      {
        for (IdxT i = 0; i < size_; ++i)
        {
          y_[static_cast<size_t>(i)]  = 0.0;
          yp_[static_cast<size_t>(i)] = 0.0;
        }

        return 0;
      }

      template <typename ScalarP, typename IdxP>
      int Ieeest<ScalarP, IdxP>::tagDifferentiable()
      {
        tag_[0]  = true;
        tag_[1]  = true;
        tag_[2]  = true;
        tag_[3]  = true;
        tag_[4]  = true;
        tag_[5]  = true;
        tag_[6]  = true;
        tag_[7]  = false;
        tag_[8]  = false;
        tag_[9]  = false;
        tag_[10] = false;
        tag_[11] = false;
        tag_[12] = false;

        return 0;
      }

      template <typename ScalarP, typename IdxP>
      __attribute__((always_inline)) inline int Ieeest<ScalarP, IdxP>::evaluateInternalResidual(
          ScalarT*                  y,
          ScalarT*                  yp,
          [[maybe_unused]] ScalarT* wb,
          ScalarT*                  ws,
          ScalarT*                  f)
      {
        ScalarT x1  = y[0];
        ScalarT x2  = y[1];
        ScalarT x3  = y[2];
        ScalarT x4  = y[3];
        ScalarT x5  = y[4];
        ScalarT x6  = y[5];
        ScalarT x7  = y[6];
        ScalarT v4  = y[7];
        ScalarT v5  = y[8];
        ScalarT v6  = y[9];
        ScalarT v7  = y[10];
        ScalarT vss = y[11];
        ScalarT vs  = y[12];

        ScalarT x1_dot = yp[0];
        ScalarT x2_dot = yp[1];
        ScalarT x3_dot = yp[2];
        ScalarT x4_dot = yp[3];
        ScalarT x5_dot = yp[4];
        ScalarT x6_dot = yp[5];
        ScalarT x7_dot = yp[6];

        ScalarT u   = ws[0];
        ScalarT vct = ws[1];

        f[0] = -x1_dot + use_notch_ * x2;
        f[1] = -x2_dot + (use_4th_order_ + use_3rd_order_) * x3
               + use_2nd_order_ * (-a0_ * x1 - a1_ * x2 + u) * safe_inv_a2_;
        f[2] = -x3_dot + use_4th_order_ * x4
               + use_3rd_order_ * (-a0_ * x1 - a1_ * x2 - a2_ * x3 + u) * safe_inv_a3_;
        f[3] = -x4_dot + use_4th_order_ * (-a0_ * x1 - a1_ * x2 - a2_ * x3 - a3_ * x4 + u) * safe_inv_a4_;

        f[4] = -x5_dot + use_T2_block_ * (v4 - x5) * safe_inv_T2_;
        f[5] = -x6_dot + use_T4_block_ * (v5 - x6) * safe_inv_T4_;
        f[6] = -x7_dot + use_T6_block_ * (v6 - x7) * safe_inv_T6_;

        f[7]  = -v4 + bypass_notch_ * u + use_notch_ * (x1 + A5_ * x2 + (use_4th_order_ + use_3rd_order_) * A6_ * x3);
        f[8]  = -v5 + bypass_T2_block_ * v4 + use_T2_block_ * (x5 + T1_ * (v4 - x5) * safe_inv_T2_);
        f[9]  = -v6 + bypass_T4_block_ * v5 + use_T4_block_ * (x6 + T3_ * (v5 - x6) * safe_inv_T4_);
        f[10] = -v7 + bypass_T6_block_ * Ks_ * v6 + use_T6_block_ * Ks_ * T5_ * (v6 - x7) * safe_inv_T6_;

        f[11] = -vss + v7 * Math::indicator_zero(Lsmin_, Lsmax_, v7)
                + Lsmin_ * Math::sigmoid(Lsmin_ - v7)
                + Lsmax_ * Math::sigmoid(v7 - Lsmax_);
        // Cutout: Vcl=Vcu=0 means disabled (passthrough). Otherwise, sigmoid window.
        f[12] = -vs + bypass_cutout_ * vss
                + use_cutout_ * vss * Math::sigmoid(vct - Vcl_) * Math::sigmoid(Vcu_ - vct);

        return 0;
      }

      template <typename ScalarP, typename IdxP>
      int Ieeest<ScalarP, IdxP>::evaluateResidual()
      {
        if (signals_.template isAttached<IeeestExternalVariables::U>())
        {
          ws_[0]         = signals_.template readExternalVariable<IeeestExternalVariables::U>();
          ws_indices_[0] = signals_.template readExternalVariableIndex<IeeestExternalVariables::U>();
        }

        ws_[1] = (Vcl_ + Vcu_) / TWO<RealT>;
        if (signals_.template isAttached<IeeestExternalVariables::VCT>())
        {
          ws_[1]         = signals_.template readExternalVariable<IeeestExternalVariables::VCT>();
          ws_indices_[1] = signals_.template readExternalVariableIndex<IeeestExternalVariables::VCT>();
        }

        evaluateInternalResidual(y_.data(), yp_.data(), wb_.data(), ws_.data(), f_.data());

        return 0;
      }

      template <typename ScalarP, typename IdxP>
      void Ieeest<ScalarP, IdxP>::initializeMonitor()
      {
        using Variable = typename ModelDataT::MonitorableVariables;
        monitor_->set(Variable::vs, [this]
                      { return y_[12]; });
      }

    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
