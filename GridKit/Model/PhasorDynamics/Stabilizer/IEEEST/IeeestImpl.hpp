#pragma once

/**
 * @file IeeestImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the IEEEST Power System Stabilizer.
 */

#include <cmath>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/Ieeest.hpp>
#include <GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/IeeestData.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Stabilizer
    {
      using Log = ::GridKit::Utilities::Logger;

      template <typename scalar_type, typename index_type>
      Ieeest<scalar_type, index_type>::Ieeest()
      {
        size_ = 12;
      }

      template <typename scalar_type, typename index_type>
      Ieeest<scalar_type, index_type>::Ieeest(const ModelDataT& data)
        : monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        initializeMonitor();
        size_ = 12;
      }

      template <typename scalar_type, typename index_type>
      Ieeest<scalar_type, index_type>::~Ieeest()
      {
      }

      template <typename scalar_type, typename index_type>
      void Ieeest<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Parameter = typename ModelDataT::Parameters;
        if (data.parameters.contains(Parameter::A1))
        {
          A1_ = std::get<RealT>(data.parameters.at(Parameter::A1));
        }
        if (data.parameters.contains(Parameter::A2))
        {
          A2_ = std::get<RealT>(data.parameters.at(Parameter::A2));
        }
        if (data.parameters.contains(Parameter::A3))
        {
          A3_ = std::get<RealT>(data.parameters.at(Parameter::A3));
        }
        if (data.parameters.contains(Parameter::A4))
        {
          A4_ = std::get<RealT>(data.parameters.at(Parameter::A4));
        }
        if (data.parameters.contains(Parameter::A5))
        {
          A5_ = std::get<RealT>(data.parameters.at(Parameter::A5));
        }
        if (data.parameters.contains(Parameter::A6))
        {
          A6_ = std::get<RealT>(data.parameters.at(Parameter::A6));
        }
        if (data.parameters.contains(Parameter::T1))
        {
          T1_ = std::get<RealT>(data.parameters.at(Parameter::T1));
        }
        if (data.parameters.contains(Parameter::T2))
        {
          T2_ = std::get<RealT>(data.parameters.at(Parameter::T2));
        }
        if (data.parameters.contains(Parameter::T3))
        {
          T3_ = std::get<RealT>(data.parameters.at(Parameter::T3));
        }
        if (data.parameters.contains(Parameter::T4))
        {
          T4_ = std::get<RealT>(data.parameters.at(Parameter::T4));
        }
        if (data.parameters.contains(Parameter::T5))
        {
          T5_ = std::get<RealT>(data.parameters.at(Parameter::T5));
        }
        if (data.parameters.contains(Parameter::T6))
        {
          T6_ = std::get<RealT>(data.parameters.at(Parameter::T6));
        }
        if (data.parameters.contains(Parameter::Ks))
        {
          Ks_ = std::get<RealT>(data.parameters.at(Parameter::Ks));
        }
        if (data.parameters.contains(Parameter::Lsmin))
        {
          Lsmin_ = std::get<RealT>(data.parameters.at(Parameter::Lsmin));
        }
        if (data.parameters.contains(Parameter::Lsmax))
        {
          Lsmax_ = std::get<RealT>(data.parameters.at(Parameter::Lsmax));
        }
        if (data.parameters.contains(Parameter::Vcl))
        {
          Vcl_ = std::get<RealT>(data.parameters.at(Parameter::Vcl));
        }
        if (data.parameters.contains(Parameter::Vcu))
        {
          Vcu_ = std::get<RealT>(data.parameters.at(Parameter::Vcu));
        }
        if (data.parameters.contains(Parameter::Tdelay))
        {
          Tdelay_ = std::get<RealT>(data.parameters.at(Parameter::Tdelay));
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

        use_T4_block_    = static_cast<RealT>(T4_ != 0.0);
        bypass_T4_block_ = 1.0 - use_T4_block_;

        use_T6_block_    = static_cast<RealT>(T6_ != 0.0);
        bypass_T6_block_ = 1.0 - use_T6_block_;
      }

      template <typename scalar_type, typename index_type>
      int Ieeest<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Ieeest<scalar_type, index_type>::allocate()
      {
        if (!this->allocated_)
        {
          this->allocateVectors(this->size_);
        }
        auto size = static_cast<size_t>(size_);

        tag_.resize(size);

        variable_indices_.resize(size);
        residual_indices_.resize(size);
        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        ws_.resize(1);
        ws_indices_.resize(1);
        ws_[0]         = 0.0;
        ws_indices_[0] = INVALID_INDEX<IdxT>;

        if (signals_.template isAssigned<IeeestInternalVariables::VSS>())
        {
          auto* y = y_.getData();
          signals_.template getSignalNode<IeeestInternalVariables::VSS>()->set(
              &y[11], &(this->getVariableIndex(11)));
        }

        this->allocated_ = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Ieeest<scalar_type, index_type>::verify() const
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

        if (a4_ == 0 && a3_ == 0 && a2_ == 0 && a1_ != 0)
        {
          Log::error() << "Ieeest: a2, a3, and a4 are all zero - no valid notch filter\n";
          ret += 1;
        }

        return ret;
      }

      template <typename scalar_type, typename index_type>
      int Ieeest<scalar_type, index_type>::initialize()
      {
        auto* y  = y_.getData();
        auto* yp = yp_.getData();

        for (IdxT i = 0; i < size_; ++i)
        {
          y[static_cast<size_t>(i)]  = 0.0;
          yp[static_cast<size_t>(i)] = 0.0;
        }

        y_.setDataUpdated();
        yp_.setDataUpdated();

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Ieeest<scalar_type, index_type>::tagDifferentiable()
      {
        tag_[0]  = true;
        tag_[1]  = true;
        tag_[2]  = true;
        tag_[3]  = true;
        tag_[4]  = (T2_ != 0.0);
        tag_[5]  = (T4_ != 0.0);
        tag_[6]  = (T6_ != 0.0);
        tag_[7]  = false;
        tag_[8]  = false;
        tag_[9]  = false;
        tag_[10] = false;
        tag_[11] = false;

        return 0;
      }

      /**
       * @brief Compute the absolute tolerance for each variable in the model
       *
       * @param rel_tol The relative tolerance which can be used to pick the
       *        absolute tolerance.
       * @tparam scalar_type Scalar data type
       * @tparam index_type Index data type
       * @return int 0 if successful, non-zero otherwise.
       *
       * This represents a "noise" level close to zero for which pure relative
       * error cannot be used.
       */
      template <typename scalar_type, typename index_type>
      int Ieeest<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline int Ieeest<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT*                  y,
          const ScalarT*                  yp,
          [[maybe_unused]] const ScalarT* wb,
          const ScalarT*                  ws,
          ScalarT*                        f)
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

        ScalarT x1_dot = yp[0];
        ScalarT x2_dot = yp[1];
        ScalarT x3_dot = yp[2];
        ScalarT x4_dot = yp[3];
        ScalarT x5_dot = yp[4];
        ScalarT x6_dot = yp[5];
        ScalarT x7_dot = yp[6];

        ScalarT u = ws[0];

        f[0] = -x1_dot + use_notch_ * x2;
        f[1] = -x2_dot + (use_4th_order_ + use_3rd_order_) * x3
               + use_2nd_order_ * (-a0_ * x1 - a1_ * x2 + u) * safe_inv_a2_;
        f[2] = -x3_dot + use_4th_order_ * x4
               + use_3rd_order_ * (-a0_ * x1 - a1_ * x2 - a2_ * x3 + u) * safe_inv_a3_;
        f[3]  = -x4_dot + use_4th_order_ * (-a0_ * x1 - a1_ * x2 - a2_ * x3 - a3_ * x4 + u) * safe_inv_a4_;
        f[4]  = -T2_ * x5_dot - x5 + v4;
        f[5]  = -T4_ * x6_dot - x6 + v5;
        f[6]  = -T6_ * x7_dot - x7 + v6;
        f[7]  = -v4 + bypass_notch_ * u + use_notch_ * (x1 + A5_ * x2 + (use_4th_order_ + use_3rd_order_) * A6_ * x3);
        f[8]  = use_T2_block_ * (-T2_ * (v5 - x5) + T1_ * (v4 - x5)) + bypass_T2_block_ * (v4 - v5);
        f[9]  = use_T4_block_ * (-T4_ * (v6 - x6) + T3_ * (v5 - x6)) + bypass_T4_block_ * (v5 - v6);
        f[10] = use_T6_block_ * (-T6_ * v7 + Ks_ * T5_ * (v6 - x7)) + bypass_T6_block_ * (Ks_ * v6 - v7);
        f[11] = -vss + Math::clamp(v7, Lsmin_, Lsmax_);

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Ieeest<scalar_type, index_type>::evaluateResidual()
      {
        if (signals_.template isAttached<IeeestExternalVariables::U>())
        {
          ws_[0]         = signals_.template readExternalVariable<IeeestExternalVariables::U>();
          ws_indices_[0] = signals_.template readExternalVariableIndex<IeeestExternalVariables::U>();
        }

        const auto* y  = y_.getData();
        const auto* yp = yp_.getData();
        auto*       f  = f_.getData();
        evaluateInternalResidual(y, yp, wb_.data(), ws_.data(), f);
        f_.setDataUpdated();

        return 0;
      }

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* Ieeest<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      template <typename scalar_type, typename index_type>
      void Ieeest<scalar_type, index_type>::initializeMonitor()
      {
        using Variable = typename ModelDataT::MonitorableVariables;
        monitor_->set(Variable::vss, [this]
                      { return y_.getData()[11]; });
      }

    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
