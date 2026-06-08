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
      void Ieeest<scalar_type, index_type>::setDerivedParameters()
      {
        a1_ = A1_ + A3_;
        a2_ = A2_ + A4_ + A1_ * A3_;
        a3_ = A1_ * A4_ + A2_ * A3_;
        a4_ = A2_ * A4_;
      }

      template <typename scalar_type, typename index_type>
      void Ieeest<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
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

        setDerivedParameters();
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
        auto size = static_cast<size_t>(size_);
        f_.resize(size);
        y_.resize(size);
        yp_.resize(size);
        tag_.resize(size);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        ws_.resize(1);
        ws_indices_.resize(1);
        ws_[0]         = 0.0;
        ws_indices_[0] = INVALID_INDEX<IdxT>;

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        if (signals_.template isAssigned<IeeestInternalVariables::VSS>())
        {
          signals_.template getSignalNode<IeeestInternalVariables::VSS>()->set(
              &y_[11], &(this->getVariableIndex(11)));
        }

        tagDifferentiable();

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

        if (a2_ == ZERO<RealT> && a3_ == ZERO<RealT> && a4_ == ZERO<RealT> && a1_ != ZERO<RealT>)
        {
          Log::error() << "Ieeest: a2, a3, and a4 are all zero - no valid notch filter\n";
          ret += 1;
        }

        return ret;
      }

      template <typename scalar_type, typename index_type>
      int Ieeest<scalar_type, index_type>::initialize()
      {
        for (IdxT i = 0; i < size_; ++i)
        {
          y_[static_cast<size_t>(i)]  = 0.0;
          yp_[static_cast<size_t>(i)] = 0.0;
        }

        ScalarT u{0.0};
        if (signals_.template isAttached<IeeestExternalVariables::U>())
        {
          u              = signals_.template readExternalVariable<IeeestExternalVariables::U>();
          ws_[0]         = u;
          ws_indices_[0] = signals_.template readExternalVariableIndex<IeeestExternalVariables::U>();
        }

        const ScalarT zero{0.0};
        const ScalarT x1 = u;
        const ScalarT x2 = zero;
        const ScalarT x3 = zero;
        const ScalarT x4 = zero;
        const ScalarT v4 = x1 + A5_ * x2 + A6_ * x3;
        const ScalarT x5 = v4;
        const ScalarT v5 = v4;
        const ScalarT x6 = v5;
        const ScalarT v6 = v5;
        const ScalarT x7 = v6;
        const ScalarT v7 = (T6_ == ZERO<RealT>) ? Ks_ * v6 : zero;

        y_[0]  = x1;
        y_[1]  = x2;
        y_[2]  = x3;
        y_[3]  = x4;
        y_[4]  = x5;
        y_[5]  = x6;
        y_[6]  = x7;
        y_[7]  = v4;
        y_[8]  = v5;
        y_[9]  = v6;
        y_[10] = v7;
        y_[11] = Math::clamp(v7, Lsmin_, Lsmax_);

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Ieeest<scalar_type, index_type>::tagDifferentiable()
      {
        tag_[0]  = (a2_ != ZERO<RealT> || a3_ != ZERO<RealT> || a4_ != ZERO<RealT>);
        tag_[1]  = tag_[0];
        tag_[2]  = (a3_ != ZERO<RealT> || a4_ != ZERO<RealT>);
        tag_[3]  = (a4_ != ZERO<RealT>);
        tag_[4]  = (T2_ != ZERO<RealT>);
        tag_[5]  = (T4_ != ZERO<RealT>);
        tag_[6]  = (T6_ != ZERO<RealT>);
        tag_[7]  = false;
        tag_[8]  = false;
        tag_[9]  = false;
        tag_[10] = false;
        tag_[11] = false;

        return 0;
      }

      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline int Ieeest<scalar_type, index_type>::evaluateInternalResidual(
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

        ScalarT x1_dot = yp[0];
        ScalarT x2_dot = yp[1];
        ScalarT x3_dot = yp[2];
        ScalarT x4_dot = yp[3];
        ScalarT x5_dot = yp[4];
        ScalarT x6_dot = yp[5];
        ScalarT x7_dot = yp[6];

        ScalarT u = ws[0];

        f[0]  = -tag_[0] * x1_dot + x2;
        f[1]  = -tag_[1] * x2_dot + x3;
        f[2]  = -tag_[2] * x3_dot + x4;
        f[3]  = -a4_ * x4_dot - x1 - a1_ * x2 - a2_ * x3 - a3_ * x4 + u;
        f[4]  = -T2_ * x5_dot - x5 + v4;
        f[5]  = -T4_ * x6_dot - x6 + v5;
        f[6]  = -T6_ * x7_dot - x7 + v6;
        f[7]  = -v4 + x1 + A5_ * x2 + A6_ * x3;
        f[8]  = tag_[4] * (-T2_ * (v5 - x5) + T1_ * (v4 - x5)) + (1 - tag_[4]) * (v4 - v5);
        f[9]  = tag_[5] * (-T4_ * (v6 - x6) + T3_ * (v5 - x6)) + (1 - tag_[5]) * (v5 - v6);
        f[10] = tag_[6] * (-T6_ * v7 + Ks_ * T5_ * (v6 - x7)) + (1 - tag_[6]) * (Ks_ * v6 - v7);
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

        evaluateInternalResidual(y_.data(), yp_.data(), wb_.data(), ws_.data(), f_.data());

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
                      { return y_[11]; });
      }

    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
