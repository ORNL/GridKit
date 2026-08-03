#pragma once

/**
 * @file IeeestImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the IEEEST Power System Stabilizer.
 */

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

      /**
       * @brief Construct a sized but unconfigured IEEEST stabilizer.
       */
      template <typename scalar_type, typename index_type>
      Ieeest<scalar_type, index_type>::Ieeest()
      {
        size_ = static_cast<IdxT>(IeeestInternalVariables::MAXIMUM);
      }

      /**
       * @brief Construct an IEEEST stabilizer from model data.
       *
       * @param[in] data Parameters and monitored-variable selections.
       */
      template <typename scalar_type, typename index_type>
      Ieeest<scalar_type, index_type>::Ieeest(const ModelDataT& data)
        : monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        initializeMonitor();
        size_ = static_cast<IdxT>(IeeestInternalVariables::MAXIMUM);
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

        setDerivedParameters();
      }

      template <typename scalar_type, typename index_type>
      void Ieeest<scalar_type, index_type>::setDerivedParameters()
      {
        a0_ = 1;
        a1_ = A1_ + A3_;
        a2_ = A2_ + A4_ + A1_ * A3_;
        a3_ = A1_ * A4_ + A2_ * A3_;
        a4_ = A2_ * A4_;

        O1_ = static_cast<RealT>(a1_ + a2_ + a3_ + a4_ > 0.0);
        O2_ = static_cast<RealT>(a2_ + a3_ + a4_ > 0.0);
        O3_ = static_cast<RealT>(a3_ + a4_ > 0.0);

        // A denominator time constant that is not positive bypasses its block.
        bypass_T2_ = static_cast<RealT>(T2_ <= 0.0);
        bypass_T4_ = static_cast<RealT>(T4_ <= 0.0);
        bypass_T6_ = static_cast<RealT>(T6_ <= 0.0);
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
        if (!allocated_)
        {
          this->allocateVectors(size_);
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

        allocated_ = true;
        return 0;
      }

      /**
       * @brief Validate the IEEEST configuration.
       *
       * Checks that the required input signal has a linked source and that
       * the filter constants are valid
       *
       * @return Number of configuration errors, zero when valid.
       */
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

        if (A1_ < 0.0 || A2_ < 0.0 || A3_ < 0.0 || A4_ < 0.0)
        {
          Log::error() << "Ieeest: filter constants A1-A4 must be nonnegative\n";
          ret += 1;
        }

        if (T2_ < 0.0 || T4_ < 0.0 || T6_ < 0.0)
        {
          Log::error() << "Ieeest: time constants T2, T4, and T6 must be nonnegative\n";
          ret += 1;
        }

        return ret;
      }

      /**
       * @brief Initialize IEEEST to steady state at the input operating point.
       *
       * The notch filter and both lead-lag blocks have unity DC gain, so every
       * block output settles at the input \f$u_0\f$ while the washout removes
       * the DC component. The initial residual therefore vanishes for any
       * \f$u_0\f$ and the stabilizer comes online at rest.
       *
       * @pre allocate() has completed and the component driving the input
       *      signal has been initialized.
       * @post All internal derivatives are zero and the initial residual is
       *       zero up to the smoothing error of the output limiter.
       */
      template <typename scalar_type, typename index_type>
      int Ieeest<scalar_type, index_type>::initialize()
      {
        const auto X1  = static_cast<size_t>(IeeestInternalVariables::X1);
        const auto X5  = static_cast<size_t>(IeeestInternalVariables::X5);
        const auto X6  = static_cast<size_t>(IeeestInternalVariables::X6);
        const auto X7  = static_cast<size_t>(IeeestInternalVariables::X7);
        const auto V4  = static_cast<size_t>(IeeestInternalVariables::V4);
        const auto V5  = static_cast<size_t>(IeeestInternalVariables::V5);
        const auto V6  = static_cast<size_t>(IeeestInternalVariables::V6);
        const auto V7  = static_cast<size_t>(IeeestInternalVariables::V7);
        const auto VSS = static_cast<size_t>(IeeestInternalVariables::VSS);

        auto* y  = y_.getData();
        auto* yp = yp_.getData();

        for (IdxT i = 0; i < size_; ++i)
        {
          y[static_cast<size_t>(i)]  = 0.0;
          yp[static_cast<size_t>(i)] = 0.0;
        }

        ScalarT u{0.0};
        if (signals_.template isAttached<IeeestExternalVariables::U>())
        {
          u              = signals_.template readExternalVariable<IeeestExternalVariables::U>();
          ws_[0]         = u;
          ws_indices_[0] = signals_.template readExternalVariableIndex<IeeestExternalVariables::U>();
        }

        // The chain states above x1 hold derivatives of the filtered signal and
        // stay at zero.
        y[X1] = u;
        y[V4] = u;
        y[X5] = u;
        y[V5] = u;
        y[X6] = u;
        y[V6] = u;
        y[X7] = u;

        // A bypassed washout is a plain gain and passes the input through.
        y[V7]  = bypass_T6_ * Ks_ * u;
        y[VSS] = Math::clamp(y[V7], Lsmin_, Lsmax_);

        y_.setDataUpdated();
        yp_.setDataUpdated();

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Ieeest<scalar_type, index_type>::tagDifferentiable()
      {
        tag_[0]  = (O1_ > 0.0);
        tag_[1]  = (O2_ > 0.0);
        tag_[2]  = (O3_ > 0.0);
        tag_[3]  = (a4_ > 0.0);
        tag_[4]  = (T2_ > 0.0);
        tag_[5]  = (T4_ > 0.0);
        tag_[6]  = (T6_ > 0.0);
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

        // Chain rows up to the notch order, then one accumulated sum row that
        // defines the rate of the highest active state. An inactive chain row
        // loses its derivative and pins the next state instead.
        f[0]  = -O1_ * x1_dot + x2;
        f[1]  = -O2_ * x2_dot + x3;
        f[2]  = -O3_ * x3_dot + x4;
        f[3]  = -a0_ * x1 - a1_ * x2 - a2_ * x3 - a3_ * x4 - a4_ * x4_dot + u;
        f[4]  = -T2_ * x5_dot - x5 + v4;
        f[5]  = -T4_ * x6_dot - x6 + v5;
        f[6]  = -T6_ * x7_dot - x7 + v6;
        f[7]  = -v4 + x1 + A5_ * x2 + A6_ * x3;
        f[8]  = -T2_ * (v5 - x5) + T1_ * (v4 - x5) + bypass_T2_ * (v4 - v5);
        f[9]  = -T4_ * (v6 - x6) + T3_ * (v5 - x6) + bypass_T4_ * (v5 - v6);
        f[10] = -T6_ * v7 + Ks_ * T5_ * (v6 - x7) + bypass_T6_ * (Ks_ * v6 - v7);
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
