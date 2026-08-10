#pragma once

/**
 * @file IeeestImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the IEEEST Power System Stabilizer.
 */

#include <algorithm>
#include <cmath>
#include <variant>

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
        size_ = static_cast<IdxT>(IeeestInternalVariables::MAXIMUM);
        setDerivedParameters();
      }

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
        const auto size = static_cast<size_t>(size_);

        tag_.assign(size, false);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        const auto signal_size = static_cast<size_t>(IeeestExternalVariables::MAXIMUM);
        ws_.resize(static_cast<IdxT>(signal_size));
        ws_.setToZero();
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        if (signals_.template isAssigned<IeeestInternalVariables::VSS>())
        {
          auto* y = y_.getData();
          signals_.template getSignalNode<IeeestInternalVariables::VSS>()->set(
              &y[static_cast<size_t>(IeeestInternalVariables::VSS)],
              &(this->getVariableIndex(static_cast<IdxT>(IeeestInternalVariables::VSS))));
        }

        allocated_ = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Ieeest<scalar_type, index_type>::verify() const
      {
        int ret = static_cast<int>(parameter_error_count_);

        auto check = [&](bool condition, const char* message)
        {
          if (!condition)
          {
            Log::error() << "Ieeest: " << message << '\n';
            ret += 1;
          }
        };

        if (!signals_.template isAttached<IeeestExternalVariables::U>())
        {
          Log::error() << "Ieeest: required input signal U is not attached\n";
          ret += 1;
        }
        else if (!signals_.template isLinked<IeeestExternalVariables::U>())
        {
          Log::error() << "Ieeest: input signal U attached with no linked source\n";
          ret += 1;
        }

        check(std::isfinite(a1_) && std::isfinite(a2_)
                  && std::isfinite(a3_) && std::isfinite(a4_),
              "expanded denominator coefficients must be finite");

        IdxT numerator_order = static_cast<IdxT>(0);
        if (A6_ != ZERO<RealT>)
        {
          numerator_order = static_cast<IdxT>(2);
        }
        else if (A5_ != ZERO<RealT>)
        {
          numerator_order = static_cast<IdxT>(1);
        }
        check(numerator_order <= order_,
              "numerator order must not exceed denominator order");
        check(Lsmin_ < Lsmax_, "Lsmin must be less than Lsmax");

        RealT       leading_coefficient = ONE<RealT>;
        const char* leading_name        = nullptr;
        switch (order_)
        {
        case 0:
          break;
        case 1:
          leading_coefficient = a1_;
          leading_name        = "a1";
          break;
        case 2:
          leading_coefficient = a2_;
          leading_name        = "a2";
          break;
        case 3:
          leading_coefficient = a3_;
          leading_name        = "a3";
          break;
        case 4:
          leading_coefficient = a4_;
          leading_name        = "a4";
          break;
        }

        if (leading_name != nullptr && leading_coefficient == ZERO<RealT>)
        {
          Log::error() << "Ieeest: active leading denominator coefficient '"
                       << leading_name << "' must be nonzero\n";
          ret += 1;
        }

        return ret;
      }

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
        const auto U   = static_cast<size_t>(IeeestExternalVariables::U);

        if (verify() != 0)
        {
          Log::error() << "Ieeest: cannot initialize with invalid configuration\n";
          return 1;
        }

        const ScalarT u = signals_.template readExternalVariable<IeeestExternalVariables::U>();
        if (!std::isfinite(static_cast<RealT>(u)))
        {
          Log::error() << "Ieeest: initial input signal U must be finite\n";
          return 1;
        }

        auto* y  = y_.getData();
        auto* yp = yp_.getData();
        auto* ws = ws_.getData();
        std::fill_n(y, static_cast<size_t>(size_), ScalarT{ZERO<RealT>});
        std::fill_n(yp, static_cast<size_t>(size_), ScalarT{ZERO<RealT>});

        ws[U]          = u;
        ws_indices_[U] = signals_.template readExternalVariableIndex<IeeestExternalVariables::U>();

        switch (order_)
        {
        case 0:
          break;
        case 1:
        case 2:
        case 3:
        case 4:
          y[X1] = u;
          break;
        }

        y[X5]  = u;
        y[X6]  = u;
        y[X7]  = u;
        y[V4]  = u;
        y[V5]  = u;
        y[V6]  = u;
        y[V7]  = ZERO<RealT>;
        y[VSS] = Math::clamp(y[V7], Lsmin_, Lsmax_);

        y_.setDataUpdated();
        yp_.setDataUpdated();

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Ieeest<scalar_type, index_type>::tagDifferentiable()
      {
        std::fill(tag_.begin(), tag_.end(), false);
        constexpr auto X7 = static_cast<size_t>(IeeestInternalVariables::X7);
        std::fill_n(tag_.begin(), X7 + 1, true);

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
        const auto X1  = static_cast<size_t>(IeeestInternalVariables::X1);
        const auto X2  = static_cast<size_t>(IeeestInternalVariables::X2);
        const auto X3  = static_cast<size_t>(IeeestInternalVariables::X3);
        const auto X4  = static_cast<size_t>(IeeestInternalVariables::X4);
        const auto X5  = static_cast<size_t>(IeeestInternalVariables::X5);
        const auto X6  = static_cast<size_t>(IeeestInternalVariables::X6);
        const auto X7  = static_cast<size_t>(IeeestInternalVariables::X7);
        const auto V4  = static_cast<size_t>(IeeestInternalVariables::V4);
        const auto V5  = static_cast<size_t>(IeeestInternalVariables::V5);
        const auto V6  = static_cast<size_t>(IeeestInternalVariables::V6);
        const auto V7  = static_cast<size_t>(IeeestInternalVariables::V7);
        const auto VSS = static_cast<size_t>(IeeestInternalVariables::VSS);
        const auto U   = static_cast<size_t>(IeeestExternalVariables::U);

        const ScalarT x1  = y[X1];
        const ScalarT x2  = y[X2];
        const ScalarT x3  = y[X3];
        const ScalarT x4  = y[X4];
        const ScalarT x5  = y[X5];
        const ScalarT x6  = y[X6];
        const ScalarT x7  = y[X7];
        const ScalarT v4  = y[V4];
        const ScalarT v5  = y[V5];
        const ScalarT v6  = y[V6];
        const ScalarT v7  = y[V7];
        const ScalarT vss = y[VSS];

        const ScalarT x1_dot = yp[X1];
        const ScalarT x2_dot = yp[X2];
        const ScalarT x3_dot = yp[X3];
        const ScalarT x4_dot = yp[X4];
        const ScalarT x5_dot = yp[X5];
        const ScalarT x6_dot = yp[X6];
        const ScalarT x7_dot = yp[X7];

        const ScalarT u = ws[U];

        switch (order_)
        {
        case 0:
          f[X1] = -x1_dot;
          f[X2] = -x2_dot;
          f[X3] = -x3_dot;
          f[X4] = -x4_dot;
          f[V4] = -v4 + u;
          break;
        case 1:
        {
          const ScalarT x1_rhs = (u - x1) / a1_;
          f[X1]                = -x1_dot + x1_rhs;
          f[X2]                = -x2_dot;
          f[X3]                = -x3_dot;
          f[X4]                = -x4_dot;
          f[V4]                = -v4 + x1 + A5_ * x1_rhs;
          break;
        }
        case 2:
        {
          const ScalarT x2_rhs = (u - x1 - a1_ * x2) / a2_;
          f[X1]                = -x1_dot + x2;
          f[X2]                = -x2_dot + x2_rhs;
          f[X3]                = -x3_dot;
          f[X4]                = -x4_dot;
          f[V4]                = -v4 + x1 + A5_ * x2 + A6_ * x2_rhs;
          break;
        }
        case 3:
        {
          const ScalarT x3_rhs = (u - x1 - a1_ * x2 - a2_ * x3) / a3_;
          f[X1]                = -x1_dot + x2;
          f[X2]                = -x2_dot + x3;
          f[X3]                = -x3_dot + x3_rhs;
          f[X4]                = -x4_dot;
          f[V4]                = -v4 + x1 + A5_ * x2 + A6_ * x3;
          break;
        }
        case 4:
        {
          const ScalarT x4_rhs =
              (u - x1 - a1_ * x2 - a2_ * x3 - a3_ * x4) / a4_;
          f[X1] = -x1_dot + x2;
          f[X2] = -x2_dot + x3;
          f[X3] = -x3_dot + x4;
          f[X4] = -x4_dot + x4_rhs;
          f[V4] = -v4 + x1 + A5_ * x2 + A6_ * x3;
          break;
        }
        }

        const ScalarT x5_rhs = (v4 - x5) / T2_;
        const ScalarT x6_rhs = (v5 - x6) / T4_;
        const ScalarT x7_rhs = (v6 - x7) / T6_;

        f[X5]  = -x5_dot + x5_rhs;
        f[X6]  = -x6_dot + x6_rhs;
        f[X7]  = -x7_dot + x7_rhs;
        f[V5]  = -v5 + x5 + T1_ * x5_rhs;
        f[V6]  = -v6 + x6 + T3_ * x6_rhs;
        f[V7]  = -v7 + Ks_ * T5_ * x7_rhs;
        f[VSS] = -vss + Math::clamp(v7, Lsmin_, Lsmax_);

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Ieeest<scalar_type, index_type>::evaluateResidual()
      {
        auto* ws = ws_.getData();
        const auto U = static_cast<size_t>(IeeestExternalVariables::U);

        ws[U]          = signals_.template readExternalVariable<IeeestExternalVariables::U>();
        ws_indices_[U] = signals_.template readExternalVariableIndex<IeeestExternalVariables::U>();

        const auto* y  = y_.getData();
        const auto* yp = yp_.getData();
        auto*       f  = f_.getData();
        evaluateInternalResidual(y, yp, nullptr, ws, f);
        f_.setDataUpdated();

        return 0;
      }

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* Ieeest<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      //
      //  Private methods
      //

      template <typename scalar_type, typename index_type>
      void Ieeest<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Params = typename ModelDataT::Parameters;

        parameter_error_count_ = 0;

        auto load_real = [&](auto key, RealT& target, const char* name) -> bool
        {
          if (!data.parameters.contains(key))
          {
            return false;
          }

          const auto& value = data.parameters.at(key);
          RealT       parsed_value{};
          if (const auto* real_value = std::get_if<RealT>(&value))
          {
            parsed_value = *real_value;
          }
          else if (const auto* index_value = std::get_if<IdxT>(&value))
          {
            parsed_value = static_cast<RealT>(*index_value);
          }
          else
          {
            Log::error() << "Ieeest: parameter '" << name << "' must be numeric\n";
            ++parameter_error_count_;
            return false;
          }

          if (!std::isfinite(parsed_value))
          {
            Log::error() << "Ieeest: parameter '" << name << "' must be finite\n";
            ++parameter_error_count_;
            return false;
          }

          target = parsed_value;
          return true;
        };

        load_real(Params::A1, A1_, "A1");
        load_real(Params::A2, A2_, "A2");
        load_real(Params::A3, A3_, "A3");
        load_real(Params::A4, A4_, "A4");
        load_real(Params::A5, A5_, "A5");
        load_real(Params::A6, A6_, "A6");
        load_real(Params::T1, T1_, "T1");
        load_real(Params::T2, T2_, "T2");
        load_real(Params::T3, T3_, "T3");
        load_real(Params::T4, T4_, "T4");
        load_real(Params::T5, T5_, "T5");
        load_real(Params::T6, T6_, "T6");
        load_real(Params::Ks, Ks_, "Ks");
        load_real(Params::Lsmin, Lsmin_, "Lsmin");
        load_real(Params::Lsmax, Lsmax_, "Lsmax");
        if (load_real(Params::Vcl, Vcl_, "Vcl") && Vcl_ != ZERO<RealT>)
        {
          Log::warning() << "Ieeest: nonzero Vcl requests a lower input cutout, "
                            "but input cutout is not implemented and Vcl is ignored\n";
        }
        if (load_real(Params::Vcu, Vcu_, "Vcu") && Vcu_ != ZERO<RealT>)
        {
          Log::warning() << "Ieeest: nonzero Vcu requests an upper input cutout, "
                            "but input cutout is not implemented and Vcu is ignored\n";
        }
        if (load_real(Params::Tdelay, Tdelay_, "Tdelay") && Tdelay_ != ZERO<RealT>)
        {
          Log::warning() << "Ieeest: nonzero Tdelay requests an input delay, "
                            "but input delay is not implemented and Tdelay is ignored\n";
        }

        setDerivedParameters();
      }

      template <typename scalar_type, typename index_type>
      void Ieeest<scalar_type, index_type>::initializeMonitor()
      {
        using Variable = typename ModelDataT::MonitorableVariables;

        constexpr auto VSS = static_cast<size_t>(IeeestInternalVariables::VSS);

        monitor_->set(Variable::vss, [this]
                      { return y_.getData()[VSS]; });
      }

      template <typename scalar_type, typename index_type>
      void Ieeest<scalar_type, index_type>::setDerivedParameters()
      {
        auto check_non_negative = [&](RealT value, const char* name)
        {
          if (value < ZERO<RealT>)
          {
            Log::error() << "Ieeest: " << name << " must be non-negative\n";
            ++parameter_error_count_;
          }
        };

        check_non_negative(T2_, "T2");
        check_non_negative(T4_, "T4");
        check_non_negative(T6_, "T6");

        if (T2_ < TIME_CONSTANT_MINIMUM
            || T4_ < TIME_CONSTANT_MINIMUM
            || T6_ < TIME_CONSTANT_MINIMUM)
        {
          Log::warning() << "Ieeest: T2, T4, and T6 below "
                         << TIME_CONSTANT_MINIMUM
                         << " s are raised to preserve Hessenberg form\n";
        }

        T2_ = std::max(T2_, TIME_CONSTANT_MINIMUM);
        T4_ = std::max(T4_, TIME_CONSTANT_MINIMUM);
        T6_ = std::max(T6_, TIME_CONSTANT_MINIMUM);

        a1_ = A1_ + A3_;
        a2_ = A2_ + A4_ + A1_ * A3_;
        a3_ = A1_ * A4_ + A2_ * A3_;
        a4_ = A2_ * A4_;

        auto factor_order = [](RealT linear, RealT quadratic) -> IdxT
        {
          if (quadratic != ZERO<RealT>)
          {
            return static_cast<IdxT>(2);
          }
          if (linear != ZERO<RealT>)
          {
            return static_cast<IdxT>(1);
          }
          return static_cast<IdxT>(0);
        };

        order_ = factor_order(A1_, A2_) + factor_order(A3_, A4_);
      }

    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
