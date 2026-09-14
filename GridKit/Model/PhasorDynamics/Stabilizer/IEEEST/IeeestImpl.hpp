#pragma once

/**
 * @file IeeestImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the IEEEST Power System Stabilizer.
 */

#include <algorithm>
#include <cmath>
#include <mutex>
#include <variant>

#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/Ieeest.hpp>
#include <GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/IeeestData.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Enum.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Stabilizer
    {
      using Log = ::GridKit::Utilities::Logger;

      template <typename scalar_type, typename index_type, size_t order>
      Ieeest<scalar_type, index_type, order>::Ieeest()
      {
        size_ = static_cast<IdxT>(Utilities::enum_size<InternalVariablesT>());
        setDerivedParameters();
      }

      template <typename scalar_type, typename index_type, size_t order>
      Ieeest<scalar_type, index_type, order>::Ieeest(const ModelDataT& data)
        : monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        initializeMonitor();
        size_ = static_cast<IdxT>(Utilities::enum_size<InternalVariablesT>());
      }

      template <typename scalar_type, typename index_type, size_t order>
      Ieeest<scalar_type, index_type, order>::~Ieeest()
      {
      }

      template <typename scalar_type, typename index_type, size_t order>
      void Ieeest<scalar_type, index_type, order>::initializeParameters(const ModelDataT& data)
      {
        using Params = typename ModelDataT::Parameters;

        parameters_valid_ = true;

        loadRealParameter(data, Params::A1, A1_);
        loadRealParameter(data, Params::A2, A2_);
        loadRealParameter(data, Params::A3, A3_);
        loadRealParameter(data, Params::A4, A4_);
        loadRealParameter(data, Params::A5, A5_);
        loadRealParameter(data, Params::A6, A6_);
        loadRealParameter(data, Params::T1, T1_);
        loadRealParameter(data, Params::T2, T2_);
        loadRealParameter(data, Params::T3, T3_);
        loadRealParameter(data, Params::T4, T4_);
        loadRealParameter(data, Params::T5, T5_);
        loadRealParameter(data, Params::T6, T6_);
        loadRealParameter(data, Params::Ks, Ks_);
        loadRealParameter(data, Params::Lsmin, Lsmin_);
        loadRealParameter(data, Params::Lsmax, Lsmax_);
        loadRealParameter(data, Params::Vcl, Vcl_);
        loadRealParameter(data, Params::Vcu, Vcu_);
        loadRealParameter(data, Params::Tdelay, Tdelay_);

        if (Vcl_ != ZERO<RealT>)
        {
          Log::warning() << "Ieeest: nonzero Vcl requests lower input cutout, which is not implemented\n";
        }
        if (Vcu_ != ZERO<RealT>)
        {
          Log::warning() << "Ieeest: nonzero Vcu requests upper input cutout, which is not implemented\n";
        }
        if (Tdelay_ != ZERO<RealT>)
        {
          Log::warning() << "Ieeest: nonzero Tdelay requests input delay, which is not implemented\n";
        }

        setDerivedParameters();
      }

      template <typename scalar_type, typename index_type, size_t order>
      int Ieeest<scalar_type, index_type, order>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      template <typename scalar_type, typename index_type, size_t order>
      int Ieeest<scalar_type, index_type, order>::allocate()
      {
        if (!allocated_)
        {
          this->allocateVectors(size_);
        }
        const auto size = static_cast<size_t>(size_);

        tag_.resize(size);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        const auto signal_size = Utilities::enum_size<IeeestExternalVariables>();
        ws_.resize(static_cast<IdxT>(signal_size));
        ws_.setToZero();
        ws_indices_.resize(signal_size);
        ws_indices_[U] = INVALID_INDEX<IdxT>;

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        if (auto output_port = ports_.out.template port<IeeestSignalOutputs::output>())
        {
          auto* y = y_.getData();
          output_port.link(
              &y[VSS],
              &(this->getVariableIndex(static_cast<IdxT>(VSS))));
        }

        allocated_ = true;
        return 0;
      }

      template <typename scalar_type, typename index_type, size_t order>
      int Ieeest<scalar_type, index_type, order>::verify() const
      {
        if (!parameters_valid_)
        {
          return 1;
        }

        const auto input = ports_.in.template port<IeeestSignalInputs::input>();
        if (!input.connected())
        {
          Log::error() << "Ieeest: required input signal U is not attached\n";
          return 1;
        }
        if (!input.linked())
        {
          Log::error() << "Ieeest: input signal U attached with no linked source\n";
          return 1;
        }
        if (!std::isfinite(a1_) || !std::isfinite(a2_)
            || !std::isfinite(a3_) || !std::isfinite(a4_))
        {
          Log::error() << "Ieeest: expanded denominator coefficients must be finite\n";
          return 1;
        }
        if (Lsmin_ >= Lsmax_)
        {
          Log::error() << "Ieeest: Lsmin must be less than Lsmax\n";
          return 1;
        }

        if (ieeestNotchOrder(A1_, A2_, A3_, A4_) != order)
        {
          Log::error() << "Ieeest: notch coefficients do not match template order " << order << "\n";
          return 1;
        }
        if constexpr (order == 0)
        {
          if (A5_ != ZERO<RealT> || A6_ != ZERO<RealT>)
          {
            Log::error() << "Ieeest: order zero requires A5 and A6 to be zero\n";
            return 1;
          }
        }
        else if constexpr (order == 1)
        {
          if (A6_ != ZERO<RealT>)
          {
            Log::error() << "Ieeest: order one requires A6 to be zero\n";
            return 1;
          }
        }

        const RealT a[] = {ONE<RealT>, a1_, a2_, a3_, a4_};
        if (a[order] == ZERO<RealT>)
        {
          Log::error() << "Ieeest: active leading coefficient a" << order
                       << " must be nonzero\n";
          return 1;
        }
        if (!std::isfinite(inv_an_))
        {
          Log::error() << "Ieeest: reciprocal of active leading coefficient a" << order
                       << " must be finite\n";
          return 1;
        }
        for (size_t i = 1; i < order; ++i)
        {
          if (!std::isfinite(a[i] * inv_an_))
          {
            Log::error() << "Ieeest: normalized denominator coefficient a" << i
                         << "/a" << order << " must be finite\n";
            return 1;
          }
        }
        if constexpr (order == 1)
        {
          if (!std::isfinite(A5_ * inv_an_))
          {
            Log::error() << "Ieeest: notch feedthrough coefficient A5/a1 must be finite\n";
            return 1;
          }
        }
        else if constexpr (order == 2)
        {
          if (!std::isfinite(A6_ * inv_an_)
              || !std::isfinite(A5_ - A6_ * (a1_ * inv_an_)))
          {
            Log::error() << "Ieeest: notch feedthrough coefficients must be finite\n";
            return 1;
          }
        }
        if (!std::isfinite(T1_ * inv_T2_) || !std::isfinite(T3_ * inv_T4_)
            || !std::isfinite(Ks_ * T5_) || !std::isfinite(Ks_ * T5_ * inv_T6_))
        {
          Log::error() << "Ieeest: lead-lag and washout coefficients must be finite\n";
          return 1;
        }
        return 0;
      }

      template <typename scalar_type, typename index_type, size_t order>
      int Ieeest<scalar_type, index_type, order>::initialize()
      {
        if (verify() != 0)
        {
          Log::error() << "Ieeest: cannot initialize with invalid configuration\n";
          return 1;
        }

        const ScalarT u = ports_.in.template port<IeeestSignalInputs::input>().readSignal();
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
        ws_indices_[U] = ports_.in.template port<IeeestSignalInputs::input>().signalVariableIndex();

        if constexpr (order >= 1)
        {
          y[static_cast<size_t>(InternalVariablesT::X1)] = u;
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

      template <typename scalar_type, typename index_type, size_t order>
      int Ieeest<scalar_type, index_type, order>::tagDifferentiable()
      {
        if constexpr (order >= 1)
        {
          tag_[static_cast<size_t>(InternalVariablesT::X1)] = true;
        }
        if constexpr (order >= 2)
        {
          tag_[static_cast<size_t>(InternalVariablesT::X2)] = true;
        }
        if constexpr (order >= 3)
        {
          tag_[static_cast<size_t>(InternalVariablesT::X3)] = true;
        }
        if constexpr (order == 4)
        {
          tag_[static_cast<size_t>(InternalVariablesT::X4)] = true;
        }
        tag_[X5]  = true;
        tag_[X6]  = true;
        tag_[X7]  = true;
        tag_[V4]  = false;
        tag_[V5]  = false;
        tag_[V6]  = false;
        tag_[V7]  = false;
        tag_[VSS] = false;

        return 0;
      }

      /**
       * @brief Compute the absolute tolerance for each variable in the model
       *
       * @param rel_tol The relative tolerance which can be used to pick the
       *        absolute tolerance.
       * @tparam scalar_type Scalar data type
       * @tparam index_type Index data type
       * @tparam order Notch-denominator degree.
       * @return int 0 if successful, non-zero otherwise.
       *
       * This represents a "noise" level close to zero for which pure relative
       * error cannot be used.
       */
      template <typename scalar_type, typename index_type, size_t order>
      int Ieeest<scalar_type, index_type, order>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      template <typename scalar_type, typename index_type, size_t order>
      [[gnu::always_inline]] inline int Ieeest<scalar_type, index_type, order>::evaluateInternalResidual(
          const ScalarT*                  y,
          const ScalarT*                  yp,
          [[maybe_unused]] const ScalarT* wb,
          const ScalarT*                  ws,
          ScalarT*                        f)
      {

        const ScalarT x5  = y[X5];
        const ScalarT x6  = y[X6];
        const ScalarT x7  = y[X7];
        const ScalarT v4  = y[V4];
        const ScalarT v5  = y[V5];
        const ScalarT v6  = y[V6];
        const ScalarT v7  = y[V7];
        const ScalarT vss = y[VSS];

        const ScalarT x5_dot = yp[X5];
        const ScalarT x6_dot = yp[X6];
        const ScalarT x7_dot = yp[X7];

        const ScalarT u = ws[U];

        // Notch filter -- order-specific realization
        if constexpr (order == 0)
        {
          f[V4] = -v4 + u;
        }
        else if constexpr (order == 1)
        {
          constexpr auto X1 = static_cast<size_t>(InternalVariablesT::X1);

          const ScalarT x1     = y[X1];
          const ScalarT x1_dot = yp[X1];

          const ScalarT x1_rhs = (u - x1) * inv_an_;

          f[X1] = -x1_dot + x1_rhs;
          f[V4] = -v4 + x1 + A5_ * x1_rhs;
        }
        else if constexpr (order == 2)
        {
          constexpr auto X1 = static_cast<size_t>(InternalVariablesT::X1);
          constexpr auto X2 = static_cast<size_t>(InternalVariablesT::X2);

          const ScalarT x1     = y[X1];
          const ScalarT x2     = y[X2];
          const ScalarT x1_dot = yp[X1];
          const ScalarT x2_dot = yp[X2];

          const ScalarT x2_rhs = (u - x1 - a1_ * x2) * inv_an_;

          f[X1] = -x1_dot + x2;
          f[X2] = -x2_dot + x2_rhs;
          f[V4] = -v4 + x1 + A5_ * x2 + A6_ * x2_rhs;
        }
        else if constexpr (order == 3)
        {
          constexpr auto X1 = static_cast<size_t>(InternalVariablesT::X1);
          constexpr auto X2 = static_cast<size_t>(InternalVariablesT::X2);
          constexpr auto X3 = static_cast<size_t>(InternalVariablesT::X3);

          const ScalarT x1     = y[X1];
          const ScalarT x2     = y[X2];
          const ScalarT x3     = y[X3];
          const ScalarT x1_dot = yp[X1];
          const ScalarT x2_dot = yp[X2];
          const ScalarT x3_dot = yp[X3];

          const ScalarT x3_rhs = (u - x1 - a1_ * x2 - a2_ * x3) * inv_an_;

          f[X1] = -x1_dot + x2;
          f[X2] = -x2_dot + x3;
          f[X3] = -x3_dot + x3_rhs;
          f[V4] = -v4 + x1 + A5_ * x2 + A6_ * x3;
        }
        else
        {
          constexpr auto X1 = static_cast<size_t>(InternalVariablesT::X1);
          constexpr auto X2 = static_cast<size_t>(InternalVariablesT::X2);
          constexpr auto X3 = static_cast<size_t>(InternalVariablesT::X3);
          constexpr auto X4 = static_cast<size_t>(InternalVariablesT::X4);

          const ScalarT x1     = y[X1];
          const ScalarT x2     = y[X2];
          const ScalarT x3     = y[X3];
          const ScalarT x4     = y[X4];
          const ScalarT x1_dot = yp[X1];
          const ScalarT x2_dot = yp[X2];
          const ScalarT x3_dot = yp[X3];
          const ScalarT x4_dot = yp[X4];

          const ScalarT x4_rhs = (u - x1 - a1_ * x2 - a2_ * x3 - a3_ * x4) * inv_an_;

          f[X1] = -x1_dot + x2;
          f[X2] = -x2_dot + x3;
          f[X3] = -x3_dot + x4;
          f[X4] = -x4_dot + x4_rhs;
          f[V4] = -v4 + x1 + A5_ * x2 + A6_ * x3;
        }

        // Lead-lags and washout -- shared across all orders
        const ScalarT x5_rhs = (v4 - x5) * inv_T2_;
        const ScalarT x6_rhs = (v5 - x6) * inv_T4_;
        const ScalarT x7_rhs = (v6 - x7) * inv_T6_;

        f[X5]  = -x5_dot + x5_rhs;
        f[X6]  = -x6_dot + x6_rhs;
        f[X7]  = -x7_dot + x7_rhs;
        f[V5]  = -v5 + x5 + T1_ * x5_rhs;
        f[V6]  = -v6 + x6 + T3_ * x6_rhs;
        f[V7]  = -v7 + Ks_ * T5_ * x7_rhs;
        f[VSS] = -vss + Math::clamp(v7, Lsmin_, Lsmax_);

        return 0;
      }

      template <typename scalar_type, typename index_type, size_t order>
      int Ieeest<scalar_type, index_type, order>::evaluateResidual()
      {
        auto* ws = ws_.getData();

        ws[U]          = ports_.in.template port<IeeestSignalInputs::input>().readSignal();
        ws_indices_[U] = ports_.in.template port<IeeestSignalInputs::input>().signalVariableIndex();

        const auto* y  = y_.getData();
        const auto* yp = yp_.getData();
        auto*       f  = f_.getData();
        evaluateInternalResidual(y, yp, nullptr, ws, f);
        f_.setDataUpdated();

        return 0;
      }

      template <typename scalar_type, typename index_type, size_t order>
      const Model::VariableMonitorBase* Ieeest<scalar_type, index_type, order>::getMonitor() const
      {
        return monitor_.get();
      }

      //
      //  Private methods
      //

      /// Load a finite numeric parameter, retaining the default when omitted.
      template <typename scalar_type, typename index_type, size_t order>
      void Ieeest<scalar_type, index_type, order>::loadRealParameter(
          const ModelDataT& data, IeeestParameters parameter, RealT& value)
      {
        const auto entry = data.parameters.find(parameter);
        if (entry == data.parameters.end())
        {
          return;
        }

        RealT parsed{};
        if (const auto* real = std::get_if<RealT>(&entry->second))
        {
          parsed = *real;
        }
        else if (const auto* integer = std::get_if<IdxT>(&entry->second))
        {
          parsed = static_cast<RealT>(*integer);
        }
        else
        {
          Log::error() << "Ieeest: parameter '" << magic_enum::enum_name(parameter)
                       << "' must be numeric\n";
          parameters_valid_ = false;
          return;
        }
        if (!std::isfinite(parsed))
        {
          Log::error() << "Ieeest: parameter '" << magic_enum::enum_name(parameter)
                       << "' must be finite\n";
          parameters_valid_ = false;
          return;
        }
        value = parsed;
      }

      template <typename scalar_type, typename index_type, size_t order>
      void Ieeest<scalar_type, index_type, order>::initializeMonitor()
      {
        using Variable = typename ModelDataT::MonitorableVariables;

        monitor_->set(Variable::vss, [this]
                      { return y_.getData()[VSS]; });
      }

      template <typename scalar_type, typename index_type, size_t order>
      void Ieeest<scalar_type, index_type, order>::setDerivedParameters()
      {
        if (T2_ < ZERO<RealT> || T4_ < ZERO<RealT> || T6_ < ZERO<RealT>)
        {
          Log::error() << "Ieeest: T2, T4, and T6 must be non-negative\n";
          parameters_valid_ = false;
          return;
        }

        if (T2_ < TIME_CONSTANT_MINIMUM
            || T4_ < TIME_CONSTANT_MINIMUM
            || T6_ < TIME_CONSTANT_MINIMUM)
        {
          static std::once_flag time_constant_warning_flag;
          std::call_once(time_constant_warning_flag, []
                         { Log::warning() << "Ieeest: T2, T4, and T6 below "
                                          << TIME_CONSTANT_MINIMUM
                                          << " s are raised to preserve Hessenberg form\n"; });
        }

        T2_ = std::max(T2_, TIME_CONSTANT_MINIMUM);
        T4_ = std::max(T4_, TIME_CONSTANT_MINIMUM);
        T6_ = std::max(T6_, TIME_CONSTANT_MINIMUM);

        a1_ = A1_ + A3_;
        a2_ = A2_ + A4_ + A1_ * A3_;
        a3_ = A1_ * A4_ + A2_ * A3_;
        a4_ = A2_ * A4_;

        // Keep parameter reciprocals outside the Enzyme kernel: differentiated
        // quotients can overflow even when the normalized coefficients are finite.
        const RealT a[] = {ONE<RealT>, a1_, a2_, a3_, a4_};
        if (a[order] != ZERO<RealT>)
        {
          inv_an_ = ONE<RealT> / a[order];
        }
        inv_T2_ = ONE<RealT> / T2_;
        inv_T4_ = ONE<RealT> / T4_;
        inv_T6_ = ONE<RealT> / T6_;
      }

    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
