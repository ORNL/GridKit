#pragma once

/**
 * @file IeeestImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the IEEEST Power System Stabilizer.
 */

#include <iostream>

#include <GridKit/Model/ConfigurationChecks.hpp>
#include <GridKit/Model/ParameterReader.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeSet.hpp>
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
        setDerivedParameters();
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

      /**
       * @brief Read the parameters out of the model data
       *
       * Every parameter is optional and keeps its documented default when
       * omitted. A non-numeric or non-finite value throws
       * std::invalid_argument. Integer JSON values are accepted for real
       * parameters.
       *
       * @param[in] data Parameters and monitored-variable selections.
       */
      template <typename scalar_type, typename index_type>
      void Ieeest<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Params = typename ModelDataT::Parameters;

        const Model::ParameterReader reader(data, "Ieeest");

        A1_     = reader.readReal(Params::A1).value_or(A1_);
        A2_     = reader.readReal(Params::A2).value_or(A2_);
        A3_     = reader.readReal(Params::A3).value_or(A3_);
        A4_     = reader.readReal(Params::A4).value_or(A4_);
        A5_     = reader.readReal(Params::A5).value_or(A5_);
        A6_     = reader.readReal(Params::A6).value_or(A6_);
        T1_     = reader.readReal(Params::T1).value_or(T1_);
        T2_     = reader.readReal(Params::T2).value_or(T2_);
        T3_     = reader.readReal(Params::T3).value_or(T3_);
        T4_     = reader.readReal(Params::T4).value_or(T4_);
        T5_     = reader.readReal(Params::T5).value_or(T5_);
        T6_     = reader.readReal(Params::T6).value_or(T6_);
        Ks_     = reader.readReal(Params::Ks).value_or(Ks_);
        Lsmin_  = reader.readReal(Params::Lsmin).value_or(Lsmin_);
        Lsmax_  = reader.readReal(Params::Lsmax).value_or(Lsmax_);
        Vcl_    = reader.readReal(Params::Vcl).value_or(Vcl_);
        Vcu_    = reader.readReal(Params::Vcu).value_or(Vcu_);
        Tdelay_ = reader.readReal(Params::Tdelay).value_or(Tdelay_);
        setDerivedParameters();
      }

      /**
       * @brief Resolve the notch coefficients and block masks
       *
       * Expands the notch denominator and precomputes the order masks, safe
       * inverse coefficients, and lead-lag bypass masks so the residual
       * stays branch-free.
       */
      template <typename scalar_type, typename index_type>
      void Ieeest<scalar_type, index_type>::setDerivedParameters()
      {
        a0_ = ONE<RealT>;
        a1_ = A1_ + A3_;
        a2_ = A2_ + A4_ + A1_ * A3_;
        a3_ = A1_ * A4_ + A2_ * A3_;
        a4_ = A2_ * A4_;

        use_notch_    = static_cast<RealT>(a2_ != ZERO<RealT> || a3_ != ZERO<RealT> || a4_ != ZERO<RealT>);
        bypass_notch_ = ONE<RealT> - use_notch_;

        use_4th_order_ = static_cast<RealT>(a4_ != ZERO<RealT>);
        use_3rd_order_ = static_cast<RealT>(a4_ == ZERO<RealT> && a3_ != ZERO<RealT>);
        use_2nd_order_ = static_cast<RealT>(a4_ == ZERO<RealT> && a3_ == ZERO<RealT> && a2_ != ZERO<RealT>);
        safe_inv_a4_   = use_4th_order_ / (a4_ + (ONE<RealT> - use_4th_order_));
        safe_inv_a3_   = use_3rd_order_ / (a3_ + (ONE<RealT> - use_3rd_order_));
        safe_inv_a2_   = use_2nd_order_ / (a2_ + (ONE<RealT> - use_2nd_order_));

        use_T2_block_    = static_cast<RealT>(T2_ != ZERO<RealT>);
        bypass_T2_block_ = ONE<RealT> - use_T2_block_;

        use_T4_block_    = static_cast<RealT>(T4_ != ZERO<RealT>);
        bypass_T4_block_ = ONE<RealT> - use_T4_block_;

        use_T6_block_    = static_cast<RealT>(T6_ != ZERO<RealT>);
        bypass_T6_block_ = ONE<RealT> - use_T6_block_;
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
        ws_.setToZero();
        ws_indices_.resize(1);
        ws_indices_[0] = INVALID_INDEX<IdxT>;

        if (auto output_port = ports_.out.template port<IeeestSignalOutputs::output>())
        {
          output_port.link(&y_.getData()[11], &(this->getVariableIndex(11)));
        }

        allocated_ = true;
        return 0;
      }

      /**
       * @brief Validate the IEEEST configuration
       *
       * Checks the output limits, the realizable notch order, the
       * unsupported input-cutout and delay features, and the required input
       * signal.
       *
       * @return The configuration checks; passed() when valid.
       */
      template <typename scalar_type, typename index_type>
      Model::ConfigurationChecks Ieeest<scalar_type, index_type>::verify() const
      {
        Model::ConfigurationChecks checks("Ieeest");

        checks.check(Lsmin_ <= Lsmax_, "Lsmin must be less than or equal to Lsmax");

        // The notch realization covers denominator orders two through four.
        const bool first_order_notch = a1_ != ZERO<RealT> && a2_ == ZERO<RealT>
                                       && a3_ == ZERO<RealT> && a4_ == ZERO<RealT>;
        checks.check(!first_order_notch, "A1 and A3 alone give a first-order notch, which is not supported");

        checks.check(Vcl_ == ZERO<RealT> && Vcu_ == ZERO<RealT>,
                     "Vcl and Vcu must be zero; input cutout is not supported");
        checks.check(Tdelay_ == ZERO<RealT>, "Tdelay must be zero; input delay is not supported");

        const auto input_port = ports_.in.template port<IeeestSignalInputs::input>();
        if (!input_port.connected())
        {
          checks.fail("required input signal is not assigned");
        }
        else if (!input_port.linked())
        {
          checks.fail("input signal attached with no linked source");
        }

        return checks;
      }

      template <typename scalar_type, typename index_type>
      int Ieeest<scalar_type, index_type>::initialize()
      {
        if (!this->verifyAndLog())
        {
          Log::error() << "Ieeest: cannot initialize with invalid configuration\n";
          return 1;
        }

        auto          input_port = ports_.in.template port<IeeestSignalInputs::input>();
        const ScalarT u          = input_port.readSignal();

        auto* y  = y_.getData();
        auto* yp = yp_.getData();

        for (IdxT i = 0; i < size_; ++i)
        {
          y[static_cast<size_t>(i)]  = 0.0;
          yp[static_cast<size_t>(i)] = 0.0;
        }

        ws_.getData()[0] = u;
        ws_indices_[0]   = input_port.signalVariableIndex();

        y[0] = use_notch_ * u;
        y[4] = u;
        y[5] = u;
        y[6] = u;
        y[7] = u;
        y[8] = u;
        y[9] = u;

        // Preserve the current T6 = 0 bypass behavior.
        y[10] = bypass_T6_block_ * Ks_ * u;
        y[11] = Math::clamp(y[10], Lsmin_, Lsmax_);

        // For DependencyTracking::Variable, set variable numbers
        if constexpr (std::is_same_v<scalar_type, DependencyTracking::Variable>)
        {
          this->initializeDependencyTrackingVariableNumbers();
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

        const ScalarT x2_rhs = (use_4th_order_ + use_3rd_order_) * x3
                               + use_2nd_order_ * (-a0_ * x1 - a1_ * x2 + u) * safe_inv_a2_;

        f[0] = -x1_dot + use_notch_ * x2;
        f[1] = -x2_dot + x2_rhs;
        f[2] = -x3_dot + use_4th_order_ * x4
               + use_3rd_order_ * (-a0_ * x1 - a1_ * x2 - a2_ * x3 + u) * safe_inv_a3_;
        f[3]  = -x4_dot + use_4th_order_ * (-a0_ * x1 - a1_ * x2 - a2_ * x3 - a3_ * x4 + u) * safe_inv_a4_;
        f[4]  = -T2_ * x5_dot - x5 + v4;
        f[5]  = -T4_ * x6_dot - x6 + v5;
        f[6]  = -T6_ * x7_dot - x7 + v6;
        f[7]  = -v4 + bypass_notch_ * u + use_notch_ * (x1 + A5_ * x2 + A6_ * x2_rhs);
        f[8]  = use_T2_block_ * (-T2_ * (v5 - x5) + T1_ * (v4 - x5)) + bypass_T2_block_ * (v4 - v5);
        f[9]  = use_T4_block_ * (-T4_ * (v6 - x6) + T3_ * (v5 - x6)) + bypass_T4_block_ * (v5 - v6);
        f[10] = use_T6_block_ * (-T6_ * v7 + Ks_ * T5_ * (v6 - x7)) + bypass_T6_block_ * (Ks_ * v6 - v7);
        f[11] = -vss + Math::clamp(v7, Lsmin_, Lsmax_);

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Ieeest<scalar_type, index_type>::evaluateResidual()
      {
        auto* ws = ws_.getData();

        if (auto input_port = ports_.in.template port<IeeestSignalInputs::input>())
        {
          ws[0]          = input_port.readSignal();
          ws_indices_[0] = input_port.signalVariableIndex();
        }

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
