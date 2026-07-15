/**
 * @file HygovImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the HYGOV governor model.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <variant>

#include <GridKit/Model/PhasorDynamics/Governor/HYGOV/Hygov.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/HYGOV/HygovData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      using Log = ::GridKit::Utilities::Logger;

      template <typename scalar_type, typename index_type>
      Hygov<scalar_type, index_type>::Hygov()
      {
        size_ = static_cast<IdxT>(HygovInternalVariables::MAXIMUM);
      }

      template <typename scalar_type, typename index_type>
      Hygov<scalar_type, index_type>::Hygov(const ModelDataT& data)
        : monitor_(std::make_unique<MonitorT>(data))
      {
        initModelParams(data);
        initializeMonitor();

        size_ = static_cast<IdxT>(HygovInternalVariables::MAXIMUM);
      }

      template <typename scalar_type, typename index_type>
      Hygov<scalar_type, index_type>::~Hygov()
      {
      }

      template <typename scalar_type, typename index_type>
      void Hygov<scalar_type, index_type>::setDerivedParameters()
      {
        va_component_base_ = Trate_ * static_cast<RealT>(1.0e6);

        Tr_  = std::max(Tr_, TIME_CONSTANT_MINIMUM);
        Tf_  = std::max(Tf_, TIME_CONSTANT_MINIMUM);
        Tg_  = std::max(Tg_, TIME_CONSTANT_MINIMUM);
        Tw_  = std::max(Tw_, TIME_CONSTANT_MINIMUM);
        Tnp_ = std::max(Tnp_, TIME_CONSTANT_MINIMUM);

        leadlag_gain_ = Tn_ / Tnp_;
      }

      template <typename scalar_type, typename index_type>
      scalar_type Hygov<scalar_type, index_type>::toComponentBase(scalar_type value) const
      {
        return value * va_system_base_ / va_component_base_;
      }

      template <typename scalar_type, typename index_type>
      scalar_type Hygov<scalar_type, index_type>::toSystemBase(scalar_type value) const
      {
        return value / toComponentBase(static_cast<ScalarT>(ONE<RealT>));
      }

      template <typename scalar_type, typename index_type>
      void Hygov<scalar_type, index_type>::initModelParams(const ModelDataT& data)
      {
        using Params = typename ModelDataT::Parameters;

        parameter_error_count_ = 0;

        Trate_ = ZERO<RealT>;
        Rperm_ = static_cast<RealT>(0.04);
        Rtemp_ = static_cast<RealT>(0.3);
        Tr_    = static_cast<RealT>(5.0);
        Tf_    = static_cast<RealT>(0.05);
        Tg_    = static_cast<RealT>(0.5);
        Velm_  = static_cast<RealT>(0.2);
        Gmax_  = ONE<RealT>;
        Gmin_  = ZERO<RealT>;
        Tw_    = ONE<RealT>;
        At_    = static_cast<RealT>(1.2);
        Dturb_ = static_cast<RealT>(0.5);
        Qnl_   = static_cast<RealT>(0.05);
        Tn_    = ZERO<RealT>;
        Tnp_   = ZERO<RealT>;
        db1_   = ZERO<RealT>;
        db2_   = ZERO<RealT>;
        Hdam_  = ONE<RealT>;
        Gv_.fill(ZERO<RealT>);
        Pgv_.fill(ZERO<RealT>);

        auto load_real = [&](auto key, RealT& target, const char* name)
        {
          if (!data.parameters.contains(key))
          {
            return;
          }

          const auto& value = data.parameters.at(key);
          if (const auto* real_value = std::get_if<RealT>(&value))
          {
            target = *real_value;
          }
          else if (const auto* index_value = std::get_if<IdxT>(&value))
          {
            target = static_cast<RealT>(*index_value);
          }
          else
          {
            Log::error() << "Hygov: parameter '" << name << "' must be numeric\n";
            ++parameter_error_count_;
          }
        };

        auto load_required_real = [&](auto key, RealT& target, const char* name)
        {
          if (!data.parameters.contains(key))
          {
            Log::error() << "Hygov: missing required parameter '" << name << "'\n";
            ++parameter_error_count_;
            return;
          }

          load_real(key, target, name);
        };

        load_required_real(Params::Trate, Trate_, "Trate");
        load_real(Params::Rperm, Rperm_, "Rperm");
        load_real(Params::Rtemp, Rtemp_, "Rtemp");
        load_real(Params::Tr, Tr_, "Tr");
        load_real(Params::Tf, Tf_, "Tf");
        load_real(Params::Tg, Tg_, "Tg");
        load_real(Params::Velm, Velm_, "Velm");
        load_real(Params::Gmax, Gmax_, "Gmax");
        load_real(Params::Gmin, Gmin_, "Gmin");
        load_real(Params::Tw, Tw_, "Tw");
        load_real(Params::At, At_, "At");
        load_real(Params::Dturb, Dturb_, "Dturb");
        load_real(Params::Qnl, Qnl_, "Qnl");
        load_real(Params::Tn, Tn_, "Tn");
        load_real(Params::Tnp, Tnp_, "Tnp");
        load_real(Params::db1, db1_, "db1");
        load_real(Params::db2, db2_, "db2");
        load_real(Params::Hdam, Hdam_, "Hdam");
        load_real(Params::Gv0, Gv_[0], "Gv0");
        load_real(Params::Gv1, Gv_[1], "Gv1");
        load_real(Params::Gv2, Gv_[2], "Gv2");
        load_real(Params::Gv3, Gv_[3], "Gv3");
        load_real(Params::Gv4, Gv_[4], "Gv4");
        load_real(Params::Gv5, Gv_[5], "Gv5");
        load_real(Params::Pgv0, Pgv_[0], "Pgv0");
        load_real(Params::Pgv1, Pgv_[1], "Pgv1");
        load_real(Params::Pgv2, Pgv_[2], "Pgv2");
        load_real(Params::Pgv3, Pgv_[3], "Pgv3");
        load_real(Params::Pgv4, Pgv_[4], "Pgv4");
        load_real(Params::Pgv5, Pgv_[5], "Pgv5");

        auto check_nonnegative = [&](RealT value, const char* name)
        {
          if (value < ZERO<RealT>)
          {
            Log::error() << "Hygov: parameter '" << name << "' must be non-negative\n";
            ++parameter_error_count_;
          }
        };

        check_nonnegative(Tr_, "Tr");
        check_nonnegative(Tf_, "Tf");
        check_nonnegative(Tg_, "Tg");
        check_nonnegative(Tw_, "Tw");
        check_nonnegative(Tnp_, "Tnp");

        const bool source_default_curve =
            std::all_of(Gv_.begin(), Gv_.end(), [](RealT value)
                        { return value == ZERO<RealT>; })
            && std::all_of(Pgv_.begin(), Pgv_.end(), [](RealT value)
                           { return value == ZERO<RealT>; });
        if (source_default_curve)
        {
          Gv_  = {ZERO<RealT>,
                  static_cast<RealT>(0.2),
                  static_cast<RealT>(0.4),
                  static_cast<RealT>(0.6),
                  static_cast<RealT>(0.8),
                  ONE<RealT>};
          Pgv_ = Gv_;
        }

        setDerivedParameters();
      }

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* Hygov<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      template <typename scalar_type, typename index_type>
      void Hygov<scalar_type, index_type>::initializeMonitor()
      {
        using Variable = typename ModelDataT::MonitorableVariables;
        auto index     = [](HygovInternalVariables variable)
        {
          return static_cast<size_t>(variable);
        };

        monitor_->set(Variable::pmech, [this, index]
                      { return y_.getData()[index(HygovInternalVariables::PMECH)]; });
        monitor_->set(Variable::filter, [this, index]
                      { return y_.getData()[index(HygovInternalVariables::XF)]; });
        monitor_->set(Variable::desiredgate, [this, index]
                      { return y_.getData()[index(HygovInternalVariables::C)]; });
        monitor_->set(Variable::gate, [this, index]
                      { return y_.getData()[index(HygovInternalVariables::G)]; });
        monitor_->set(Variable::flow, [this, index]
                      { return y_.getData()[index(HygovInternalVariables::Q)]; });
        monitor_->set(Variable::head, [this, index]
                      { return y_.getData()[index(HygovInternalVariables::H)]; });
      }

      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::allocate()
      {
        size_     = static_cast<IdxT>(HygovInternalVariables::MAXIMUM);
        auto size = static_cast<size_t>(size_);

        if (!allocated_)
        {
          this->allocateVectors(size_);
        }

        tag_.assign(size, false);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        wb_.clear();

        auto signal_size = static_cast<size_t>(HygovExternalVariables::MAXIMUM);
        ws_.assign(signal_size, ScalarT{0});
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        if (signals_.template isAssigned<HygovInternalVariables::PMECH>())
        {
          auto* y = y_.getData();
          signals_.template getSignalNode<HygovInternalVariables::PMECH>()->set(
              &y[static_cast<size_t>(HygovInternalVariables::PMECH)],
              &(this->getVariableIndex(static_cast<IdxT>(HygovInternalVariables::PMECH))));
        }

        allocated_ = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::verify() const
      {
        int ret = parameter_error_count_;

        auto check = [&](bool condition, const char* message)
        {
          if (!condition)
          {
            Log::error() << "Hygov: " << message << '\n';
            ret += 1;
          }
        };

        check(Trate_ > ZERO<RealT>, "Trate must be positive");
        check(Rtemp_ != ZERO<RealT>, "Rtemp must be nonzero");
        check(Tr_ >= ZERO<RealT>, "Tr must be non-negative");
        check(Tf_ >= ZERO<RealT>, "Tf must be non-negative");
        check(Tg_ >= ZERO<RealT>, "Tg must be non-negative");
        check(Tw_ >= ZERO<RealT>, "Tw must be non-negative");
        check(Tn_ >= ZERO<RealT>, "Tn must be non-negative");
        check(Tnp_ >= ZERO<RealT>, "Tnp must be non-negative");
        check(Velm_ >= ZERO<RealT>, "Velm must be non-negative");
        check(Gmin_ <= Gmax_, "Gmin must be less than or equal to Gmax");
        check(At_ > ZERO<RealT>, "At must be positive");
        check(Dturb_ >= ZERO<RealT>, "Dturb must be non-negative");
        check(db1_ >= ZERO<RealT>, "db1 must be non-negative");
        check(Hdam_ > ZERO<RealT>, "Hdam must be positive");

        for (size_t i = 1; i < Gv_.size(); ++i)
        {
          check(Gv_[i - 1] < Gv_[i], "Gv points must be strictly increasing");
          check(Pgv_[i - 1] <= Pgv_[i], "Pgv points must be non-decreasing");
        }

        if (signals_.template isAttached<HygovExternalVariables::OMEGA>()
            && !signals_.template isLinked<HygovExternalVariables::OMEGA>())
        {
          Log::error() << "Hygov: omega signal attached with no linked source\n";
          ret += 1;
        }

        if (signals_.template isAttached<HygovExternalVariables::PREF>()
            && !signals_.template isLinked<HygovExternalVariables::PREF>())
        {
          Log::error() << "Hygov: pref signal attached with no linked source\n";
          ret += 1;
        }

        if (signals_.template isAttached<HygovExternalVariables::PAUX>()
            && !signals_.template isLinked<HygovExternalVariables::PAUX>())
        {
          Log::error() << "Hygov: paux signal attached with no linked source\n";
          ret += 1;
        }

        return ret;
      }

      template <typename scalar_type, typename index_type>
      scalar_type Hygov<scalar_type, index_type>::gatePower(scalar_type gate) const
      {
        return ScalarT{Pgv_[0]}
               + Math::linseg(gate, Gv_[0], Gv_[1], Pgv_[1] - Pgv_[0])
               + Math::linseg(gate, Gv_[1], Gv_[2], Pgv_[2] - Pgv_[1])
               + Math::linseg(gate, Gv_[2], Gv_[3], Pgv_[3] - Pgv_[2])
               + Math::linseg(gate, Gv_[3], Gv_[4], Pgv_[4] - Pgv_[3])
               + Math::linseg(gate, Gv_[4], Gv_[5], Pgv_[5] - Pgv_[4]);
      }

      template <typename scalar_type, typename index_type>
      typename Hygov<scalar_type, index_type>::RealT
      Hygov<scalar_type, index_type>::invertGatePower(
          typename Hygov<scalar_type, index_type>::RealT pgv) const
      {
        static constexpr RealT tol = static_cast<RealT>(1.0e-10);

        if (std::abs(pgv - Pgv_[0]) <= tol)
        {
          return Gv_[0];
        }

        for (size_t i = 0; i < 5; ++i)
        {
          if (Pgv_[i + 1] <= Pgv_[i])
          {
            continue;
          }

          if (Pgv_[i] - tol <= pgv && pgv <= Pgv_[i + 1] + tol)
          {
            const RealT fraction = (pgv - Pgv_[i]) / (Pgv_[i + 1] - Pgv_[i]);
            return Gv_[i] + fraction * (Gv_[i + 1] - Gv_[i]);
          }
        }

        return std::numeric_limits<RealT>::quiet_NaN();
      }

      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::initialize()
      {
        if (parameter_error_count_ > 0 || verify() > 0)
        {
          Log::error() << "Hygov: cannot initialize with invalid configuration\n";
          return 1;
        }

        const auto XN      = static_cast<size_t>(HygovInternalVariables::XN);
        const auto XF      = static_cast<size_t>(HygovInternalVariables::XF);
        const auto C       = static_cast<size_t>(HygovInternalVariables::C);
        const auto G       = static_cast<size_t>(HygovInternalVariables::G);
        const auto Q       = static_cast<size_t>(HygovInternalVariables::Q);
        const auto OMEGADB = static_cast<size_t>(HygovInternalVariables::OMEGADB);
        const auto EF      = static_cast<size_t>(HygovInternalVariables::EF);
        const auto FC      = static_cast<size_t>(HygovInternalVariables::FC);
        const auto RC      = static_cast<size_t>(HygovInternalVariables::RC);
        const auto PGV     = static_cast<size_t>(HygovInternalVariables::PGV);
        const auto H       = static_cast<size_t>(HygovInternalVariables::H);
        const auto PMECH   = static_cast<size_t>(HygovInternalVariables::PMECH);

        auto* y  = y_.getData();
        auto* yp = yp_.getData();

        ScalarT omega0{ZERO<RealT>};
        if (signals_.template isAttached<HygovExternalVariables::OMEGA>())
        {
          omega0 = signals_.template readExternalVariable<HygovExternalVariables::OMEGA>();
        }

        paux_set_ = ScalarT{ZERO<RealT>};
        if (signals_.template isAttached<HygovExternalVariables::PAUX>())
        {
          paux_set_ = signals_.template readExternalVariable<HygovExternalVariables::PAUX>();
        }

        const ScalarT paux0  = toComponentBase(paux_set_);
        const ScalarT pmech0 = toComponentBase(y[PMECH]);
        y[H]                 = Hdam_;
        y[Q]                 = Qnl_ + pmech0 / (At_ * y[H]);
        y[PGV]               = y[Q] / std::sqrt(y[H]);

        const RealT gate0 = invertGatePower(static_cast<RealT>(y[PGV]));
        if (std::isnan(gate0))
        {
          Log::error() << "Hygov: initial Pgv is outside the invertible gate curve\n";
          return 1;
        }

        y[G] = gate0;
        y[C] = y[G];

        if (y[C] < Gmin_ || y[C] > Gmax_)
        {
          Log::error() << "Hygov: initialized gate is outside Gmin/Gmax\n";
          return 1;
        }

        y[OMEGADB] = Math::deadband1(omega0, -db1_, db1_);
        y[XN]      = y[OMEGADB];
        y[XF]      = ZERO<RealT>;
        y[EF]      = ZERO<RealT>;
        y[FC]      = ZERO<RealT>;
        y[RC]      = ZERO<RealT>;
        y[PMECH]   = toSystemBase(pmech0);

        const ScalarT yomega = y[XN] + leadlag_gain_ * (y[OMEGADB] - y[XN]);
        pref_set_            = toSystemBase(y[EF] - paux0 + yomega + Rperm_ * y[C]);
        if (signals_.template isAttached<HygovExternalVariables::PREF>())
        {
          signals_.template writeExternalVariable<HygovExternalVariables::PREF>(pref_set_);
        }

        for (IdxT i = 0; i < size_; ++i)
        {
          yp[i] = ZERO<RealT>;
        }

        y_.setDataUpdated();
        yp_.setDataUpdated();
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::tagDifferentiable()
      {
        std::fill(tag_.begin(), tag_.end(), false);
        tag_[static_cast<size_t>(HygovInternalVariables::XN)] = true;
        tag_[static_cast<size_t>(HygovInternalVariables::XF)] = true;
        tag_[static_cast<size_t>(HygovInternalVariables::C)]  = true;
        tag_[static_cast<size_t>(HygovInternalVariables::G)]  = true;
        tag_[static_cast<size_t>(HygovInternalVariables::Q)]  = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline int
      Hygov<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT*                  y,
          const ScalarT*                  yp,
          [[maybe_unused]] const ScalarT* wb,
          const ScalarT*                  ws,
          ScalarT*                        f)
      {
        const auto XN      = static_cast<size_t>(HygovInternalVariables::XN);
        const auto XF      = static_cast<size_t>(HygovInternalVariables::XF);
        const auto C       = static_cast<size_t>(HygovInternalVariables::C);
        const auto G       = static_cast<size_t>(HygovInternalVariables::G);
        const auto Q       = static_cast<size_t>(HygovInternalVariables::Q);
        const auto OMEGADB = static_cast<size_t>(HygovInternalVariables::OMEGADB);
        const auto EF      = static_cast<size_t>(HygovInternalVariables::EF);
        const auto FC      = static_cast<size_t>(HygovInternalVariables::FC);
        const auto RC      = static_cast<size_t>(HygovInternalVariables::RC);
        const auto PGV     = static_cast<size_t>(HygovInternalVariables::PGV);
        const auto H       = static_cast<size_t>(HygovInternalVariables::H);
        const auto PMECH   = static_cast<size_t>(HygovInternalVariables::PMECH);

        const auto OMEGA = static_cast<size_t>(HygovExternalVariables::OMEGA);
        const auto PREF  = static_cast<size_t>(HygovExternalVariables::PREF);
        const auto PAUX  = static_cast<size_t>(HygovExternalVariables::PAUX);

        const ScalarT xn      = y[XN];
        const ScalarT xf      = y[XF];
        const ScalarT c       = y[C];
        const ScalarT g       = y[G];
        const ScalarT q       = y[Q];
        const ScalarT omegadb = y[OMEGADB];
        const ScalarT ef      = y[EF];
        const ScalarT fc      = y[FC];
        const ScalarT rc      = y[RC];
        const ScalarT pgv     = y[PGV];
        const ScalarT head    = y[H];
        const ScalarT pmech   = y[PMECH];

        const ScalarT omega = ws[OMEGA];
        const ScalarT pref  = toComponentBase(ws[PREF]);
        const ScalarT paux  = toComponentBase(ws[PAUX]);

        const ScalarT yomega = xn + leadlag_gain_ * (omegadb - xn);

        f[XN]      = -yp[XN] + (omegadb - xn) / Tnp_;
        f[XF]      = -yp[XF] + (ef - xf) / Tf_;
        f[C]       = -yp[C] + Math::antiwindup(c, rc, Gmin_, Gmax_);
        f[G]       = -yp[G] + (c - g) / Tg_;
        f[Q]       = -yp[Q] + (Hdam_ - head) / Tw_;
        f[OMEGADB] = -omegadb + Math::deadband1(omega, -db1_, db1_);
        f[EF]      = -ef + pref + paux - yomega - Rperm_ * c;
        f[FC]      = -fc + (xf / Tr_ + (ef - xf) / Tf_) / Rtemp_;
        f[RC]      = -rc + Math::clamp(fc, -Velm_, Velm_);
        f[PGV]     = -pgv + gatePower(g);
        f[H]       = -q * q + head * pgv * pgv;
        f[PMECH]   = -toComponentBase(pmech) + At_ * head * (q - Qnl_) - Dturb_ * omega * g;

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::evaluateResidual()
      {
        const auto OMEGA = static_cast<size_t>(HygovExternalVariables::OMEGA);
        const auto PREF  = static_cast<size_t>(HygovExternalVariables::PREF);
        const auto PAUX  = static_cast<size_t>(HygovExternalVariables::PAUX);

        ws_[OMEGA] = ZERO<RealT>;
        ws_[PREF]  = pref_set_;
        ws_[PAUX]  = paux_set_;
        std::fill(ws_indices_.begin(), ws_indices_.end(), INVALID_INDEX<IdxT>);

        if (signals_.template isAttached<HygovExternalVariables::OMEGA>())
        {
          ws_[OMEGA] = signals_.template readExternalVariable<HygovExternalVariables::OMEGA>();
          ws_indices_[OMEGA] =
              signals_.template readExternalVariableIndex<HygovExternalVariables::OMEGA>();
        }

        if (signals_.template isAttached<HygovExternalVariables::PREF>())
        {
          ws_[PREF] = signals_.template readExternalVariable<HygovExternalVariables::PREF>();
          ws_indices_[PREF] =
              signals_.template readExternalVariableIndex<HygovExternalVariables::PREF>();
        }

        if (signals_.template isAttached<HygovExternalVariables::PAUX>())
        {
          ws_[PAUX] = signals_.template readExternalVariable<HygovExternalVariables::PAUX>();
          ws_indices_[PAUX] =
              signals_.template readExternalVariableIndex<HygovExternalVariables::PAUX>();
        }

        const auto* y  = y_.getData();
        const auto* yp = yp_.getData();
        auto*       f  = f_.getData();
        evaluateInternalResidual(y, yp, wb_.data(), ws_.data(), f);

        f_.setDataUpdated();
        return 0;
      }
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
