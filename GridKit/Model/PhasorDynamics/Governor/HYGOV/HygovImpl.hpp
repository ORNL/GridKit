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

      template <class ScalarT, typename IdxT>
      Hygov<ScalarT, IdxT>::Hygov()
      {
        size_ = static_cast<IdxT>(HygovInternalVariables::MAXIMUM);
      }

      template <class ScalarT, typename IdxT>
      Hygov<ScalarT, IdxT>::Hygov(const model_data_type& data)
        : monitor_(std::make_unique<MonitorT>(data))
      {
        initModelParams(data);
        initializeMonitor();
        size_ = static_cast<IdxT>(HygovInternalVariables::MAXIMUM);
      }

      template <class ScalarT, typename IdxT>
      Hygov<ScalarT, IdxT>::~Hygov()
      {
      }

      template <class ScalarT, typename IdxT>
      ScalarT Hygov<ScalarT, IdxT>::toComponentBase(ScalarT value) const
      {
        return value * pmech_mva_base_ / Trate_;
      }

      template <class ScalarT, typename IdxT>
      ScalarT Hygov<ScalarT, IdxT>::toPmechBase(ScalarT value) const
      {
        return value * Trate_ / pmech_mva_base_;
      }

      template <class ScalarT, typename IdxT>
      void Hygov<ScalarT, IdxT>::initModelParams(const model_data_type& data)
      {
        using Params = typename model_data_type::Parameters;

        parameter_error_count_ = 0;

        auto loadRequiredReal = [&](auto key, RealT& target, const char* name)
        {
          if (!data.parameters.contains(key))
          {
            Log::error() << "Hygov: missing required parameter '" << name << "'\n";
            ++parameter_error_count_;
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

        auto loadOptionalReal = [&](auto key, RealT& target, const char* name)
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

        loadRequiredReal(Params::Trate, Trate_, "Trate");
        loadRequiredReal(Params::mva, pmech_mva_base_, "mva");
        loadRequiredReal(Params::Rperm, Rperm_, "Rperm");
        loadRequiredReal(Params::Rtemp, Rtemp_, "Rtemp");
        loadRequiredReal(Params::Tr, Tr_, "Tr");
        loadRequiredReal(Params::Tf, Tf_, "Tf");
        loadRequiredReal(Params::Tg, Tg_, "Tg");
        loadRequiredReal(Params::Velm, Velm_, "Velm");
        loadRequiredReal(Params::Gmax, Gmax_, "Gmax");
        loadRequiredReal(Params::Gmin, Gmin_, "Gmin");
        loadRequiredReal(Params::Tw, Tw_, "Tw");
        loadRequiredReal(Params::At, At_, "At");
        loadRequiredReal(Params::Dturb, Dturb_, "Dturb");
        loadRequiredReal(Params::Qnl, Qnl_, "Qnl");
        loadRequiredReal(Params::Tn, Tn_, "Tn");
        loadRequiredReal(Params::Tnp, Tnp_, "Tnp");
        loadRequiredReal(Params::db1, db1_, "db1");
        loadOptionalReal(Params::db2, db2_, "db2");
        loadRequiredReal(Params::Hdam, Hdam_, "Hdam");

        loadRequiredReal(Params::Gv0, Gv_[0], "Gv0");
        loadRequiredReal(Params::Gv1, Gv_[1], "Gv1");
        loadRequiredReal(Params::Gv2, Gv_[2], "Gv2");
        loadRequiredReal(Params::Gv3, Gv_[3], "Gv3");
        loadRequiredReal(Params::Gv4, Gv_[4], "Gv4");
        loadRequiredReal(Params::Gv5, Gv_[5], "Gv5");
        loadRequiredReal(Params::Pgv0, Pgv_[0], "Pgv0");
        loadRequiredReal(Params::Pgv1, Pgv_[1], "Pgv1");
        loadRequiredReal(Params::Pgv2, Pgv_[2], "Pgv2");
        loadRequiredReal(Params::Pgv3, Pgv_[3], "Pgv3");
        loadRequiredReal(Params::Pgv4, Pgv_[4], "Pgv4");
        loadRequiredReal(Params::Pgv5, Pgv_[5], "Pgv5");

        const bool source_default_curve =
            std::all_of(Gv_.begin(), Gv_.end(), [](RealT value)
                        { return value == ZERO<RealT>; })
            && std::all_of(Pgv_.begin(), Pgv_.end(), [](RealT value)
                           { return value == ZERO<RealT>; });
        if (source_default_curve)
        {
          Gv_  = {ZERO<RealT>, static_cast<RealT>(0.2), static_cast<RealT>(0.4), static_cast<RealT>(0.6), static_cast<RealT>(0.8), ONE<RealT>};
          Pgv_ = Gv_;
        }

        leadlag_gain_ = (Tnp_ > ZERO<RealT>) ? Tn_ / Tnp_ : ZERO<RealT>;
      }

      template <class ScalarT, typename IdxT>
      const Model::VariableMonitorBase* Hygov<ScalarT, IdxT>::getMonitor() const
      {
        return monitor_.get();
      }

      template <class ScalarT, typename IdxT>
      void Hygov<ScalarT, IdxT>::initializeMonitor()
      {
        using Variable = typename model_data_type::MonitorableVariables;
        auto index     = [](HygovInternalVariables variable)
        {
          return static_cast<size_t>(variable);
        };

        monitor_->set(Variable::pmech, [this, index]
                      { return y_[index(HygovInternalVariables::PMECH)]; });
        monitor_->set(Variable::leadlag, [this, index]
                      { return y_[index(HygovInternalVariables::XN)]; });
        monitor_->set(Variable::filter, [this, index]
                      { return y_[index(HygovInternalVariables::XF)]; });
        monitor_->set(Variable::desiredgate, [this, index]
                      { return y_[index(HygovInternalVariables::C)]; });
        monitor_->set(Variable::gate, [this, index]
                      { return y_[index(HygovInternalVariables::G)]; });
        monitor_->set(Variable::flow, [this, index]
                      { return y_[index(HygovInternalVariables::Q)]; });
        monitor_->set(Variable::head, [this, index]
                      { return y_[index(HygovInternalVariables::H)]; });
        monitor_->set(Variable::pgv, [this, index]
                      { return y_[index(HygovInternalVariables::PGV)]; });
      }

      template <class ScalarT, typename IdxT>
      int Hygov<ScalarT, IdxT>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      template <class ScalarT, typename IdxT>
      int Hygov<ScalarT, IdxT>::allocate()
      {
        size_     = static_cast<IdxT>(HygovInternalVariables::MAXIMUM);
        auto size = static_cast<size_t>(size_);

        f_.assign(size, ScalarT{0});
        y_.assign(size, ScalarT{0});
        yp_.assign(size, ScalarT{0});
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
          signals_.template getSignalNode<HygovInternalVariables::PMECH>()->set(
              &y_[static_cast<size_t>(HygovInternalVariables::PMECH)],
              &(this->getVariableIndex(static_cast<IdxT>(HygovInternalVariables::PMECH))));
        }

        return 0;
      }

      template <class ScalarT, typename IdxT>
      int Hygov<ScalarT, IdxT>::verify() const
      {
        int ret = static_cast<int>(parameter_error_count_);

        auto check = [&](bool condition, const char* message)
        {
          if (!condition)
          {
            Log::error() << "Hygov: " << message << '\n';
            ret += 1;
          }
        };

        check(Trate_ > ZERO<RealT>, "Trate must be positive");
        check(pmech_mva_base_ > ZERO<RealT>, "mva must be positive");
        check(Rperm_ > ZERO<RealT>, "Rperm must be positive");
        check(Rtemp_ > ZERO<RealT>, "Rtemp must be positive");
        check(Tr_ > ZERO<RealT>, "Tr must be positive");
        check(Tf_ > ZERO<RealT>, "Tf must be positive");
        check(Tg_ >= ZERO<RealT>, "Tg must be non-negative");
        check(Tw_ > ZERO<RealT>, "Tw must be positive");
        check(Tn_ >= ZERO<RealT>, "Tn must be non-negative");
        check(Tnp_ >= ZERO<RealT>, "Tnp must be non-negative");
        check(Tnp_ > ZERO<RealT> || Tn_ == ZERO<RealT>, "Tn must be zero when Tnp is zero");
        check(Velm_ > ZERO<RealT>, "Velm must be positive");
        check(Gmin_ <= Gmax_, "Gmin must be less than or equal to Gmax");
        check(At_ > ZERO<RealT>, "At must be positive");
        check(Dturb_ >= ZERO<RealT>, "Dturb must be non-negative");
        check(Qnl_ >= ZERO<RealT>, "Qnl must be non-negative");
        check(db1_ >= ZERO<RealT>, "db1 must be non-negative");
        check(db2_ == ZERO<RealT>, "nonzero db2 is not supported");
        check(Hdam_ > ZERO<RealT>, "Hdam must be positive");

        for (size_t i = 1; i < Gv_.size(); ++i)
        {
          check(Gv_[i - 1] < Gv_[i], "Gv points must be strictly increasing");
          check(Pgv_[i - 1] <= Pgv_[i], "Pgv points must be non-decreasing");
        }
        check(Pgv_[0] >= ZERO<RealT>, "Pgv0 must be non-negative");

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

      template <class ScalarT, typename IdxT>
      ScalarT Hygov<ScalarT, IdxT>::gatePower(ScalarT gate) const
      {
        return ScalarT{Pgv_[0]}
               + Math::linseg(gate, Gv_[0], Gv_[1], Pgv_[1] - Pgv_[0])
               + Math::linseg(gate, Gv_[1], Gv_[2], Pgv_[2] - Pgv_[1])
               + Math::linseg(gate, Gv_[2], Gv_[3], Pgv_[3] - Pgv_[2])
               + Math::linseg(gate, Gv_[3], Gv_[4], Pgv_[4] - Pgv_[3])
               + Math::linseg(gate, Gv_[4], Gv_[5], Pgv_[5] - Pgv_[4]);
      }

      template <class ScalarT, typename IdxT>
      typename Hygov<ScalarT, IdxT>::RealT Hygov<ScalarT, IdxT>::invertGatePower(
          typename Hygov<ScalarT, IdxT>::RealT pgv) const
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

      template <class ScalarT, typename IdxT>
      int Hygov<ScalarT, IdxT>::initialize()
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
        const auto YOMEGA  = static_cast<size_t>(HygovInternalVariables::YOMEGA);
        const auto EF      = static_cast<size_t>(HygovInternalVariables::EF);
        const auto FC      = static_cast<size_t>(HygovInternalVariables::FC);
        const auto RC      = static_cast<size_t>(HygovInternalVariables::RC);
        const auto PGV     = static_cast<size_t>(HygovInternalVariables::PGV);
        const auto PGVSAFE = static_cast<size_t>(HygovInternalVariables::PGVSAFE);
        const auto H       = static_cast<size_t>(HygovInternalVariables::H);
        const auto PMECH   = static_cast<size_t>(HygovInternalVariables::PMECH);

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

        const ScalarT pmech0 = toComponentBase(y_[PMECH]);
        y_[H]                = Hdam_;
        y_[Q]                = Qnl_ + pmech0 / (At_ * y_[H]);
        y_[PGV]              = y_[Q] / std::sqrt(y_[H]);
        y_[PGVSAFE]          = Math::max(y_[PGV], static_cast<RealT>(0.01));

        const RealT gate0 = invertGatePower(static_cast<RealT>(y_[PGV]));
        if (std::isnan(gate0))
        {
          Log::error() << "Hygov: initial Pgv is outside the invertible gate curve\n";
          return 1;
        }

        y_[G] = gate0;
        y_[C] = y_[G];

        if (y_[C] < Gmin_ || y_[C] > Gmax_)
        {
          Log::error() << "Hygov: initialized gate is outside Gmin/Gmax\n";
          return 1;
        }

        y_[OMEGADB] = Math::deadband1(omega0, -db1_, db1_);
        y_[XN]      = y_[OMEGADB];
        y_[YOMEGA]  = y_[OMEGADB];
        y_[XF]      = ZERO<RealT>;
        y_[EF]      = ZERO<RealT>;
        y_[FC]      = ZERO<RealT>;
        y_[RC]      = ZERO<RealT>;
        y_[PMECH]   = toPmechBase(pmech0);

        pref_set_ = y_[EF] - paux_set_ + y_[YOMEGA] + Rperm_ * y_[C];
        if (signals_.template isAttached<HygovExternalVariables::PREF>())
        {
          const ScalarT pref0    = signals_.template readExternalVariable<HygovExternalVariables::PREF>();
          const RealT   pref_err = static_cast<RealT>(pref0 - pref_set_);
          if (std::abs(pref_err) > static_cast<RealT>(1.0e-10))
          {
            Log::error() << "Hygov: pref initial condition is not steady state\n";
            return 1;
          }
          pref_set_ = pref0;
        }

        std::fill(yp_.begin(), yp_.end(), ZERO<RealT>);
        return 0;
      }

      template <class ScalarT, typename IdxT>
      int Hygov<ScalarT, IdxT>::tagDifferentiable()
      {
        std::fill(tag_.begin(), tag_.end(), false);
        tag_[static_cast<size_t>(HygovInternalVariables::XN)] = (Tnp_ > ZERO<RealT>);
        tag_[static_cast<size_t>(HygovInternalVariables::XF)] = true;
        tag_[static_cast<size_t>(HygovInternalVariables::C)]  = true;
        tag_[static_cast<size_t>(HygovInternalVariables::G)]  = (Tg_ > ZERO<RealT>);
        tag_[static_cast<size_t>(HygovInternalVariables::Q)]  = true;
        return 0;
      }

      template <class ScalarT, typename IdxT>
      __attribute__((always_inline)) inline int Hygov<ScalarT, IdxT>::evaluateInternalResidual(
          ScalarT*                  y,
          ScalarT*                  yp,
          [[maybe_unused]] ScalarT* wb,
          ScalarT*                  ws,
          ScalarT*                  f)
      {
        const auto XN      = static_cast<size_t>(HygovInternalVariables::XN);
        const auto XF      = static_cast<size_t>(HygovInternalVariables::XF);
        const auto C       = static_cast<size_t>(HygovInternalVariables::C);
        const auto G       = static_cast<size_t>(HygovInternalVariables::G);
        const auto Q       = static_cast<size_t>(HygovInternalVariables::Q);
        const auto OMEGADB = static_cast<size_t>(HygovInternalVariables::OMEGADB);
        const auto YOMEGA  = static_cast<size_t>(HygovInternalVariables::YOMEGA);
        const auto EF      = static_cast<size_t>(HygovInternalVariables::EF);
        const auto FC      = static_cast<size_t>(HygovInternalVariables::FC);
        const auto RC      = static_cast<size_t>(HygovInternalVariables::RC);
        const auto PGV     = static_cast<size_t>(HygovInternalVariables::PGV);
        const auto PGVSAFE = static_cast<size_t>(HygovInternalVariables::PGVSAFE);
        const auto H       = static_cast<size_t>(HygovInternalVariables::H);
        const auto PMECH   = static_cast<size_t>(HygovInternalVariables::PMECH);

        const auto OMEGA = static_cast<size_t>(HygovExternalVariables::OMEGA);
        const auto PREF  = static_cast<size_t>(HygovExternalVariables::PREF);
        const auto PAUX  = static_cast<size_t>(HygovExternalVariables::PAUX);

        const ScalarT xn       = y[XN];
        const ScalarT xf       = y[XF];
        const ScalarT c        = y[C];
        const ScalarT g        = y[G];
        const ScalarT q        = y[Q];
        const ScalarT omegadb  = y[OMEGADB];
        const ScalarT yomega   = y[YOMEGA];
        const ScalarT ef       = y[EF];
        const ScalarT fc       = y[FC];
        const ScalarT rc       = y[RC];
        const ScalarT pgv      = y[PGV];
        const ScalarT pgv_safe = y[PGVSAFE];
        const ScalarT head     = y[H];
        const ScalarT pmech    = y[PMECH];

        const ScalarT omega = ws[OMEGA];
        const ScalarT pref  = ws[PREF];
        const ScalarT paux  = ws[PAUX];

        f[XN] = -Tnp_ * yp[XN] - xn + omegadb;
        f[XF] = -Tf_ * yp[XF] - xf + ef;
        f[FC] = -Rtemp_ * Tr_ * fc + xf + Tr_ * yp[XF];
        f[C]  = -yp[C] + Math::antiwindup(c, rc, Gmin_, Gmax_);
        f[G]  = -Tg_ * yp[G] - g + c;
        f[Q]  = -Tw_ * yp[Q] + Hdam_ - head;

        f[OMEGADB]                 = -omegadb + Math::deadband1(omega, -db1_, db1_);
        f[YOMEGA]                  = -yomega + xn + leadlag_gain_ * (omegadb - xn);
        f[EF]                      = -ef + pref + paux - yomega - Rperm_ * c;
        f[RC]                      = -rc + Math::clamp(fc, -Velm_, Velm_);
        f[PGV]                     = -pgv + gatePower(g);
        f[PGVSAFE]                 = -pgv_safe + Math::max(pgv, static_cast<RealT>(0.01));
        f[H]                       = -head + (q / pgv_safe) * (q / pgv_safe);
        const ScalarT pmech_output = At_ * head * (q - Qnl_) - Dturb_ * omega * g;
        f[PMECH]                   = -pmech + toPmechBase(pmech_output);

        return 0;
      }

      template <class ScalarT, typename IdxT>
      int Hygov<ScalarT, IdxT>::evaluateResidual()
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

        evaluateInternalResidual(y_.data(), yp_.data(), wb_.data(), ws_.data(), f_.data());
        return 0;
      }
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
