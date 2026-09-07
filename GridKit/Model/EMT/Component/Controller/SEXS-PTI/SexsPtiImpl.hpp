#pragma once

#include <GridKit/Model/EMT/ComponentInitialization.hpp>

/**
 * @file SexsPtiImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the SEXS-PTI exciter model.
 */

#include <cmath>
#include <iostream>

#include <GridKit/Model/EMT/Component/Controller/SEXS-PTI/SexsPti.hpp>
#include <GridKit/Model/EMT/Component/Controller/SEXS-PTI/SexsPtiData.hpp>
#include <GridKit/Model/EMT/Signal/Signal.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Controller
    {
      using Log = ::GridKit::Utilities::Logger;

      template <typename scalar_type, typename index_type>
      SexsPti<scalar_type, index_type>::SexsPti()
        : SexsPti(ModelDataT{})
      {
        size_ = 4;
      }

      template <typename scalar_type, typename index_type>
      SexsPti<scalar_type, index_type>::SexsPti(const ModelDataT& data)
        : monitor_(std::make_unique<MonitorT>(data))
      {
        initModelParams(data);
        initializeMonitor();
        size_ = 4;
      }

      template <typename scalar_type, typename index_type>
      SexsPti<scalar_type, index_type>::~SexsPti()
      {
      }

      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::allocate()
      {
        if (!allocated_)
        {
          this->allocateVectors(size_);
        }
        auto size = static_cast<size_t>(size_);

        tag_.resize(size);

        variable_indices_.resize(size);
        residual_indices_.resize(size);

        // Default variable and residual index mapping to local index
        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        this->allocateExternalVectors(static_cast<IdxT>(SexsPtiExternalVariables::MAXIMUM), 0);
        signals_.registerExternalVariableSignals(*this);
        uel_on_ = signals_.template isAttached<SexsPtiExternalVariables::VUEL>() ? ONE<RealT> : ZERO<RealT>;
        oel_on_ = signals_.template isAttached<SexsPtiExternalVariables::VOEL>() ? ONE<RealT> : ZERO<RealT>;

        if (signals_.template isAssigned<SexsPtiInternalVariables::EFD>())
        {
          auto* y = y_.getData();
          signals_.template getSignal<SexsPtiInternalVariables::EFD>()->set(
              &y[1], &(this->getVariableIndex(1)));
        }

        allocated_ = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::verify() const
      {
        int ret = missing_param_count_;

        if (!signals_.template isAttached<SexsPtiExternalVariables::VA>()
            || !signals_.template isAttached<SexsPtiExternalVariables::VB>()
            || !signals_.template isAttached<SexsPtiExternalVariables::VC>())
        {
          Log::error() << "SexsPti: terminal bus port must be attached\n";
          ret += 1;
        }
        for (RealT value : {V_, Tr_, Ta_, Tb_, Te_, K_, Efdmin_, Efdmax_})
        {
          if (!std::isfinite(value))
          {
            Log::error() << "SexsPti: parameters must be finite\n";
            ++ret;
          }
        }
        if (V_ <= ZERO<RealT> || Tr_ < ZERO<RealT>)
        {
          Log::error() << "SexsPti: V must be positive and Tr non-negative\n";
          ++ret;
        }
        if (Ta_ < 0.0)
        {
          Log::error() << "SexsPti: Ta must be non-negative\n";
          ret += 1;
        }
        if (Tb_ <= 0.0)
        {
          Log::error() << "SexsPti: Tb must be positive\n";
          ret += 1;
        }
        if (Te_ <= 0.0)
        {
          Log::error() << "SexsPti: Te must be positive\n";
          ret += 1;
        }
        if (K_ <= 0.0)
        {
          Log::error() << "SexsPti: K must be positive\n";
          ret += 1;
        }
        if (Efdmin_ >= Efdmax_)
        {
          Log::error() << "SexsPti: Efdmin must be less than Efdmax\n";
          ret += 1;
        }

        if (!signals_.template isAssigned<SexsPtiInternalVariables::EFD>())
        {
          Log::error() << "SexsPti: required EFD signal is not assigned\n";
          ret += 1;
        }

        auto check_attached_signal =
            [&]<SexsPtiExternalVariables variable>(const char* name)
        {
          if (signals_.template isAttached<variable>()
              && !signals_.template isLinked<variable>())
          {
            Log::error() << "SexsPti: " << name << " signal attached with no linked source\n";
            ret += 1;
          }
        };

        check_attached_signal.template operator()<SexsPtiExternalVariables::VA>("bus phase a");
        check_attached_signal.template operator()<SexsPtiExternalVariables::VB>("bus phase b");
        check_attached_signal.template operator()<SexsPtiExternalVariables::VC>("bus phase c");
        check_attached_signal.template operator()<SexsPtiExternalVariables::VREF>("vref");
        check_attached_signal.template operator()<SexsPtiExternalVariables::VS>("vs");
        check_attached_signal.template operator()<SexsPtiExternalVariables::VUEL>("vuel");
        check_attached_signal.template operator()<SexsPtiExternalVariables::VOEL>("voel");

        return ret;
      }

      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::initialize(const std::map<Outputs, RealT>& outputs)
      {
        this->validateOutputValues(outputs);
        if (!allocated_ || verify() != 0)
          return 1;
        gatherExternalVariables();
        ScalarT efd0{0.0};
        auto*   y  = y_.getData();
        auto*   yp = yp_.getData();

        if (signals_.template isAssigned<SexsPtiInternalVariables::EFD>())
        {
          efd0 = y[1];
        }

        efd0 = this->outputValue(outputs, Outputs::efd, static_cast<RealT>(efd0));

        // Setpoint members provide the defaults for unattached signals.
        auto read_signal = [&]<SexsPtiExternalVariables variable>(const ScalarT& default_value) -> ScalarT
        {
          if (signals_.template isAttached<variable>())
          {
            return signals_.template readExternalVariable<variable>();
          }
          return default_value;
        };

        const ScalarT vs   = read_signal.template operator()<SexsPtiExternalVariables::VS>(vs_set_);
        const ScalarT vuel = read_signal.template operator()<SexsPtiExternalVariables::VUEL>(vuel_set_);
        const ScalarT voel = read_signal.template operator()<SexsPtiExternalVariables::VOEL>(voel_set_);

        uel_on_ = ZERO<RealT>;
        if (signals_.template isAttached<SexsPtiExternalVariables::VUEL>())
        {
          uel_on_ = ONE<RealT>;
        }

        oel_on_ = ZERO<RealT>;
        if (signals_.template isAttached<SexsPtiExternalVariables::VOEL>())
        {
          oel_on_ = ONE<RealT>;
        }

        ScalarT Ec   = voltageMagnitude(y_ext_.data() + 4);
        ScalarT vtr  = efd0 / K_;
        ScalarT vr   = (Ta_ - Tb_) * vtr;
        ScalarT vref = Ec + vtr - vs - uel_on_ * vuel - oel_on_ * voel;

        for (const ScalarT& value : {efd0, Ec, vtr, vr, vref})
        {
          if (!std::isfinite(static_cast<RealT>(value)))
            return 1;
        }
        if (static_cast<RealT>(efd0) < Efdmin_ || static_cast<RealT>(efd0) > Efdmax_)
        {
          Log::error() << "SexsPti: initial EFD is outside its limits\n";
          return 1;
        }
        y[3] = Ec;
        y[0] = vr;
        y[1] = efd0;
        y[2] = vtr;

        for (IdxT i = 0; i < size_; ++i)
        {
          yp[static_cast<size_t>(i)] = 0.0;
        }

        vref_set_ = vref;
        vs_set_   = vs;
        vuel_set_ = vuel;
        voel_set_ = voel;

        if (signals_.template isAttached<SexsPtiExternalVariables::VREF>())
        {
          signals_.template writeExternalVariable<SexsPtiExternalVariables::VREF>(vref_set_);
        }

        y_.setDataUpdated();
        yp_.setDataUpdated();

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
      int SexsPti<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline int SexsPti<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT*                  y,
          const ScalarT*                  yp,
          const ScalarT*                  external,
          [[maybe_unused]] const ScalarT* external_dot,
          ScalarT*                        f)
      {
        const auto VREF = static_cast<size_t>(SexsPtiExternalVariables::VREF);
        const auto VS   = static_cast<size_t>(SexsPtiExternalVariables::VS);
        const auto VUEL = static_cast<size_t>(SexsPtiExternalVariables::VUEL);
        const auto VOEL = static_cast<size_t>(SexsPtiExternalVariables::VOEL);

        ScalarT vr      = y[0];
        ScalarT efd     = y[1];
        ScalarT vtr     = y[2];
        ScalarT vr_dot  = yp[0];
        ScalarT efd_dot = yp[1];

        ScalarT Ec   = y[3];
        ScalarT vref = external[VREF];
        ScalarT vs   = external[VS];
        ScalarT vuel = external[VUEL];
        ScalarT voel = external[VOEL];

        ScalarT func = (-efd + (K_ / Tb_) * (-vr + Ta_ * vtr)) / Te_;

        f[0] = -vr_dot + (-vr + Ta_ * vtr) / Tb_ - vtr;
        f[1] = -efd_dot + Math::antiwindup(efd, func, Efdmin_, Efdmax_);
        f[2] = -vtr - Ec + vref + vs + uel_on_ * vuel + oel_on_ * voel;
        f[3] = -Tr_ * yp[3] + voltageMagnitude(external + 4) - y[3];

        return 0;
      }

      template <typename scalar_type, typename index_type>
      void SexsPti<scalar_type, index_type>::gatherExternalVariables()
      {
        y_ext_[0] = vref_set_;
        y_ext_[1] = vs_set_;
        y_ext_[2] = vuel_set_;
        y_ext_[3] = voel_set_;
        Component<scalar_type, index_type>::gatherExternalVariables();
      }

      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::evaluateInternalResidual()
      {
        gatherExternalVariables();
        evaluateInternalResidual(y_.getData(), yp_.getData(), y_ext_.data(), yp_ext_.data(), f_.getData());
        f_.setDataUpdated();
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::evaluateResidual()
      {
        return evaluateInternalResidual();
      }

      template <typename scalar_type, typename index_type>
      void SexsPti<scalar_type, index_type>::initModelParams(const ModelDataT& data)
      {
        using Params = typename ModelDataT::Parameters;

        missing_param_count_ = 0;

        auto load = [&](auto param, RealT& member, const char* name)
        {
          if (data.parameters.contains(param))
          {
            member = std::get<RealT>(data.parameters.at(param));
          }
          else
          {
            Log::error() << "SexsPti: missing required parameter '" << name << "'\n";
            ++missing_param_count_;
          }
        };

        load(Params::V, V_, "V");
        if (data.parameters.contains(Params::Tr))
          Tr_ = std::get<RealT>(data.parameters.at(Params::Tr));
        load(Params::Ta, Ta_, "Ta");
        load(Params::Tb, Tb_, "Tb");
        load(Params::Te, Te_, "Te");
        load(Params::K, K_, "K");
        load(Params::Efdmax, Efdmax_, "Efdmax");
        load(Params::Efdmin, Efdmin_, "Efdmin");
      }

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* SexsPti<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      template <typename scalar_type, typename index_type>
      void SexsPti<scalar_type, index_type>::initializeMonitor()
      {
        using Variable = typename ModelDataT::MonitorableVariables;
        monitor_->set(Variable::vts, [this]
                      { return y_.getData()[3]; });
        monitor_->set(Variable::vr, [this]
                      { return y_.getData()[0]; });
        monitor_->set(Variable::vtr, [this]
                      { return y_.getData()[2]; });
        monitor_->set(Variable::efd, [this]
                      { return y_.getData()[1]; });
      }

    } // namespace Controller
  } // namespace EMT
} // namespace GridKit
