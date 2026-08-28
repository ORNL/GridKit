#pragma once

/**
 * @file SexsPtiImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the SEXS-PTI exciter model.
 */

#include <iostream>

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/SEXS-PTI/SexsPti.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/SEXS-PTI/SexsPtiData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      using Log = ::GridKit::Utilities::Logger;

      template <typename scalar_type, typename index_type>
      SexsPti<scalar_type, index_type>::SexsPti(BusT* bus)
        : bus_(bus)
      {
        size_ = 3;
      }

      template <typename scalar_type, typename index_type>
      SexsPti<scalar_type, index_type>::SexsPti(BusT*             bus,
                                                const ModelDataT& data)
        : bus_(bus),
          monitor_(std::make_unique<MonitorT>(data))
      {
        initModelParams(data);
        initializeMonitor();
        size_ = 3;
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

        wb_.resize(2);

        const auto signal_size = static_cast<size_t>(SexsPtiExternalVariables::MAXIMUM);
        ws_.assign(signal_size, ScalarT{0});
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        if (signals_.template isAssigned<SexsPtiInternalVariables::EFD>())
        {
          auto* y = y_.getData();
          signals_.template getSignalNode<SexsPtiInternalVariables::EFD>()->set(
              &y[1], &(this->getVariableIndex(1)));
        }

        allocated_ = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::verify() const
      {
        int ret = missing_param_count_;

        if (bus_ == nullptr)
        {
          Log::error() << "SexsPti: bus pointer is null\n";
          ret += 1;
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

        check_attached_signal.template operator()<SexsPtiExternalVariables::VREF>("vref");
        check_attached_signal.template operator()<SexsPtiExternalVariables::VS>("vs");
        check_attached_signal.template operator()<SexsPtiExternalVariables::VUEL>("vuel");
        check_attached_signal.template operator()<SexsPtiExternalVariables::VOEL>("voel");

        return ret;
      }

      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::initialize()
      {
        ScalarT efd0{0.0};
        auto*   y  = y_.getData();
        auto*   yp = yp_.getData();

        if (signals_.template isAssigned<SexsPtiInternalVariables::EFD>())
        {
          efd0 = y[1];
        }

        auto read_attached = [&]<SexsPtiExternalVariables variable>() -> ScalarT
        {
          if (signals_.template isAttached<variable>())
          {
            return signals_.template readExternalVariable<variable>();
          }
          return ScalarT{0};
        };

        const ScalarT vs   = read_attached.template operator()<SexsPtiExternalVariables::VS>();
        const ScalarT vuel = read_attached.template operator()<SexsPtiExternalVariables::VUEL>();
        const ScalarT voel = read_attached.template operator()<SexsPtiExternalVariables::VOEL>();

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

        ScalarT vreal = bus_->Vr();
        ScalarT vimag = bus_->Vi();
        ScalarT Ec    = std::sqrt(vreal * vreal + vimag * vimag);
        ScalarT vtr   = efd0 / K_;
        ScalarT vr    = (Ta_ - Tb_) * vtr;
        ScalarT vref  = Ec + vtr - vs - uel_on_ * vuel - oel_on_ * voel;

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

      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::tagDifferentiable()
      {
        tag_[0] = true;
        tag_[1] = true;
        tag_[2] = false;

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
          const ScalarT* y,
          const ScalarT* yp,
          const ScalarT* wb,
          const ScalarT* ws,
          ScalarT*       f)
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

        ScalarT Ec   = std::sqrt(wb[0] * wb[0] + wb[1] * wb[1]);
        ScalarT vref = ws[VREF];
        ScalarT vs   = ws[VS];
        ScalarT vuel = ws[VUEL];
        ScalarT voel = ws[VOEL];

        ScalarT func = (-efd + (K_ / Tb_) * (-vr + Ta_ * vtr)) / Te_;

        f[0] = -vr_dot + (-vr + Ta_ * vtr) / Tb_ - vtr;
        f[1] = -efd_dot + Math::antiwindup(efd, func, Efdmin_, Efdmax_);
        f[2] = -vtr - Ec + vref + vs + uel_on_ * vuel + oel_on_ * voel;

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::evaluateResidual()
      {
        // Attached signals are read live; unattached ones keep the latched value.
        auto read_signal = [&]<SexsPtiExternalVariables variable>(const ScalarT& latched)
        {
          const auto index   = static_cast<size_t>(variable);
          ws_[index]         = latched;
          ws_indices_[index] = INVALID_INDEX<IdxT>;
          if (signals_.template isAttached<variable>())
          {
            ws_[index]         = signals_.template readExternalVariable<variable>();
            ws_indices_[index] = signals_.template readExternalVariableIndex<variable>();
          }
        };

        read_signal.template operator()<SexsPtiExternalVariables::VREF>(vref_set_);
        read_signal.template operator()<SexsPtiExternalVariables::VS>(vs_set_);
        read_signal.template operator()<SexsPtiExternalVariables::VUEL>(vuel_set_);
        read_signal.template operator()<SexsPtiExternalVariables::VOEL>(voel_set_);

        wb_[0] = bus_->Vr();
        wb_[1] = bus_->Vi();

        const auto* y  = y_.getData();
        const auto* yp = yp_.getData();
        auto*       f  = f_.getData();
        evaluateInternalResidual(y, yp, wb_.data(), ws_.data(), f);

        f_.setDataUpdated();

        return 0;
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
        monitor_->set(Variable::efd, [this]
                      { return y_.getData()[1]; });
      }

    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
