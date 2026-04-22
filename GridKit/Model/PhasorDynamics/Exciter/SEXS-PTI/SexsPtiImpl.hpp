#pragma once

/**
 * @file SexsPtiImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the SEXS-PTI exciter model.
 */

#include <cmath>
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

      template <class ScalarT, typename IdxT>
      SexsPti<ScalarT, IdxT>::SexsPti(bus_type* bus)
        : bus_(bus)
      {
        size_ = 3;
      }

      template <class ScalarT, typename IdxT>
      SexsPti<ScalarT, IdxT>::SexsPti(bus_type*              bus,
                                      const model_data_type& data)
        : bus_(bus),
          monitor_(std::make_unique<MonitorT>(data))
      {
        initModelParams(data);
        initializeMonitor();
        size_ = 3;
      }

      template <class ScalarT, typename IdxT>
      SexsPti<ScalarT, IdxT>::~SexsPti()
      {
      }

      template <class ScalarT, typename IdxT>
      int SexsPti<ScalarT, IdxT>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      template <class ScalarT, typename IdxT>
      int SexsPti<ScalarT, IdxT>::allocate()
      {
        auto size = static_cast<size_t>(size_);
        f_.resize(size);
        y_.resize(size);
        yp_.resize(size);
        tag_.resize(size);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        wb_.resize(2);

        ws_.resize(1);
        ws_indices_.resize(1);
        ws_[0]         = 0.0;
        ws_indices_[0] = INVALID_INDEX<IdxT>;

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        if (signals_.template isAssigned<SexsPtiInternalVariables::EFD>())
        {
          signals_.template getSignalNode<SexsPtiInternalVariables::EFD>()->set(
              &y_[1], &(this->getVariableIndex(1)));
        }

        return 0;
      }

      template <class ScalarT, typename IdxT>
      int SexsPti<ScalarT, IdxT>::verify() const
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

        if (signals_.template isAttached<SexsPtiExternalVariables::VS>())
        {
          if (!signals_.template isLinked<SexsPtiExternalVariables::VS>())
          {
            Log::error() << "SexsPti: VS signal attached with no linked source\n";
            ret += 1;
          }
        }

        return ret;
      }

      template <class ScalarT, typename IdxT>
      int SexsPti<ScalarT, IdxT>::initialize()
      {
        ScalarT efd0{0.0};
        if (signals_.template isAssigned<SexsPtiInternalVariables::EFD>())
        {
          efd0 = y_[1];
        }

        ScalarT vreal = bus_->Vr();
        ScalarT vimag = bus_->Vi();
        ScalarT Ec    = std::sqrt(vreal * vreal + vimag * vimag);
        ScalarT vtr   = efd0 / K_;
        ScalarT vr    = (Ta_ - Tb_) * vtr;

        vref_ = Ec + vtr;

        y_[0] = vr;
        y_[1] = efd0;
        y_[2] = vtr;

        for (IdxT i = 0; i < size_; ++i)
        {
          yp_[static_cast<size_t>(i)] = 0.0;
        }

        return 0;
      }

      template <class ScalarT, typename IdxT>
      int SexsPti<ScalarT, IdxT>::tagDifferentiable()
      {
        tag_[0] = true;
        tag_[1] = true;
        tag_[2] = false;

        return 0;
      }

      template <class ScalarT, typename IdxT>
      __attribute__((always_inline)) inline int SexsPti<ScalarT, IdxT>::evaluateInternalResidual(
          ScalarT* y,
          ScalarT* yp,
          ScalarT* wb,
          ScalarT* ws,
          ScalarT* f)
      {
        ScalarT vr      = y[0];
        ScalarT efd     = y[1];
        ScalarT vtr     = y[2];
        ScalarT vr_dot  = yp[0];
        ScalarT efd_dot = yp[1];

        ScalarT Ec = std::sqrt(wb[0] * wb[0] + wb[1] * wb[1]);
        ScalarT vs = ws[0];

        ScalarT g           = -efd + (K_ / Tb_) * (-vr + Ta_ * vtr);
        ScalarT g_direction = g / std::sqrt(ONE<RealT> + g * g);
        ScalarT efd_ind     = Math::indicator(Efdmin_, Efdmax_, efd, g_direction);

        f[0] = -vr_dot + (-vr + Ta_ * vtr) / Tb_ - vtr;
        f[1] = -efd_dot + efd_ind * g / Te_;
        f[2] = -vtr - Ec + vref_ + vs + vOEL_ + vUEL_;

        return 0;
      }

      template <class ScalarT, typename IdxT>
      int SexsPti<ScalarT, IdxT>::evaluateResidual()
      {
        ws_[0]         = 0.0;
        ws_indices_[0] = INVALID_INDEX<IdxT>;
        if (signals_.template isAttached<SexsPtiExternalVariables::VS>())
        {
          ws_[0]         = signals_.template readExternalVariable<SexsPtiExternalVariables::VS>();
          ws_indices_[0] = signals_.template readExternalVariableIndex<SexsPtiExternalVariables::VS>();
        }

        wb_[0] = bus_->Vr();
        wb_[1] = bus_->Vi();

        evaluateInternalResidual(y_.data(), yp_.data(), wb_.data(), ws_.data(), f_.data());

        return 0;
      }

      template <class ScalarT, typename IdxT>
      void SexsPti<ScalarT, IdxT>::initModelParams(const model_data_type& data)
      {
        using Params = typename model_data_type::Parameters;

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

      template <class ScalarT, typename IdxT>
      const Model::VariableMonitorBase* SexsPti<ScalarT, IdxT>::getMonitor() const
      {
        return monitor_.get();
      }

      template <class ScalarT, typename IdxT>
      void SexsPti<ScalarT, IdxT>::initializeMonitor()
      {
        using Variable = typename model_data_type::MonitorableVariables;
        monitor_->set(Variable::efd, [this]
                      { return y_[1]; });
      }

    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
