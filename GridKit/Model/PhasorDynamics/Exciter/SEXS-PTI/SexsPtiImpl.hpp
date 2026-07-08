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
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeSet.hpp>
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
      SexsPti<scalar_type, index_type>::SexsPti(BusT* bus, const ModelDataT& data)
        : bus_(bus),
          monitor_(std::make_unique<MonitorT>(data))
      {
        initModelParams(data);
        initializeMonitor();
        size_ = 3;
      }

      template <typename scalar_type, typename index_type>
      SexsPti<scalar_type, index_type>::SexsPti(
          BusT*             bus,
          const ModelDataT& data,
          SignalNodeSetT&   signal_nodes)
        : bus_(bus),
          ports_(data, signal_nodes),
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
        auto size = static_cast<size_t>(size_);
        f_.resize(size);
        y_.resize(size);
        yp_.resize(size);
        tag_.resize(size);
        abs_tol_.resize(size);
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

        if (auto efd_port = ports_.out[SexsPtiSignalOutputs::efd])
        {
          efd_port.link(&y_[1], &(this->getVariableIndex(1)));
        }

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

        auto efd_port = ports_.out[SexsPtiSignalOutputs::efd];
        if (!efd_port.connected() || !efd_port.linked())
        {
          Log::error() << "SexsPti: required EFD signal is not assigned\n";
          ret += 1;
        }

        auto vs_port = ports_.in[SexsPtiSignalInputs::vs];
        if (vs_port.connected() && !vs_port.linked())
        {
          Log::error() << "SexsPti: VS signal attached with no linked source\n";
          ret += 1;
        }

        return ret;
      }

      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::initialize()
      {
        ScalarT efd0{0.0};
        if (ports_.out[SexsPtiSignalOutputs::efd])
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
        std::fill(abs_tol_.begin(), abs_tol_.end(), rel_tol);
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
        ScalarT vr      = y[0];
        ScalarT efd     = y[1];
        ScalarT vtr     = y[2];
        ScalarT vr_dot  = yp[0];
        ScalarT efd_dot = yp[1];

        ScalarT Ec = std::sqrt(wb[0] * wb[0] + wb[1] * wb[1]);
        ScalarT vs = ws[0];

        ScalarT func = (-efd + (K_ / Tb_) * (-vr + Ta_ * vtr)) / Te_;

        f[0] = -vr_dot + (-vr + Ta_ * vtr) / Tb_ - vtr;
        f[1] = -efd_dot + Math::antiwindup(efd, func, Efdmin_, Efdmax_);
        f[2] = -vtr - Ec + vref_ + vs + vOEL_ + vUEL_;

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::evaluateResidual()
      {
        ws_[0]         = 0.0;
        ws_indices_[0] = INVALID_INDEX<IdxT>;
        if (auto vs_port = ports_.in[SexsPtiSignalInputs::vs])
        {
          ws_[0]         = vs_port.readSignal();
          ws_indices_[0] = vs_port.signalVariableIndex();
        }

        wb_[0] = bus_->Vr();
        wb_[1] = bus_->Vi();

        evaluateInternalResidual(y_.data(), yp_.data(), wb_.data(), ws_.data(), f_.data());

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
                      { return y_[1]; });
      }

    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
