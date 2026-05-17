/**
 * @file RegcaImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the REGCA phasor-dynamics converter model.
 */

#pragma once

#include <algorithm>
#include <variant>

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/Regca.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/RegcaData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Converter
    {
      using Log = ::GridKit::Utilities::Logger;

      template <class ScalarT, typename IdxT>
      Regca<ScalarT, IdxT>::Regca(bus_type* bus)
        : bus_(bus)
      {
        size_ = static_cast<IdxT>(RegcaInternalVariables::MAXIMUM);
      }

      template <class ScalarT, typename IdxT>
      Regca<ScalarT, IdxT>::Regca(bus_type* bus, const model_data_type& data)
        : bus_(bus),
          monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        initializeMonitor();
        size_ = static_cast<IdxT>(RegcaInternalVariables::MAXIMUM);
      }

      template <class ScalarT, typename IdxT>
      Regca<ScalarT, IdxT>::~Regca()
      {
      }

      template <class ScalarT, typename IdxT>
      void Regca<ScalarT, IdxT>::initializeParameters(const model_data_type& data)
      {
        using Params = typename model_data_type::Parameters;
        using Ports  = typename model_data_type::Ports;

        auto loadReal = [&](auto key, RealT& target)
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
        };

        auto loadBool = [&](auto key, bool& target)
        {
          if (!data.parameters.contains(key))
          {
            return;
          }

          const auto& value = data.parameters.at(key);
          if (const auto* bool_value = std::get_if<bool>(&value))
          {
            target = *bool_value;
          }
          else if (const auto* index_value = std::get_if<IdxT>(&value))
          {
            target = (*index_value != 0);
          }
        };

        loadReal(Params::P0, P0_);
        loadReal(Params::Q0, Q0_);
        loadReal(Params::Sconv, Sconv_);
        loadReal(Params::Tg, Tg_);
        loadReal(Params::TM, TM_);
        loadReal(Params::Rqmax, Rqmax_);
        loadReal(Params::Rqmin, Rqmin_);
        loadReal(Params::Rpmax, Rpmax_);
        loadBool(Params::sL, sL_);
        loadReal(Params::IL1, IL1_);
        loadReal(Params::VL0, VL0_);
        loadReal(Params::VL1, VL1_);
        loadReal(Params::VA0, VA0_);
        loadReal(Params::VA1, VA1_);
        loadReal(Params::Vhvmax, Vhvmax_);

        if (data.ports.contains(Ports::bus))
        {
          bus_id_ = data.ports.at(Ports::bus);
        }
      }

      template <class ScalarT, typename IdxT>
      const Model::VariableMonitorBase* Regca<ScalarT, IdxT>::getMonitor() const
      {
        return monitor_.get();
      }

      template <class ScalarT, typename IdxT>
      void Regca<ScalarT, IdxT>::initializeMonitor()
      {
        using Variable = typename model_data_type::MonitorableVariables;
        auto index     = [](RegcaInternalVariables variable)
        {
          return static_cast<size_t>(variable);
        };

        monitor_->set(Variable::ir, [this, index]
                      { return y_[index(RegcaInternalVariables::IR)]; });
        monitor_->set(Variable::ii, [this, index]
                      { return y_[index(RegcaInternalVariables::II)]; });
        monitor_->set(Variable::p, []
                      { return ScalarT{0}; });
        monitor_->set(Variable::q, []
                      { return ScalarT{0}; });
        monitor_->set(Variable::vt, [this, index]
                      { return y_[index(RegcaInternalVariables::VT)]; });
        monitor_->set(Variable::vm, [this, index]
                      { return y_[index(RegcaInternalVariables::VM)]; });
        monitor_->set(Variable::ip, [this, index]
                      { return y_[index(RegcaInternalVariables::IP)]; });
        monitor_->set(Variable::iq, [this, index]
                      { return y_[index(RegcaInternalVariables::IQ)]; });
        monitor_->set(Variable::iqextra, [this, index]
                      { return y_[index(RegcaInternalVariables::IQEXTRA)]; });
        monitor_->set(Variable::il, [this, index]
                      { return y_[index(RegcaInternalVariables::IL)]; });
        monitor_->set(Variable::lp, [this, index]
                      { return y_[index(RegcaInternalVariables::LP)]; });
        monitor_->set(Variable::up, [this, index]
                      { return y_[index(RegcaInternalVariables::UP)]; });
      }

      template <class ScalarT, typename IdxT>
      int Regca<ScalarT, IdxT>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      template <class ScalarT, typename IdxT>
      int Regca<ScalarT, IdxT>::allocate()
      {
        size_     = static_cast<IdxT>(RegcaInternalVariables::MAXIMUM);
        auto size = static_cast<size_t>(size_);

        f_.assign(size, ScalarT{0});
        y_.assign(size, ScalarT{0});
        yp_.assign(size, ScalarT{0});
        tag_.assign(size, false);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        wb_.assign(2, ScalarT{0});
        h_.assign(2, ScalarT{0});

        auto signal_size = static_cast<size_t>(RegcaExternalVariables::MAXIMUM);
        ws_.assign(signal_size, ScalarT{0});
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        return 0;
      }

      template <class ScalarT, typename IdxT>
      int Regca<ScalarT, IdxT>::verify() const
      {
        int ret = 0;

        if (bus_ == nullptr)
        {
          Log::error() << "Regca: bus pointer is null\n";
          ret += 1;
        }

        if (signals_.template isAttached<RegcaExternalVariables::IPCMD>())
        {
          if (!signals_.template isLinked<RegcaExternalVariables::IPCMD>())
          {
            Log::error() << "Regca: ipcmd signal attached with no linked source\n";
            ret += 1;
          }
        }

        if (signals_.template isAttached<RegcaExternalVariables::IQCMD>())
        {
          if (!signals_.template isLinked<RegcaExternalVariables::IQCMD>())
          {
            Log::error() << "Regca: iqcmd signal attached with no linked source\n";
            ret += 1;
          }
        }

        return ret;
      }

      template <class ScalarT, typename IdxT>
      int Regca<ScalarT, IdxT>::initialize()
      {
        std::fill(y_.begin(), y_.end(), ScalarT{0});
        std::fill(yp_.begin(), yp_.end(), ScalarT{0});
        return 0;
      }

      template <class ScalarT, typename IdxT>
      int Regca<ScalarT, IdxT>::tagDifferentiable()
      {
        std::fill(tag_.begin(), tag_.end(), false);
        tag_[static_cast<size_t>(RegcaInternalVariables::VM)] = true;
        tag_[static_cast<size_t>(RegcaInternalVariables::IQ)] = true;
        tag_[static_cast<size_t>(RegcaInternalVariables::IP)] = true;
        return 0;
      }

      template <class ScalarT, typename IdxT>
      __attribute__((always_inline)) inline int Regca<ScalarT, IdxT>::evaluateInternalResidual(
          [[maybe_unused]] ScalarT* y,
          [[maybe_unused]] ScalarT* yp,
          [[maybe_unused]] ScalarT* wb,
          [[maybe_unused]] ScalarT* ws,
          ScalarT*                  f)
      {
        for (IdxT i = 0; i < size_; ++i)
        {
          f[static_cast<size_t>(i)] = ScalarT{0};
        }

        return 0;
      }

      template <class ScalarT, typename IdxT>
      __attribute__((always_inline)) inline int Regca<ScalarT, IdxT>::evaluateBusResidual(
          [[maybe_unused]] ScalarT* y,
          [[maybe_unused]] ScalarT* yp,
          [[maybe_unused]] ScalarT* wb,
          ScalarT*                  h)
      {
        h[0] = ScalarT{0};
        h[1] = ScalarT{0};
        return 0;
      }

      template <class ScalarT, typename IdxT>
      int Regca<ScalarT, IdxT>::evaluateResidual()
      {
        std::fill(ws_.begin(), ws_.end(), ScalarT{0});
        std::fill(ws_indices_.begin(), ws_indices_.end(), INVALID_INDEX<IdxT>);

        if (signals_.template isAttached<RegcaExternalVariables::IPCMD>())
        {
          const auto index = static_cast<size_t>(RegcaExternalVariables::IPCMD);
          ws_[index]       = signals_.template readExternalVariable<RegcaExternalVariables::IPCMD>();
          ws_indices_[index] =
              signals_.template readExternalVariableIndex<RegcaExternalVariables::IPCMD>();
        }

        if (signals_.template isAttached<RegcaExternalVariables::IQCMD>())
        {
          const auto index = static_cast<size_t>(RegcaExternalVariables::IQCMD);
          ws_[index]       = signals_.template readExternalVariable<RegcaExternalVariables::IQCMD>();
          ws_indices_[index] =
              signals_.template readExternalVariableIndex<RegcaExternalVariables::IQCMD>();
        }

        wb_[0] = Vr();
        wb_[1] = Vi();

        evaluateInternalResidual(y_.data(), yp_.data(), wb_.data(), ws_.data(), f_.data());
        evaluateBusResidual(y_.data(), yp_.data(), wb_.data(), h_.data());

        Ir() += h_[0];
        Ii() += h_[1];

        return 0;
      }
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
