#pragma once

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename ScalarT, typename IdxT>
    BusBase<ScalarT, IdxT>::BusBase(const BusData<RealT, IdxT>& data)
      : bus_id_(data.bus_id),
        monitor_(std::make_unique<MonitorT>("Bus_" + data.name, data.monitored_variables))
    {
      using Variable = typename BusData<RealT, IdxT>::MonitorableVariables;
      monitor_->set(Variable::Vr, [this]
                    { return Vr(); });
      monitor_->set(Variable::Vi, [this]
                    { return Vi(); });
      monitor_->set(Variable::Vm, [this]
                    { return std::sqrt(Vr() * Vr() + Vi() * Vi()); });
      monitor_->set(Variable::Va, [this]
                    { return std::atan2(Vi(), Vr()); });
    }

    template <typename ScalarT, typename IdxT>
    BusBase<ScalarT, IdxT>::~BusBase()
    {
      if (J_rows_buffer_ != nullptr)
      {
        delete[] J_rows_buffer_;
        delete[] J_cols_buffer_;
        delete[] J_vals_buffer_;
        J_rows_buffer_ = nullptr;
        J_cols_buffer_ = nullptr;
        J_vals_buffer_ = nullptr;
      }
    }

    template <typename ScalarT, typename IdxT>
    inline const Model::VariableMonitorBase* BusBase<ScalarT, IdxT>::getMonitor() const
    {
      return monitor_.get();
    }
  } // namespace PhasorDynamics
} // namespace GridKit
