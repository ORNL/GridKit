#pragma once

#include <algorithm>
#include <bitset>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include <rapidjson/reader.h>

#include <Model/PhasorDynamics/Branch/BranchData.hpp>
#include <Model/PhasorDynamics/Bus/BusData.hpp>
#include <Model/PhasorDynamics/BusFault/BusFaultData.hpp>
#include <Model/PhasorDynamics/Load/LoadData.hpp>
#include <Model/PhasorDynamics/SynchronousMachine/GENROUwS/GenrouData.hpp>
#include <Model/PhasorDynamics/SystemModelData.hpp>

namespace GridKit
{
  using namespace rapidjson;

  /// RapidJSON handler for the Grid Dynamics case format
  template <typename RealT = double, typename IdxT = size_t>
  struct GridDynamicsFormatHandler
    : public BaseReaderHandler<UTF8<>,
                               GridDynamicsFormatHandler<RealT,
                                                         IdxT>>
  {
    using SystemModelData =
        GridKit::PhasorDynamics::SystemModelData<RealT, IdxT>;
    using BranchDataT   = PhasorDynamics::BranchData<RealT, IdxT>;
    using BusDataT      = PhasorDynamics::BusData<RealT, IdxT>;
    using BusFaultDataT = PhasorDynamics::BusFaultData<RealT, IdxT>;
    using GenrouDataT   = PhasorDynamics::GenrouData<RealT, IdxT>;
    using LoadDataT     = PhasorDynamics::LoadData<RealT, IdxT>;

    /// Default constructor for this handler. Initializes the
    /// class such that it expects the outer object of the
    /// case format
    GridDynamicsFormatHandler() : state_(State::ExpectOuterObject)
    {
    }

    bool Uint(unsigned u)
    {
      switch (state_)
      {
      case State::ExpectFormatVersion:
        // technically we should work with this and change parsing of the
        // rest of the file based on this versioning but that's not
        // necessary just yet
        system_model.format_version = u;
        state_                      = State::ExpectHeaderKeyOrObjectClose;
        break;
      case State::ExpectFormatRevision:
        // likewise here
        system_model.format_revision = u;
        state_                       = State::ExpectHeaderKeyOrObjectClose;
        break;
      case State::ExpectBusNumber:
        bus_data.bus_id = u;
        state_          = State::ExpectBusKeyOrObjectClose;
        break;
      default:
        return false;
      }
      return true;
    }

    bool String(const char* str, SizeType length, bool)
    {
      auto s = std::string(str, length);
      switch (state_)
      {
      case State::ExpectCaseName:
        system_model.case_name = s;
        state_                 = State::ExpectHeaderKeyOrObjectClose;
        break;
      case State::ExpectCaseDateTime:
        // this should be validated to be iso 8601 -- see the comment in the
        // system model data structure
        system_model.case_date_time = s;
        state_                      = State::ExpectHeaderKeyOrObjectClose;
        break;
      case State::ExpectCaseDescription:
        system_model.case_description = s;
        state_                        = State::ExpectHeaderKeyOrObjectClose;
        break;
      case State::ExpectCaseComments:
        system_model.case_comments = s;
        state_                     = State::ExpectHeaderKeyOrObjectClose;
        break;
      case State::ExpectBusClass:
        if (s == "bus")
        {
          bus_data.bus_type = BusDataT::BusType::Default;
        }
        else if (s == "infinite_bus")
        {
          bus_data.bus_type = BusDataT::BusType::Slack;
        }
        else
        {
          return false;
        }
        state_ = State::ExpectBusKeyOrObjectClose;
        break;
      case State::ExpectBusName:
        bus_data.name = s;
        state_        = State::ExpectBusKeyOrObjectClose;
        break;
      case State::ExpectMonitoredBusVariableOrArrayClose:
        if (s == "Vr")
        {
          bus_data.monitored_variables.set(BusDataT::MonitorableVariables::Vr);
        }
        else if (s == "Vi")
        {
          bus_data.monitored_variables.set(BusDataT::MonitorableVariables::Vi);
        }
        else if (s == "Vm")
        {
          bus_data.monitored_variables.set(BusDataT::MonitorableVariables::Vm);
        }
        else if (s == "Va")
        {
          bus_data.monitored_variables.set(BusDataT::MonitorableVariables::Va);
        }
        else
        {
          return false;
        }
      default:
        return false;
      }
      return true;
    }

    bool Double(double d)
    {
      switch (state_)
      {
      case State::ExpectFreqBase:
        system_model.freq_base = d;
        state_                 = State::ExpectHeaderKeyOrObjectClose;
        break;
      case State::ExpectVaBase:
        system_model.va_base = d;
        state_               = State::ExpectHeaderKeyOrObjectClose;
        break;
      case State::ExpectBusInitialParameterVr:
        bus_data.Vr0 = d;
        state_       = State::ExpectBusInitialParameterKeyOrObjectClose;
        break;
      case State::ExpectBusInitialParameterVi:
        bus_data.Vi0 = d;
        state_       = State::ExpectBusInitialParameterKeyOrObjectClose;
        break;
      case State::ExpectBusVBase:
        // TODO: set this. there doesn't seem to be a field on the bus data
        //       structure for this
        state_ = State::ExpectBusKeyOrObjectClose;
        break;
      case State::ExpectBusFreqBase:
        bus_data.freq_base = d;
        state_             = State::ExpectBusKeyOrObjectClose;
        break;
      case State::ExpectBusVaBase:
        bus_data.va_base = d;
        state_           = State::ExpectBusKeyOrObjectClose;
        break;
      default:
        return false;
      }
      return true;
    }

    bool Key(const char* str, SizeType length, bool)
    {
      auto key = std::string_view(str, length);
      switch (state_)
      {
      case State::ExpectInnerKey:
        if (key == "header")
        {
          state_ = State::ExpectHeader;
        }
        else if (key == "buses")
        {
          state_ = State::ExpectBuses;
        }
        else if (key == "devices")
        {
          state_ = State::ExpectDevices;
        }
        else
        {
          return false;
        }
        break;
      case State::ExpectHeaderKeyOrObjectClose:
        if (key == "format_version")
        {
          state_ = State::ExpectFormatVersion;
        }
        else if (key == "format_revision")
        {
          state_ = State::ExpectFormatRevision;
        }
        else if (key == "case_name")
        {
          state_ = State::ExpectCaseName;
        }
        else if (key == "case_date_time")
        {
          state_ = State::ExpectCaseDateTime;
        }
        else if (key == "case_description")
        {
          state_ = State::ExpectCaseDescription;
        }
        else if (key == "case_comments")
        {
          state_ = State::ExpectCaseComments;
        }
        else if (key == "freq_base")
        {
          state_ = State::ExpectFreqBase;
        }
        else if (key == "va_base")
        {
          state_ = State::ExpectVaBase;
        }
        else
        {
          return false;
        }
        break;
      case State::ExpectBusKeyOrObjectClose:
        if (key == "number")
        {
          state_ = State::ExpectBusNumber;
        }
        else if (key == "class")
        {
          state_ = State::ExpectBusClass;
        }
        else if (key == "name")
        {
          state_ = State::ExpectBusName;
        }
        else if (key == "init")
        {
          state_ = State::ExpectBusInitialParameters;
        }
        else if (key == "v_base")
        {
          state_ = State::ExpectBusVBase;
        }
        else if (key == "mon")
        {
          state_ = State::ExpectBusMonitor;
        }
        else if (key == "freq_base")
        {
          state_ = State::ExpectBusFreqBase;
        }
        else if (key == "va_base")
        {
          state_ = State::ExpectBusVaBase;
        }
        else if (key == "extension")
        {
          state_ = State::ExpectBusExtension;
        }
        else
        {
          return false;
        }
        break;
      case State::ExpectDeviceKeyOrObjectClose:
        if (key == "class")
        {
          state_ = State::ExpectDeviceClass;
        }
        else if (key == "ports")
        {
          state_ = State::ExpectDevicePorts;
        }
        else if (key == "id")
        {
          state_ = State::ExpectDeviceId;
        }
        else if (key == "params")
        {
          state_ = State::ExpectDeviceInitialParameters;
        }
        else if (key == "mon")
        {
          state_ = State::ExpectDeviceMonitor;
        }
        else if (key == "va_base")
        {
          state_ = State::ExpectDeviceVaBase;
        }
        else if (key == "freq_base")
        {
          state_ = State::ExpectDeviceFreqBase;
        }
        else if (key == "extension")
        {
          state_ = State::ExpectDeviceExtension;
        }
        else
        {
          return false;
        }
        break;
      case State::ExpectBusInitialParameterKeyOrObjectClose:
        // NOTE: these are hardcoded for now because there are so few of them.
        if (key == "Vr")
        {
          state_ = State::ExpectBusInitialParameterVr;
        }
        else if (key == "Vi")
        {
          state_ = State::ExpectBusInitialParameterVi;
        }
        else
        {
          return false;
        }
        break;
      default:
        return false;
      }
      return true;
    }

    bool StartObject()
    {
      // TODO: we should probably erase prior state information here when
      //       a new object is entered that needs to store certain
      //       state information
      switch (state_)
      {
      case State::ExpectHeader:
        state_ = State::ExpectHeaderKeyOrObjectClose;
        break;
      case State::ExpectOuterObject:
        state_ = State::ExpectInnerKey;
        break;
      case State::ExpectBusOrArrayClose:
        state_   = State::ExpectBusKeyOrObjectClose;
        bus_data = BusDataT();
        break;
      case State::ExpectDeviceOrArrayClose:
        state_ = State::ExpectDeviceKeyOrObjectClose;
        device_input_parameters.clear();
        device_ports.clear();
        device_id.clear();
        device_va_base.reset();
        device_freq_base.reset();
        monitored_device_variables.reset();
        break;
      case State::ExpectBusInitialParameters:
        state_ = State::ExpectBusInitialParameterKeyOrObjectClose;
        break;
      case State::ExpectBusExtension:
        state_ = State::IgnoreUntilExtensionObjectClose;
        break;
      default:
        return false;
      }

      return true;
    }

    bool EndObject(size_t length)
    {
      // TODO: validate the length of the object
      switch (state_)
      {
      case State::ExpectBusKeyOrObjectClose:
        system_model.bus.push_back(bus_data);
        state_ = State::ExpectBusOrArrayClose;
        break;
      case State::ExpectDeviceKeyOrObjectClose:
        state_ = State::ExpectDeviceOrArrayClose;
        break;
      case State::ExpectBusInitialParameterKeyOrObjectClose:
        state_ = State::ExpectBusKeyOrObjectClose;
        break;
      case State::IgnoreUntilExtensionObjectClose:
        state_ = State::ExpectBusKeyOrObjectClose;
        break;
      default:
        return false;
      }

      return true;
    }

    bool StartArray()
    {
      switch (state_)
      {
      case State::ExpectBuses:
        state_ = State::ExpectBusOrArrayClose;
        break;
      case State::ExpectDevices:
        state_ = State::ExpectDeviceOrArrayClose;
        break;
      case State::ExpectBusMonitor:
        state_ = State::ExpectMonitoredBusVariableOrArrayClose;
        break;
      case State::ExpectDeviceMonitor:
        state_ = State::ExpectMonitoredDeviceVariableOrArrayClose;
        break;
      default:
        return false;
      }

      return true;
    }

    bool EndArray(size_t length)
    {
      // TODO: validate the length of the array
      switch (state_)
      {
      case State::ExpectBusOrArrayClose:
      case State::ExpectDeviceOrArrayClose:
        state_ = State::ExpectInnerKey;
        break;
      case State::ExpectMonitoredBusVariableOrArrayClose:
        state_ = State::ExpectBusKeyOrObjectClose;
        break;
      case State::ExpectMonitoredDeviceVariableOrArrayClose:
        state_ = State::ExpectDeviceKeyOrObjectClose;
        break;
      default:
        return false;
      }

      return true;
    }

    /// Default handler for JSON elements
    bool Default()
    {
      return false;
    }

    /// Bus data structure used to hold information during construction
    BusDataT bus_data;

    /// The internal `SystemModelData` being constructed
    SystemModelData system_model;

    /// Enumeration used to indicate the kind of parameter being
    /// set for a bus
    enum class BusParameter : size_t
    {
      Vr,
      Vi,
      Maximum,
    };

    /// Enumeration used to indicate the kind of port being set
    /// for a device
    enum class DevicePort : size_t
    {
      Bus1,
      Bus2,
      Bus,
      ExciterSignal,
      GovernorSignal,
      ControlSignal,
      Maximum,
    };

    /// Enumeration used to indicate the kind of parameter
    /// being set for a device
    enum class DeviceParameter : size_t
    {
      R,
      X,
      G,
      B,
      Pz,
      Qz,
      Pi,
      Qi,
      Pp,
      Qp,
      P0,
      Q0,
      H,
      D,
      Ra,
      Tdop,
      Tdopp,
      Tqopp,
      Tqop,
      Xd,
      Xdp,
      Xdpp,
      Xq,
      Xqp,
      Xqpp,
      Xl,
      S10,
      S12,
      state0,
    };

    /// Enumeration used to assign indices of variables able to be monitored
    /// in the bitset
    enum class MonitorableDeviceVariables : size_t
    {
      Ir1,
      Ii1,
      Im1,
      P1,
      Q1,
      Ir2,
      Ii2,
      Im2,
      P2,
      Q2,
      State,
      P,
      Q,
      Delta,
      Omega,
      Maximum,
    };

    /// Storage for input parameters to devices
    std::vector<std::pair<DeviceParameter, RealT>>
        device_input_parameters;

    /// Storage for ports of devices
    std::vector<std::pair<DevicePort, IdxT>>
        device_ports;

    /// Storage for the ID of devices
    std::string device_id;

    /// Storage for the system power base override of devices
    std::optional<RealT> device_va_base;

    /// Storage for the system frequency base override
    /// of devices
    std::optional<RealT> device_freq_base;

    /// Bitfield for tracking the device variables being monitored
    std::bitset<static_cast<std::underlying_type_t<MonitorableDeviceVariables>>(MonitorableDeviceVariables::Maximum) - 1> monitored_device_variables;

    /// Enumeration used to track the current state of
    /// the parser
    enum class State : size_t
    {
      ExpectOuterObject,
      ExpectInnerKey,
      ExpectHeader,
      ExpectHeaderKeyOrObjectClose,
      ExpectFormatVersion,
      ExpectFormatRevision,
      ExpectCaseName,
      ExpectCaseDateTime,
      ExpectCaseDescription,
      ExpectCaseComments,
      ExpectFreqBase,
      ExpectVaBase,
      ExpectBuses,
      ExpectBusOrArrayClose,
      ExpectBusKeyOrObjectClose,
      ExpectBusNumber,
      ExpectBusClass,
      ExpectBusName,
      ExpectBusInitialParameters,
      ExpectBusInitialParameterKeyOrObjectClose,
      ExpectBusInitialParameterVr,
      ExpectBusInitialParameterVi,
      ExpectBusVBase,
      ExpectBusMonitor,
      ExpectMonitoredBusVariableOrArrayClose,
      ExpectBusFreqBase,
      ExpectBusVaBase,
      ExpectBusExtension,
      ExpectDevices,
      ExpectDeviceOrArrayClose,
      ExpectDeviceKeyOrObjectClose,
      ExpectDeviceClass,
      ExpectDevicePorts,
      ExpectDeviceId,
      ExpectDeviceInitialParameters,
      ExpectDeviceMonitor,
      ExpectMonitoredDeviceVariableOrArrayClose,
      ExpectDeviceVaBase,
      ExpectDeviceFreqBase,
      ExpectDeviceExtension,
      IgnoreUntilExtensionObjectClose,
    } state_;
  };
} // namespace GridKit
