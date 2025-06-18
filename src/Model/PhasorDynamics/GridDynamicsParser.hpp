#pragma once

#include <algorithm>
#include <bitset>
#include <iostream>
#include <optional>
#include <string>
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
  class GridDynamicsFormatHandler
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

    bool Key(const char* str, SizeType length, bool)
    {
      auto key = std::string(str, length);
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
        state_ = State::ExpectBusKeyOrObjectClose;
        current_structure_.emplace(BusDataT());
        monitoring_parameter_bitfield_.reset();
        input_parameter_and_port_bitfield_.reset();
        break;
      case State::ExpectDeviceOrArrayClose:
        state_ = State::ExpectDeviceKeyOrObjectClose;
        monitoring_parameter_bitfield_.reset();
        input_parameter_and_port_bitfield_.reset();
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
        state_ = State::ExpectBusOrArrayClose;
        break;
      case State::ExpectDeviceKeyOrObjectClose:
        state_ = State::ExpectDeviceOrArrayClose;
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
        // TODO
        break;
      case State::ExpectDeviceMonitor:
        // TODO
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

    /// Type-safe union of the values that are constructed
    std::variant<BusDataT,
                 BranchDataT,
                 BusFaultDataT,
                 GenrouDataT,
                 LoadDataT>
        current_structure_;

    /// The internal [`SystemModelData`] being constructed
    SystemModelData system_model_;

    /// Enumeration used to indicate the kind of parameter being
    /// set for a bus
    enum class BusParameter : size_t
    {
      Vr,
      Vi,
      Va,
      Vb,
      Vc,
      X,
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

    /// Storage for input parameters to devices
    std::vector<std::pair<DeviceParameter, RealT>>
        device_input_parameters_;

    /// Storage for ports of devices
    std::vector<std::pair<DevicePort, IdxT>>
        device_ports_;

    /// Storage for the ID of devices
    std::string device_id_;

    /// Storage for the system power base override of devices
    std::optional<RealT> device_va_base_;

    /// Storage for the system frequency base override
    /// of devices
    std::optional<RealT> device_freq_base_;

    /// Bitfield for ensuring no device ports are repeated nor
    /// any bus or device input parameters
    std::bitset<std::max(static_cast<std::underlying_type_t<DevicePort>>(DevicePort::Maximum) + static_cast<std::underlying_type_t<DeviceParameter>>(DeviceParameter::Maximum) - 2,
                         static_cast<std::underlying_type_t<BusParameter>>(BusParameter::Maximum) - 1)>
        input_parameter_and_port_bitfield_;

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
      ExpectBusVBase,
      ExpectBusMonitor,
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
      ExpectDeviceVaBase,
      ExpectDeviceFreqBase,
      ExpectDeviceExtension,
    } state_;
  };
} // namespace GridKit
