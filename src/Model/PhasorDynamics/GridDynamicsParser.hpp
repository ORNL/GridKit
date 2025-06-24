#pragma once

#include <algorithm>
#include <bitset>
#include <cstdint>
#include <iostream>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
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
    GridDynamicsFormatHandler()
      : state(State::ExpectOuterObject)
    {
    }

    bool Bool(bool b)
    {
      switch (state)
      {
      case State::ExpectDeviceInitialParameterMapping:
        if (staged_device_parameter != DeviceParameter::state0)
        {
          return false;
        }
        device_input_parameters.insert({staged_device_parameter, b});
        state = State::ExpectDeviceInitialParameterKeyOrObjectClose;
        break;
      default:
        return false;
      }
      return true;
    }

    bool Uint64(std::uint64_t i)
    {
      switch (state)
      {
      case State::ExpectBusNumber:
        bus_data.bus_id = static_cast<IdxT>(i);
        state           = State::ExpectBusKeyOrObjectClose;
        break;
      case State::ExpectDevicePortMapping:
        device_ports.insert({staged_device_port,
                             static_cast<IdxT>(i)});
        state = State::ExpectDevicePortKeyOrObjectClose;
        break;
      case State::ExpectDeviceInitialParameterMapping:
        if (staged_device_parameter != DeviceParameter::UnitId)
        {
          return false;
        }
        device_input_parameters.insert({staged_device_parameter, static_cast<IdxT>(i)});
        state = State::ExpectDeviceInitialParameterKeyOrObjectClose;
        break;
      default:
        return false;
      }
      return true;
    }

    bool Uint(unsigned i)
    {
      switch (state)
      {
      case State::ExpectFormatVersion:
        // technically we should work with this and change parsing of the
        // rest of the file based on this versioning but that's not
        // necessary just yet
        system_model.format_version = i;
        state                       = State::ExpectHeaderKeyOrObjectClose;
        break;
      case State::ExpectFormatRevision:
        // likewise here
        system_model.format_revision = i;
        state                        = State::ExpectHeaderKeyOrObjectClose;
        break;
      case State::ExpectBusNumber:
        bus_data.bus_id = static_cast<IdxT>(i);
        state           = State::ExpectBusKeyOrObjectClose;
        break;
      case State::ExpectDevicePortMapping:
        device_ports.insert({staged_device_port,
                             static_cast<IdxT>(i)});
        state = State::ExpectDevicePortKeyOrObjectClose;
        break;
      case State::ExpectDeviceInitialParameterMapping:
        if (staged_device_parameter != DeviceParameter::UnitId)
        {
          return false;
        }
        device_input_parameters.insert({staged_device_parameter, static_cast<IdxT>(i)});
        state = State::ExpectDeviceInitialParameterKeyOrObjectClose;
        break;
      default:
        return false;
      }
      return true;
    }

    bool String(const char* str, SizeType length, bool)
    {
      auto s = std::string(str, length);
      switch (state)
      {
      case State::ExpectCaseName:
        system_model.case_name = s;
        state                  = State::ExpectHeaderKeyOrObjectClose;
        break;
      case State::ExpectCaseDateTime:
        // this should be validated to be iso 8601 -- see the comment in the
        // system model data structure
        system_model.case_date_time = s;
        state                       = State::ExpectHeaderKeyOrObjectClose;
        break;
      case State::ExpectCaseDescription:
        system_model.case_description = s;
        state                         = State::ExpectHeaderKeyOrObjectClose;
        break;
      case State::ExpectCaseComments:
        system_model.case_comments = s;
        state                      = State::ExpectHeaderKeyOrObjectClose;
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
        state = State::ExpectBusKeyOrObjectClose;
        break;
      case State::ExpectBusName:
        bus_data.name = s;
        state         = State::ExpectBusKeyOrObjectClose;
        break;
      case State::ExpectMonitoredBusVariableOrArrayClose:
        if (s == "Vr")
        {
          bus_data.monitored_variables.set(static_cast<size_t>(BusDataT::MonitorableVariables::Vr));
        }
        else if (s == "Vi")
        {
          bus_data.monitored_variables.set(static_cast<size_t>(BusDataT::MonitorableVariables::Vi));
        }
        else if (s == "Vm")
        {
          bus_data.monitored_variables.set(static_cast<size_t>(BusDataT::MonitorableVariables::Vm));
        }
        else if (s == "Va")
        {
          bus_data.monitored_variables.set(static_cast<size_t>(BusDataT::MonitorableVariables::Va));
        }
        else
        {
          return false;
        }
        break;
      case State::ExpectDeviceClass:
        device_class = s;
        state        = State::ExpectDeviceKeyOrObjectClose;
        break;
      case State::ExpectDeviceId:
        device_id = s;
        state     = State::ExpectDeviceKeyOrObjectClose;
        break;
      case State::ExpectMonitoredDeviceVariableOrArrayClose:
        if (s == "ir1")
        {
          monitored_device_variables.set(static_cast<size_t>(MonitorableDeviceVariables::Ir1));
        }
        else if (s == "ii1")
        {
          monitored_device_variables.set(static_cast<size_t>(MonitorableDeviceVariables::Ii1));
        }
        else if (s == "im1")
        {
          monitored_device_variables.set(static_cast<size_t>(MonitorableDeviceVariables::Im1));
        }
        else if (s == "p1")
        {
          monitored_device_variables.set(static_cast<size_t>(MonitorableDeviceVariables::P1));
        }
        else if (s == "q1")
        {
          monitored_device_variables.set(static_cast<size_t>(MonitorableDeviceVariables::Q1));
        }
        else if (s == "ir2")
        {
          monitored_device_variables.set(static_cast<size_t>(MonitorableDeviceVariables::Ir2));
        }
        else if (s == "ii2")
        {
          monitored_device_variables.set(static_cast<size_t>(MonitorableDeviceVariables::Ii2));
        }
        else if (s == "im2")
        {
          monitored_device_variables.set(static_cast<size_t>(MonitorableDeviceVariables::Im2));
        }
        else if (s == "p2")
        {
          monitored_device_variables.set(static_cast<size_t>(MonitorableDeviceVariables::P2));
        }
        else if (s == "q2")
        {
          monitored_device_variables.set(static_cast<size_t>(MonitorableDeviceVariables::Q2));
        }
        else if (s == "state")
        {
          monitored_device_variables.set(static_cast<size_t>(MonitorableDeviceVariables::State));
        }
        else if (s == "p")
        {
          monitored_device_variables.set(static_cast<size_t>(MonitorableDeviceVariables::P));
        }
        else if (s == "q")
        {
          monitored_device_variables.set(static_cast<size_t>(MonitorableDeviceVariables::Q));
        }
        else if (s == "ir")
        {
          monitored_device_variables.set(static_cast<size_t>(MonitorableDeviceVariables::Ir));
        }
        else if (s == "ii")
        {
          monitored_device_variables.set(static_cast<size_t>(MonitorableDeviceVariables::Ii));
        }
        else if (s == "delta")
        {
          monitored_device_variables.set(static_cast<size_t>(MonitorableDeviceVariables::Delta));
        }
        else if (s == "omega")
        {
          monitored_device_variables.set(static_cast<size_t>(MonitorableDeviceVariables::Omega));
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

    bool Double(double d)
    {
      RealT r = static_cast<RealT>(d);
      switch (state)
      {
      case State::ExpectFreqBase:
        system_model.freq_base = r;
        state                  = State::ExpectHeaderKeyOrObjectClose;
        break;
      case State::ExpectVaBase:
        system_model.va_base = r;
        state                = State::ExpectHeaderKeyOrObjectClose;
        break;
      case State::ExpectBusInitialParameterVr:
        bus_data.Vr0 = r;
        state        = State::ExpectBusInitialParameterKeyOrObjectClose;
        break;
      case State::ExpectBusInitialParameterVi:
        bus_data.Vi0 = r;
        state        = State::ExpectBusInitialParameterKeyOrObjectClose;
        break;
      case State::ExpectBusVBase:
        bus_data.v_base = r;
        state           = State::ExpectBusKeyOrObjectClose;
        break;
      case State::ExpectBusFreqBase:
        bus_data.freq_base = r;
        state              = State::ExpectBusKeyOrObjectClose;
        break;
      case State::ExpectBusVaBase:
        bus_data.va_base = r;
        state            = State::ExpectBusKeyOrObjectClose;
        break;
      case State::ExpectDeviceInitialParameterMapping:
        if (staged_device_parameter == DeviceParameter::UnitId
            || staged_device_parameter == DeviceParameter::state0)
        {
          return false;
        }
        device_input_parameters.insert({staged_device_parameter, r});
        state = State::ExpectDeviceInitialParameterKeyOrObjectClose;
        break;
      case State::ExpectDeviceVaBase:
        device_va_base = r;
        state          = State::ExpectDeviceKeyOrObjectClose;
        break;
      case State::ExpectDeviceFreqBase:
        device_freq_base = r;
        state            = State::ExpectDeviceKeyOrObjectClose;
        break;
      default:
        return false;
      }
      return true;
    }

    bool Key(const char* str, SizeType length, bool)
    {
      auto key = std::string_view(str, length);
      switch (state)
      {
      case State::ExpectInnerKey:
        if (key == "header")
        {
          state = State::ExpectHeader;
        }
        else if (key == "buses")
        {
          state = State::ExpectBuses;
        }
        else if (key == "devices")
        {
          state = State::ExpectDevices;
        }
        else
        {
          return false;
        }
        break;
      case State::ExpectHeaderKeyOrObjectClose:
        if (key == "format_version")
        {
          state = State::ExpectFormatVersion;
        }
        else if (key == "format_revision")
        {
          state = State::ExpectFormatRevision;
        }
        else if (key == "case_name")
        {
          state = State::ExpectCaseName;
        }
        else if (key == "case_date_time")
        {
          state = State::ExpectCaseDateTime;
        }
        else if (key == "case_description")
        {
          state = State::ExpectCaseDescription;
        }
        else if (key == "case_comments")
        {
          state = State::ExpectCaseComments;
        }
        else if (key == "freq_base")
        {
          state = State::ExpectFreqBase;
        }
        else if (key == "va_base")
        {
          state = State::ExpectVaBase;
        }
        else
        {
          return false;
        }
        break;
      case State::ExpectBusKeyOrObjectClose:
        if (key == "number")
        {
          state = State::ExpectBusNumber;
        }
        else if (key == "class")
        {
          state = State::ExpectBusClass;
        }
        else if (key == "name")
        {
          state = State::ExpectBusName;
        }
        else if (key == "init")
        {
          state = State::ExpectBusInitialParameters;
        }
        else if (key == "v_base")
        {
          state = State::ExpectBusVBase;
        }
        else if (key == "mon")
        {
          state = State::ExpectBusMonitor;
        }
        else if (key == "freq_base")
        {
          state = State::ExpectBusFreqBase;
        }
        else if (key == "va_base")
        {
          state = State::ExpectBusVaBase;
        }
        else if (key == "extension")
        {
          state = State::ExpectBusExtension;
        }
        else
        {
          return false;
        }
        break;
      case State::ExpectDeviceKeyOrObjectClose:
        if (key == "class")
        {
          state = State::ExpectDeviceClass;
        }
        else if (key == "ports")
        {
          state = State::ExpectDevicePorts;
        }
        else if (key == "id")
        {
          state = State::ExpectDeviceId;
        }
        else if (key == "params")
        {
          state = State::ExpectDeviceInitialParameters;
        }
        else if (key == "mon")
        {
          state = State::ExpectDeviceMonitor;
        }
        else if (key == "va_base")
        {
          state = State::ExpectDeviceVaBase;
        }
        else if (key == "freq_base")
        {
          state = State::ExpectDeviceFreqBase;
        }
        else if (key == "extension")
        {
          state = State::ExpectDeviceExtension;
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
          state = State::ExpectBusInitialParameterVr;
        }
        else if (key == "Vi")
        {
          state = State::ExpectBusInitialParameterVi;
        }
        else
        {
          return false;
        }
        break;
      case State::ExpectDevicePortKeyOrObjectClose:
        if (key == "bus1")
        {
          staged_device_port = DevicePort::Bus1;
        }
        else if (key == "bus2")
        {
          staged_device_port = DevicePort::Bus2;
        }
        else if (key == "bus")
        {
          staged_device_port = DevicePort::Bus;
        }
        else if (key == "exciter_signal")
        {
          staged_device_port = DevicePort::ExciterSignal;
        }
        else if (key == "governor_signal")
        {
          staged_device_port = DevicePort::GovernorSignal;
        }
        else if (key == "control_signal")
        {
          staged_device_port = DevicePort::ControlSignal;
        }
        else
        {
          return false;
        }
        state = State::ExpectDevicePortMapping;
        break;
      case State::ExpectDeviceInitialParameterKeyOrObjectClose:
        if (key == "R")
        {
          staged_device_parameter = DeviceParameter::R;
        }
        else if (key == "X")
        {
          staged_device_parameter = DeviceParameter::X;
        }
        else if (key == "G")
        {
          staged_device_parameter = DeviceParameter::G;
        }
        else if (key == "B")
        {
          staged_device_parameter = DeviceParameter::B;
        }
        else if (key == "Pz")
        {
          staged_device_parameter = DeviceParameter::Pz;
        }
        else if (key == "Qz")
        {
          staged_device_parameter = DeviceParameter::Qz;
        }
        else if (key == "Pi")
        {
          staged_device_parameter = DeviceParameter::Pi;
        }
        else if (key == "Qi")
        {
          staged_device_parameter = DeviceParameter::Qi;
        }
        else if (key == "Pp")
        {
          staged_device_parameter = DeviceParameter::Pp;
        }
        else if (key == "Qp")
        {
          staged_device_parameter = DeviceParameter::Qp;
        }
        else if (key == "unit_id")
        {
          staged_device_parameter = DeviceParameter::UnitId;
        }
        else if (key == "p0")
        {
          staged_device_parameter = DeviceParameter::P0;
        }
        else if (key == "q0")
        {
          staged_device_parameter = DeviceParameter::Q0;
        }
        else if (key == "H")
        {
          staged_device_parameter = DeviceParameter::H;
        }
        else if (key == "D")
        {
          staged_device_parameter = DeviceParameter::D;
        }
        else if (key == "Ra")
        {
          staged_device_parameter = DeviceParameter::Ra;
        }
        else if (key == "Tdop")
        {
          staged_device_parameter = DeviceParameter::Tdop;
        }
        else if (key == "Tdopp")
        {
          staged_device_parameter = DeviceParameter::Tdopp;
        }
        else if (key == "Tqopp")
        {
          staged_device_parameter = DeviceParameter::Tqopp;
        }
        else if (key == "Tqop")
        {
          staged_device_parameter = DeviceParameter::Tqop;
        }
        else if (key == "Xd")
        {
          staged_device_parameter = DeviceParameter::Xd;
        }
        else if (key == "Xdp")
        {
          staged_device_parameter = DeviceParameter::Xdp;
        }
        else if (key == "Xdpp")
        {
          staged_device_parameter = DeviceParameter::Xdpp;
        }
        else if (key == "Xq")
        {
          staged_device_parameter = DeviceParameter::Xq;
        }
        else if (key == "Xqp")
        {
          staged_device_parameter = DeviceParameter::Xqp;
        }
        else if (key == "Xqpp")
        {
          staged_device_parameter = DeviceParameter::Xqpp;
        }
        else if (key == "Xl")
        {
          staged_device_parameter = DeviceParameter::Xl;
        }
        else if (key == "S10")
        {
          staged_device_parameter = DeviceParameter::S10;
        }
        else if (key == "S12")
        {
          staged_device_parameter = DeviceParameter::S12;
        }
        else if (key == "state0")
        {
          staged_device_parameter = DeviceParameter::state0;
        }
        else
        {
          return false;
        }
        state = State::ExpectDeviceInitialParameterMapping;
        break;
      default:
        return false;
      }
      return true;
    }

    bool StartObject()
    {
      switch (state)
      {
      case State::ExpectHeader:
        state = State::ExpectHeaderKeyOrObjectClose;
        break;
      case State::ExpectOuterObject:
        state = State::ExpectInnerKey;
        break;
      case State::ExpectBusOrArrayClose:
        state    = State::ExpectBusKeyOrObjectClose;
        bus_data = BusDataT();
        break;
      case State::ExpectDeviceOrArrayClose:
        state = State::ExpectDeviceKeyOrObjectClose;
        device_class.clear();
        device_input_parameters.clear();
        device_ports.clear();
        device_id.clear();
        device_va_base.reset();
        device_freq_base.reset();
        monitored_device_variables.reset();
        break;
      case State::ExpectBusInitialParameters:
        state = State::ExpectBusInitialParameterKeyOrObjectClose;
        break;
      case State::ExpectBusExtension:
        state = State::IgnoreUntilBusExtensionObjectClose;
        break;
      case State::ExpectDevicePorts:
        state = State::ExpectDevicePortKeyOrObjectClose;
        break;
      case State::ExpectDeviceInitialParameters:
        state = State::ExpectDeviceInitialParameterKeyOrObjectClose;
        break;
      case State::ExpectDeviceExtension:
        state = State::IgnoreUntilDeviceExtensionObjectClose;
        break;
      default:
        return false;
      }

      return true;
    }

    bool EndObject(size_t)
    {
      switch (state)
      {
      case State::ExpectBusKeyOrObjectClose:
        system_model.bus.push_back(bus_data);
        state = State::ExpectBusOrArrayClose;
        break;
      case State::ExpectDeviceKeyOrObjectClose:
        if (device_class == "branch")
        {
          auto branch_data = BranchDataT{};

          if (auto R = device_input_parameters.find(DeviceParameter::R);
              R != device_input_parameters.end())
          {
            branch_data.R = std::get<RealT>(R->second);
          }
          else
          {
            return false;
          }

          if (auto X = device_input_parameters.find(DeviceParameter::X);
              X != device_input_parameters.end())
          {
            branch_data.X = std::get<RealT>(X->second);
          }
          else
          {
            return false;
          }

          if (auto G = device_input_parameters.find(DeviceParameter::G);
              G != device_input_parameters.end())
          {
            branch_data.G = std::get<RealT>(G->second);
          }
          else
          {
            return false;
          }

          if (auto B = device_input_parameters.find(DeviceParameter::B);
              B != device_input_parameters.end())
          {
            branch_data.B = std::get<RealT>(B->second);
          }
          else
          {
            return false;
          }

          if (auto bus1 = device_ports.find(DevicePort::Bus1);
              bus1 != device_ports.end())
          {
            branch_data.bus1_id = bus1->second;
          }
          else
          {
            return false;
          }

          if (auto bus2 = device_ports.find(DevicePort::Bus2);
              bus2 != device_ports.end())
          {
            branch_data.bus2_id = bus2->second;
          }
          else
          {
            return false;
          }

          branch_data.freq_base             = device_freq_base;
          branch_data.va_base               = device_va_base;
          branch_data.disambiguation_string = device_id;

          branch_data.monitored_variables[static_cast<size_t>(BranchDataT::MonitorableVariables::Ir1)] = monitored_device_variables[static_cast<size_t>(MonitorableDeviceVariables::Ir1)];
          branch_data.monitored_variables[static_cast<size_t>(BranchDataT::MonitorableVariables::Ii1)] = monitored_device_variables[static_cast<size_t>(MonitorableDeviceVariables::Ii1)];
          branch_data.monitored_variables[static_cast<size_t>(BranchDataT::MonitorableVariables::Im1)] = monitored_device_variables[static_cast<size_t>(MonitorableDeviceVariables::Im1)];
          branch_data.monitored_variables[static_cast<size_t>(BranchDataT::MonitorableVariables::P1)]  = monitored_device_variables[static_cast<size_t>(MonitorableDeviceVariables::P1)];
          branch_data.monitored_variables[static_cast<size_t>(BranchDataT::MonitorableVariables::Q1)]  = monitored_device_variables[static_cast<size_t>(MonitorableDeviceVariables::Q1)];
          branch_data.monitored_variables[static_cast<size_t>(BranchDataT::MonitorableVariables::Ir2)] = monitored_device_variables[static_cast<size_t>(MonitorableDeviceVariables::Ir2)];
          branch_data.monitored_variables[static_cast<size_t>(BranchDataT::MonitorableVariables::Ii2)] = monitored_device_variables[static_cast<size_t>(MonitorableDeviceVariables::Ii2)];
          branch_data.monitored_variables[static_cast<size_t>(BranchDataT::MonitorableVariables::Im2)] = monitored_device_variables[static_cast<size_t>(MonitorableDeviceVariables::Im2)];
          branch_data.monitored_variables[static_cast<size_t>(BranchDataT::MonitorableVariables::P2)]  = monitored_device_variables[static_cast<size_t>(MonitorableDeviceVariables::P2)];
          branch_data.monitored_variables[static_cast<size_t>(BranchDataT::MonitorableVariables::Q2)]  = monitored_device_variables[static_cast<size_t>(MonitorableDeviceVariables::Q2)];

          system_model.branch.push_back(branch_data);
        }
        else if (device_class == "static_load")
        {
          // TODO: there doesn't appear to be an implementation
          //       of this in gridkit
          return false;
        }
        else if (device_class == "GENROU")
        {
          auto genrou_data = GenrouDataT{};

          if (auto unit_id = device_input_parameters.find(DeviceParameter::UnitId);
              unit_id != device_input_parameters.end())
          {
            genrou_data.unit_id = std::get<IdxT>(unit_id->second);
          }
          else
          {
            return false;
          }

          if (auto p0 = device_input_parameters.find(DeviceParameter::P0);
              p0 != device_input_parameters.end())
          {
            genrou_data.p0 = std::get<RealT>(p0->second);
          }
          else
          {
            return false;
          }

          if (auto q0 = device_input_parameters.find(DeviceParameter::Q0);
              q0 != device_input_parameters.end())
          {
            genrou_data.q0 = std::get<RealT>(q0->second);
          }
          else
          {
            return false;
          }

          if (auto H = device_input_parameters.find(DeviceParameter::H);
              H != device_input_parameters.end())
          {
            genrou_data.H = std::get<RealT>(H->second);
          }
          else
          {
            return false;
          }

          if (auto D = device_input_parameters.find(DeviceParameter::D);
              D != device_input_parameters.end())
          {
            genrou_data.D = std::get<RealT>(D->second);
          }
          else
          {
            return false;
          }

          if (auto Ra = device_input_parameters.find(DeviceParameter::Ra);
              Ra != device_input_parameters.end())
          {
            genrou_data.Ra = std::get<RealT>(Ra->second);
          }
          else
          {
            return false;
          }

          if (auto Tdop = device_input_parameters.find(DeviceParameter::Tdop);
              Tdop != device_input_parameters.end())
          {
            genrou_data.Tdop = std::get<RealT>(Tdop->second);
          }
          else
          {
            return false;
          }

          if (auto Tdopp = device_input_parameters.find(DeviceParameter::Tdopp);
              Tdopp != device_input_parameters.end())
          {
            genrou_data.Tdopp = std::get<RealT>(Tdopp->second);
          }
          else
          {
            return false;
          }

          if (auto Tqop = device_input_parameters.find(DeviceParameter::Tqop);
              Tqop != device_input_parameters.end())
          {
            genrou_data.Tqop = std::get<RealT>(Tqop->second);
          }
          else
          {
            return false;
          }

          if (auto Tqopp = device_input_parameters.find(DeviceParameter::Tqopp);
              Tqopp != device_input_parameters.end())
          {
            genrou_data.Tqopp = std::get<RealT>(Tqopp->second);
          }
          else
          {
            return false;
          }

          if (auto Xd = device_input_parameters.find(DeviceParameter::Xd);
              Xd != device_input_parameters.end())
          {
            genrou_data.Xd = std::get<RealT>(Xd->second);
          }
          else
          {
            return false;
          }

          if (auto Xdp = device_input_parameters.find(DeviceParameter::Xdp);
              Xdp != device_input_parameters.end())
          {
            genrou_data.Xdp = std::get<RealT>(Xdp->second);
          }
          else
          {
            return false;
          }

          if (auto Xdpp = device_input_parameters.find(DeviceParameter::Xdpp);
              Xdpp != device_input_parameters.end())
          {
            genrou_data.Xdpp = std::get<RealT>(Xdpp->second);
          }
          else
          {
            return false;
          }

          if (auto Xq = device_input_parameters.find(DeviceParameter::Xq);
              Xq != device_input_parameters.end())
          {
            genrou_data.Xq = std::get<RealT>(Xq->second);
          }
          else
          {
            return false;
          }

          if (auto Xqp = device_input_parameters.find(DeviceParameter::Xqp);
              Xqp != device_input_parameters.end())
          {
            genrou_data.Xqp = std::get<RealT>(Xqp->second);
          }
          else
          {
            return false;
          }

          if (auto Xqpp = device_input_parameters.find(DeviceParameter::Xqpp);
              Xqpp != device_input_parameters.end())
          {
            genrou_data.Xqpp = std::get<RealT>(Xqpp->second);
          }
          else
          {
            return false;
          }

          if (auto Xl = device_input_parameters.find(DeviceParameter::Xl);
              Xl != device_input_parameters.end())
          {
            genrou_data.Xl = std::get<RealT>(Xl->second);
          }
          else
          {
            return false;
          }

          if (auto S10 = device_input_parameters.find(DeviceParameter::S10);
              S10 != device_input_parameters.end())
          {
            genrou_data.S10 = std::get<RealT>(S10->second);
          }
          else
          {
            return false;
          }

          if (auto S12 = device_input_parameters.find(DeviceParameter::S12);
              S12 != device_input_parameters.end())
          {
            genrou_data.S12 = std::get<RealT>(S12->second);
          }
          else
          {
            return false;
          }

          if (auto bus = device_ports.find(DevicePort::Bus);
              bus != device_ports.end())
          {
            genrou_data.bus_id = bus->second;
          }
          else
          {
            return false;
          }

          genrou_data.freq_base             = device_freq_base;
          genrou_data.va_base               = device_va_base;
          genrou_data.disambiguation_string = device_id;

          genrou_data.monitored_variables[static_cast<size_t>(GenrouDataT::MonitorableVariables::Ir)]    = monitored_device_variables[static_cast<size_t>(MonitorableDeviceVariables::Ir)];
          genrou_data.monitored_variables[static_cast<size_t>(GenrouDataT::MonitorableVariables::Ii)]    = monitored_device_variables[static_cast<size_t>(MonitorableDeviceVariables::Ii)];
          genrou_data.monitored_variables[static_cast<size_t>(GenrouDataT::MonitorableVariables::P)]     = monitored_device_variables[static_cast<size_t>(MonitorableDeviceVariables::P)];
          genrou_data.monitored_variables[static_cast<size_t>(GenrouDataT::MonitorableVariables::Q)]     = monitored_device_variables[static_cast<size_t>(MonitorableDeviceVariables::Q)];
          genrou_data.monitored_variables[static_cast<size_t>(GenrouDataT::MonitorableVariables::Delta)] = monitored_device_variables[static_cast<size_t>(MonitorableDeviceVariables::Delta)];
          genrou_data.monitored_variables[static_cast<size_t>(GenrouDataT::MonitorableVariables::Omega)] = monitored_device_variables[static_cast<size_t>(MonitorableDeviceVariables::Omega)];

          system_model.genrou.push_back(genrou_data);
        }
        else if (device_class == "bus_fault")
        {
          auto bus_fault_data = BusFaultDataT{};

          if (auto R = device_input_parameters.find(DeviceParameter::R);
              R != device_input_parameters.end())
          {
            bus_fault_data.R = std::get<RealT>(R->second);
          }
          else
          {
            return false;
          }

          if (auto X = device_input_parameters.find(DeviceParameter::X);
              X != device_input_parameters.end())
          {
            bus_fault_data.X = std::get<RealT>(X->second);
          }
          else
          {
            return false;
          }

          if (auto status = device_input_parameters.find(DeviceParameter::state0);
              status != device_input_parameters.end())
          {
            bus_fault_data.status = std::get<bool>(status->second);
          }
          else
          {
            return false;
          }

          if (auto bus = device_ports.find(DevicePort::Bus);
              bus != device_ports.end())
          {
            bus_fault_data.bus_id = bus->second;
          }
          else
          {
            return false;
          }

          bus_fault_data.freq_base             = device_freq_base;
          bus_fault_data.va_base               = device_va_base;
          bus_fault_data.disambiguation_string = device_id;

          bus_fault_data.monitored_variables[static_cast<size_t>(BusFaultDataT::MonitorableVariables::State)] = monitored_device_variables[static_cast<size_t>(MonitorableDeviceVariables::State)];
          bus_fault_data.monitored_variables[static_cast<size_t>(BusFaultDataT::MonitorableVariables::Ir)]    = monitored_device_variables[static_cast<size_t>(MonitorableDeviceVariables::Ir)];
          bus_fault_data.monitored_variables[static_cast<size_t>(BusFaultDataT::MonitorableVariables::Ii)]    = monitored_device_variables[static_cast<size_t>(MonitorableDeviceVariables::Ii)];

          system_model.bus_fault.push_back(bus_fault_data);
        }
        else
        {
          return false;
        }
        state = State::ExpectDeviceOrArrayClose;
        break;
      case State::ExpectBusInitialParameterKeyOrObjectClose:
        state = State::ExpectBusKeyOrObjectClose;
        break;
      case State::IgnoreUntilBusExtensionObjectClose:
        state = State::ExpectBusKeyOrObjectClose;
        break;
      case State::ExpectDevicePortKeyOrObjectClose:
        state = State::ExpectDeviceKeyOrObjectClose;
        break;
      case State::IgnoreUntilDeviceExtensionObjectClose:
        state = State::ExpectDeviceKeyOrObjectClose;
        break;
      default:
        return false;
      }

      return true;
    }

    bool StartArray()
    {
      switch (state)
      {
      case State::ExpectBuses:
        state = State::ExpectBusOrArrayClose;
        break;
      case State::ExpectDevices:
        state = State::ExpectDeviceOrArrayClose;
        break;
      case State::ExpectBusMonitor:
        state = State::ExpectMonitoredBusVariableOrArrayClose;
        break;
      case State::ExpectDeviceMonitor:
        state = State::ExpectMonitoredDeviceVariableOrArrayClose;
        break;
      default:
        return false;
      }

      return true;
    }

    bool EndArray(size_t)
    {
      switch (state)
      {
      case State::ExpectBusOrArrayClose:
      case State::ExpectDeviceOrArrayClose:
        state = State::ExpectInnerKey;
        break;
      case State::ExpectMonitoredBusVariableOrArrayClose:
        state = State::ExpectBusKeyOrObjectClose;
        break;
      case State::ExpectMonitoredDeviceVariableOrArrayClose:
        state = State::ExpectDeviceKeyOrObjectClose;
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
      UnitId,
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
      Ir,
      Ii,
      Delta,
      Omega,
      Maximum,
    };

    /// Device class
    std::string device_class;

    /// Mapping of device input parameters to values
    std::map<DeviceParameter, std::variant<RealT, IdxT, bool>> device_input_parameters;

    /// Place for holding the device initial parameter to map to a value
    DeviceParameter staged_device_parameter;

    /// Mapping of ports to bus indices
    std::map<DevicePort, IdxT> device_ports;

    /// Place for holding the device port to map to a bus index
    DevicePort staged_device_port;

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
      IgnoreUntilBusExtensionObjectClose,
      ExpectDevices,
      ExpectDeviceOrArrayClose,
      ExpectDeviceKeyOrObjectClose,
      ExpectDeviceClass,
      ExpectDevicePorts,
      ExpectDevicePortKeyOrObjectClose,
      ExpectDevicePortMapping,
      ExpectDeviceId,
      ExpectDeviceInitialParameters,
      ExpectDeviceInitialParameterKeyOrObjectClose,
      ExpectDeviceInitialParameterMapping,
      ExpectDeviceMonitor,
      ExpectMonitoredDeviceVariableOrArrayClose,
      ExpectDeviceVaBase,
      ExpectDeviceFreqBase,
      ExpectDeviceExtension,
      IgnoreUntilDeviceExtensionObjectClose,
    } state;
  };
} // namespace GridKit
