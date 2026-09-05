#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include <GridKit/Model/EMT/Component/Bus/BusData.hpp>
#include <GridKit/Model/EMT/Component/Controller/TGOV1/Tgov1Data.hpp>
#include <GridKit/Model/EMT/Component/Line/LineLumped/LineLumpedData.hpp>
#include <GridKit/Model/EMT/Component/Load/LoadZ/LoadZData.hpp>
#include <GridKit/Model/EMT/Component/Source/DependentVoltageSource/DependentVoltageSourceData.hpp>
#include <GridKit/Model/EMT/Component/Source/Machine/MachineData.hpp>
#include <GridKit/Model/EMT/Component/Source/VoltageSource/VoltageSourceData.hpp>
#include <GridKit/Model/EMT/Component/Switch/SwitchData.hpp>
#include <GridKit/Model/EMT/Signal/SignalData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Data for one compositional EMT scope.
     *
     * A container has the same model body as an EMT case: local signals and
     * devices, including other containers. Inputs and outputs name the local
     * endpoints that form the container boundary.
     */
    template <typename real_type = double, typename index_type = size_t>
    struct ContainerData
    {
      using RealT = real_type;
      using IdxT  = index_type;

      using BusDataT                    = BusData<RealT, IdxT>;
      using DependentVoltageSourceDataT = DependentVoltageSourceData<RealT, IdxT>;
      using LineLumpedDataT             = LineLumpedData<RealT, IdxT>;
      using LoadZDataT                  = LoadZData<RealT, IdxT>;
      using MachineDataT                = MachineData<RealT, IdxT>;
      using SignalDataT                 = SignalData<RealT, IdxT>;
      using SwitchDataT                 = SwitchData<RealT, IdxT>;
      using Tgov1DataT                  = Controller::Tgov1Data<RealT, IdxT>;
      using VoltageSourceDataT          = VoltageSourceData<RealT, IdxT>;

      /// Identifier within the parent scope. Empty only for the root.
      std::string id;

      /// Public input name to its source endpoint in the parent scope.
      std::map<std::string, std::string> inputs;

      /// Public output name to an internal scalar signal or electrical bus.
      std::map<std::string, std::string> outputs;

      std::vector<SignalDataT>   signal;    ///< Signals local to this scope
      std::vector<ContainerData> container; ///< Child scopes

      std::vector<BusDataT>                    bus;
      std::vector<DependentVoltageSourceDataT> dependent_voltage_source;
      std::vector<LineLumpedDataT>             line_lumped;
      std::vector<LoadZDataT>                  loadz;
      std::vector<MachineDataT>                machine;
      std::vector<SwitchDataT>                 sw;
      std::vector<Tgov1DataT>                  gov;
      std::vector<VoltageSourceDataT>          voltage_source;
    };
  } // namespace EMT
} // namespace GridKit
