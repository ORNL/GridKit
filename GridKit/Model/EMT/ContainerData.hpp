#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include <GridKit/Model/EMT/Component/Bus/BusData.hpp>
#include <GridKit/Model/EMT/Component/Controller/DCLink/DcLinkData.hpp>
#include <GridKit/Model/EMT/Component/Controller/GASTPTI/GastPtiData.hpp>
#include <GridKit/Model/EMT/Component/Controller/IEEEST/IeeestData.hpp>
#include <GridKit/Model/EMT/Component/Controller/IEEET1/Ieeet1Data.hpp>
#include <GridKit/Model/EMT/Component/Controller/InnerCurrentControl/InnerCurrentControlData.hpp>
#include <GridKit/Model/EMT/Component/Controller/OuterPowerControl/OuterPowerControlData.hpp>
#include <GridKit/Model/EMT/Component/Controller/OuterVoltageControl/OuterVoltageControlData.hpp>
#include <GridKit/Model/EMT/Component/Controller/PWM/PwmData.hpp>
#include <GridKit/Model/EMT/Component/Controller/SEXS-PTI/SexsPtiData.hpp>
#include <GridKit/Model/EMT/Component/Controller/TGOV1/Tgov1Data.hpp>
#include <GridKit/Model/EMT/Component/Filter/FilterData.hpp>
#include <GridKit/Model/EMT/Component/Line/LineDistributed/LineDistributedData.hpp>
#include <GridKit/Model/EMT/Component/Line/LineLumped/LineLumpedData.hpp>
#include <GridKit/Model/EMT/Component/Load/LoadZ/LoadZData.hpp>
#include <GridKit/Model/EMT/Component/Source/DependentVoltageSource/DependentVoltageSourceData.hpp>
#include <GridKit/Model/EMT/Component/Source/Machine/MachineData.hpp>
#include <GridKit/Model/EMT/Component/Source/REGFMA/RegfmaData.hpp>
#include <GridKit/Model/EMT/Component/Source/VoltageSource/VoltageSourceData.hpp>
#include <GridKit/Model/EMT/Component/Switch/SwitchData.hpp>
#include <GridKit/Model/EMT/Component/Transformer/TransformerData.hpp>
#include <GridKit/Model/EMT/Operators/Converter/ConverterData.hpp>
#include <GridKit/Model/EMT/Operators/Modulation/ModulationData.hpp>
#include <GridKit/Model/EMT/Operators/Reference/PLL/PllData.hpp>
#include <GridKit/Model/EMT/Operators/Reference/Park/ParkData.hpp>
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

      using InnerCurrentControlDataT    = Controller::InnerCurrentControlData<RealT, IdxT>;
      using OuterVoltageControlDataT    = Controller::OuterVoltageControlData<RealT, IdxT>;
      using OuterPowerControlDataT      = Controller::OuterPowerControlData<RealT, IdxT>;
      using ParkDataT                   = ParkData<RealT, IdxT>;
      using PllDataT                    = PllData<RealT, IdxT>;
      using FilterDataT                 = FilterData<RealT, IdxT>;
      using ModulationDataT             = ModulationData<RealT, IdxT>;
      using PwmDataT                    = Controller::PwmData<RealT, IdxT>;
      using DcLinkDataT                 = Controller::DcLinkData<RealT, IdxT>;
      using ConverterDataT              = ConverterData<RealT, IdxT>;
      using BusDataT                    = BusData<RealT, IdxT>;
      using DependentVoltageSourceDataT = DependentVoltageSourceData<RealT, IdxT>;
      using LineLumpedDataT             = LineLumpedData<RealT, IdxT>;
      using LoadZDataT                  = LoadZData<RealT, IdxT>;
      using MachineDataT                = MachineData<RealT, IdxT>;
      using RegfmaDataT                 = RegfmaData<RealT, IdxT>;
      using SignalDataT                 = SignalData<RealT, IdxT>;
      using SwitchDataT                 = SwitchData<RealT, IdxT>;
      using TransformerDataT            = TransformerData<RealT, IdxT>;
      using IeeestDataT                 = Controller::IeeestData<RealT, IdxT>;
      using GastPtiDataT                = Controller::GastPtiData<RealT, IdxT>;
      using Tgov1DataT                  = Controller::Tgov1Data<RealT, IdxT>;
      using SexsPtiDataT                = Controller::SexsPtiData<RealT, IdxT>;
      using Ieeet1DataT                 = Controller::Ieeet1Data<RealT, IdxT>;
      using VoltageSourceDataT          = VoltageSourceData<RealT, IdxT>;

      /// Identifier within the parent scope. Empty only for the root.
      std::string id;

      /// Public input name to its source endpoint in the parent scope.
      std::map<std::string, std::string> inputs;

      /// Public output name to an internal scalar signal or electrical bus.
      std::map<std::string, std::string> outputs;

      std::vector<SignalDataT>   signal;    ///< Signals local to this scope
      std::vector<ContainerData> container; ///< Child scopes

      std::vector<InnerCurrentControlDataT>         inner_current_control;
      std::vector<OuterVoltageControlDataT>         outer_voltage_control;
      std::vector<OuterPowerControlDataT>           outer_power_control;
      std::vector<ParkDataT>                        park;
      std::vector<PllDataT>                         pll;
      std::vector<FilterDataT>                      filter;
      std::vector<ModulationDataT>                  modulation;
      std::vector<PwmDataT>                         pwm;
      std::vector<DcLinkDataT>                      dc_link;
      std::vector<ConverterDataT>                   converter;
      std::vector<BusDataT>                         bus;
      std::vector<DependentVoltageSourceDataT>      dependent_voltage_source;
      std::vector<LineLumpedDataT>                  line_lumped;
      std::vector<LineDistributedData<RealT, IdxT>> line_distributed;
      std::vector<LoadZDataT>                       loadz;
      std::vector<MachineDataT>                     machine;
      std::vector<RegfmaDataT>                      regfma;
      std::vector<SwitchDataT>                      sw;
      std::vector<TransformerDataT>                 transformer;
      std::vector<IeeestDataT>                      ieeest;
      std::vector<GastPtiDataT>                     gastpti;
      std::vector<Tgov1DataT>                       gov;
      std::vector<SexsPtiDataT>                     sexs_pti;
      std::vector<Ieeet1DataT>                      exciter;
      std::vector<VoltageSourceDataT>               voltage_source;
    };
  } // namespace EMT
} // namespace GridKit
