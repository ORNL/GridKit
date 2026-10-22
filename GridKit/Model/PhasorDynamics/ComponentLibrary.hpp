#pragma once

#include <cstddef>
#include <type_traits>

#include <GridKit/Model/PhasorDynamics/Branch/Branch.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusSignalVoltageIn/BusSignalVoltageIn.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusSignalVoltageOut/BusSignalVoltageOut.hpp>
#include <GridKit/Model/PhasorDynamics/Controller/REECB/Reecb.hpp>
#include <GridKit/Model/PhasorDynamics/Controller/REPCA/Repca.hpp>
#include <GridKit/Model/PhasorDynamics/Converter/REGCA/Regca.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/ESDC1A/Esdc1a.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/IEEET1/Ieeet1.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/SEXS-PTI/SexsPti.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/GASTPTI/GastPti.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/HYGOV/Hygov.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/Tgov1/Tgov1.hpp>
#include <GridKit/Model/PhasorDynamics/Load/LoadZ/LoadZ.hpp>
#include <GridKit/Model/PhasorDynamics/Load/LoadZIP/LoadZIP.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SignalSource/ConstantSignalSource.hpp>
#include <GridKit/Model/PhasorDynamics/SignalSource/FunctionSignalSource.hpp>
#include <GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/Ieeest.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENROU/Genrou.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENSAL/Gensal.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GenClassical/GenClassical.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Component models with a manifest, in `SystemModel` construction order
    template <typename ScalarT = double, typename IdxT = size_t>
    struct ComponentLibrary
    {
      /// Calls `visit` with each model type and its data type, in construction order
      template <typename Visitor>
      static constexpr void forEach(Visitor&& visit)
      {
        using RealT = typename ScalarTraits<ScalarT>::RealT;
        visit(std::type_identity<Genrou<ScalarT, IdxT>>{},
              std::type_identity<GenrouData<RealT, IdxT>>{});
        visit(std::type_identity<Exciter::Ieeet1<ScalarT, IdxT>>{},
              std::type_identity<Exciter::Ieeet1Data<RealT, IdxT>>{});
      }
    };
  } // namespace PhasorDynamics
} // namespace GridKit
