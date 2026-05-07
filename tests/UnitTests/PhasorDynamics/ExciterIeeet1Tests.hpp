#pragma once

#include <cmath>
#include <iostream>
#include <limits>

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/IEEET1/Ieeet1.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/IEEET1/Ieeet1Data.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class ExciterIeeet1Tests
    {
    public:
      using RealT = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      ExciterIeeet1Tests()  = default;
      ~ExciterIeeet1Tests() = default;

      TestOutcome constructor()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(3.0, 4.0);
        auto                               data    = makeTestData();
        auto*                              exciter = new PhasorDynamics::Exciter::Ieeet1<ScalarT, IdxT>(&bus, data);

        success *= (exciter != nullptr);
        success *= (exciter->size() == 9);
        success *= (exciter->getMonitor() != nullptr);

        delete exciter;

        return success.report(__func__);
      }

    private:
      auto makeTestData() -> PhasorDynamics::Exciter::Ieeet1Data<RealT, IdxT>
      {
        using Params = PhasorDynamics::Exciter::Ieeet1Parameters;

        PhasorDynamics::Exciter::Ieeet1Data<RealT, IdxT> data;
        data.device_class          = "exciter";
        data.disambiguation_string = "ieeet1_test";
        data.monitored_variables.insert(PhasorDynamics::Exciter::Ieeet1MonitorableVariables::efd);

        data.parameters[Params::Tr]      = 0.001;
        data.parameters[Params::Ka]      = 50.0;
        data.parameters[Params::Ta]      = 0.04;
        data.parameters[Params::Ke]      = -0.06;
        data.parameters[Params::Te]      = 0.6;
        data.parameters[Params::Kf]      = 0.09;
        data.parameters[Params::Tf]      = 1.46;
        data.parameters[Params::Vrmin]   = -1.0;
        data.parameters[Params::Vrmax]   = 1.0;
        data.parameters[Params::E1]      = 2.8;
        data.parameters[Params::E2]      = 3.373;
        data.parameters[Params::Se1]     = 0.04;
        data.parameters[Params::Se2]     = 0.33;
        data.parameters[Params::Ispdlim] = 0.0;

        return data;
      }
    };
  } // namespace Testing
} // namespace GridKit
