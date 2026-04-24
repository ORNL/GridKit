#pragma once

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/SEXS-PTI/SexsPti.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/SEXS-PTI/SexsPtiData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class ExciterSexsPtiTests
    {
    public:
      using RealT = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      ExciterSexsPtiTests()  = default;
      ~ExciterSexsPtiTests() = default;

      TestOutcome constructor()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(3.0, 4.0);
        auto                               data    = makeTestData();
        auto*                              exciter = new PhasorDynamics::Exciter::SexsPti<ScalarT, IdxT>(&bus, data);

        success *= (exciter != nullptr);
        success *= (exciter->size() == 3);
        success *= (exciter->getMonitor() != nullptr);

        delete exciter;

        return success.report(__func__);
      }

      TestOutcome zeroInitialResidual()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(3.0, 4.0);
        bus.allocate();
        bus.initialize();

        PhasorDynamics::SignalNode<ScalarT, IdxT> efd_node;
        ScalarT                                   efd_value{0.0};
        IdxT                                      efd_index = INVALID_INDEX<IdxT>;
        efd_node.set(&efd_value, &efd_index);

        auto                                            data = makeTestData();
        PhasorDynamics::Exciter::SexsPti<ScalarT, IdxT> exciter(&bus, data);
        exciter.getSignals().template assignSignalNode<PhasorDynamics::Exciter::SexsPtiInternalVariables::EFD>(&efd_node);

        exciter.allocate();
        efd_node.init(1.2);
        success *= (exciter.verify() == 0);
        exciter.initialize();
        exciter.evaluateResidual();

        success *= efd_node.linked();
        success *= (efd_node.getVariableIndex() == 1);
        success *= isEqual(efd_node.read(), static_cast<ScalarT>(1.2), static_cast<ScalarT>(1.0e-12));
        success *= isEqual(exciter.y()[0], static_cast<ScalarT>(-0.048), static_cast<ScalarT>(1.0e-12));
        success *= isEqual(exciter.y()[1], static_cast<ScalarT>(1.2), static_cast<ScalarT>(1.0e-12));
        success *= isEqual(exciter.y()[2], static_cast<ScalarT>(0.12), static_cast<ScalarT>(1.0e-12));

        const auto& f = exciter.getResidual();
        for (size_t i = 0; i < f.size(); ++i)
        {
          if (!isEqual(f[i], static_cast<ScalarT>(0.0), static_cast<ScalarT>(1.0e-10)))
          {
            std::cout << "Non-zero SEXS-PTI residual at index " << i << ": " << f[i] << "\n";
            success = false;
          }
        }

        return success.report(__func__);
      }

      TestOutcome vsSignal()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(3.0, 4.0);
        bus.allocate();
        bus.initialize();

        PhasorDynamics::SignalNode<ScalarT, IdxT> efd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> vs_node;
        ScalarT                                   efd_value{0.0};
        ScalarT                                   vs_value{0.25};
        IdxT                                      efd_index = INVALID_INDEX<IdxT>;
        IdxT                                      vs_index  = 4;
        efd_node.set(&efd_value, &efd_index);
        vs_node.set(&vs_value, &vs_index);

        auto                                            data = makeTestData();
        PhasorDynamics::Exciter::SexsPti<ScalarT, IdxT> exciter(&bus, data);
        exciter.getSignals().template assignSignalNode<PhasorDynamics::Exciter::SexsPtiInternalVariables::EFD>(&efd_node);
        exciter.getSignals().template attachSignalNode<PhasorDynamics::Exciter::SexsPtiExternalVariables::VS>(&vs_node);

        exciter.allocate();
        efd_node.init(1.2);
        success *= (exciter.verify() == 0);
        exciter.initialize();
        exciter.evaluateResidual();

        success *= isEqual(exciter.getResidual()[2], vs_value, static_cast<ScalarT>(1.0e-12));

        return success.report(__func__);
      }

      TestOutcome antiWindupLimiter()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(3.0, 4.0);
        bus.allocate();
        bus.initialize();

        auto                                            data = makeTestData();
        PhasorDynamics::Exciter::SexsPti<ScalarT, IdxT> exciter(&bus, data);
        exciter.allocate();
        exciter.initialize();

        exciter.y()[0]  = -1.0;
        exciter.y()[1]  = 10.0;
        exciter.y()[2]  = 1.0;
        exciter.yp()[0] = 0.0;
        exciter.yp()[1] = 0.0;
        exciter.evaluateResidual();
        auto blocked = std::abs(exciter.getResidual()[1]);

        exciter.y()[0] = 1.0;
        exciter.y()[1] = 10.0;
        exciter.y()[2] = 0.0;
        exciter.evaluateResidual();
        auto returning = std::abs(exciter.getResidual()[1]);

        success *= (blocked < static_cast<ScalarT>(1.0e-8));
        success *= (returning > static_cast<ScalarT>(1.0));

        // Regression guard: the limiter direction follows the SEXS g term.
        // A small positive g should be blocked when Efd is above Efdmax.
        exciter.y()[0] = -0.2575;
        exciter.y()[1] = 5.05;
        exciter.y()[2] = 0.0;
        exciter.evaluateResidual();
        auto small_push_above_limit  = std::abs(exciter.getResidual()[1]);
        success                     *= (small_push_above_limit < static_cast<ScalarT>(1.0e-5));

        return success.report(__func__);
      }

      TestOutcome parameterValidation()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);

        auto missing = makeTestData();
        missing.parameters.erase(PhasorDynamics::Exciter::SexsPtiParameters::K);
        PhasorDynamics::Exciter::SexsPti<ScalarT, IdxT> missing_model(&bus, missing);
        success *= (missing_model.verify() > 0);

        auto invalid                                                       = makeTestData();
        invalid.parameters[PhasorDynamics::Exciter::SexsPtiParameters::Tb] = 0.0;
        PhasorDynamics::Exciter::SexsPti<ScalarT, IdxT> invalid_model(&bus, invalid);
        success *= (invalid_model.verify() > 0);

        return success.report(__func__);
      }

      TestOutcome systemAssembly()
      {
        using namespace PhasorDynamics;
        using namespace PhasorDynamics::Exciter;

        TestStatus success = true;

        SystemModelData<RealT, IdxT> data;
        data.bus.resize(1);
        data.bus[0].bus_id   = 1;
        data.bus[0].bus_type = BusData<RealT, IdxT>::BusType::DEFAULT;
        data.bus[0].Vr0      = 1.0;
        data.bus[0].Vi0      = 0.0;

        data.signal.resize(1);
        data.signal[0].signal_id = 3;
        data.signal[0].name      = "EFD";

        auto exciter                     = makeTestData();
        exciter.ports[SexsPtiPorts::bus] = 1;
        exciter.ports[SexsPtiPorts::efd] = 3;
        data.sexspti.push_back(exciter);

        SystemModel<ScalarT, IdxT> system(data);
        success *= (system.allocate() == 0);
        success *= (system.initialize() == 0);
        success *= (system.evaluateResidual() == 0);
        success *= (system.size() == 5);

        auto missing_efd = data;
        missing_efd.sexspti[0].ports.erase(SexsPtiPorts::efd);
        SystemModel<ScalarT, IdxT> missing_efd_system(missing_efd);
        success *= (missing_efd_system.verify() > 0);

        return success.report(__func__);
      }

      TestOutcome systemSignalFreshness()
      {
        using namespace PhasorDynamics;
        using namespace PhasorDynamics::Exciter;

        TestStatus success = true;

        SystemModelData<RealT, IdxT> data;
        data.bus.resize(1);
        data.bus[0].bus_id   = 1;
        data.bus[0].bus_type = BusData<RealT, IdxT>::BusType::DEFAULT;
        data.bus[0].Vr0      = 1.0;
        data.bus[0].Vi0      = 0.0;

        data.signal.resize(2);
        data.signal[0].signal_id = 10;
        data.signal[0].name      = "source";
        data.signal[1].signal_id = 20;
        data.signal[1].name      = "consumer_efd";

        auto consumer                     = makeTestData();
        consumer.disambiguation_string    = "consumer";
        consumer.ports[SexsPtiPorts::bus] = 1;
        consumer.ports[SexsPtiPorts::efd] = 20;
        consumer.ports[SexsPtiPorts::vs]  = 10;
        data.sexspti.push_back(consumer);

        auto source                     = makeTestData();
        source.disambiguation_string    = "source";
        source.ports[SexsPtiPorts::bus] = 1;
        source.ports[SexsPtiPorts::efd] = 10;
        data.sexspti.push_back(source);

        SystemModel<ScalarT, IdxT> system(data);
        success *= (system.allocate() == 0);
        success *= (system.initialize() == 0);

        constexpr IdxT consumer_vtr_residual  = 4;
        constexpr IdxT source_efd_variable    = 6;
        system.y()[source_efd_variable]       = 0.75;
        success                              *= (system.evaluateResidual() == 0);
        success                              *= isEqual(system.getResidual()[consumer_vtr_residual],
                           static_cast<ScalarT>(0.75),
                           static_cast<ScalarT>(1.0e-12));

        return success.report(__func__);
      }

    private:
      auto makeTestData() -> PhasorDynamics::Exciter::SexsPtiData<RealT, IdxT>
      {
        using Params = PhasorDynamics::Exciter::SexsPtiParameters;

        PhasorDynamics::Exciter::SexsPtiData<RealT, IdxT> data;
        data.device_class          = "exciter";
        data.disambiguation_string = "sexspti_test";
        data.monitored_variables.insert(PhasorDynamics::Exciter::SexsPtiMonitorableVariables::efd);

        data.parameters[Params::Ta]     = 0.1;
        data.parameters[Params::Tb]     = 0.5;
        data.parameters[Params::Te]     = 0.8;
        data.parameters[Params::K]      = 10.0;
        data.parameters[Params::Efdmax] = 5.0;
        data.parameters[Params::Efdmin] = -5.0;

        return data;
      }
    };

  } // namespace Testing
} // namespace GridKit
