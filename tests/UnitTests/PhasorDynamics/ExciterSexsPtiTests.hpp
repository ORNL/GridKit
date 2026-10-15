#pragma once

#include <iostream>
#include <limits>
#include <stdexcept>

#include <GridKit/Definitions.hpp>
#include <GridKit/Model/ConfigurationChecks.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/SEXS-PTI/SexsPti.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/SEXS-PTI/SexsPtiData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    using Log = GridKit::Utilities::Logger;

    template <class ScalarT, typename IdxT>
    class ExciterSexsPtiTests
    {
    public:
      using RealT  = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;
      using Data   = PhasorDynamics::Exciter::SexsPtiData<RealT, IdxT>;
      using Params = PhasorDynamics::Exciter::SexsPtiParameters;

      ExciterSexsPtiTests()  = default;
      ~ExciterSexsPtiTests() = default;

      // Init and saturated-limiter checks are exact algebraic identities
      // (no iteration, sigmoid underflowed to 0/1 at the depths tested here).
      static constexpr ScalarT kTol = 100 * std::numeric_limits<RealT>::epsilon();

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
        efd_node.link(&efd_value, &efd_index);

        using namespace GridKit::PhasorDynamics::Exciter;

        auto data    = makeTestData();
        auto exciter = SexsPti<ScalarT, IdxT>(&bus, data);
        exciter.getPorts().out.template port<SexsPtiSignalOutputs::efd>().connect(&efd_node);

        exciter.allocate();
        efd_node.init(1.2);
        success *= (exciter.verify().passed());
        success *= (exciter.initialize() == 0);
        exciter.evaluateResidual();

        success       *= efd_node.linked();
        success       *= (efd_node.getVariableIndex() == 1);
        success       *= isEqual(efd_node.read(), static_cast<ScalarT>(1.2), kTol);
        const auto* y  = exciter.y().getData();
        success       *= isEqual(y[0], static_cast<ScalarT>(-0.048), kTol);
        success       *= isEqual(y[1], static_cast<ScalarT>(1.2), kTol);
        success       *= isEqual(y[2], static_cast<ScalarT>(0.12), kTol);

        const auto& f      = exciter.getResidual();
        const auto* f_data = f.getData();
        for (size_t i = 0; i < f.getSize(); ++i)
        {
          if (!isEqual(f_data[i], static_cast<ScalarT>(0.0), kTol))
          {
            std::cout << "Non-zero SEXS-PTI residual at index " << i << ": " << f_data[i] << "\n";
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
        efd_node.link(&efd_value, &efd_index);
        vs_node.link(&vs_value, &vs_index);

        using namespace GridKit::PhasorDynamics::Exciter;

        auto data    = makeTestData();
        auto exciter = SexsPti<ScalarT, IdxT>(&bus, data);
        exciter.getPorts().out.template port<SexsPtiSignalOutputs::efd>().connect(&efd_node);
        exciter.getPorts().in.template port<SexsPtiSignalInputs::vs>().connect(&vs_node);

        exciter.allocate();
        efd_node.init(1.2);
        success *= (exciter.verify().passed());
        success *= (exciter.initialize() == 0);
        exciter.evaluateResidual();

        // vref absorbs the attached vs at initialization; later changes pass through
        const auto* f  = exciter.getResidual().getData();
        success       *= isEqual(f[2], static_cast<ScalarT>(0.0), kTol);

        vs_value = 0.35;
        exciter.evaluateResidual();
        success *= isEqual(f[2], static_cast<ScalarT>(0.1), kTol);

        return success.report(__func__);
      }

      TestOutcome vrefAndLimiterPorts()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(3.0, 4.0);
        bus.allocate();
        bus.initialize();

        PhasorDynamics::SignalNode<ScalarT, IdxT> efd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> vref_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> vuel_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> voel_node;
        ScalarT                                   efd_value{0.0};
        ScalarT                                   vref_value{0.0};
        ScalarT                                   vuel_value{-0.4};
        ScalarT                                   voel_value{0.2};
        IdxT                                      efd_index  = INVALID_INDEX<IdxT>;
        IdxT                                      vref_index = 5;
        IdxT                                      vuel_index = 6;
        IdxT                                      voel_index = 7;
        efd_node.link(&efd_value, &efd_index);
        vref_node.link(&vref_value, &vref_index);
        vuel_node.link(&vuel_value, &vuel_index);
        voel_node.link(&voel_value, &voel_index);

        auto                                            data = makeTestData();
        PhasorDynamics::Exciter::SexsPti<ScalarT, IdxT> exciter(&bus, data);
        exciter.getPorts().out.template port<PhasorDynamics::Exciter::SexsPtiSignalOutputs::efd>().connect(&efd_node);
        exciter.getPorts().in.template port<PhasorDynamics::Exciter::SexsPtiSignalInputs::vref>().connect(&vref_node);
        exciter.getPorts().in.template port<PhasorDynamics::Exciter::SexsPtiSignalInputs::vuel>().connect(&vuel_node);
        exciter.getPorts().in.template port<PhasorDynamics::Exciter::SexsPtiSignalInputs::voel>().connect(&voel_node);

        exciter.allocate();
        efd_node.init(1.2);
        success *= (exciter.initialize() == 0);
        exciter.evaluateResidual();

        // vref absorbs both limiter inputs; later changes enter with unit gain.
        const auto* f  = exciter.getResidual().getData();
        success       *= isEqual(vref_node.read(), static_cast<ScalarT>(5.32), kTol);
        success       *= isEqual(f[2], static_cast<ScalarT>(0.0), kTol);

        vuel_value = -0.3;
        exciter.evaluateResidual();
        success *= isEqual(f[2], static_cast<ScalarT>(0.1), kTol);

        vuel_value = -0.4;
        voel_value = 0.3;
        exciter.evaluateResidual();
        success *= isEqual(f[2], static_cast<ScalarT>(0.1), kTol);

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
        PhasorDynamics::SignalNode<ScalarT, IdxT>       efd_node;
        exciter.getPorts().out.template port<PhasorDynamics::Exciter::SexsPtiSignalOutputs::efd>().connect(&efd_node);
        exciter.allocate();
        efd_node.init(0.0);
        success *= (exciter.initialize() == 0);

        auto* y  = exciter.y().getData();
        auto* yp = exciter.yp().getData();

        yp[0] = 0.0;
        yp[1] = 0.0;

        // Windup: Efd = 10 is far above Efdmax = 5 and the pre-limit
        // derivative f = +15 drives further past the limit. The indicator
        // saturates to 0, so residual[1] = 0 exactly.
        y[0] = -1.0;
        y[1] = 10.0;
        y[2] = 1.0;
        exciter.y().setDataUpdated();
        exciter.yp().setDataUpdated();
        exciter.evaluateResidual();
        const auto* f  = exciter.getResidual().getData();
        success       *= isEqual(f[1], static_cast<ScalarT>(0.0), kTol);

        // Release: same over-limit Efd, but f = -37.5 < 0 restores toward
        // the interior. The indicator saturates to 1, so residual[1] = f.
        y[0] = 1.0;
        y[1] = 10.0;
        y[2] = 0.0;
        exciter.y().setDataUpdated();
        exciter.evaluateResidual();
        success *= isEqual(f[1], static_cast<ScalarT>(-37.5), kTol);

        // Mirror (windup below Efdmin): Efd = -10 with f = -15 drives further
        // past the lower limit. Indicator saturates to 0, residual[1] = 0.
        y[0] = 1.0;
        y[1] = -10.0;
        y[2] = -1.0;
        exciter.y().setDataUpdated();
        exciter.evaluateResidual();
        success *= isEqual(f[1], static_cast<ScalarT>(0.0), kTol);

        // Mirror (release above Efdmin): Efd = -10 with f = +37.5 pulls back
        // toward the interior. Indicator saturates to 1, residual[1] = f.
        y[0] = -1.0;
        y[1] = -10.0;
        y[2] = 0.0;
        exciter.y().setDataUpdated();
        exciter.evaluateResidual();
        success *= isEqual(f[1], static_cast<ScalarT>(37.5), kTol);

        // Regression guard: Efd barely above Efdmax with a small positive f.
        // Here the sigmoid is not fully saturated, so the residual is small
        // but not exactly zero; it should still be orders of magnitude below
        // f = 0.125 itself.
        y[0] = -0.2575;
        y[1] = 5.05;
        y[2] = 0.0;
        exciter.y().setDataUpdated();
        exciter.evaluateResidual();
        success *= (std::abs(f[1]) < static_cast<ScalarT>(0.1));

        return success.report(__func__);
      }

      /// Every verify() error class, including parameter types and
      /// finiteness, negative floored lags, parameter relationships, bus
      /// ownership, the required efd output, and signal linkage.
      TestOutcome validation()
      {
        TestStatus success = true;

        const auto previous_verbosity = Log::verbosity();
        // Suppress expected errors and warnings from the invalid cases below.
        // Use EVERYTHING to inspect those diagnostics.
        Log::setVerbosity(Log::Verbosity::NONE);

        // SEXS-PTI has no required parameters: an empty parameter set is the
        // documented-default model.
        success *= verifyData(makeTestData()).passed();
        success *= verifyData(Data{}).passed();

        success *= invalidParameterCase(Params::K, 0.0);
        success *= invalidParameterCase(Params::Efdmin, 5.0);
        success *= invalidParameterCase(Params::Tb, -0.1);
        success *= invalidParameterCase(Params::Te, -0.1);

        for (const Params parameter : {Params::Ta,
                                       Params::Tb,
                                       Params::Te,
                                       Params::K,
                                       Params::Efdmax,
                                       Params::Efdmin})
        {
          success *= invalidParameterCase(parameter, std::numeric_limits<RealT>::quiet_NaN());
          success *= invalidParameterCase(parameter, std::numeric_limits<RealT>::infinity());
        }

        // Integer JSON values are accepted for real parameters; booleans are
        // not numeric.
        auto integer_real                  = makeTestData();
        integer_real.parameters[Params::K] = static_cast<IdxT>(10);

        success *= verifyData(integer_real).passed();
        success *= invalidParameterCase(Params::K, true);

        // The lead time constant enters only a numerator, so any sign is valid.
        auto negative_lead                   = makeTestData();
        negative_lead.parameters[Params::Ta] = -0.1;

        success *= verifyData(negative_lead).passed();

        // A model without the required efd output assignment is rejected.
        PhasorDynamics::Bus<ScalarT, IdxT>              bus(1.0, 0.0);
        PhasorDynamics::Exciter::SexsPti<ScalarT, IdxT> unassigned(&bus, makeTestData());
        success *= !unassigned.verify().passed();

        PhasorDynamics::SignalNode<ScalarT, IdxT>       busless_efd;
        PhasorDynamics::Exciter::SexsPti<ScalarT, IdxT> busless(nullptr, makeTestData());
        busless.getPorts().out.template port<PhasorDynamics::Exciter::SexsPtiSignalOutputs::efd>().connect(&busless_efd);
        success *= !busless.verify().passed();

        success *= unlinkedSignalRejected<PhasorDynamics::Exciter::SexsPtiSignalInputs::vref>();
        success *= unlinkedSignalRejected<PhasorDynamics::Exciter::SexsPtiSignalInputs::vs>();
        success *= unlinkedSignalRejected<PhasorDynamics::Exciter::SexsPtiSignalInputs::vuel>();
        success *= unlinkedSignalRejected<PhasorDynamics::Exciter::SexsPtiSignalInputs::voel>();

        // Zero lags use the documented numerical floor, so the model still
        // initializes and no residual row divides by zero. The floor scales
        // the efd row by 1/(Tb Te), so rounding is not compared to zero.
        auto zero_time                   = makeTestData();
        zero_time.parameters[Params::Tb] = 0.0;
        zero_time.parameters[Params::Te] = 0.0;

        PhasorDynamics::Bus<ScalarT, IdxT>              floored_bus(3.0, 4.0);
        PhasorDynamics::Exciter::SexsPti<ScalarT, IdxT> floored(&floored_bus, zero_time);
        PhasorDynamics::SignalNode<ScalarT, IdxT>       floored_efd;
        floored.getPorts().out.template port<PhasorDynamics::Exciter::SexsPtiSignalOutputs::efd>().connect(&floored_efd);
        success *= (floored_bus.allocate() == 0);
        success *= (floored_bus.initialize() == 0);
        success *= (floored.allocate() == 0);
        floored_efd.init(1.2);
        success *= (floored.initialize() == 0);
        success *= (floored.evaluateResidual() == 0);

        const auto& f      = floored.getResidual();
        const auto* f_data = f.getData();
        for (IdxT i = 0; i < f.getSize(); ++i)
        {
          success *= std::isfinite(f_data[i]);
        }

        Log::setVerbosity(previous_verbosity);
        return success.report(__func__);
      }

      TestOutcome systemAssembly()
      {
        using namespace PhasorDynamics;
        using namespace PhasorDynamics::Exciter;
        using Data         = SexsPtiData<RealT, IdxT>;
        using SignalOutput = typename Data::SignalOutputs;
        using Buses        = typename Data::Buses;

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

        auto exciter                              = makeTestData();
        exciter.buses[Buses::bus]                 = 1;
        exciter.signal_outputs[SignalOutput::efd] = 3;
        data.sexspti.push_back(exciter);

        SystemModel<ScalarT, IdxT> system(data);
        success *= (system.allocate() == 0);
        success *= (system.initialize() == 0);
        success *= (system.evaluateResidual() == 0);
        success *= (system.size() == 5);

        auto missing_efd = data;
        missing_efd.sexspti[0].signal_outputs.erase(SignalOutput::efd);
        SystemModel<ScalarT, IdxT> missing_efd_system(missing_efd);
        const auto                 previous_verbosity = Log::verbosity();
        // Suppress the expected missing-EFD configuration error below.
        // Use EVERYTHING to inspect the diagnostic.
        Log::setVerbosity(Log::Verbosity::NONE);
        success *= (!missing_efd_system.verify().passed());
        Log::setVerbosity(previous_verbosity);

        return success.report(__func__);
      }

      TestOutcome systemSignalFreshness()
      {
        using namespace PhasorDynamics;
        using namespace PhasorDynamics::Exciter;
        using Data         = SexsPtiData<RealT, IdxT>;
        using SignalInput  = typename Data::SignalInputs;
        using SignalOutput = typename Data::SignalOutputs;
        using Buses        = typename Data::Buses;

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

        auto consumer                              = makeTestData();
        consumer.disambiguation_string             = "consumer";
        consumer.buses[Buses::bus]                 = 1;
        consumer.signal_outputs[SignalOutput::efd] = 20;
        consumer.signal_inputs[SignalInput::vs]    = 10;
        data.sexspti.push_back(consumer);

        auto source                              = makeTestData();
        source.disambiguation_string             = "source";
        source.buses[Buses::bus]                 = 1;
        source.signal_outputs[SignalOutput::efd] = 10;
        data.sexspti.push_back(source);

        SystemModel<ScalarT, IdxT> system(data);
        success *= (system.allocate() == 0);
        success *= (system.initialize() == 0);

        constexpr IdxT consumer_vtr_residual = 4;
        constexpr IdxT source_efd_variable   = 6;
        auto*          y                     = system.y().getData();
        y[source_efd_variable]               = 0.75;
        system.y().setDataUpdated();
        success       *= (system.evaluateResidual() == 0);
        const auto* f  = system.getResidual().getData();
        success       *= isEqual(f[consumer_vtr_residual],
                           static_cast<ScalarT>(0.75),
                           kTol);

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /**
       * @brief Checks Jacobian evaluation.
       */
      TestOutcome jacobian()
      {
        TestStatus success = true;

        auto tol = 10 * std::numeric_limits<RealT>::epsilon();

        // Jacobian via DependencyTracking
        std::vector<DependencyTracking::Variable::DependencyMap> dependency_tracking_jacobian = DependencyTrackingJacobian();

        // Jacobian via Enzyme
        std::vector<DependencyTracking::Variable::DependencyMap> enzyme_jacobian = EnzymeJacobian();

        /// Compare DependencyTracking dependencies to Enzyme's
        for (size_t i = 0; i < dependency_tracking_jacobian.size(); ++i)
        {
          success *= (GridKit::Testing::isEqual(dependency_tracking_jacobian[i], enzyme_jacobian[i], tol));
        }

        return success.report(__func__);
      }

    private:
      std::vector<DependencyTracking::Variable::DependencyMap> DependencyTrackingJacobian()
      {
        auto data = makeTestData();

        DependencyTracking::Variable                                         Vr1{3.0};
        DependencyTracking::Variable                                         Vi1{4.0};
        PhasorDynamics::Bus<DependencyTracking::Variable, IdxT>              bus(Vr1, Vi1);
        PhasorDynamics::Exciter::SexsPti<DependencyTracking::Variable, IdxT> exciter(&bus, data);
        PhasorDynamics::SignalNode<DependencyTracking::Variable, IdxT>       efd_node;
        exciter.getPorts().out.template port<PhasorDynamics::Exciter::SexsPtiSignalOutputs::efd>().connect(&efd_node);

        bus.allocate();
        exciter.allocate();

        for (size_t i = 0; i < bus.size(); ++i)
        {
          bus.setVariableIndex(i, i + exciter.size()); // Reset bus variable indices
          bus.setResidualIndex(i, i + exciter.size()); // Reset bus residual indices
        }

        bus.initialize();
        efd_node.init(DependencyTracking::Variable{1.2});
        exciter.initialize();

        exciter.updateTime(0.0, 1.0); // Set alpha to 1.0 to verify d/dy' term

        bus.evaluateResidual();
        exciter.evaluateResidual();

        bus.evaluateJacobian();
        exciter.evaluateJacobian();
        auto* model_jacobian = exciter.getCsrJacobian();
        auto& output_stream  = Log::misc();
        output_stream << "Sparse Csr Matrix: SexsPti DependencyTracking Jacobian\n";
        model_jacobian->print(output_stream);

        return GridKit::Testing::MapFromCsr(model_jacobian);
      }

      std::vector<DependencyTracking::Variable::DependencyMap> EnzymeJacobian()
      {
        auto data = makeTestData();

        PhasorDynamics::Bus<ScalarT, IdxT>              bus(3.0, 4.0);
        PhasorDynamics::Exciter::SexsPti<ScalarT, IdxT> exciter(&bus, data);
        PhasorDynamics::SignalNode<ScalarT, IdxT>       efd_node;
        exciter.getPorts().out.template port<PhasorDynamics::Exciter::SexsPtiSignalOutputs::efd>().connect(&efd_node);

        bus.allocate();
        exciter.allocate();

        for (size_t i = 0; i < bus.size(); ++i)
        {
          bus.setVariableIndex(i, i + exciter.size()); // Reset bus variable indices
          bus.setResidualIndex(i, i + exciter.size()); // Reset bus residual indices
        }

        bus.initialize();
        efd_node.init(1.2);
        exciter.initialize();

        exciter.updateTime(0.0, 1.0); // Set alpha to 1.0 to verify d/dy' term

        bus.evaluateResidual();
        exciter.evaluateResidual();

        bus.evaluateJacobian();
        exciter.evaluateJacobian();
        exciter.constructCsr();
        auto* model_jacobian = exciter.getCsrJacobian();
        auto& output_stream  = Log::misc();
        output_stream << "Sparse Csr Matrix: SexsPti Enzyme Jacobian\n";
        model_jacobian->print(output_stream);

        return GridKit::Testing::MapFromCsr(model_jacobian);
      }
#endif

    private:
      /// Verify a model built from @p data with its required efd output assigned.
      Model::ConfigurationChecks verifyData(const Data& data) const
      {
        PhasorDynamics::Bus<ScalarT, IdxT>              bus(1.0, 0.0);
        PhasorDynamics::SignalNode<ScalarT, IdxT>       efd_node;
        PhasorDynamics::Exciter::SexsPti<ScalarT, IdxT> exciter(&bus, data);
        exciter.getPorts().out.template port<PhasorDynamics::Exciter::SexsPtiSignalOutputs::efd>().connect(&efd_node);
        return exciter.verify();
      }

      template <typename ValueT>
      bool invalidParameterCase(Params parameter, ValueT value) const
      {
        auto data                  = makeTestData();
        data.parameters[parameter] = value;
        return configurationRejected(
            [&]()
            {
              return verifyData(data);
            });
      }

      template <PhasorDynamics::Exciter::SexsPtiSignalInputs input>
      bool unlinkedSignalRejected() const
      {
        PhasorDynamics::Bus<ScalarT, IdxT>              bus(1.0, 0.0);
        PhasorDynamics::SignalNode<ScalarT, IdxT>       efd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT>       unlinked_node;
        PhasorDynamics::Exciter::SexsPti<ScalarT, IdxT> exciter(&bus, makeTestData());
        exciter.getPorts().out.template port<PhasorDynamics::Exciter::SexsPtiSignalOutputs::efd>().connect(&efd_node);
        exciter.getPorts().in.template port<input>().connect(&unlinked_node);
        return !exciter.verify().passed();
      }

      auto makeTestData() const -> Data
      {
        Data data;
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
