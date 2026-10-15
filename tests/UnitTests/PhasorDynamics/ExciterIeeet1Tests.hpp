#pragma once

#include <iostream>
#include <limits>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/ConfigurationChecks.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/IEEET1/Ieeet1.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/IEEET1/Ieeet1Data.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    using Log = ::GridKit::Utilities::Logger;

    template <class ScalarT, typename IdxT>
    class ExciterIeeet1Tests
    {
    public:
      using RealT  = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;
      using Data   = PhasorDynamics::Exciter::Ieeet1Data<RealT, IdxT>;
      using Params = PhasorDynamics::Exciter::Ieeet1Parameters;

      ExciterIeeet1Tests()  = default;
      ~ExciterIeeet1Tests() = default;

      static constexpr RealT kTol =
          static_cast<RealT>(4.0) * std::numeric_limits<RealT>::epsilon();

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

      TestOutcome zeroInitialResidual()
      {
        TestStatus success = true;

        auto data = makeTestData();

        PhasorDynamics::Bus<ScalarT, IdxT>             bus(3.0, 4.0);
        PhasorDynamics::Exciter::Ieeet1<ScalarT, IdxT> exciter(&bus, data);
        PhasorDynamics::SignalNode<ScalarT, IdxT>      efd_node;
        exciter.getPorts().out.template port<PhasorDynamics::Exciter::Ieeet1SignalOutputs::efd>().connect(&efd_node);

        success *= (bus.allocate() == 0);
        success *= (exciter.allocate() == 0);

        success *= (bus.initialize() == 0);
        efd_node.init(1.2);
        success *= (exciter.initialize() == 0);

        success *= (exciter.evaluateResidual() == 0);

        const auto& residual      = exciter.getResidual();
        const auto* residual_data = residual.getData();
        for (size_t i = 0; i < residual.getSize(); ++i)
        {
          if (!isEqual(residual_data[i], static_cast<ScalarT>(0.0), kTol))
          {
            std::cout << "Non-zero Ieeet1 residual at index " << i << ": " << residual_data[i] << "\n";
            success = false;
          }
        }

        return success.report(__func__);
      }

      TestOutcome zeroTimeConstantsAndDisabledSaturation()
      {
        TestStatus success = true;

        auto data                    = makeTestData();
        data.parameters[Params::Tr]  = 0.0;
        data.parameters[Params::Ta]  = 0.0;
        data.parameters[Params::Ke]  = 0.0;
        data.parameters[Params::Te]  = 0.0;
        data.parameters[Params::Kf]  = 0.0;
        data.parameters[Params::Tf]  = 0.0;
        data.parameters[Params::E1]  = 0.0;
        data.parameters[Params::E2]  = 0.0;
        data.parameters[Params::Se1] = 0.0;
        data.parameters[Params::Se2] = 0.0;

        PhasorDynamics::Bus<ScalarT, IdxT>             bus(3.0, 4.0);
        PhasorDynamics::Exciter::Ieeet1<ScalarT, IdxT> exciter(&bus, data);
        PhasorDynamics::SignalNode<ScalarT, IdxT>      efd_node;
        ScalarT                                        efd_value{0.0};
        IdxT                                           efd_index = INVALID_INDEX<IdxT>;

        efd_node.link(&efd_value, &efd_index);
        exciter.getPorts().out.template port<PhasorDynamics::Exciter::Ieeet1SignalOutputs::efd>().connect(&efd_node);

        bus.allocate();
        exciter.allocate();

        bus.initialize();
        efd_node.init(1.2);
        success *= (exciter.initialize() == 0);
        exciter.tagDifferentiable();

        success *= (exciter.tag()[0]);
        success *= (exciter.tag()[1]);
        success *= (exciter.tag()[2]);
        success *= (exciter.tag()[3]);

        auto*       y  = exciter.y().getData();
        auto*       yp = exciter.yp().getData();
        const auto* f  = exciter.getResidual().getData();

        success *= isEqual(y[2], static_cast<ScalarT>(1.2));
        success *= isEqual(y[6], static_cast<ScalarT>(0.0));
        success *= isEqual(y[7], static_cast<ScalarT>(1.2));
        success *= isEqual(y[8], static_cast<ScalarT>(0.0));

        for (IdxT i = 0; i < exciter.y().getSize(); ++i)
        {
          success *= std::isfinite(y[i]);
        }

        exciter.evaluateResidual();
        for (IdxT i = 0; i < exciter.getResidual().getSize(); ++i)
        {
          success *= std::isfinite(f[i]);
          success *= isEqual(f[i], static_cast<ScalarT>(0.0), 100.0 * std::numeric_limits<RealT>::epsilon());
        }

        y[2] = 4.0;
        exciter.y().setDataUpdated();
        exciter.evaluateResidual();
        success *= isEqual(f[8], static_cast<ScalarT>(0.0));
        y[2]     = 1.2;

        yp[0] = 123.0;
        exciter.y().setDataUpdated();
        exciter.yp().setDataUpdated();
        exciter.evaluateResidual();
        success *= isEqual(f[0], static_cast<ScalarT>(-123.0));
        yp[0]    = 0.0;

        y[0] = 4.0;
        exciter.y().setDataUpdated();
        exciter.yp().setDataUpdated();
        exciter.evaluateResidual();
        success *= isEqual(f[0], static_cast<ScalarT>(1.0e3));

        y[0] = 5.0;
        y[4] = 0.02;
        exciter.y().setDataUpdated();
        exciter.evaluateResidual();
        success *= isEqual(f[1], static_cast<ScalarT>(900.0));

        y[4] = 0.0;
        y[1] = 1.0;
        exciter.y().setDataUpdated();
        exciter.evaluateResidual();
        success *= isEqual(f[2], static_cast<ScalarT>(900.0));

        y[1] = 0.0;
        y[5] = 1.0;
        exciter.y().setDataUpdated();
        exciter.evaluateResidual();
        success *= isEqual(f[3], static_cast<ScalarT>(1.0e3));

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

        // IEEET1 has no required parameters: an empty parameter set is the
        // documented-default model.
        success *= verifyData(makeTestData()).passed();
        success *= verifyData(Data{}).passed();

        success *= invalidParameterCase(Params::Ka, 0.0);
        success *= invalidParameterCase(Params::Vrmin, 2.0);
        success *= invalidParameterCase(Params::Ispdlim, 0.5);
        success *= invalidParameterCase(Params::Tr, -0.1);
        success *= invalidParameterCase(Params::Ta, -0.1);
        success *= invalidParameterCase(Params::Te, -0.1);
        success *= invalidParameterCase(Params::Tf, -0.1);
        success *= invalidParameterCase(Params::Se1, -0.1);
        success *= invalidParameterCase(Params::E1, -1.0);
        success *= invalidParameterCase(Params::Se2, 0.01);

        for (const Params parameter : {Params::Tr,
                                       Params::Ka,
                                       Params::Ta,
                                       Params::Ke,
                                       Params::Te,
                                       Params::Kf,
                                       Params::Tf,
                                       Params::Vrmin,
                                       Params::Vrmax,
                                       Params::E1,
                                       Params::E2,
                                       Params::Se1,
                                       Params::Se2,
                                       Params::Ispdlim})
        {
          success *= invalidParameterCase(parameter, std::numeric_limits<RealT>::quiet_NaN());
          success *= invalidParameterCase(parameter, std::numeric_limits<RealT>::infinity());
        }

        // Integer JSON values are accepted for real parameters; booleans are
        // not numeric.
        auto integer_real                   = makeTestData();
        integer_real.parameters[Params::Ka] = static_cast<IdxT>(50);

        success *= verifyData(integer_real).passed();
        success *= invalidParameterCase(Params::Ka, true);

        // A model without the required efd output assignment is rejected.
        PhasorDynamics::Bus<ScalarT, IdxT>             bus(1.0, 0.0);
        PhasorDynamics::Exciter::Ieeet1<ScalarT, IdxT> unassigned(&bus, makeTestData());
        success *= !unassigned.verify().passed();

        PhasorDynamics::SignalNode<ScalarT, IdxT>      busless_efd;
        PhasorDynamics::Exciter::Ieeet1<ScalarT, IdxT> busless(nullptr, makeTestData());
        busless.getPorts().out.template port<PhasorDynamics::Exciter::Ieeet1SignalOutputs::efd>().connect(&busless_efd);
        success *= !busless.verify().passed();

        success *= unlinkedSignalRejected<PhasorDynamics::Exciter::Ieeet1SignalInputs::speed>();
        success *= unlinkedSignalRejected<PhasorDynamics::Exciter::Ieeet1SignalInputs::vref>();
        success *= unlinkedSignalRejected<PhasorDynamics::Exciter::Ieeet1SignalInputs::vs>();
        success *= unlinkedSignalRejected<PhasorDynamics::Exciter::Ieeet1SignalInputs::vuel>();
        success *= unlinkedSignalRejected<PhasorDynamics::Exciter::Ieeet1SignalInputs::voel>();

        Log::setVerbosity(previous_verbosity);
        return success.report(__func__);
      }

      TestOutcome vrefAndLimiterPorts()
      {
        TestStatus success = true;

        auto data = makeTestData();

        PhasorDynamics::Bus<ScalarT, IdxT>             bus(3.0, 4.0);
        PhasorDynamics::Exciter::Ieeet1<ScalarT, IdxT> exciter(&bus, data);
        PhasorDynamics::SignalNode<ScalarT, IdxT>      efd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT>      vref_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT>      vuel_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT>      voel_node;
        ScalarT                                        efd_value{0.0};
        ScalarT                                        vref_value{0.0};
        ScalarT                                        vuel_value{-0.4};
        ScalarT                                        voel_value{0.2};
        IdxT                                           efd_index  = INVALID_INDEX<IdxT>;
        IdxT                                           vref_index = 11;
        IdxT                                           vuel_index = 12;
        IdxT                                           voel_index = 13;

        efd_node.link(&efd_value, &efd_index);
        vref_node.link(&vref_value, &vref_index);
        vuel_node.link(&vuel_value, &vuel_index);
        voel_node.link(&voel_value, &voel_index);
        exciter.getPorts().out.template port<PhasorDynamics::Exciter::Ieeet1SignalOutputs::efd>().connect(&efd_node);
        exciter.getPorts().in.template port<PhasorDynamics::Exciter::Ieeet1SignalInputs::vref>().connect(&vref_node);
        exciter.getPorts().in.template port<PhasorDynamics::Exciter::Ieeet1SignalInputs::vuel>().connect(&vuel_node);
        exciter.getPorts().in.template port<PhasorDynamics::Exciter::Ieeet1SignalInputs::voel>().connect(&voel_node);

        bus.allocate();
        exciter.allocate();

        bus.initialize();
        efd_node.init(1.2);
        success *= (exciter.initialize() == 0);
        exciter.evaluateResidual();

        const auto* y = exciter.y().getData();
        const auto* f = exciter.getResidual().getData();

        // vref absorbs both limiter inputs; later changes enter with unit gain.
        success *= isEqual(vref_node.read(),
                           y[0] + y[4] + y[5] - vuel_value - voel_value,
                           kTol);
        success *= isEqual(f[4], static_cast<ScalarT>(0.0), kTol);

        vuel_value = -0.3;
        exciter.evaluateResidual();
        success *= isEqual(f[4], static_cast<ScalarT>(0.1), kTol);

        vuel_value = -0.4;
        voel_value = 0.3;
        exciter.evaluateResidual();
        success *= isEqual(f[4], static_cast<ScalarT>(0.1), kTol);

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

        // Compare DependencyTracking dependencies to Enzyme's
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

        DependencyTracking::Variable                                        Vr1{3.0};
        DependencyTracking::Variable                                        Vi1{4.0};
        PhasorDynamics::Bus<DependencyTracking::Variable, IdxT>             bus(Vr1, Vi1);
        PhasorDynamics::Exciter::Ieeet1<DependencyTracking::Variable, IdxT> exciter(&bus, data);
        PhasorDynamics::SignalNode<DependencyTracking::Variable, IdxT>      efd_node;
        exciter.getPorts().out.template port<PhasorDynamics::Exciter::Ieeet1SignalOutputs::efd>().connect(&efd_node);

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
        output_stream << "Sparse Csr Matrix: Ieeet1 DependencyTracking Jacobian\n";
        model_jacobian->print(output_stream);

        return GridKit::Testing::MapFromCsr(model_jacobian);
      }

      std::vector<DependencyTracking::Variable::DependencyMap> EnzymeJacobian()
      {
        auto data = makeTestData();

        PhasorDynamics::Bus<ScalarT, IdxT>             bus(3.0, 4.0);
        PhasorDynamics::Exciter::Ieeet1<ScalarT, IdxT> exciter(&bus, data);
        PhasorDynamics::SignalNode<ScalarT, IdxT>      efd_node;
        exciter.getPorts().out.template port<PhasorDynamics::Exciter::Ieeet1SignalOutputs::efd>().connect(&efd_node);

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
        output_stream << "Sparse Csr Matrix: Ieeet1 Enzyme Jacobian\n";
        model_jacobian->print(output_stream);

        return GridKit::Testing::MapFromCsr(model_jacobian);
      }
#endif

    private:
      /// Verify a model built from @p data with its required efd output assigned.
      Model::ConfigurationChecks verifyData(const Data& data) const
      {
        PhasorDynamics::Bus<ScalarT, IdxT>             bus(1.0, 0.0);
        PhasorDynamics::SignalNode<ScalarT, IdxT>      efd_node;
        PhasorDynamics::Exciter::Ieeet1<ScalarT, IdxT> exciter(&bus, data);
        exciter.getPorts().out.template port<PhasorDynamics::Exciter::Ieeet1SignalOutputs::efd>().connect(&efd_node);
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

      template <PhasorDynamics::Exciter::Ieeet1SignalInputs input>
      bool unlinkedSignalRejected() const
      {
        PhasorDynamics::Bus<ScalarT, IdxT>             bus(1.0, 0.0);
        PhasorDynamics::SignalNode<ScalarT, IdxT>      efd_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT>      unlinked_node;
        PhasorDynamics::Exciter::Ieeet1<ScalarT, IdxT> exciter(&bus, makeTestData());
        exciter.getPorts().out.template port<PhasorDynamics::Exciter::Ieeet1SignalOutputs::efd>().connect(&efd_node);
        exciter.getPorts().in.template port<input>().connect(&unlinked_node);
        return !exciter.verify().passed();
      }

      auto makeTestData() const -> Data
      {
        Data data;
        data.device_class          = "exciter";
        data.disambiguation_string = "ieeet1_test";
        data.monitored_variables.insert(PhasorDynamics::Exciter::Ieeet1MonitorableVariables::efd);

        data.parameters[Params::Tr]      = 0.0;
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
