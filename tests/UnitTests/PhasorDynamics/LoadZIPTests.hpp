#pragma once

#include <iostream>
#include <sstream>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <GridKit/Model/PhasorDynamics/Load/LoadZIP/LoadZIP.hpp>
#include <GridKit/Model/VariableMonitorController.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class LoadZIPTests
    {
    public:
      using RealT = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;
      using DataT = PhasorDynamics::LoadZIPData<RealT, IdxT>;

      LoadZIPTests()  = default;
      ~LoadZIPTests() = default;

      TestOutcome constructor()
      {
        TestStatus success = true;

        auto* bus = new PhasorDynamics::Bus<ScalarT, IdxT>(1.0, 0.0);

        PhasorDynamics::Component<ScalarT, IdxT>* load =
            new PhasorDynamics::LoadZIP<ScalarT, IdxT>(bus);

        success *= (load != nullptr);

        if (load)
        {
          success *= (load->size() == 0);
          success *= (load->allocate() == 0);
          success *= (load->y().getSize() == 0);
          success *= (load->yp().getSize() == 0);
          success *= (load->getResidual().getSize() == 0);
          success *= (load->absoluteTolerance().getSize() == 0);
          delete load;
        }

        auto                                      data = makeData();
        PhasorDynamics::LoadZIP<ScalarT, IdxT>    monitored_load(bus, data);
        PhasorDynamics::Component<ScalarT, IdxT>& monitored_component  = monitored_load;
        success                                                       *= (monitored_component.getMonitor() != nullptr);

        delete bus;

        return success.report(__func__);
      }

      TestOutcome initialization()
      {
        TestStatus success = true;

        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus(0.3, 0.4);
        PhasorDynamics::LoadZIP<ScalarT, IdxT>     load(&bus, 2.0, 0.5, 0.2, 0.4);

        bus.allocate();
        load.allocate();

        success *= bus.initialize() == 0;
        success *= load.initialize() == 0;

        success *= (load.size() == 0);
        success *= (load.y().getSize() == 0);
        success *= (load.yp().getSize() == 0);
        success *= (load.getResidual().getSize() == 0);
        success *= (load.absoluteTolerance().getSize() == 0);

        return success.report(__func__);
      }

      /// The ZIP anchor is the bus voltage sampled at initialization, so the
      /// load consumes exactly Pnom and Qnom at every initialized voltage.
      TestOutcome dispatchAtInitializedVoltage()
      {
        TestStatus success = true;

        const RealT Pnom{2.0};
        const RealT Qnom{0.5};

        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus(1.2, 0.9);
        PhasorDynamics::LoadZIP<ScalarT, IdxT>     load(&bus, Pnom, Qnom, 0.2, 0.4);

        bus.allocate();
        load.allocate();

        bus.initialize();
        success *= load.initialize() == 0;

        bus.evaluateResidual();
        load.evaluateResidual();
        ScalarT p  = bus.Vr() * bus.Ir() + bus.Vi() * bus.Ii();
        ScalarT q  = bus.Vi() * bus.Ir() - bus.Vr() * bus.Ii();
        success   *= isEqual(p, static_cast<ScalarT>(-Pnom), tol_);
        success   *= isEqual(q, static_cast<ScalarT>(-Qnom), tol_);

        // Reinitializing at a different voltage anchors the same dispatch
        // there.
        bus.Vr() = 0.6;
        bus.Vi() = 0.8;
        bus.y().setDataUpdated();
        success *= load.initialize() == 0;

        bus.evaluateResidual();
        load.evaluateResidual();
        p        = bus.Vr() * bus.Ir() + bus.Vi() * bus.Ii();
        q        = bus.Vi() * bus.Ir() - bus.Vr() * bus.Ii();
        success *= isEqual(p, static_cast<ScalarT>(-Pnom), tol_);
        success *= isEqual(q, static_cast<ScalarT>(-Qnom), tol_);

        // A zero initialization voltage has no anchor and must fail.
        bus.Vr() = 0.0;
        bus.Vi() = 0.0;
        bus.y().setDataUpdated();
        success *= load.initialize() != 0;

        return success.report(__func__);
      }

      TestOutcome residual()
      {
        TestStatus success = true;

        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus(0.3, 0.4);
        PhasorDynamics::LoadZIP<ScalarT, IdxT>     load(&bus, 2.0, 0.5, 0.2, 0.4);

        bus.allocate();
        load.allocate();

        bus.initialize();
        load.initialize();

        bus.Vr() = 0.9;
        bus.Vi() = 1.2;
        bus.y().setDataUpdated();

        bus.evaluateResidual();
        load.evaluateResidual();

        success *= isEqual(bus.Ir(), static_cast<ScalarT>(-368.0 / 75.0), tol_);
        success *= isEqual(bus.Ii(), static_cast<ScalarT>(-299.0 / 75.0), tol_);

        return success.report(__func__);
      }

      TestOutcome monitor()
      {
        TestStatus success = true;

        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus(0.3, 0.4);
        auto                                       data = makeData();
        PhasorDynamics::LoadZIP<ScalarT, IdxT>     load(&bus, data);

        bus.allocate();
        load.allocate();

        bus.initialize();
        load.initialize();

        bus.Vr() = 0.9;
        bus.Vi() = 1.2;

        RealT                                     time = 0.0;
        Model::VariableMonitorController<ScalarT> controller(time);
        PhasorDynamics::Component<ScalarT, IdxT>& component = load;
        controller.addMonitor(component.getMonitor());

        std::stringstream os;
        controller.addSink({Model::VariableMonitorFormat::CSV}, os);
        controller.print();

        auto values = Tokenizer<RealT>(os.str(), ',')();
        if (values.size() == 6)
        {
          const RealT ir  = -368.0 / 75.0;
          const RealT ii  = -299.0 / 75.0;
          success        *= isEqual(values[1], ir, tol_);
          success        *= isEqual(values[2], ii, tol_);
          success        *= isEqual(values[3], std::sqrt(ir * ir + ii * ii), tol_);
          success        *= isEqual(values[4], -9.2, tol_);
          success        *= isEqual(values[5], -2.3, tol_);
        }
        else
        {
          success = false;
        }

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      TestOutcome jacobian()
      {
        TestStatus success = true;

        const RealT Pnom{2.0};
        const RealT Qnom{0.5};
        const RealT alphaI{4.0};
        const RealT alphaP{2.0};

        auto                                                           dependency_tracking_jacobian = DependencyTrackingJacobian(Pnom, Qnom, alphaI, alphaP);
        auto                                                           enzyme_jacobian              = EnzymeJacobian(Pnom, Qnom, alphaI, alphaP);
        const std::vector<DependencyTracking::Variable::DependencyMap> expected                     = {
            {{0, 22.72}, {1, 38.96}},
            {{0, 26.96}, {1, 25.28}}};

        for (size_t i = 0; i < expected.size(); ++i)
        {
          success *= GridKit::Testing::isEqual(dependency_tracking_jacobian[i], expected[i], tol_);
          success *= GridKit::Testing::isEqual(enzyme_jacobian[i], expected[i], tol_);
        }

        return success.report(__func__);
      }

    private:
      std::vector<DependencyTracking::Variable::DependencyMap> DependencyTrackingJacobian(
          const RealT Pnom, const RealT Qnom, const RealT alphaI, const RealT alphaP)
      {
        DependencyTracking::Variable Vr{0.3};
        DependencyTracking::Variable Vi{0.4};

        PhasorDynamics::Bus<DependencyTracking::Variable, IdxT>     bus(Vr, Vi);
        PhasorDynamics::LoadZIP<DependencyTracking::Variable, IdxT> load(&bus, Pnom, Qnom, alphaI, alphaP);

        bus.allocate();
        load.allocate();

        bus.initialize();
        load.initialize();

        auto* bus_y = bus.y().getData();
        for (size_t i = 0; i < bus.size(); ++i)
        {
          bus_y[i].setVariableNumber(i);
        }
        bus.y().setDataUpdated();

        bus.evaluateResidual();
        load.evaluateResidual(); ///< Computes the residual and the Jacobian values by tracking
                                 ///< the dependencies
        const auto& residual = bus.getResidual();

        std::vector<DependencyTracking::Variable::DependencyMap> dependencies(residual.getSize());
        for (IdxT i = 0; i < residual.getSize(); ++i)
        {
          dependencies[i] = residual.getData()[i].getDependencies();
        }

        return dependencies;
      }

      std::vector<DependencyTracking::Variable::DependencyMap> EnzymeJacobian(
          const RealT Pnom, const RealT Qnom, const RealT alphaI, const RealT alphaP)
      {
        ScalarT Vr{0.3};
        ScalarT Vi{0.4};

        PhasorDynamics::Bus<ScalarT, IdxT>     bus(Vr, Vi);
        PhasorDynamics::LoadZIP<ScalarT, IdxT> load(&bus, Pnom, Qnom, alphaI, alphaP);

        bus.allocate();
        load.allocate();

        bus.initialize();
        load.initialize();

        load.updateTime(0.0, 1.0);

        for (size_t i = 0; i < bus.size(); ++i)
        {
          bus.setVariableIndex(i, i);
          bus.setResidualIndex(i, i);
        }

        bus.evaluateResidual();
        load.evaluateResidual();

        bus.evaluateJacobian();
        load.evaluateJacobian();
        load.constructCsr();
        GridKit::LinearAlgebra::CsrMatrix<ScalarT, IdxT>* model_jacobian = load.getCsrJacobian();
        std::cout << "Sparse Csr Matrix: LoadZIP Jacobian\n";
        model_jacobian->print();

        return GridKit::Testing::MapFromCsr(model_jacobian);
      }
#endif

    private:
      static constexpr RealT tol_ = 1.0e-10;

      auto makeData() -> DataT
      {
        using Params   = PhasorDynamics::LoadZIPParameters;
        using Variable = PhasorDynamics::LoadZIPMonitorableVariables;

        DataT data;
        data.device_class          = "LoadZIP";
        data.disambiguation_string = "loadzip_test";

        data.parameters[Params::Pnom]   = static_cast<RealT>(2.0);
        data.parameters[Params::Qnom]   = static_cast<RealT>(0.5);
        data.parameters[Params::alphaI] = static_cast<RealT>(0.2);
        data.parameters[Params::alphaP] = static_cast<RealT>(0.4);

        data.monitored_variables.insert(Variable::ir);
        data.monitored_variables.insert(Variable::ii);
        data.monitored_variables.insert(Variable::im);
        data.monitored_variables.insert(Variable::p);
        data.monitored_variables.insert(Variable::q);

        return data;
      }
    }; // class LoadZIPTests

  } // namespace Testing
} // namespace GridKit
