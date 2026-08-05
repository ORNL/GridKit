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

        bus.initialize();
        load.initialize();

        const auto* y   = load.y().getData();
        const auto* yp  = load.yp().getData();
        success        *= isEqual(y[0], static_cast<ScalarT>(-3.2), tol_);
        success        *= isEqual(y[1], static_cast<ScalarT>(-2.6), tol_);
        success        *= isEqual(yp[0], static_cast<ScalarT>(0.0), tol_);
        success        *= isEqual(yp[1], static_cast<ScalarT>(0.0), tol_);

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

        const auto* y  = load.y().getData();
        ScalarT     p  = bus.Vr() * y[0] + bus.Vi() * y[1];
        ScalarT     q  = bus.Vi() * y[0] - bus.Vr() * y[1];
        success       *= isEqual(p, static_cast<ScalarT>(-Pnom), tol_);
        success       *= isEqual(q, static_cast<ScalarT>(-Qnom), tol_);

        // Reinitializing at a different voltage anchors the same dispatch
        // there.
        bus.Vr() = 0.6;
        bus.Vi() = 0.8;
        bus.y().setDataUpdated();
        success *= load.initialize() == 0;

        p        = bus.Vr() * y[0] + bus.Vi() * y[1];
        q        = bus.Vi() * y[0] - bus.Vr() * y[1];
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

        const auto* f  = load.getResidual().getData();
        success       *= isEqual(f[0], static_cast<ScalarT>(128.0 / 75.0), tol_);
        success       *= isEqual(f[1], static_cast<ScalarT>(104.0 / 75.0), tol_);
        success       *= isEqual(bus.Ir(), static_cast<ScalarT>(-3.2), tol_);
        success       *= isEqual(bus.Ii(), static_cast<ScalarT>(-2.6), tol_);

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
          success *= isEqual(values[1], -3.2, tol_);
          success *= isEqual(values[2], -2.6, tol_);
          success *= isEqual(values[3], std::sqrt(17.0), tol_);
          success *= isEqual(values[4], -2.0, tol_);
          success *= isEqual(values[5], -0.5, tol_);
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

        auto dependency_tracking_jacobian = DependencyTrackingJacobian(Pnom, Qnom, alphaI, alphaP);
        auto enzyme_jacobian              = EnzymeJacobian(Pnom, Qnom, alphaI, alphaP);

        for (size_t i = 0; i < dependency_tracking_jacobian.size(); ++i)
        {
          success *= GridKit::Testing::isEqual(dependency_tracking_jacobian[i], enzyme_jacobian[i]);
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

        auto* load_y = load.y().getData();
        for (size_t i = 0; i < load.size(); ++i)
        {
          load_y[i].setVariableNumber(i);
        }
        load.y().setDataUpdated();
        auto* bus_y = bus.y().getData();
        for (size_t i = 0; i < bus.size(); ++i)
        {
          bus_y[i].setVariableNumber(i + load.size());
        }
        bus.y().setDataUpdated();

        bus.evaluateResidual();
        load.evaluateResidual(); ///< Computes the residual and the Jacobian values by tracking
                                 ///< the dependencies
        auto&                                     residual_y_view = load.getResidual();
        std::vector<DependencyTracking::Variable> residual_y(residual_y_view.getData(), residual_y_view.getData() + residual_y_view.getSize());

        bus.initialize();
        load.initialize();

        auto* load_yp = load.yp().getData();
        for (size_t i = 0; i < load.size(); ++i)
        {
          load_yp[i].setVariableNumber(i);
        }
        load.yp().setDataUpdated();

        bus.evaluateResidual();
        load.evaluateResidual(); ///< Computes the residual and the Jacobian values by tracking
                                 ///< the dependencies
        auto&                                     residual_yp_view = load.getResidual();
        std::vector<DependencyTracking::Variable> residual_yp(residual_yp_view.getData(), residual_yp_view.getData() + residual_yp_view.getSize());

        std::vector<DependencyTracking::Variable::DependencyMap> dependencies(residual_y.size());
        for (IdxT i = 0; i < residual_y.size(); ++i)
        {
          auto dependency_y  = residual_y[i].getDependencies();
          auto dependency_yp = residual_yp[i].getDependencies();

          for (const auto& pair_y : dependency_y)
          {
            auto it_yp = dependency_yp.find(pair_y.first);
            if (it_yp != dependency_yp.end())
            {
              dependencies[i].insert(std::make_pair(pair_y.first, pair_y.second + it_yp->second));
            }
            else
            {
              dependencies[i].insert(std::make_pair(pair_y.first, pair_y.second));
            }
          }

          for (const auto& pair_yp : dependency_yp)
          {
            if (!dependency_y.contains(pair_yp.first))
            {
              dependencies[i].insert(std::make_pair(pair_yp.first, pair_yp.second));
            }
          }
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
          bus.setVariableIndex(i, i + load.size());
          bus.setResidualIndex(i, i + load.size());
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
