#pragma once

#include <iomanip>
#include <iostream>
#include <sstream>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <GridKit/Model/PhasorDynamics/Load/LoadZ/LoadZ.hpp>
#include <GridKit/Model/VariableMonitorController.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <typename scalar_type, typename index_type>
    class LoadZTests
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;
      using DataT   = PhasorDynamics::LoadZData<RealT, IdxT>;

      LoadZTests()  = default;
      ~LoadZTests() = default;

      TestOutcome constructor()
      {
        TestStatus success = true;

        auto* bus = new PhasorDynamics::Bus<ScalarT, IdxT>(1.0, 0.0);

        PhasorDynamics::Component<ScalarT, IdxT>* load =
            new PhasorDynamics::LoadZ<ScalarT, IdxT>(bus);

        success *= (load != nullptr);

        if (load)
        {
          delete load;
        }

        auto                                      data = makeData();
        PhasorDynamics::LoadZ<ScalarT, IdxT>      monitored_load(bus, data);
        PhasorDynamics::Component<ScalarT, IdxT>& monitored_component  = monitored_load;
        success                                                       *= (monitored_component.getMonitor() != nullptr);

        delete bus;

        return success.report(__func__);
      }

      TestOutcome residual()
      {
        TestStatus success = true;

        RealT R{2.0}; ///< Load resistance
        RealT X{4.0}; ///< Load reactance

        ScalarT Vr{10.0}; ///< Bus real voltage
        ScalarT Vi{20.0}; ///< Bus imaginary voltage

        const ScalarT Ir{-5.0}; ///< Solution real current
        const ScalarT Ii{0.0};  ///< Solution imaginary current

        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus(Vr, Vi);
        PhasorDynamics::LoadZ<ScalarT, IdxT>       load(&bus, R, X);

        bus.allocate();
        load.allocate();

        bus.initialize();
        load.initialize();

        bus.evaluateResidual();
        load.evaluateResidual();

        success *= isEqual(bus.Ir(), Ir);
        success *= isEqual(bus.Ii(), Ii);

        return success.report(__func__);
      }

      TestOutcome jacobian()
      {
        TestStatus success = true;

        RealT R{2.0}; ///< Load resistance
        RealT X{4.0}; ///< Load reactance

        DependencyTracking::Variable Vr{10.0}; ///< Bus real voltage
        DependencyTracking::Variable Vi{20.0}; ///< Bus imaginary voltage

        PhasorDynamics::Bus<DependencyTracking::Variable, IdxT>   bus(Vr, Vi);
        PhasorDynamics::LoadZ<DependencyTracking::Variable, IdxT> load(&bus, R, X);

        bus.allocate();
        load.allocate();

        bus.initialize();
        load.initialize();

        auto* load_y = load.y().getData();
        for (size_t i = 0; i < load.size(); ++i)
        {
          load_y[i].setVariableNumber(i); ///< load independent variables
        }
        load.y().setDataUpdated();
        auto* bus_y = bus.y().getData();
        for (size_t i = 0; i < bus.size(); ++i)
        {
          bus_y[i].setVariableNumber(i + load.size()); // Bus independent variables
        }
        bus.y().setDataUpdated();

        bus.evaluateResidual();
        load.evaluateResidual(); ///< Computes the residual and the Jacobian values by tracking
                                 ///< the dependencies

        auto&                                                    residuals     = load.getResidual();
        const auto*                                              residual_data = residuals.getData();
        std::vector<DependencyTracking::Variable::DependencyMap> ref           = analyticalJacobian(R, X);

        /// Compare dependencies computed automatically to the ones computed analytically
        for (size_t i = 0; i < residuals.getSize(); ++i)
        {
          DependencyTracking::Variable                       res           = residual_data[i];
          const DependencyTracking::Variable::DependencyMap& dependencies  = res.getDependencies();
          success                                                         *= (GridKit::Testing::isEqual(dependencies, ref[i]));
        }

        return success.report(__func__);
      }

      TestOutcome monitor()
      {
        TestStatus success = true;

        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus(10.0, 20.0);
        auto                                       data = makeData();
        PhasorDynamics::LoadZ<ScalarT, IdxT>       load(&bus, data);

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
        if (values.size() == 3)
        {
          success *= isEqual(values[1], static_cast<RealT>(-50.0), tol_);
          success *= isEqual(values[2], static_cast<RealT>(-100.0), tol_);
        }
        else
        {
          success = false;
        }

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      TestOutcome enzymeJacobian()
      {
        TestStatus success = true;

        RealT R{2.0}; ///< Load resistance
        RealT X{4.0}; ///< Load reactance

        ScalarT Vr{10.0}; ///< Bus real voltage
        ScalarT Vi{20.0}; ///< Bus imaginary voltage

        PhasorDynamics::Bus<ScalarT, IdxT>   bus(Vr, Vi);
        PhasorDynamics::LoadZ<ScalarT, IdxT> load(&bus, R, X);

        bus.allocate();
        load.allocate();

        bus.initialize();
        load.initialize();

        for (size_t i = 0; i < bus.size(); ++i)
        {
          bus.setVariableIndex(i, i + load.size()); // Reset bus variable indices
          bus.setResidualIndex(i, i + load.size()); // Reset bus residual indices
        }

        bus.evaluateJacobian();
        load.evaluateJacobian();
        load.constructCsr();
        GridKit::LinearAlgebra::CsrMatrix<ScalarT, IdxT>* model_jacobian = load.getCsrJacobian();
        std::cout << "Sparse Csr Matrix: Load Jacobian\n";
        model_jacobian->print();

        /// Compare model Jacobian wih dependencies computed analytically
        std::vector<DependencyTracking::Variable::DependencyMap> ref                = analyticalJacobian(R, X);
        std::vector<DependencyTracking::Variable::DependencyMap> model_dependencies = GridKit::Testing::MapFromCsr(model_jacobian);
        for (size_t i = 0; i < ref.size(); ++i)
        {
          success *= (GridKit::Testing::isEqual(model_dependencies[i], ref[i]));
        }

        return success.report(__func__);
      }
#endif

    private:
      static constexpr RealT tol_ = 1.0e-10;

      auto makeData() -> DataT
      {
        using Params   = PhasorDynamics::LoadZParameters;
        using Variable = PhasorDynamics::LoadZMonitorableVariables;

        DataT data;
        data.device_class          = "LoadZ";
        data.disambiguation_string = "loadz_test";

        data.parameters[Params::R] = static_cast<RealT>(2.0);
        data.parameters[Params::X] = static_cast<RealT>(4.0);

        data.monitored_variables.insert(Variable::p);
        data.monitored_variables.insert(Variable::q);

        return data;
      }

      std::vector<DependencyTracking::Variable::DependencyMap> analyticalJacobian(const RealT R,
                                                                                  const RealT X)
      {
        const RealT b = -X / (R * R + X * X);
        const RealT g = R / (R * R + X * X);

        std::vector<DependencyTracking::Variable::DependencyMap> dependencies(2);
        dependencies[0] = {{0, 1.0}, {2, g}, {3, -b}};
        dependencies[1] = {{1, 1.0}, {2, b}, {3, g}};

        return dependencies;
      }
    };

  } // namespace Testing
} // namespace GridKit
