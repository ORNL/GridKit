#pragma once

#include <iomanip>
#include <iostream>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <GridKit/Model/PhasorDynamics/Load/Load.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCOO.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class LoadTests
    {
    public:
      using RealT = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      LoadTests()  = default;
      ~LoadTests() = default;

      TestOutcome constructor()
      {
        TestStatus success = true;

        auto* bus = new PhasorDynamics::Bus<ScalarT, IdxT>(1.0, 0.0);

        PhasorDynamics::Component<ScalarT, IdxT>* load =
            new PhasorDynamics::Load<ScalarT, IdxT>(bus);

        success *= (load != nullptr);

        if (load)
        {
          delete load;
        }
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
        PhasorDynamics::Load<ScalarT, IdxT>        load(&bus, R, X);

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

        PhasorDynamics::Bus<DependencyTracking::Variable, IdxT>  bus(Vr, Vi);
        PhasorDynamics::Load<DependencyTracking::Variable, IdxT> load(&bus, R, X);

        bus.allocate();
        load.allocate();

        bus.initialize();
        load.initialize();

        for (size_t i = 0; i < load.size(); ++i)
        {
          load.y()[i].setVariableNumber(i); ///< load independent variables
        }
        for (size_t i = 0; i < bus.size(); ++i)
        {
          bus.y()[i].setVariableNumber(i + load.size()); // Bus independent variables
        }

        bus.evaluateResidual();
        load.evaluateResidual(); ///< Computes the residual and the Jacobian values by tracking
                                 ///< the dependencies

        std::vector<DependencyTracking::Variable>                residuals = load.getResidual();
        std::vector<DependencyTracking::Variable::DependencyMap> ref       = analyticalJacobian(R, X);

        /// Compare dependencies computed automatically to the ones computed analytically
        for (size_t i = 0; i < residuals.size(); ++i)
        {
          DependencyTracking::Variable                       res           = residuals[i];
          const DependencyTracking::Variable::DependencyMap& dependencies  = res.getDependencies();
          success                                                         *= (GridKit::Testing::isEqual(dependencies, ref[i]));
        }

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      TestOutcome enzyme_jacobian()
      {
        TestStatus success = true;

        RealT R{2.0}; ///< Load resistance
        RealT X{4.0}; ///< Load reactance

        ScalarT Vr{10.0}; ///< Bus real voltage
        ScalarT Vi{20.0}; ///< Bus imaginary voltage

        PhasorDynamics::Bus<ScalarT, IdxT>  bus(Vr, Vi);
        PhasorDynamics::Load<ScalarT, IdxT> load(&bus, R, X);

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

        GridKit::LinearAlgebra::COO_Matrix<ScalarT, IdxT> model_jacobian = load.getJacobian();
        model_jacobian.printMatrix("Model Jacobian");

        /// Compare model Jacobian wih dependencies computed analytically
        std::vector<DependencyTracking::Variable::DependencyMap> ref                = analyticalJacobian(R, X);
        std::vector<DependencyTracking::Variable::DependencyMap> model_dependencies = GridKit::Testing::MapFromCOO(model_jacobian);
        for (size_t i = 0; i < ref.size(); ++i)
        {
          success *= (GridKit::Testing::isEqual(model_dependencies[i], ref[i]));
        }

        return success.report(__func__);
      }
#endif

    private:
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
