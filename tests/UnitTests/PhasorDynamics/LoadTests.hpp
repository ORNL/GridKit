#pragma once

#include <iomanip>
#include <iostream>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <GridKit/Model/PhasorDynamics/Load/Load.hpp>
#include <GridKit/Utilities/MapFromCOO.hpp>
#include <GridKit/Utilities/TestHelpers.hpp>
#include <GridKit/Utilities/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class LoadTests
    {
    public:
      using real_type = typename PhasorDynamics::Component<ScalarT, IdxT>::real_type;

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

        real_type R{2.0}; ///< Load resistance
        real_type X{4.0}; ///< Load reactance

        ScalarT Vr{10.0}; ///< Bus real voltage
        ScalarT Vi{20.0}; ///< Bus imaginary voltage

        const ScalarT Ir{-5.0}; ///< Solution real current
        const ScalarT Ii{0.0};  ///< Solution imaginary current

        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus(Vr, Vi);
        PhasorDynamics::Load<ScalarT, IdxT>        load(&bus, R, X);
        bus.allocate();
        load.allocate();
        load.evaluateResidual();

        success *= isEqual(bus.Ir(), Ir);
        success *= isEqual(bus.Ii(), Ii);

        return success.report(__func__);
      }

      TestOutcome jacobian()
      {
        TestStatus success = true;

        real_type R{2.0}; ///< Load resistance
        real_type X{4.0}; ///< Load reactance

        DependencyTracking::Variable Vr{10.0}; ///< Bus real voltage
        DependencyTracking::Variable Vi{20.0}; ///< Bus imaginary voltage

        Vr.setVariableNumber(0); ///< Independent variables: first
        Vi.setVariableNumber(1); ///< Independent variables: second

        PhasorDynamics::BusInfinite<DependencyTracking::Variable, IdxT> bus(Vr, Vi);
        PhasorDynamics::Load<DependencyTracking::Variable, IdxT>        load(&bus, R, X);
        bus.allocate();
        load.allocate();
        load.evaluateResidual(); ///< Computes the residual and the Jacobian values by tracking
                                 ///< the dependencies

        std::vector<DependencyTracking::Variable>                residuals{bus.Ir(), bus.Ii()};
        std::vector<DependencyTracking::Variable::DependencyMap> ref = analyticalJacobian(R, X);

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

        real_type R{2.0}; ///< Load resistance
        real_type X{4.0}; ///< Load reactance

        ScalarT Vr{10.0}; ///< Bus real voltage
        ScalarT Vi{20.0}; ///< Bus imaginary voltage

        PhasorDynamics::Bus<ScalarT, IdxT>  bus(Vr, Vi);
        PhasorDynamics::Load<ScalarT, IdxT> load(&bus, R, X);
        bus.allocate();
        load.allocate();
        load.evaluateJacobian();
        GridKit::LinearAlgebra::COO_Matrix<ScalarT, IdxT> model_jacobian = bus.getJacobian();
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
      std::vector<DependencyTracking::Variable::DependencyMap> analyticalJacobian(const real_type R,
                                                                                  const real_type X)
      {
        const real_type b = -X / (R * R + X * X);
        const real_type g = R / (R * R + X * X);

        real_type dIr_dVr = -g;
        real_type dIr_dVi = b;

        real_type dIi_dVr = -b;
        real_type dIi_dVi = -g;

        std::vector<DependencyTracking::Variable::DependencyMap> dependencies(2);
        dependencies[0] = {{0, dIr_dVr}, {1, dIr_dVi}};
        dependencies[1] = {{0, dIi_dVr}, {1, dIi_dVi}};

        return dependencies;
      }
    };

  } // namespace Testing
} // namespace GridKit
