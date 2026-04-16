#include <iomanip>
#include <iostream>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/PhasorDynamics/Branch/Branch.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {

    template <class ScalarT, typename IdxT>
    class BranchTests
    {
    private:
      using RealT = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

    public:
      BranchTests()  = default;
      ~BranchTests() = default;

      TestOutcome constructor()
      {
        TestStatus success = true;

        auto* bus1 = new PhasorDynamics::Bus<ScalarT, IdxT>(1.0, 0.0);
        auto* bus2 = new PhasorDynamics::Bus<ScalarT, IdxT>(1.0, 0.1);

        PhasorDynamics::Component<ScalarT, IdxT>* branch =
            new PhasorDynamics::Branch<ScalarT, IdxT>(bus1, bus2);

        success *= (branch != nullptr);

        if (branch)
        {
          delete branch;
        }
        delete bus1;
        delete bus2;

        return success.report(__func__);
      }

      TestOutcome residual()
      {
        TestStatus success = true;

        RealT R{2.0}; ///< Branch series resistance
        RealT X{4.0}; ///< Branch series reactance
        RealT G{0.2}; ///< Branch shunt conductance
        RealT B{1.2}; ///< Branch shunt charging

        ScalarT Vr1{10.0}; ///< Bus-1 real voltage
        ScalarT Vi1{20.0}; ///< Bus-1 imaginary voltage
        ScalarT Vr2{30.0}; ///< Bus-2 real voltage
        ScalarT Vi2{40.0}; ///< Bus-2 imaginary voltage

        const ScalarT Ir1{17.0};  ///< Solution: real current entering bus-1
        const ScalarT Ii1{-10.0}; ///< Solution: imaginary current entering bus-1
        const ScalarT Ir2{15.0};  ///< Solution: real current entering bus-2
        const ScalarT Ii2{-20.0}; ///< Solution: imaginary current entering bus-2

        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus1(Vr1, Vi1);
        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus2(Vr2, Vi2);

        PhasorDynamics::Branch<ScalarT, IdxT> branch(&bus1, &bus2, R, X, G, B);
        branch.allocate();
        branch.evaluateResidual();

        success *= isEqual(bus1.i(0), Ir1);
        success *= isEqual(bus1.i(1), Ii1);
        success *= isEqual(bus2.i(0), Ir2);
        success *= isEqual(bus2.i(1), Ii2);

        return success.report(__func__);
      }

      TestOutcome jacobian()
      {
        TestStatus success = true;

        RealT R{2.0}; ///< Branch series resistance
        RealT X{4.0}; ///< Branch series reactance
        RealT G{0.2}; ///< Branch shunt conductance
        RealT B{1.2}; ///< Branch shunt charging

        DependencyTracking::Variable Vr1{10.0}; ///< Bus-1 real voltage
        DependencyTracking::Variable Vi1{20.0}; ///< Bus-1 imaginary voltage
        DependencyTracking::Variable Vr2{30.0}; ///< Bus-2 real voltage
        DependencyTracking::Variable Vi2{40.0}; ///< Bus-2 imaginary voltage

        Vr1.setVariableNumber(0); ///< Independent variables: first
        Vi1.setVariableNumber(1); ///< Independent variables: second
        Vr2.setVariableNumber(2); ///< Independent variables: third
        Vi2.setVariableNumber(3); ///< Independent variables: fourth

        PhasorDynamics::BusInfinite<DependencyTracking::Variable, IdxT> bus1(Vr1, Vi1);
        PhasorDynamics::BusInfinite<DependencyTracking::Variable, IdxT> bus2(Vr2, Vi2);

        PhasorDynamics::Branch<DependencyTracking::Variable, IdxT> branch(&bus1, &bus2, R, X, G, B);
        branch.allocate();
        branch.evaluateResidual(); ///< Computes the residual and the Jacobian values by tracking
                                   ///< the dependencies

        std::vector<DependencyTracking::Variable>                residuals{bus1.i(0), bus1.i(1), bus2.i(0), bus2.i(1)};
        std::vector<DependencyTracking::Variable::DependencyMap> ref = analyticalJacobian(R, X, G, B);

        /// Compare dependencies computed automatically to the ones computed analytically
        for (size_t i = 0; i < residuals.size(); ++i)
        {
          DependencyTracking::Variable                       res           = residuals[i];
          const DependencyTracking::Variable::DependencyMap& dependencies  = res.getDependencies();
          success                                                         *= (GridKit::Testing::isEqual(dependencies, ref[i]));
        }

        return success.report(__func__);
      }

      TestOutcome accessors()
      {
        TestStatus success = true;

        const RealT zero{0.0};

        RealT R{2.0}; ///< Branch series resistance
        RealT X{4.0}; ///< Branch series reactance
        RealT G{0.2}; ///< Branch shunt conductance
        RealT B{1.2}; ///< Branch shunt charging

        ScalarT Vr1{-1.0}; ///< Bus-1 real voltage
        ScalarT Vi1{-1.0}; ///< Bus-1 imaginary voltage
        ScalarT Vr2{1.0};  ///< Bus-2 real voltage
        ScalarT Vi2{1.0};  ///< Bus-2 imaginary voltage

        const ScalarT res_R{1.0};  ///< Solution: real current entering bus-1
        const ScalarT res_X{0.5};  ///< Solution: imaginary current entering bus-1
        const ScalarT res_G{3.0};  ///< Solution: real current entering bus-2
        const ScalarT res_B{-4.0}; ///< Solution: imaginary current entering bus-2

        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus1(Vr1, Vi1);
        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus2(Vr2, Vi2);

        PhasorDynamics::Branch<ScalarT, IdxT> branch(&bus1,
                                                     &bus2,
                                                     zero,
                                                     zero,
                                                     zero,
                                                     zero);
        branch.allocate();

        // Test setting branch series resistance
        branch.setR(R);
        bus1.evaluateResidual(); // <- set Ir1 to zero
        branch.evaluateResidual();
        success *= isEqual(bus1.i(0), res_R);
        branch.setR(zero);

        // Test setting branch series reactance
        branch.setX(X);
        bus1.evaluateResidual(); // <- set Ir1 to zero
        branch.evaluateResidual();
        success *= isEqual(bus1.i(0), res_X);
        branch.setX(zero);
        return success.report(__func__);

        // Test setting branch shunt conductance
        branch.setG(G);
        bus1.evaluateResidual(); // <- set Ir1 to zero
        branch.evaluateResidual();
        success *= isEqual(bus1.i(0), res_G);
        branch.setG(zero);

        // Test setting branch shunt charging
        branch.setB(B);
        bus1.evaluateResidual(); // <- set Ir1 to zero
        branch.evaluateResidual();
        success *= isEqual(bus1.i(0), res_B);
        branch.setB(zero);

        return success.report(__func__);
      }

    private:
      std::vector<DependencyTracking::Variable::DependencyMap> analyticalJacobian(const RealT R,
                                                                                  const RealT X,
                                                                                  const RealT G,
                                                                                  const RealT B)
      {
        const RealT b = -X / (R * R + X * X);
        const RealT g = R / (R * R + X * X);

        RealT dIr1_dVr1 = -(g + 0.5 * G);
        RealT dIr1_dVi1 = (b + 0.5 * B);
        RealT dIr1_dVr2 = g;
        RealT dIr1_dVi2 = -b;

        RealT dIi1_dVr1 = -(b + 0.5 * B);
        RealT dIi1_dVi1 = -(g + 0.5 * G);
        RealT dIi1_dVr2 = b;
        RealT dIi1_dVi2 = g;

        RealT dIr2_dVr1 = g;
        RealT dIr2_dVi1 = -b;
        RealT dIr2_dVr2 = -(g + 0.5 * G);
        RealT dIr2_dVi2 = (b + 0.5 * B);

        RealT dIi2_dVr1 = b;
        RealT dIi2_dVi1 = g;
        RealT dIi2_dVr2 = -(b + 0.5 * B);
        RealT dIi2_dVi2 = -(g + 0.5 * G);

        std::vector<DependencyTracking::Variable::DependencyMap> dependencies(4);
        dependencies[0] = {{0, dIr1_dVr1}, {1, dIr1_dVi1}, {2, dIr1_dVr2}, {3, dIr1_dVi2}};
        dependencies[1] = {{0, dIi1_dVr1}, {1, dIi1_dVi1}, {2, dIi1_dVr2}, {3, dIi1_dVi2}};
        dependencies[2] = {{0, dIr2_dVr1}, {1, dIr2_dVi1}, {2, dIr2_dVr2}, {3, dIr2_dVi2}};
        dependencies[3] = {{0, dIi2_dVr1}, {1, dIi2_dVi1}, {2, dIi2_dVr2}, {3, dIi2_dVi2}};

        return dependencies;
      }
    }; // class BranchTest

  } // namespace Testing
} // namespace GridKit
