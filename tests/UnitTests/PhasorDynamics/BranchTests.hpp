#include <iomanip>
#include <iostream>

#include <LinearAlgebra/SparsityPattern/Variable.hpp>
#include <Model/PhasorDynamics/Branch/Branch.hpp>
#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <Utilities/TestHelpers.hpp>
#include <Utilities/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {

    template <class ScalarT, typename IdxT>
    class BranchTests
    {
    private:
      using real_type = typename PhasorDynamics::Component<ScalarT, IdxT>::real_type;

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

        real_type R{2.0}; ///< Branch series resistance
        real_type X{4.0}; ///< Branch series reactance
        real_type G{0.2}; ///< Branch shunt conductance
        real_type B{1.2}; ///< Branch shunt charging

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
        branch.evaluateResidual();

        success *= isEqual(bus1.Ir(), Ir1);
        success *= isEqual(bus1.Ii(), Ii1);
        success *= isEqual(bus2.Ir(), Ir2);
        success *= isEqual(bus2.Ii(), Ii2);

        return success.report(__func__);
      }

      TestOutcome jacobian()
      {
        TestStatus success = true;

        real_type R{2.0}; ///< Branch series resistance
        real_type X{4.0}; ///< Branch series reactance
        real_type G{0.2}; ///< Branch shunt conductance
        real_type B{1.2}; ///< Branch shunt charging

        Sparse::Variable Vr1{10.0}; ///< Bus-1 real voltage
        Sparse::Variable Vi1{20.0}; ///< Bus-1 imaginary voltage
        Sparse::Variable Vr2{30.0}; ///< Bus-2 real voltage
        Sparse::Variable Vi2{40.0}; ///< Bus-2 imaginary voltage

        const Sparse::Variable Ir1{17.0};  ///< Solution: real current entering bus-1
        const Sparse::Variable Ii1{-10.0}; ///< Solution: imaginary current entering bus-1
        const Sparse::Variable Ir2{15.0};  ///< Solution: real current entering bus-2
        const Sparse::Variable Ii2{-20.0}; ///< Solution: imaginary current entering bus-2

        PhasorDynamics::BusInfinite<Sparse::Variable, IdxT> bus1(Vr1, Vi1);
        PhasorDynamics::BusInfinite<Sparse::Variable, IdxT> bus2(Vr2, Vi2);

        PhasorDynamics::Branch<Sparse::Variable, IdxT> branch(&bus1, &bus2, R, X, G, B);
        branch.evaluateResidual();

        {
          Sparse::Variable res = bus1.Ir();
          const Sparse::Variable::DependencyMap& dependencies =
              res.getDependencies();

          std::cout << "Dependency size: " << dependencies.size() << "\n";
        }

        return success.report(__func__);
      }

      TestOutcome accessors()
      {
        TestStatus success = true;

        const real_type zero{0.0};

        real_type R{2.0}; ///< Branch series resistance
        real_type X{4.0}; ///< Branch series reactance
        real_type G{0.2}; ///< Branch shunt conductance
        real_type B{1.2}; ///< Branch shunt charging

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

        // Test setting branch series resistance
        branch.setR(R);
        bus1.evaluateResidual(); // <- set Ir1 to zero
        branch.evaluateResidual();
        success *= isEqual(bus1.Ir(), res_R);
        branch.setR(zero);

        // Test setting branch series reactance
        branch.setX(X);
        bus1.evaluateResidual(); // <- set Ir1 to zero
        branch.evaluateResidual();
        success *= isEqual(bus1.Ir(), res_X);
        branch.setX(zero);
        return success.report(__func__);

        // Test setting branch shunt conductance
        branch.setG(G);
        bus1.evaluateResidual(); // <- set Ir1 to zero
        branch.evaluateResidual();
        success *= isEqual(bus1.Ir(), res_G);
        branch.setG(zero);

        // Test setting branch shunt charging
        branch.setB(B);
        bus1.evaluateResidual(); // <- set Ir1 to zero
        branch.evaluateResidual();
        success *= isEqual(bus1.Ir(), res_B);
        branch.setB(zero);

        return success.report(__func__);
      }
    }; // class BranchTest

  } // namespace Testing
} // namespace GridKit
