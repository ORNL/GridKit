#include <iomanip>
#include <iostream>

#include <AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <Definitions.hpp>
#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <Model/PhasorDynamics/SynchronousMachine/GENROUwS/Genrou.hpp>
#include <Utilities/MapFromCOO.hpp>
#include <Utilities/TestHelpers.hpp>
#include <Utilities/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {

    template <class ScalarT, typename IdxT>
    class GenrouTests
    {
    private:
      using real_type = typename PhasorDynamics::Component<ScalarT, IdxT>::real_type;

    public:
      GenrouTests()  = default;
      ~GenrouTests() = default;

      TestOutcome constructor()
      {
        TestStatus success = true;

        auto* bus = new PhasorDynamics::Bus<ScalarT, IdxT>(1.0, 0.0);

        PhasorDynamics::Component<ScalarT, IdxT>* machine =
            new PhasorDynamics::Genrou<ScalarT, IdxT>(bus, 1);

        success *= (machine != nullptr);

        if (machine)
        {
          delete machine;
        }
        delete bus;

        return success.report(__func__);
      }

      /**
       * @brief Checks residual evaluation.
       *
       * The test instantiates and initializes Genrou model. Properly
       * initialized model should have residual equal to zero within machine
       * precision.
       *
       * @return TestOutcome - wheter test was successful
       */
      TestOutcome residual()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT>    bus(1.0, 0.0);
        PhasorDynamics::Genrou<ScalarT, IdxT> gen(&bus,
                                                  1,
                                                  1,
                                                  0.05013,
                                                  3,
                                                  0,
                                                  0,
                                                  7,
                                                  0.04,
                                                  0.05,
                                                  0.75,
                                                  2.1,
                                                  0.2,
                                                  0.18,
                                                  0.5,
                                                  0.5,
                                                  0.18,
                                                  0.15,
                                                  0,
                                                  0);

        bus.allocate();
        bus.initialize();
        bus.evaluateResidual();

        gen.allocate();
        gen.initialize();
        gen.evaluateResidual();

        // Require results to be within machine precision
        auto tol = 10 * std::numeric_limits<real_type>::epsilon();

        const std::vector<ScalarT>& f = gen.getResidual();
        for (const auto& f_val : f)
        {
          if (!isEqual(f_val, 0.0, tol))
            success = false;
        }

        return success.report(__func__);
      }

      TestOutcome accessors()
      {
        TestStatus success = true;
        success.skipTest();

        return success.report(__func__);
      }

      // #ifdef GRIDKIT_ENABLE_ENZYME
      //       /**
      //        * @brief Checks Jacobian evaluation.
      //        */
      //       TestOutcome jacobian()
      //       {
      //         TestStatus success = true;
      //
      //         auto tol = 10 * std::numeric_limits<real_type>::epsilon();
      //
      //         // Jacobian via DependencyTracking
      //         std::vector<DependencyTracking::Variable> dependency_tracking_residuals = DependencyTrackingJacobian();
      //
      //         // Jacobian via Enzyme
      //         std::vector<DependencyTracking::Variable::DependencyMap> enzyme_jacobian = EnzymeJacobian();
      //
      //         /// Compare DependencyTracking dependencies to Enzyme's
      //         for (size_t i = 0; i < dependency_tracking_residuals.size(); ++i)
      //         {
      //           DependencyTracking::Variable                       res           = dependency_tracking_residuals[i];
      //           const DependencyTracking::Variable::DependencyMap& dependencies  = res.getDependencies();
      //           success                                                         *= (GridKit::Testing::isEqual(dependencies, enzyme_jacobian[i], tol));
      //         }
      //
      //         return success.report(__func__);
      //       }
      //
      //     private:
      //       std::vector<DependencyTracking::Variable> DependencyTrackingJacobian()
      //       {
      //         DependencyTracking::Variable                               Vr1{1.0}; ///< Bus real voltage
      //         DependencyTracking::Variable                               Vi1{0.0}; ///< Bus imaginary voltage
      //         PhasorDynamics::Bus<DependencyTracking::Variable, IdxT>    bus(Vr1, Vi1);
      //         PhasorDynamics::Genrou<DependencyTracking::Variable, IdxT> gen(&bus,
      //                                                                        1,
      //                                                                        1,
      //                                                                        0.05013,
      //                                                                        3,
      //                                                                        0,
      //                                                                        0,
      //                                                                        7,
      //                                                                        0.04,
      //                                                                        0.05,
      //                                                                        0.75,
      //                                                                        2.1,
      //                                                                        0.2,
      //                                                                        0.18,
      //                                                                        0.5,
      //                                                                        0.5,
      //                                                                        0.18,
      //                                                                        0.15,
      //                                                                        0,
      //                                                                        0);
      //
      //         bus.allocate();
      //         bus.initialize();
      //
      //         gen.allocate();
      //         gen.initialize();
      //
      //         for (size_t i = 0; i < gen.size(); ++i)
      //         {
      //           gen.y()[i].setVariableNumber(i); ///< Independent variables
      //         }
      //
      //         bus.evaluateResidual();
      //         gen.evaluateResidual(); ///< Computes the residual and the Jacobian values by tracking
      //                                 ///< the dependencies
      //         std::vector<DependencyTracking::Variable> residual = gen.getResidual();
      //
      //         /// Print the dependencies
      //         for (size_t i = 0; i < residual.size(); ++i)
      //         {
      //           std::cout << i << "th residual: ";
      //           (residual[i]).print(std::cout);
      //           std::cout << "\n";
      //         }
      //
      //         return residual;
      //       }
      //
      //       std::vector<DependencyTracking::Variable::DependencyMap> EnzymeJacobian()
      //       {
      //         ScalarT                               Vr1{1.0}; ///< Bus real voltage
      //         ScalarT                               Vi1{0.0}; ///< Bus imaginary voltage
      //         PhasorDynamics::Bus<ScalarT, IdxT>    bus(Vr1, Vi1);
      //         PhasorDynamics::Genrou<ScalarT, IdxT> gen(&bus,
      //                                                   1,
      //                                                   1,
      //                                                   0.05013,
      //                                                   3,
      //                                                   0,
      //                                                   0,
      //                                                   7,
      //                                                   0.04,
      //                                                   0.05,
      //                                                   0.75,
      //                                                   2.1,
      //                                                   0.2,
      //                                                   0.18,
      //                                                   0.5,
      //                                                   0.5,
      //                                                   0.18,
      //                                                   0.15,
      //                                                   0,
      //                                                   0);
      //
      //         bus.allocate();
      //         bus.initialize();
      //         bus.evaluateResidual();
      //
      //         gen.allocate();
      //         gen.initialize();
      //         gen.evaluateResidual();
      //
      //         gen.evaluateJacobian();
      //         GridKit::LinearAlgebra::COO_Matrix<ScalarT, IdxT> model_jacobian = gen.getJacobian();
      //         model_jacobian.printMatrix("Model Jacobian");
      //
      //         return GridKit::Testing::MapFromCOO(model_jacobian);
      //       }
      // #endif
    }; // class GenrouTest

  } // namespace Testing
} // namespace GridKit
