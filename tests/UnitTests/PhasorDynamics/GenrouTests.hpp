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
      using real_type               = typename PhasorDynamics::Component<ScalarT, IdxT>::real_type;
      static constexpr ScalarT tol_ = 10 * std::numeric_limits<ScalarT>::epsilon(); // added this: was not originally there

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

      // A test to verify that the hard coded answers match those given by the residual functions
      // Hard code parameters, differential, and algebraic terms
      TestOutcome hard_coded_residual()
      {
        TestStatus success = true;

        // GenRou generator parameters
        real_type p0{1};
        real_type q0{.05013};
        real_type H{3};
        real_type D{.5};
        real_type Ra{.1};
        real_type Tdop{7};
        real_type Tdopp{.04};
        real_type Tqopp{.05};
        real_type Tqop{.75};
        real_type Xd{2.1};
        real_type Xdp{.2};
        real_type Xdpp{.5};
        real_type Xq{.18};
        real_type Xqp{.3};
        real_type Xqpp{.5};
        real_type Xl{.5};
        real_type S10{.1};
        real_type S12{.2};

        ScalarT Vr1{1.0}; ///< Bus real voltage
        ScalarT Vi1{0};   ///< Bus imaginary voltage

        PhasorDynamics::Bus<ScalarT, IdxT>    bus(Vr1, Vi1);
        PhasorDynamics::Genrou<ScalarT, IdxT> gen(&bus, 1, p0, q0, H, D, Ra, Tdop, Tdopp, Tqopp, Tqop, Xd, Xdp, Xdpp, Xq, Xqp, Xqpp, Xl, S10, S12);

        // Answer key is available only in double precision.
        // Therefore, only double precision tests are done at this time.
        const std::vector<ScalarT> res_answer = {
            -2 * M_PI * 60.0,
            -static_cast<ScalarT>(10.) / static_cast<ScalarT>(9.),
            -static_cast<ScalarT>(223.) / static_cast<ScalarT>(525.),
            -54.75,
            -9.6,
            static_cast<ScalarT>(892.) / static_cast<ScalarT>(375.),
            0.21,
            -0.07,
            -0.19223748416156686,
            1.8896749891587163,
            1.4,
            0.31,
            2.211,
            0.85,
            1.2,
            static_cast<ScalarT>(64.) / static_cast<ScalarT>(65.),
            -static_cast<ScalarT>(237.) / static_cast<ScalarT>(130.),
            -static_cast<ScalarT>(141.) / static_cast<ScalarT>(130.),
            -static_cast<ScalarT>(241.) / static_cast<ScalarT>(260.)};

        bus.allocate();
        bus.initialize();

        // Allocate but not initialize generator model
        gen.allocate();

        // Set variable values matching the answer key
        gen.y()[0]  = M_PI; // delta
        gen.y()[1]  = 2.0;  // omega
        gen.y()[2]  = 2.0;  // Eqp
        gen.y()[3]  = .1;   // psidp
        gen.y()[4]  = .01;  // psiqp
        gen.y()[5]  = .6;   // Edp
        gen.y()[6]  = .2;   // psiqp
        gen.y()[7]  = .03;  // psidpp
        gen.y()[8]  = .01;  // psipp
        gen.y()[9]  = 2;    // ksat
        gen.y()[10] = .8;   // vd
        gen.y()[11] = .4;   // vq
        gen.y()[12] = 2;    // telec
        gen.y()[13] = 1.1;  // id
        gen.y()[14] = .3;   // iq
        gen.y()[15] = .9;   // ir
        gen.y()[16] = .25;  // ii
        gen.y()[17] = .3;   // inr
        gen.y()[18] = .15;  // ini

        // Set derivative values matching the answer key
        gen.yp()[0] = 2 * M_PI * 60.0; // delta_dot
        gen.yp()[1] = -1.5;            // omega_dot
        gen.yp()[2] = 1;               // Eqp_dot
        gen.yp()[3] = 1;               // psidp_dot
        gen.yp()[4] = 1;               // psiqp_dot
        gen.yp()[5] = 1;               // Edp_dot

        gen.evaluateResidual();
        std::vector<ScalarT>& residual = gen.getResidual();

        for (size_t i = 0; i < res_answer.size(); ++i)
        {
          if (!isEqual(residual[i], res_answer[i], tol_))
          {
            std::cout << "Incorrect result for residual " << i << ": "
                      << residual[i] << " != " << res_answer[i] << "\n";
            success = false;
            break;
          }
        }

        return success.report(__func__);
      }

      TestOutcome accessors()
      {
        TestStatus success = true;
        success.skipTest();

        return success.report(__func__);
      }


#if 0 // Disabled GenrouEnzyme because of an intermittent build issue
#ifdef GRIDKIT_ENABLE_ENZYME
      /**
       * @brief Checks Jacobian evaluation.
       */
      TestOutcome jacobian()
      {
        TestStatus success = true;

        auto tol = 10 * std::numeric_limits<real_type>::epsilon();

        // Jacobian via DependencyTracking
        std::vector<DependencyTracking::Variable> dependency_tracking_residuals = DependencyTrackingJacobian();

        // Jacobian via Enzyme
        std::vector<DependencyTracking::Variable::DependencyMap> enzyme_jacobian = EnzymeJacobian();

        /// Compare DependencyTracking dependencies to Enzyme's
        for (size_t i = 0; i < dependency_tracking_residuals.size(); ++i)
        {
          DependencyTracking::Variable                       res           = dependency_tracking_residuals[i];
          const DependencyTracking::Variable::DependencyMap& dependencies  = res.getDependencies();
          success                                                         *= (GridKit::Testing::isEqual(dependencies, enzyme_jacobian[i], tol));
        }

        return success.report(__func__);
      }

    private:
      std::vector<DependencyTracking::Variable> DependencyTrackingJacobian()
      {
        DependencyTracking::Variable                               Vr1{1.0}; ///< Bus real voltage
        DependencyTracking::Variable                               Vi1{0.0}; ///< Bus imaginary voltage
        PhasorDynamics::Bus<DependencyTracking::Variable, IdxT>    bus(Vr1, Vi1);
        PhasorDynamics::Genrou<DependencyTracking::Variable, IdxT> gen(&bus,
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

        gen.allocate();
        gen.initialize();

        for (size_t i = 0; i < gen.size(); ++i)
        {
          gen.y()[i].setVariableNumber(i); ///< Independent variables
        }

        bus.evaluateResidual();
        gen.evaluateResidual(); ///< Computes the residual and the Jacobian values by tracking
                                ///< the dependencies
        std::vector<DependencyTracking::Variable> residual = gen.getResidual();

        /// Print the dependencies
        for (size_t i = 0; i < residual.size(); ++i)
        {
          std::cout << i << "th residual: ";
          (residual[i]).print(std::cout);
          std::cout << "\n";
        }

        return residual;
      }

      std::vector<DependencyTracking::Variable::DependencyMap> EnzymeJacobian()
      {
        ScalarT                               Vr1{1.0}; ///< Bus real voltage
        ScalarT                               Vi1{0.0}; ///< Bus imaginary voltage
        PhasorDynamics::Bus<ScalarT, IdxT>    bus(Vr1, Vi1);
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

        gen.evaluateJacobian();
        GridKit::LinearAlgebra::COO_Matrix<ScalarT, IdxT> model_jacobian = gen.getJacobian();
        model_jacobian.printMatrix("Model Jacobian");

        return GridKit::Testing::MapFromCOO(model_jacobian);
      }
#endif
#endif
    }; // class GenrouTest

  } // namespace Testing
} // namespace GridKit
