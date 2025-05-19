#include <iomanip>
#include <iostream>
#include <limits>

#include <Model/PhasorDynamics/SynchronousMachine/GenClassical/GenClassical.cpp>

#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <Model/PhasorDynamics/SynchronousMachine/GenClassical/GenClassical.hpp>
#include <Utilities/TestHelpers.hpp>
#include <Utilities/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {

    template <class ScalarT, typename IdxT>
    class GenClassicalTests
    {
    private:
      using real_type = typename PhasorDynamics::Component<ScalarT, IdxT>::real_type;

    public:
      GenClassicalTests()  = default;
      ~GenClassicalTests() = default;

      TestOutcome constructor()
      {
        TestStatus success = true;

        auto* bus = new PhasorDynamics::Bus<ScalarT, IdxT>(1.0, 0.0);

        PhasorDynamics::Component<ScalarT, IdxT>* machine =
            new PhasorDynamics::GenClassical<ScalarT, IdxT>(bus, 1);

        success *= (machine != nullptr);

        if (machine)
        {
          delete machine;
        }
        delete bus;

        return success.report(__func__);
      }

      /**
       * A test case to verify residual values
       */
      TestOutcome residual()
      {
        TestStatus success = true;

        // classical generator parameters
        real_type H{0.5};
        real_type D{-1.0};
        real_type Ra{0.5};
        real_type Xdp{0.5};
        real_type pmech{1.0};
        real_type ep{2.0};

        ScalarT Vr1{1.0}; ///< Bus-1 real voltage
        ScalarT Vi1{1.0}; ///< Bus-1 imaginary voltage

        const ScalarT res0{0.0};                                          /// first residual
        const ScalarT res1{0.0};                                          /// second residual
        const ScalarT res2{0.0};                                          /// third residual
        const ScalarT res3{0.0};                                          /// fourth residual
        const ScalarT res4{0.0};                                          /// fifth residual
        const ScalarT res5{0.0};                                          /// fifth residual
        const ScalarT res6{0.0};                                          /// fifth residual
        const ScalarT tol = 7 * (std::numeric_limits<double>::epsilon()); // tolerance for comparing results

        PhasorDynamics::Bus<ScalarT, IdxT>          bus(Vr1, Vi1);
        PhasorDynamics::GenClassical<ScalarT, IdxT> gen(&bus, 1, 1.0, 1.0, H, D, Ra, Xdp);
        bus.allocate();
        bus.initialize();
        gen.allocate();

        gen.y()[0] = M_PI; // delta
        gen.y()[1] = 1.0;  // omega
        gen.y()[2] = 2.0;  // telec
        gen.y()[1] = 1.0;  // omega
        gen.y()[2] = 2.0;  // telec
        gen.y()[3] = -2.0; // ir
        gen.y()[4] = -4.0; // ii
        gen.y()[5] = 1.0;  // pmech
        gen.y()[6] = 2.0;  // Ep
        gen.y()[5] = 1.0;  // pmech
        gen.y()[6] = 2.0;  // Ep

        gen.yp()[0] = 2 * M_PI * 60.0; // delta_dot
        gen.yp()[1] = -1.0;            // omega_dot
        gen.yp()[2] = 0.0;             // telec
        gen.yp()[3] = 0.0;             // ir
        gen.yp()[4] = 0.0;             // ii
        gen.yp()[5] = 0.0;             // pmech
        gen.yp()[6] = 0.0;             // Ep
        gen.yp()[0] = 2 * M_PI * 60.0; // delta_dot
        gen.yp()[1] = -1.0;            // omega_dot
        gen.yp()[2] = 0.0;             // telec
        gen.yp()[3] = 0.0;             // ir
        gen.yp()[4] = 0.0;             // ii
        gen.yp()[5] = 0.0;             // pmech
        gen.yp()[6] = 0.0;             // Ep

        gen.evaluateResidual();

        std::vector<ScalarT> residual = gen.getResidual();

        success *= isEqual(residual[0], res0, tol);
        success *= isEqual(residual[1], res1, tol);
        success *= isEqual(residual[2], res2, tol);
        success *= isEqual(residual[3], res3, tol);
        success *= isEqual(residual[4], res4, tol);
        // success *= isEqual(residual[5], res5, tol);
        // success *= isEqual(residual[6], res6, tol);

        return success.report(__func__);
      }

      /**
       *
       * Verifies correctness of the system initialization
       */
      TestOutcome initial()
      {
        TestStatus success = true;

        // classical generator parameters
        real_type p0{3.0};
        real_type q0{-1.0};
        real_type H{1.0};
        real_type D{1.0};
        real_type Ra{0.6};
        real_type Xdp{0.2};

        ScalarT Vr1{1.0}; ///< Bus-1 real voltage
        ScalarT Vi1{1.0}; ///< Bus-1 imaginary voltage

        const ScalarT delta{M_PI / 4.0};   /// first residual
        const ScalarT omega{0.0};          /// second residual
        const ScalarT Te{6.0};             /// third residual
        const ScalarT ir{1.0};             /// fourth residual
        const ScalarT ii{2.0};             /// fifth residual
        const ScalarT pmech{6.0};          /// fifth residual
        const ScalarT Ep{2.0 * sqrt(2.0)}; /// fifth residual

        const ScalarT tol = 5.0 * (std::numeric_limits<double>::epsilon()); // tolerance for comparing result

        PhasorDynamics::Bus<ScalarT, IdxT>          bus(Vr1, Vi1);
        PhasorDynamics::GenClassical<ScalarT, IdxT> gen(&bus, 1, p0, q0, H, D, Ra, Xdp);
        bus.allocate();
        bus.initialize();
        gen.allocate();
        gen.initialize();

        success *= isEqual(gen.y()[0], delta, tol);
        success *= isEqual(gen.y()[1], omega, tol);
        success *= isEqual(gen.y()[2], Te, tol);
        success *= isEqual(gen.y()[3], ir, tol);
        success *= isEqual(gen.y()[4], ii, tol);
        success *= isEqual(gen.y()[5], pmech, tol);
        success *= isEqual(gen.y()[6], Ep, tol);

        success *= isEqual(gen.yp()[0], 0.0, tol);
        success *= isEqual(gen.yp()[1], 0.0, tol);
        success *= isEqual(gen.yp()[2], 0.0, tol);
        success *= isEqual(gen.yp()[3], 0.0, tol);
        success *= isEqual(gen.yp()[4], 0.0, tol);
        success *= isEqual(gen.yp()[5], 0.0, tol);
        success *= isEqual(gen.yp()[6], 0.0, tol);

        return success.report(__func__);
      }

      /*
       *Verifies the residual evaluates to zero for the initial conditions
       */
      TestOutcome zeroInitialResidual()
      {
        TestStatus success = true;

        // classical generator parameters
        real_type p0{3.0};
        real_type q0{-1.0};
        real_type H{1.0};
        real_type D{1.0};
        real_type Ra{0.6};
        real_type Xdp{0.2};

        ScalarT Vr1{1.0}; ///< Bus-1 real voltage
        ScalarT Vi1{1.0}; ///< Bus-1 imaginary voltage

        const ScalarT delta{M_PI / 4.0};   /// first residual
        const ScalarT omega{0.0};          /// second residual
        const ScalarT Te{6.0};             /// third residual
        const ScalarT ir{1.0};             /// fourth residual
        const ScalarT ii{2.0};             /// fifth residual
        const ScalarT pmech{6.0};          /// fifth residual
        const ScalarT Ep{2.0 * sqrt(2.0)}; /// fifth residual

        const ScalarT tol = 5.0 * (std::numeric_limits<double>::epsilon()); // tolerance for comparing result

        PhasorDynamics::Bus<ScalarT, IdxT>          bus(Vr1, Vi1);
        PhasorDynamics::GenClassical<ScalarT, IdxT> gen(&bus, 1, p0, q0, H, D, Ra, Xdp);
        bus.allocate();
        bus.initialize();
        gen.allocate();
        gen.initialize();
        gen.evaluateResidual();
        std::vector<double> res = gen.getResidual();

        success *= isEqual(res[0], 0.0, tol);
        success *= isEqual(res[1], 0.0, tol);
        success *= isEqual(res[2], 0.0, tol);
        success *= isEqual(res[3], 0.0, tol);
        success *= isEqual(res[4], 0.0, tol);
        success *= isEqual(res[5], 0.0, tol);
        success *= isEqual(res[6], 0.0, tol);

        return success.report(__func__);
      }

    }; // class BranchTest

  } // namespace Testing
} // namespace GridKit
