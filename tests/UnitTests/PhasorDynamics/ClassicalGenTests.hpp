#include <iomanip>
#include <iostream>

#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <Model/PhasorDynamics/SynchronousMachine/ClassicalGenerator/ClassicalGen.hpp>
#include <Model/PhasorDynamics/SynchronousMachine/ClassicalGenerator/ClassicalGen.cpp>
#include <Utilities/TestHelpers.hpp>
#include <Utilities/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {

    template <class ScalarT, typename IdxT>
    class ClassicalGenTests
    {
    private:
      using real_type = typename PhasorDynamics::Component<ScalarT, IdxT>::real_type;

    public:
      ClassicalGenTests()  = default;
      ~ClassicalGenTests() = default;

      TestOutcome constructor()
      {
        TestStatus success = true;

        auto* bus = new PhasorDynamics::Bus<ScalarT, IdxT>(1.0, 0.0);

        PhasorDynamics::Component<ScalarT, IdxT>* machine =
            new PhasorDynamics::ClassicalGen<ScalarT, IdxT>(bus, 1);

        success *= (machine != nullptr);

        if (machine)
        {
          delete machine;
        }
        delete bus;

        return success.report(__func__);
      }

      TestOutcome residual()
      {
        TestStatus success = true;
        
        // classical generator parameters
        real_type H{0.1}; 
        real_type D{2.35};
        real_type Ra{1.5}; 
        real_type Xdp{4.5}; 
        real_type pmech{5.0}; 
        real_type ep{2.5}; 

        ScalarT Vr1{2.0}; ///< Bus-1 real voltage
        ScalarT Vi1{1.5}; ///< Bus-1 imaginary voltage

        const ScalarT res0{-1128.973355292326};  /// first residual
        const ScalarT res1{27.5625000000000}; /// second residual
        const ScalarT res2{4.102511525891203};  /// third residual
        const ScalarT res3{8.164018441425924}; /// fourth residual
        const ScalarT res4{2.089603682931281}; /// fifth residual
        const ScalarT res5{5.0}; /// fifth residual
        const ScalarT res6{2.5}; /// fifth residual

        const ScalarT tol{0.0000001}; //tolerance for comparing result

        PhasorDynamics::Bus<ScalarT, IdxT> bus(Vr1, Vi1);
        PhasorDynamics::ClassicalGen<ScalarT, IdxT> gen(&bus, 1, 1, 1, H, D, Ra, Xdp);
        bus.allocate();
        bus.initialize();
        gen.allocate(); 

        

        gen.y()[0] = 1.0; //delta
        gen.y()[1] = 3.0; //omega
        gen.y()[2] = 4.0; //telec
        gen.y()[3] = 8.0; //ir
        gen.y()[4] = 2.0; //ii
        gen.y()[5] = 5.0; //pmech
        gen.y()[6] = 2.5; //Ep

        gen.yp()[0] = 2.0; //delta_dot
        gen.yp()[1] = 5.0; //omega_dot
        gen.yp()[2] = 0.0; //telec
        gen.yp()[3] = 0.0; //ir
        gen.yp()[4] = 0.0; //ii
        gen.yp()[5] = 0.0; //pmech
        gen.yp()[6] = 0.0; //Ep

        

        gen.evaluateResidual();

        std::vector<ScalarT> residual = gen.getResidual();

        success *= isEqual(residual[0], res0, tol);
        success *= isEqual(residual[1], res1, tol);
        success *= isEqual(residual[2], res2, tol);
        success *= isEqual(residual[3], res3, tol);
        success *= isEqual(residual[4], res4, tol);
        success *= isEqual(residual[5], res5, tol);
        success *= isEqual(residual[6], res6, tol);

        return success.report(__func__);
      }

      TestOutcome initial()
      {
        TestStatus success = true;
        
        // classical generator parameters
        real_type p0{1.12}; 
        real_type q0{0.35};
        real_type H{0.85}; 
        real_type D{2.77};
        real_type Ra{1.55}; 
        real_type Xdp{2.15}; 

        ScalarT Vr1{2.0}; ///< Bus-1 real voltage
        ScalarT Vi1{1.5}; ///< Bus-1 imaginary voltage

        const ScalarT delta{0.2562082853611203};  /// first residual
        const ScalarT omega{0.0}; /// second residual
        const ScalarT Te{1.461471200000000};  /// third residual
        const ScalarT ir{0.4424000000000000}; /// fourth residual
        const ScalarT ii{0.1568000000000000}; /// fifth residual
        const ScalarT pmech{1.461471200000000}; /// fifth residual
        const ScalarT Ep{3.124841691990172}; /// fifth residual

        const ScalarT tol{0.0000001}; //tolerance for comparing result

        PhasorDynamics::Bus<ScalarT, IdxT> bus(Vr1, Vi1);
        PhasorDynamics::ClassicalGen<ScalarT, IdxT> gen(&bus, 1, p0, q0, H, D, Ra, Xdp);
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


      
    }; // class BranchTest

  } // namespace Testing
} // namespace GridKit
