#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <limits>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/PhasorDynamics/Branch/Branch.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace Testing
  {
    using Log = ::GridKit::Utilities::Logger;

    template <class ScalarT, typename IdxT>
    class BranchTests
    {
    private:
      using RealT = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

    public:
      BranchTests()  = default;
      ~BranchTests() = default;

      /// Wire a branch to the voltage signal nodes of its two terminal buses.
      template <typename T>
      static void wireBranch(PhasorDynamics::BusBase<T, IdxT>&    bus1,
                             PhasorDynamics::BusBase<T, IdxT>&    bus2,
                             PhasorDynamics::Branch<T, IdxT>&     branch,
                             PhasorDynamics::SignalNode<T, IdxT>& bus1_vr,
                             PhasorDynamics::SignalNode<T, IdxT>& bus1_vi,
                             PhasorDynamics::SignalNode<T, IdxT>& bus2_vr,
                             PhasorDynamics::SignalNode<T, IdxT>& bus2_vi)
      {
        using PhasorDynamics::BranchExternalVariables;
        using PhasorDynamics::BusInternalVariables;

        bus1.getSignals().template assignSignalNode<BusInternalVariables::VR>(&bus1_vr);
        bus1.getSignals().template assignSignalNode<BusInternalVariables::VI>(&bus1_vi);
        bus2.getSignals().template assignSignalNode<BusInternalVariables::VR>(&bus2_vr);
        bus2.getSignals().template assignSignalNode<BusInternalVariables::VI>(&bus2_vi);
        branch.getSignals().template attachSignalNode<BranchExternalVariables::VR1>(&bus1_vr);
        branch.getSignals().template attachSignalNode<BranchExternalVariables::VI1>(&bus1_vi);
        branch.getSignals().template attachSignalNode<BranchExternalVariables::VR2>(&bus2_vr);
        branch.getSignals().template attachSignalNode<BranchExternalVariables::VI2>(&bus2_vi);
      }

      TestOutcome constructor()
      {
        // Verifies Branch construction through the component interface.
        TestStatus success = true;

        PhasorDynamics::Component<ScalarT, IdxT>* branch =
            new PhasorDynamics::Branch<ScalarT, IdxT>();

        success *= (branch != nullptr);

        if (branch)
        {
          delete branch;
        }

        return success.report(__func__);
      }

      TestOutcome residual()
      {
        // Verifies nominal pi-branch current residual contributions.
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

        PhasorDynamics::Bus<ScalarT, IdxT>        bus1(Vr1, Vi1);
        PhasorDynamics::Bus<ScalarT, IdxT>        bus2(Vr2, Vi2);
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus1_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus1_vi;
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus2_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus2_vi;
        PhasorDynamics::Branch<ScalarT, IdxT>     branch(R, X, G, B);
        wireBranch(bus1, bus2, branch, bus1_vr, bus1_vi, bus2_vr, bus2_vi);

        bus1.allocate();
        bus1.initialize();
        bus1.evaluateResidual();
        bus2.allocate();
        bus2.initialize();
        bus2.evaluateResidual();

        branch.allocate();
        branch.evaluateResidual();

        success *= isEqual(branch.getExternalResidual()[0], Ir1);
        success *= isEqual(branch.getExternalResidual()[1], Ii1);
        success *= isEqual(branch.getExternalResidual()[2], Ir2);
        success *= isEqual(branch.getExternalResidual()[3], Ii2);

        return success.report(__func__);
      }

      TestOutcome offNominalResidual()
      {
        // Verifies tap and phase-shift current residual contributions.
        TestStatus success = true;

        RealT R{2.0};
        RealT X{4.0};
        RealT G{0.2};
        RealT B{1.2};
        RealT tap{1.25};
        RealT phase{0.3};

        ScalarT Vr1{10.0};
        ScalarT Vi1{20.0};
        ScalarT Vr2{30.0};
        ScalarT Vi2{40.0};

        const ScalarT Ir1{12.719793434963478};
        const ScalarT Ii1{-4.047960563981182};
        const ScalarT Ir2{13.821345956502421};
        const ScalarT Ii2{-21.182080826645354};

        PhasorDynamics::Bus<ScalarT, IdxT>        bus1(Vr1, Vi1);
        PhasorDynamics::Bus<ScalarT, IdxT>        bus2(Vr2, Vi2);
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus1_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus1_vi;
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus2_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus2_vi;
        PhasorDynamics::Branch<ScalarT, IdxT>     branch(R, X, G, B, tap, phase);
        wireBranch(bus1, bus2, branch, bus1_vr, bus1_vi, bus2_vr, bus2_vi);

        bus1.allocate();
        bus1.initialize();
        bus1.evaluateResidual();
        bus2.allocate();
        bus2.initialize();
        bus2.evaluateResidual();

        branch.allocate();
        branch.evaluateResidual();

        success *= isEqual(branch.getExternalResidual()[0], Ir1);
        success *= isEqual(branch.getExternalResidual()[1], Ii1);
        success *= isEqual(branch.getExternalResidual()[2], Ir2);
        success *= isEqual(branch.getExternalResidual()[3], Ii2);

        return success.report(__func__);
      }

      TestOutcome jacobian()
      {
        // Verifies nominal branch residual derivatives.
        TestStatus success = true;

        RealT R{2.0}; ///< Branch series resistance
        RealT X{4.0}; ///< Branch series reactance
        RealT G{0.2}; ///< Branch shunt conductance
        RealT B{1.2}; ///< Branch shunt charging

        DependencyTracking::Variable Vr1{10.0}; ///< Bus-1 real voltage
        DependencyTracking::Variable Vi1{20.0}; ///< Bus-1 imaginary voltage
        DependencyTracking::Variable Vr2{30.0}; ///< Bus-2 real voltage
        DependencyTracking::Variable Vi2{40.0}; ///< Bus-2 imaginary voltage

        PhasorDynamics::Bus<DependencyTracking::Variable, IdxT>        bus1(Vr1, Vi1);
        PhasorDynamics::Bus<DependencyTracking::Variable, IdxT>        bus2(Vr2, Vi2);
        PhasorDynamics::SignalNode<DependencyTracking::Variable, IdxT> bus1_vr;
        PhasorDynamics::SignalNode<DependencyTracking::Variable, IdxT> bus1_vi;
        PhasorDynamics::SignalNode<DependencyTracking::Variable, IdxT> bus2_vr;
        PhasorDynamics::SignalNode<DependencyTracking::Variable, IdxT> bus2_vi;
        PhasorDynamics::Branch<DependencyTracking::Variable, IdxT>     branch(R, X, G, B);
        wireBranch(bus1, bus2, branch, bus1_vr, bus1_vi, bus2_vr, bus2_vi);

        bus1.allocate();
        bus1.initialize();
        bus1.evaluateResidual();
        bus2.allocate();
        bus2.initialize();
        bus2.evaluateResidual();
        bus1.Vr().setVariableNumber(0);
        bus1.Vi().setVariableNumber(1);
        bus2.Vr().setVariableNumber(2);
        bus2.Vi().setVariableNumber(3);

        branch.allocate();
        branch.evaluateResidual(); ///< Computes the residual and the Jacobian values by tracking
                                   ///< the dependencies

        const auto&                                              f_ext = branch.getExternalResidual();
        std::vector<DependencyTracking::Variable>                residuals{f_ext[0], f_ext[1], f_ext[2], f_ext[3]};
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

      TestOutcome offNominalJacobian()
      {
        TestStatus success = true;

        RealT R{2.0};
        RealT X{4.0};
        RealT G{0.2};
        RealT B{1.2};
        RealT tap{1.25};
        RealT phase{0.3};

        DependencyTracking::Variable Vr1{10.0};
        DependencyTracking::Variable Vi1{20.0};
        DependencyTracking::Variable Vr2{30.0};
        DependencyTracking::Variable Vi2{40.0};

        Vr1.setVariableNumber(0);
        Vi1.setVariableNumber(1);
        Vr2.setVariableNumber(2);
        Vi2.setVariableNumber(3);

        PhasorDynamics::BusInfinite<DependencyTracking::Variable, IdxT> bus1(Vr1, Vi1);
        PhasorDynamics::BusInfinite<DependencyTracking::Variable, IdxT> bus2(Vr2, Vi2);
        PhasorDynamics::SignalNode<DependencyTracking::Variable, IdxT>  bus1_vr;
        PhasorDynamics::SignalNode<DependencyTracking::Variable, IdxT>  bus1_vi;
        PhasorDynamics::SignalNode<DependencyTracking::Variable, IdxT>  bus2_vr;
        PhasorDynamics::SignalNode<DependencyTracking::Variable, IdxT>  bus2_vi;

        PhasorDynamics::Branch<DependencyTracking::Variable, IdxT> branch(R,
                                                                          X,
                                                                          G,
                                                                          B,
                                                                          tap,
                                                                          phase);
        wireBranch(bus1, bus2, branch, bus1_vr, bus1_vi, bus2_vr, bus2_vi);

        bus1.allocate();
        bus2.allocate();
        branch.allocate();
        branch.evaluateResidual();

        const auto&                                              f_ext = branch.getExternalResidual();
        std::vector<DependencyTracking::Variable>                residuals{f_ext[0], f_ext[1], f_ext[2], f_ext[3]};
        std::vector<DependencyTracking::Variable::DependencyMap> ref = analyticalJacobian(R, X, G, B, tap, phase);

        for (size_t i = 0; i < residuals.size(); ++i)
        {
          DependencyTracking::Variable                       res           = residuals[i];
          const DependencyTracking::Variable::DependencyMap& dependencies  = res.getDependencies();
          success                                                         *= (GridKit::Testing::isEqual(dependencies, ref[i]));
        }

        return success.report(__func__);
      }

      TestOutcome parameterSetters()
      {
        // Verifies parameter setters refresh derived admittance values.
        TestStatus success = true;

        const RealT R{2.0};
        const RealT X{4.0};
        const RealT G{0.2};
        const RealT B{1.2};
        const RealT tap{1.25};
        const RealT phase{0.3};

        const ScalarT Vr1{10.0};
        const ScalarT Vi1{20.0};
        const ScalarT Vr2{30.0};
        const ScalarT Vi2{40.0};

        PhasorDynamics::Bus<ScalarT, IdxT>        ref_bus1(Vr1, Vi1);
        PhasorDynamics::Bus<ScalarT, IdxT>        ref_bus2(Vr2, Vi2);
        PhasorDynamics::SignalNode<ScalarT, IdxT> ref_bus1_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ref_bus1_vi;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ref_bus2_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ref_bus2_vi;
        PhasorDynamics::Branch<ScalarT, IdxT>     ref_branch(R, X, G, B, tap, phase);
        wireBranch(ref_bus1, ref_bus2, ref_branch, ref_bus1_vr, ref_bus1_vi, ref_bus2_vr, ref_bus2_vi);
        ref_bus1.allocate();
        ref_bus1.initialize();
        ref_bus1.evaluateResidual();
        ref_bus2.allocate();
        ref_bus2.initialize();
        ref_bus2.evaluateResidual();
        ref_branch.allocate();

        PhasorDynamics::Bus<ScalarT, IdxT>        test_bus1(Vr1, Vi1);
        PhasorDynamics::Bus<ScalarT, IdxT>        test_bus2(Vr2, Vi2);
        PhasorDynamics::SignalNode<ScalarT, IdxT> test_bus1_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT> test_bus1_vi;
        PhasorDynamics::SignalNode<ScalarT, IdxT> test_bus2_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT> test_bus2_vi;
        PhasorDynamics::Branch<ScalarT, IdxT>     test_branch(1.0, 1.0, 0.0, 0.0);
        wireBranch(test_bus1, test_bus2, test_branch, test_bus1_vr, test_bus1_vi, test_bus2_vr, test_bus2_vi);
        test_bus1.allocate();
        test_bus1.initialize();
        test_bus1.evaluateResidual();
        test_bus2.allocate();
        test_bus2.initialize();
        test_bus2.evaluateResidual();
        test_branch.allocate();

        test_branch.setR(R);
        test_branch.setX(X);
        test_branch.setG(G);
        test_branch.setB(B);
        test_branch.setTap(tap);
        test_branch.setPhase(phase);

        ref_branch.evaluateResidual();
        test_branch.evaluateResidual();

        success *= isEqual(test_branch.getExternalResidual()[0], ref_branch.getExternalResidual()[0]);
        success *= isEqual(test_branch.getExternalResidual()[1], ref_branch.getExternalResidual()[1]);
        success *= isEqual(test_branch.getExternalResidual()[2], ref_branch.getExternalResidual()[2]);
        success *= isEqual(test_branch.getExternalResidual()[3], ref_branch.getExternalResidual()[3]);

        return success.report(__func__);
      }

      TestOutcome parameterValidation()
      {
        // Verifies invalid branch parameters are rejected.
        TestStatus success = true;

        Log::setVerbosity(Log::Verbosity::EVERYTHING);
        Log::misc() << "Testing that invalid branch parameters are rejected. "
                    << "Logged errors are are expected.\n";
        Log::setVerbosity(Log::Verbosity::WARNINGS);

        PhasorDynamics::Bus<ScalarT, IdxT>        bus1(1.0, 0.0);
        PhasorDynamics::Bus<ScalarT, IdxT>        bus2(1.0, 0.0);
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus1_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus1_vi;
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus2_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT> bus2_vi;

        PhasorDynamics::Branch<ScalarT, IdxT> valid_branch(0.0, 0.1, 0.0, 0.0);
        wireBranch(bus1, bus2, valid_branch, bus1_vr, bus1_vi, bus2_vr, bus2_vi);
        success *= (valid_branch.verify() == 0);

        PhasorDynamics::Branch<ScalarT, IdxT> unwired_branch(0.0, 0.1, 0.0, 0.0);
        success *= (unwired_branch.verify() != 0);

        PhasorDynamics::Branch<ScalarT, IdxT> zero_impedance_branch(0.0, 0.0, 0.0, 0.0);
        wireBranch(bus1, bus2, zero_impedance_branch, bus1_vr, bus1_vi, bus2_vr, bus2_vi);
        success *= (zero_impedance_branch.verify() != 0);

        PhasorDynamics::Branch<ScalarT, IdxT> zero_tap_branch(0.0, 0.1, 0.0, 0.0, 0.0, 0.0);
        wireBranch(bus1, bus2, zero_tap_branch, bus1_vr, bus1_vi, bus2_vr, bus2_vi);
        success *= (zero_tap_branch.verify() != 0);

        PhasorDynamics::Branch<ScalarT, IdxT> negative_tap_branch(0.0, 0.1, 0.0, 0.0, -1.0, 0.0);
        wireBranch(bus1, bus2, negative_tap_branch, bus1_vr, bus1_vi, bus2_vr, bus2_vi);
        success *= (negative_tap_branch.verify() != 0);

        const RealT                           nan = std::numeric_limits<RealT>::quiet_NaN();
        PhasorDynamics::Branch<ScalarT, IdxT> nonfinite_branch(nan, 0.1, 0.0, 0.0);
        wireBranch(bus1, bus2, nonfinite_branch, bus1_vr, bus1_vi, bus2_vr, bus2_vi);
        success *= (nonfinite_branch.verify() != 0);

        return success.report(__func__);
      }

      TestOutcome dataConstructorDefaults()
      {
        // Verifies omitted data parameters use tap and phase defaults.
        TestStatus success = true;

        const RealT R{2.0};
        const RealT X{4.0};
        const RealT G{0.2};
        const RealT B{1.2};

        const ScalarT Vr1{10.0};
        const ScalarT Vi1{20.0};
        const ScalarT Vr2{30.0};
        const ScalarT Vi2{40.0};

        using Data      = typename PhasorDynamics::Branch<ScalarT, IdxT>::ModelDataT;
        using Parameter = typename Data::Parameters;
        using Buses     = typename Data::Buses;

        Data data;
        data.buses[Buses::bus1]       = 1;
        data.buses[Buses::bus2]       = 2;
        data.parameters[Parameter::R] = R;
        data.parameters[Parameter::X] = X;
        data.parameters[Parameter::G] = G;
        data.parameters[Parameter::B] = B;

        PhasorDynamics::Bus<ScalarT, IdxT>        data_bus1(Vr1, Vi1);
        PhasorDynamics::Bus<ScalarT, IdxT>        data_bus2(Vr2, Vi2);
        PhasorDynamics::SignalNode<ScalarT, IdxT> data_bus1_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT> data_bus1_vi;
        PhasorDynamics::SignalNode<ScalarT, IdxT> data_bus2_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT> data_bus2_vi;

        PhasorDynamics::Bus<ScalarT, IdxT>        ref_bus1(Vr1, Vi1);
        PhasorDynamics::Bus<ScalarT, IdxT>        ref_bus2(Vr2, Vi2);
        PhasorDynamics::SignalNode<ScalarT, IdxT> ref_bus1_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ref_bus1_vi;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ref_bus2_vr;
        PhasorDynamics::SignalNode<ScalarT, IdxT> ref_bus2_vi;

        PhasorDynamics::Branch<ScalarT, IdxT> data_branch(data);
        PhasorDynamics::Branch<ScalarT, IdxT> ref_branch(R, X, G, B, 1.0, 0.0);
        wireBranch(data_bus1, data_bus2, data_branch, data_bus1_vr, data_bus1_vi, data_bus2_vr, data_bus2_vi);
        wireBranch(ref_bus1, ref_bus2, ref_branch, ref_bus1_vr, ref_bus1_vi, ref_bus2_vr, ref_bus2_vi);

        data_bus1.allocate();
        data_bus1.initialize();
        data_bus1.evaluateResidual();
        data_bus2.allocate();
        data_bus2.initialize();
        data_bus2.evaluateResidual();

        ref_bus1.allocate();
        ref_bus1.initialize();
        ref_bus1.evaluateResidual();
        ref_bus2.allocate();
        ref_bus2.initialize();
        ref_bus2.evaluateResidual();

        data_branch.allocate();
        ref_branch.allocate();

        data_branch.evaluateResidual();
        ref_branch.evaluateResidual();

        success *= isEqual(data_branch.getExternalResidual()[0], ref_branch.getExternalResidual()[0]);
        success *= isEqual(data_branch.getExternalResidual()[1], ref_branch.getExternalResidual()[1]);
        success *= isEqual(data_branch.getExternalResidual()[2], ref_branch.getExternalResidual()[2]);
        success *= isEqual(data_branch.getExternalResidual()[3], ref_branch.getExternalResidual()[3]);

        return success.report(__func__);
      }

    private:
      std::vector<DependencyTracking::Variable::DependencyMap> analyticalJacobian(const RealT R,
                                                                                  const RealT X,
                                                                                  const RealT G,
                                                                                  const RealT B,
                                                                                  const RealT tap   = 1.0,
                                                                                  const RealT phase = 0.0)
      {
        const RealT denom   = R * R + X * X;
        const RealT g       = R / denom;
        const RealT b       = -X / denom;
        const RealT inv_tap = RealT{1.0} / tap;

        const std::complex<RealT> ybr{g, b};
        const std::complex<RealT> ysh{G, B};
        const std::complex<RealT> rotation{std::cos(phase), std::sin(phase)};
        const std::complex<RealT> ydiag = -(ybr + RealT{0.5} * ysh);

        const std::complex<RealT> y11 = ydiag * inv_tap * inv_tap;
        const std::complex<RealT> y12 = ybr * rotation * inv_tap;
        const std::complex<RealT> y21 = ybr * std::conj(rotation) * inv_tap;
        const std::complex<RealT> y22 = ydiag;

        std::vector<DependencyTracking::Variable::DependencyMap> dependencies(4);
        dependencies[0] = {{0, y11.real()}, {1, -y11.imag()}, {2, y12.real()}, {3, -y12.imag()}};
        dependencies[1] = {{0, y11.imag()}, {1, y11.real()}, {2, y12.imag()}, {3, y12.real()}};
        dependencies[2] = {{0, y21.real()}, {1, -y21.imag()}, {2, y22.real()}, {3, -y22.imag()}};
        dependencies[3] = {{0, y21.imag()}, {1, y21.real()}, {2, y22.imag()}, {3, y22.real()}};

        return dependencies;
      }
    }; // class BranchTest

  } // namespace Testing
} // namespace GridKit
