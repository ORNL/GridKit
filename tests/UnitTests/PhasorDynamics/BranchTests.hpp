#include <complex>
#include <iomanip>
#include <iostream>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/PhasorDynamics/Branch/Branch.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
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

      TestOutcome constructor()
      {
        // Verifies Branch construction through the component interface.
        TestStatus success = true;

        auto* bus1 = new PhasorDynamics::Bus<ScalarT, IdxT>(1.0, 0.0);
        auto* bus2 = new PhasorDynamics::Bus<ScalarT, IdxT>(1.0, 0.1);

        PhasorDynamics::Component<ScalarT, IdxT>* branch =
            new PhasorDynamics::Branch<ScalarT, IdxT>(bus1, bus2);

        success *= (branch != nullptr);
        success *= branch->size()
                   == static_cast<IdxT>(PhasorDynamics::BranchInternalVariables::MAXIMUM);

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

        PhasorDynamics::Bus<ScalarT, IdxT> bus1(Vr1, Vi1);
        PhasorDynamics::Bus<ScalarT, IdxT> bus2(Vr2, Vi2);
        bus1.allocate();
        bus1.initialize();
        bus1.evaluateResidual();
        bus2.allocate();
        bus2.initialize();
        bus2.evaluateResidual();

        PhasorDynamics::Branch<ScalarT, IdxT> branch(&bus1, &bus2, R, X, G, B);
        branch.allocate();
        branch.initialize();
        branch.tagDifferentiable();
        branch.setAbsoluteTolerance(static_cast<RealT>(1e-6));
        branch.evaluateResidual();

        const auto* y          = branch.y().getData();
        const auto* yp         = branch.yp().getData();
        const auto* f          = branch.getResidual().getData();
        const auto* tolerances = branch.absoluteTolerance().getData();
        const auto  size       = static_cast<std::size_t>(branch.size());

        success *= isEqual(y[0], Ir1);
        success *= isEqual(y[1], Ii1);
        success *= isEqual(y[2], Ir2);
        success *= isEqual(y[3], Ii2);
        success *= isEqual(bus1.Ir(), Ir1);
        success *= isEqual(bus1.Ii(), Ii1);
        success *= isEqual(bus2.Ir(), Ir2);
        success *= isEqual(bus2.Ii(), Ii2);
        for (std::size_t i = 0; i < size; ++i)
        {
          success *= isEqual(yp[i], ScalarT{0.0});
          success *= isEqual(f[i], ScalarT{0.0}, static_cast<RealT>(100.0) * std::numeric_limits<RealT>::epsilon());
          success *= !branch.tag()[i];
          success *= isEqual(tolerances[i], static_cast<ScalarT>(1e-6));
        }

        // The buses consume the current variables while the Branch equations
        // continue to enforce I = YV.
        constexpr ScalarT current_perturbation[] = {1.0, 2.0, 3.0, 4.0};
        auto*             mutable_y              = branch.y().getData();
        for (std::size_t i = 0; i < size; ++i)
        {
          mutable_y[i] += current_perturbation[i];
        }
        branch.y().setDataUpdated();
        bus1.evaluateResidual();
        bus2.evaluateResidual();
        branch.evaluateResidual();

        success *= isEqual(bus1.Ir(), Ir1 + current_perturbation[0]);
        success *= isEqual(bus1.Ii(), Ii1 + current_perturbation[1]);
        success *= isEqual(bus2.Ir(), Ir2 + current_perturbation[2]);
        success *= isEqual(bus2.Ii(), Ii2 + current_perturbation[3]);
        for (std::size_t i = 0; i < size; ++i)
        {
          success *= isEqual(f[i], -current_perturbation[i], static_cast<RealT>(100.0) * std::numeric_limits<RealT>::epsilon());
        }

        return success.report(__func__);
      }

      TestOutcome offNominalResidual()
      {
        // Verifies tap, phase-shift, line-shunt, and bus-1 magnetizing-shunt contributions.
        TestStatus success = true;

        RealT R{2.0};
        RealT X{4.0};
        RealT G{0.4};
        RealT B{0.8};
        RealT Gmag{0.2};
        RealT Bmag{1.2};
        RealT tap{1.25};
        RealT phase{0.3};

        ScalarT Vr1{10.0};
        ScalarT Vi1{20.0};
        ScalarT Vr2{30.0};
        ScalarT Vi2{40.0};

        const ScalarT Ir1{33.679793434963472};
        const ScalarT Ii1{-22.927960563981181};
        const ScalarT Ir2{2.821345956502423};
        const ScalarT Ii2{-19.182080826645358};

        PhasorDynamics::Bus<ScalarT, IdxT> bus1(Vr1, Vi1);
        PhasorDynamics::Bus<ScalarT, IdxT> bus2(Vr2, Vi2);
        bus1.allocate();
        bus1.initialize();
        bus1.evaluateResidual();
        bus2.allocate();
        bus2.initialize();
        bus2.evaluateResidual();

        using Data      = typename PhasorDynamics::Branch<ScalarT, IdxT>::ModelDataT;
        using Parameter = typename Data::Parameters;

        Data data;
        data.parameters[Parameter::R]     = R;
        data.parameters[Parameter::X]     = X;
        data.parameters[Parameter::G]     = G;
        data.parameters[Parameter::B]     = B;
        data.parameters[Parameter::Gmag]  = Gmag;
        data.parameters[Parameter::Bmag]  = Bmag;
        data.parameters[Parameter::tap]   = tap;
        data.parameters[Parameter::phase] = phase;

        PhasorDynamics::Branch<ScalarT, IdxT> branch(&bus1, &bus2, data);
        branch.allocate();
        branch.initialize();
        branch.evaluateResidual();

        success *= isEqual(bus1.Ir(), Ir1);
        success *= isEqual(bus1.Ii(), Ii1);
        success *= isEqual(bus2.Ir(), Ir2);
        success *= isEqual(bus2.Ii(), Ii2);

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

        PhasorDynamics::Bus<DependencyTracking::Variable, IdxT> bus1(Vr1, Vi1);
        PhasorDynamics::Bus<DependencyTracking::Variable, IdxT> bus2(Vr2, Vi2);
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

        PhasorDynamics::Branch<DependencyTracking::Variable, IdxT> branch(&bus1, &bus2, R, X, G, B);
        branch.allocate();
        branch.initialize();
        for (size_t i = 0; i < static_cast<size_t>(branch.size()); ++i)
        {
          branch.y().getData()[i].setVariableNumber(i + 4);
        }
        branch.y().setDataUpdated();
        branch.evaluateResidual(); ///< Computes the residual and the Jacobian values by tracking
                                   ///< the dependencies

        const auto*                               branch_residual = branch.getResidual().getData();
        std::vector<DependencyTracking::Variable> residuals{
            branch_residual[0],
            branch_residual[1],
            branch_residual[2],
            branch_residual[3],
            bus1.Ir(),
            bus1.Ii(),
            bus2.Ir(),
            bus2.Ii()};
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

        PhasorDynamics::Branch<DependencyTracking::Variable, IdxT> branch(&bus1,
                                                                          &bus2,
                                                                          R,
                                                                          X,
                                                                          G,
                                                                          B,
                                                                          tap,
                                                                          phase);
        branch.allocate();
        branch.initialize();
        for (size_t i = 0; i < static_cast<size_t>(branch.size()); ++i)
        {
          branch.y().getData()[i].setVariableNumber(i + 4);
        }
        branch.y().setDataUpdated();
        branch.evaluateResidual();

        const auto*                               branch_residual = branch.getResidual().getData();
        std::vector<DependencyTracking::Variable> residuals{
            branch_residual[0],
            branch_residual[1],
            branch_residual[2],
            branch_residual[3]};
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

        PhasorDynamics::Bus<ScalarT, IdxT> ref_bus1(Vr1, Vi1);
        PhasorDynamics::Bus<ScalarT, IdxT> ref_bus2(Vr2, Vi2);
        ref_bus1.allocate();
        ref_bus1.initialize();
        ref_bus1.evaluateResidual();
        ref_bus2.allocate();
        ref_bus2.initialize();
        ref_bus2.evaluateResidual();
        PhasorDynamics::Branch<ScalarT, IdxT> ref_branch(&ref_bus1, &ref_bus2, R, X, G, B, tap, phase);
        ref_branch.allocate();

        PhasorDynamics::Bus<ScalarT, IdxT> test_bus1(Vr1, Vi1);
        PhasorDynamics::Bus<ScalarT, IdxT> test_bus2(Vr2, Vi2);
        test_bus1.allocate();
        test_bus1.initialize();
        test_bus1.evaluateResidual();
        test_bus2.allocate();
        test_bus2.initialize();
        test_bus2.evaluateResidual();
        PhasorDynamics::Branch<ScalarT, IdxT> test_branch(&test_bus1, &test_bus2, 1.0, 1.0, 0.0, 0.0);
        test_branch.allocate();

        test_branch.setR(R);
        test_branch.setX(X);
        test_branch.setG(G);
        test_branch.setB(B);
        test_branch.setTap(tap);
        test_branch.setPhase(phase);

        ref_branch.initialize();
        test_branch.initialize();
        ref_branch.evaluateResidual();
        test_branch.evaluateResidual();

        success *= isEqual(test_bus1.Ir(), ref_bus1.Ir());
        success *= isEqual(test_bus1.Ii(), ref_bus1.Ii());
        success *= isEqual(test_bus2.Ir(), ref_bus2.Ir());
        success *= isEqual(test_bus2.Ii(), ref_bus2.Ii());

        return success.report(__func__);
      }

      TestOutcome parameterValidation()
      {
        // Verifies invalid branch parameters are rejected.
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus1(1.0, 0.0);
        PhasorDynamics::Bus<ScalarT, IdxT> bus2(1.0, 0.0);

        PhasorDynamics::Branch<ScalarT, IdxT> valid_branch(&bus1, &bus2, 0.0, 0.1, 0.0, 0.0);
        success *= (valid_branch.verify() == 0);

        const auto previous_verbosity = Log::verbosity();
        // Suppress expected errors from the invalid branch configurations below.
        // Use EVERYTHING to inspect those diagnostics.
        Log::setVerbosity(Log::Verbosity::NONE);

        PhasorDynamics::Branch<ScalarT, IdxT> zero_impedance_branch(&bus1, &bus2, 0.0, 0.0, 0.0, 0.0);
        success *= (zero_impedance_branch.verify() != 0);

        PhasorDynamics::Branch<ScalarT, IdxT> zero_tap_branch(&bus1, &bus2, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0);
        success *= (zero_tap_branch.verify() != 0);

        PhasorDynamics::Branch<ScalarT, IdxT> negative_tap_branch(&bus1, &bus2, 0.0, 0.1, 0.0, 0.0, -1.0, 0.0);
        success *= (negative_tap_branch.verify() != 0);

        const RealT                           nan = std::numeric_limits<RealT>::quiet_NaN();
        PhasorDynamics::Branch<ScalarT, IdxT> nonfinite_branch(&bus1, &bus2, nan, 0.1, 0.0, 0.0);
        success *= (nonfinite_branch.verify() != 0);

        Log::setVerbosity(previous_verbosity);
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

        PhasorDynamics::Bus<ScalarT, IdxT> data_bus1(Vr1, Vi1);
        PhasorDynamics::Bus<ScalarT, IdxT> data_bus2(Vr2, Vi2);
        data_bus1.allocate();
        data_bus1.initialize();
        data_bus1.evaluateResidual();
        data_bus2.allocate();
        data_bus2.initialize();
        data_bus2.evaluateResidual();

        PhasorDynamics::Bus<ScalarT, IdxT> ref_bus1(Vr1, Vi1);
        PhasorDynamics::Bus<ScalarT, IdxT> ref_bus2(Vr2, Vi2);
        ref_bus1.allocate();
        ref_bus1.initialize();
        ref_bus1.evaluateResidual();
        ref_bus2.allocate();
        ref_bus2.initialize();
        ref_bus2.evaluateResidual();

        PhasorDynamics::Branch<ScalarT, IdxT> data_branch(&data_bus1, &data_bus2, data);
        PhasorDynamics::Branch<ScalarT, IdxT> ref_branch(&ref_bus1, &ref_bus2, R, X, G, B, 1.0, 0.0);

        data_branch.allocate();
        ref_branch.allocate();
        data_branch.initialize();
        ref_branch.initialize();
        data_branch.evaluateResidual();
        ref_branch.evaluateResidual();

        success *= isEqual(data_bus1.Ir(), ref_bus1.Ir());
        success *= isEqual(data_bus1.Ii(), ref_bus1.Ii());
        success *= isEqual(data_bus2.Ir(), ref_bus2.Ir());
        success *= isEqual(data_bus2.Ii(), ref_bus2.Ii());

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
        const std::complex<RealT> ydiag = -ybr;

        const std::complex<RealT> y11 = ydiag * inv_tap * inv_tap - RealT{0.5} * ysh;
        const std::complex<RealT> y12 = ybr * rotation * inv_tap;
        const std::complex<RealT> y21 = ybr * std::conj(rotation) * inv_tap;
        const std::complex<RealT> y22 = ydiag - RealT{0.5} * ysh;

        std::vector<DependencyTracking::Variable::DependencyMap> dependencies(8);
        dependencies[0] = {{0, y11.real()}, {1, -y11.imag()}, {2, y12.real()}, {3, -y12.imag()}, {4, -1.0}};
        dependencies[1] = {{0, y11.imag()}, {1, y11.real()}, {2, y12.imag()}, {3, y12.real()}, {5, -1.0}};
        dependencies[2] = {{0, y21.real()}, {1, -y21.imag()}, {2, y22.real()}, {3, -y22.imag()}, {6, -1.0}};
        dependencies[3] = {{0, y21.imag()}, {1, y21.real()}, {2, y22.imag()}, {3, y22.real()}, {7, -1.0}};
        dependencies[4] = {{4, 1.0}};
        dependencies[5] = {{5, 1.0}};
        dependencies[6] = {{6, 1.0}};
        dependencies[7] = {{7, 1.0}};

        return dependencies;
      }
    }; // class BranchTest

  } // namespace Testing
} // namespace GridKit
