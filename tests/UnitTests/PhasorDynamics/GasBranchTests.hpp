/**
 * @file GasBranchTests.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Tests for the gas branch (pipe) model.
 */
#include <limits>
#include <numbers>
#include <type_traits>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Constants.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Gas/GasBranch/GasBranch.hpp>
#include <GridKit/Model/PhasorDynamics/Gas/GasBus/GasBus.hpp>
#include <GridKit/Model/PhasorDynamics/Gas/GasBusInfinite/GasBusInfinite.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class GasBranchTests
    {
    private:
      using RealT           = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;
      using GasBusDataT     = PhasorDynamics::Gas::GasBusData<RealT, IdxT>;
      using GasBranchDataT  = PhasorDynamics::Gas::GasBranchData<RealT, IdxT>;
      using GasBusT         = PhasorDynamics::Gas::GasBus<ScalarT, IdxT>;
      using GasBusInfiniteT = PhasorDynamics::Gas::GasBusInfinite<ScalarT, IdxT>;
      using GasBranchT      = PhasorDynamics::Gas::GasBranch<ScalarT, IdxT>;
      using DependencyMap   = DependencyTracking::Variable::DependencyMap;

      static constexpr RealT tol_   = 10 * std::numeric_limits<RealT>::epsilon();
      static constexpr RealT pi_    = std::numbers::pi_v<RealT>;
      static constexpr RealT alpha_ = 2.0;

      static GasBusDataT busData(RealT pressure)
      {
        GasBusDataT data;
        data.p0 = pressure;
        return data;
      }

      /// A 2 km pipe of 0.25 m diameter, so A = pi/64 m^2, with a = 500 m/s
      static GasBranchDataT branchData(RealT dxmax)
      {
        using Parameter = typename GasBranchDataT::Parameters;

        GasBranchDataT data;
        data.device_class                  = "GasBranch";
        data.disambiguation_string         = "GP";
        data.parameters[Parameter::L]      = RealT{2000.0};
        data.parameters[Parameter::D]      = RealT{0.25};
        data.parameters[Parameter::lambda] = RealT{0.01};
        data.parameters[Parameter::a]      = RealT{500.0};
        data.parameters[Parameter::dxmax]  = dxmax;
        return data;
      }

      /// The same pipe with the GasTranSim simple CNGA gas law
      static GasBranchDataT cngaData(RealT dxmax)
      {
        using Parameter = typename GasBranchDataT::Parameters;

        auto data                      = branchData(dxmax);
        data.parameters[Parameter::b1] = RealT{1.00300865};
        data.parameters[Parameter::b2] = RealT{2.96848838e-8};
        return data;
      }

      /// Allocate and initialize, numbering the bus variables after the pipe's
      template <typename T>
      static void connect(PhasorDynamics::Gas::GasBusBase<T, IdxT>& bus1,
                          PhasorDynamics::Gas::GasBusBase<T, IdxT>& bus2,
                          PhasorDynamics::Gas::GasBranch<T, IdxT>&  branch)
      {
        bus1.allocate();
        bus2.allocate();
        branch.allocate();

        IdxT next = branch.size();
        for (auto* bus : {&bus1, &bus2})
        {
          if (bus->size() > 0)
          {
            bus->setVariableIndex(0, next);
            bus->setResidualIndex(0, next);
            ++next;
          }
        }

        bus1.initialize();
        bus2.initialize();
        branch.initialize();
      }

      static void assign(double& variable, double value)
      {
        variable = value;
      }

      static void assign(DependencyTracking::Variable& variable, double value)
      {
        variable.setValue(value);
      }

      /// A state away from steady flow for N = 2, with fluxes in both directions
      template <typename T>
      static void setState(PhasorDynamics::Gas::GasBus<T, IdxT>&    bus1,
                           PhasorDynamics::Gas::GasBus<T, IdxT>&    bus2,
                           PhasorDynamics::Gas::GasBranch<T, IdxT>& branch)
      {
        assign(bus1.yp().getData()[0], 2000.0);
        assign(bus2.yp().getData()[0], -3000.0);

        auto* y  = branch.y().getData();
        auto* yp = branch.yp().getData();
        assign(y[0], 4.0e6);
        assign(y[1], 100.0);
        assign(y[2], -70.0);
        assign(yp[0], 10.0);
        assign(yp[1], 1.0);
        assign(yp[2], 2.0);
      }

    public:
      GasBranchTests()  = default;
      ~GasBranchTests() = default;

      /// The segment count is the length over the maximum segment length, rounded up
      TestOutcome constructor()
      {
        TestStatus success = true;

        GasBusT bus1(busData(5.0e6));
        GasBusT bus2(busData(4.0e6));

        GasBranchT three_segments(&bus1, &bus2, branchData(800.0));
        GasBranchT one_segment(&bus1, &bus2, branchData(5000.0));

        success *= (three_segments.size() == 5);
        success *= (one_segment.size() == 1);

        return success.report(__func__);
      }

      TestOutcome parameterValidation()
      {
        TestStatus success = true;

        using Parameter = typename GasBranchDataT::Parameters;

        GasBusT bus1(busData(5.0e6));
        GasBusT bus2(busData(4.0e6));

        success *= GasBranchT(&bus1, &bus2, branchData(1000.0)).verify().passed();
        success *= !GasBranchT(nullptr, &bus2, branchData(1000.0)).verify().passed();

        for (auto parameter : {Parameter::L, Parameter::D, Parameter::lambda, Parameter::a, Parameter::b1, Parameter::dxmax})
        {
          auto data                   = branchData(1000.0);
          data.parameters[parameter]  = RealT{0.0};
          success                    *= !GasBranchT(&bus1, &bus2, data).verify().passed();
        }

        auto negative_b2                       = branchData(1000.0);
        negative_b2.parameters[Parameter::b2]  = RealT{-1.0e-8};
        success                               *= !GasBranchT(&bus1, &bus2, negative_b2).verify().passed();

        return success.report(__func__);
      }

      /**
       * @brief Residual away from steady flow, derived by hand.
       *
       * Nodes hold 6, 4, and 3 MPa, so the mean segment densities are 20 and
       * 14 kg/m^3 and the friction terms are 10 and -7.
       */
      TestOutcome residual()
      {
        TestStatus success = true;

        GasBusT    bus1(busData(6.0e6));
        GasBusT    bus2(busData(3.0e6));
        GasBranchT branch(&bus1, &bus2, branchData(1000.0));
        connect(bus1, bus2, branch);
        setState(bus1, bus2, branch);

        bus1.evaluateResidual();
        bus2.evaluateResidual();
        branch.evaluateResidual();

        const auto* f  = branch.getResidual().getData();
        success       *= isEqual(f[0], ScalarT{42490.0}, tol_); // -10 + 250 (100 + 70)
        success       *= isEqual(f[1], ScalarT{1989.0}, tol_);  // -1 + 2000 - 10
        success       *= isEqual(f[2], ScalarT{1005.0}, tol_);  // -2 + 1000 + 7

        // Flux and half-segment storage: -(pi/64)(100 + 4) and (pi/64)(-70 + 6)
        success *= isEqual(bus1.q(), ScalarT{-13.0 * pi_ / 8.0}, tol_);
        success *= isEqual(bus2.q(), ScalarT{-pi_}, tol_);

        return success.report(__func__);
      }

      /// Initialization gives the steady flow, with p^2 falling linearly along the pipe
      TestOutcome initialization()
      {
        TestStatus success = true;

        GasBusT    bus1(busData(6.0e6));
        GasBusT    bus2(busData(4.0e6));
        GasBranchT branch(&bus1, &bus2, branchData(500.0));
        connect(bus1, bus2, branch);

        // phi^2 = D (p_0^2 - p_N^2) / (lambda a^2 L) = 0.25 (2e13) / 5e6
        const auto* y  = branch.y().getData();
        const auto* yp = branch.yp().getData();
        for (IdxT k = 1; k <= 4; ++k)
        {
          success *= isEqual(y[k + 2], ScalarT{1000.0}, tol_);
        }
        for (IdxT j = 1; j < 4; ++j)
        {
          success *= isEqual(y[j - 1] * y[j - 1], ScalarT{3.6e13 - 0.5e13 * static_cast<RealT>(j)}, tol_);
        }
        for (IdxT i = 0; i < branch.size(); ++i)
        {
          success *= isEqual(yp[i], ScalarT{0.0}, tol_);
        }

        return success.report(__func__);
      }

      /// The initialized CNGA pipe is at rest, and both terminals carry the same mass flow
      TestOutcome initializationResidual()
      {
        TestStatus success = true;

        GasBusT    bus1(busData(6.0e6));
        GasBusT    bus2(busData(4.0e6));
        GasBranchT branch(&bus1, &bus2, cngaData(500.0));
        connect(bus1, bus2, branch);

        bus1.evaluateResidual();
        bus2.evaluateResidual();
        branch.evaluateResidual();

        // Momentum rows cancel terms of order (p_0 - p_N)/dx = 4000 Pa/m
        const auto* f = branch.getResidual().getData();
        for (IdxT i = 0; i < branch.size(); ++i)
        {
          success *= isEqual(f[i], ScalarT{0.0}, RealT{4000.0} * tol_);
        }

        success *= isEqual(bus1.q(), -bus2.q(), tol_);

        return success.report(__func__);
      }

      /// With the node balance met, the linepack changes only by the terminal mass flows
      TestOutcome massConservation()
      {
        TestStatus success = true;

        GasBusT    bus1(busData(6.0e6));
        GasBusT    bus2(busData(3.0e6));
        GasBranchT branch(&bus1, &bus2, branchData(1000.0));
        connect(bus1, bus2, branch);
        setState(bus1, bus2, branch);
        assign(branch.yp().getData()[0], 42500.0);

        bus1.evaluateResidual();
        bus2.evaluateResidual();
        branch.evaluateResidual();

        // Linepack A dx / a^2 [(p_0 + p_2)/2 + p_1] at the rates above
        const RealT area          = pi_ / 64.0;
        const RealT linepack_rate = area * 1000.0 / 250000.0 * ((2000.0 - 3000.0) / 2.0 + 42500.0);

        success *= isEqual(branch.getResidual().getData()[0], ScalarT{0.0}, tol_);
        success *= isEqual(bus1.q() + bus2.q(), ScalarT{-linepack_rate}, tol_);

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// The hand-written CNGA Jacobian matches DependencyTracking, bus rows included
      TestOutcome jacobian()
      {
        TestStatus success = true;

        const auto reference    = jacobianRows<DependencyTracking::Variable>();
        const auto hand_written = jacobianRows<ScalarT>();

        success *= (reference.size() == hand_written.size());
        for (size_t i = 0; i < reference.size(); ++i)
        {
          success *= isEqual(reference[i], hand_written[i], tol_);
        }

        return success.report(__func__);
      }

      /// A fixed-pressure bus 1 has no row and no column in the pipe's Jacobian
      TestOutcome fixedPressureBus()
      {
        TestStatus success = true;

        GasBusInfiniteT bus1(busData(6.0e6));
        GasBusT         bus2(busData(3.0e6));
        GasBranchT      branch(&bus1, &bus2, branchData(1000.0));
        connect(bus1, bus2, branch);

        branch.updateTime(0.0, alpha_);
        bus2.evaluateResidual();
        branch.evaluateResidual();
        branch.evaluateJacobian();

        // 6N + 1 entries less the two of the bus-1 row and the bus-1 column of segment 1
        success *= (branch.getCooJacobian()->getNnz() == 10);

        return success.report(__func__);
      }

    private:
      /// Jacobian rows of the pipe followed by the rows of buses 1 and 2
      template <typename T>
      std::vector<DependencyMap> jacobianRows()
      {
        PhasorDynamics::Gas::GasBus<T, IdxT>    bus1(busData(6.0e6));
        PhasorDynamics::Gas::GasBus<T, IdxT>    bus2(busData(3.0e6));
        PhasorDynamics::Gas::GasBranch<T, IdxT> branch(&bus1, &bus2, cngaData(1000.0));
        connect(bus1, bus2, branch);
        setState(bus1, bus2, branch);

        branch.updateTime(0.0, alpha_);
        bus1.evaluateResidual();
        bus2.evaluateResidual();
        branch.evaluateResidual();
        branch.evaluateJacobian();

        if constexpr (std::is_same_v<T, DependencyTracking::Variable>)
        {
          auto rows = MapFromCsr(branch.getCsrJacobian());
          rows.push_back(withDerivatives(bus1.q().getDependencies()));
          rows.push_back(withDerivatives(bus2.q().getDependencies()));
          return rows;
        }
        else
        {
          auto*                      jacobian = branch.getCooJacobian();
          std::vector<DependencyMap> rows(static_cast<size_t>(branch.size()) + 2);
          for (IdxT i = 0; i < jacobian->getNnz(); ++i)
          {
            rows[jacobian->getRowData()[i]][jacobian->getColData()[i]] += jacobian->getValues()[i];
          }
          return rows;
        }
      }

      /// Fold derivative dependencies (odd numbers) into their variables, scaled by alpha
      static DependencyMap withDerivatives(const DependencyMap& dependencies)
      {
        DependencyMap row;
        for (const auto& [number, value] : dependencies)
        {
          RealT scale = ONE<RealT>;
          if (number % 2 == 1)
          {
            scale = alpha_;
          }
          row[number / 2] += scale * value;
        }
        return row;
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
