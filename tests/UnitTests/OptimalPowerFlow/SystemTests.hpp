#pragma once

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <IpIpoptApplication.hpp>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/OptimalPowerFlow/SystemModel.hpp>
#include <GridKit/Model/OptimalPowerFlow/SystemModelData.hpp>
#include <GridKit/Solver/Optimization/SteadyState.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    namespace OPF = GridKit::OptimalPowerFlow;

    /**
     * @brief System assembly, derivatives, and an Ipopt solve on a three-bus
     * network
     *
     * Generator 1 is cheaper, but the rating of branch 1-2 forces generator
     * 3 to cover part of the load at bus 2.
     */
    class SystemTests
    {
      using RealT     = double;
      using IdxT      = size_t;
      using VariableT = DependencyTracking::Variable;
      using SystemT   = OPF::SystemModel<RealT, IdxT>;
      using TrackedT  = OPF::SystemModel<VariableT, IdxT>;

    public:
      /// Balance and magnitude rows at every bus, and flow rows at the rated branch only
      TestOutcome allocation()
      {
        TestStatus success = true;

        SystemT system(data());
        success *= system.allocate() == 0;

        // Three buses and two generators with two variables each
        success *= system.size() == 10;
        success *= system.sizeConstraints() == 11;

        const RealT* g_lower = system.gLower().getData();
        const RealT* g_upper = system.gUpper().getData();
        for (const IdxT i : BALANCE_ROWS)
        {
          success *= isEqual(g_lower[i], ZERO<RealT>) && isEqual(g_upper[i], ZERO<RealT>);
        }
        for (const IdxT i : MAGNITUDE_ROWS)
        {
          success *= isEqual(g_lower[i], VMIN * VMIN) && isEqual(g_upper[i], VMAX * VMAX);
        }
        for (const IdxT i : FLOW_ROWS)
        {
          success *= std::isinf(g_lower[i]) && isEqual(g_upper[i], SMAX * SMAX);
        }

        // The imaginary voltage of reference bus 1 is fixed at zero
        success *= isEqual(system.xLower().getData()[1], ZERO<RealT>) && isEqual(system.xUpper().getData()[1], ZERO<RealT>);

        return success.report(__func__);
      }

      /// The reference bus keeps zero imaginary voltage, and the others move
      TestOutcome reference()
      {
        TestStatus success = true;

        SystemT system(data());
        system.allocate();
        success *= optimize(system);

        const RealT* x = system.x().getData();

        success *= isEqual(x[1], ZERO<RealT>, BALANCE_TOL);
        success *= std::abs(x[3]) > BALANCE_TOL && std::abs(x[5]) > BALANCE_TOL;

        return success.report(__func__);
      }

      /// Enzyme system Jacobian equals the dependency-tracking Jacobian
      TestOutcome jacobian()
      {
        TestStatus success = true;

        SystemT system(data());
        system.allocate();
        setVariables(system);
        system.evaluateJacobian();
        const auto enzyme = MapFromCsr(system.getCsrJacobian());

        TrackedT tracked(data());
        tracked.allocate();
        setVariables(tracked);
        tracked.evaluateConstraints();
        tracked.evaluateJacobian();
        const auto dependencies = MapFromCsr(tracked.getCsrJacobian());

        success *= enzyme.size() == dependencies.size();
        for (IdxT i = 0; i < enzyme.size(); ++i)
        {
          success *= isEqual(enzyme[i], dependencies[i], DERIVATIVE_TOL);
        }

        return success.report(__func__);
      }

      /// Enzyme objective gradient equals the dependency-tracking gradient
      TestOutcome gradient()
      {
        TestStatus success = true;

        SystemT system(data());
        system.allocate();
        setVariables(system);
        system.evaluateGradient();

        TrackedT tracked(data());
        tracked.allocate();
        setVariables(tracked);
        tracked.evaluateObjective();
        tracked.evaluateGradient();

        const RealT* enzyme       = system.gradient().getData();
        const RealT* dependencies = tracked.gradient().getData();
        for (IdxT j = 0; j < system.size(); ++j)
        {
          success *= isEqual(enzyme[j], dependencies[j], DERIVATIVE_TOL);
        }

        return success.report(__func__);
      }

      /// The Hessian keeps its pattern and is linear in the factors
      TestOutcome hessian()
      {
        TestStatus success = true;

        SystemT system(data());
        system.allocate();
        setVariables(system);

        const IdxT               nnz = system.getCsrHessian()->getNnz();
        const std::vector<RealT> zeros(system.sizeConstraints(), ZERO<RealT>);
        std::vector<RealT>       lambda(system.sizeConstraints());
        for (IdxT i = 0; i < lambda.size(); ++i)
        {
          lambda[i] = LAMBDA_SCALE * static_cast<RealT>(i + 1);
        }

        const std::vector<RealT> objective   = hessianValues(system, SIGMA, zeros);
        const std::vector<RealT> constraints = hessianValues(system, ZERO<RealT>, lambda);
        const std::vector<RealT> lagrangian  = hessianValues(system, SIGMA, lambda);

        success *= system.getCsrHessian()->getNnz() == nnz;
        success *= objective.size() == nnz && constraints.size() == nnz && lagrangian.size() == nnz;
        for (IdxT k = 0; k < nnz; ++k)
        {
          success *= isEqual(lagrangian[k], objective[k] + constraints[k], DERIVATIVE_TOL);
        }

        return success.report(__func__);
      }

      /// Ipopt meets the balance with the flow limit binding
      TestOutcome solve()
      {
        TestStatus success = true;

        SystemT system(data());
        system.allocate();
        success *= optimize(system);

        system.evaluateConstraints();
        const RealT* g = system.g().getData();
        for (const IdxT i : BALANCE_ROWS)
        {
          success *= isEqual(g[i], ZERO<RealT>, BALANCE_TOL);
        }
        success *= isEqual(g[FLOW_ROWS[0]], SMAX * SMAX, BALANCE_TOL);

        const RealT* x = system.x().getData();

        success *= x[6] > x[8];
        success *= x[8] > ZERO<RealT>;

        return success.report(__func__);
      }

      /// Data with an island without a reference bus, or with a repeated bus, is rejected
      TestOutcome configuration()
      {
        TestStatus success = true;

        auto no_reference             = data();
        no_reference.bus[0].reference = false;

        success *= configurationRejected([&]()
                                         { return SystemT(no_reference).verify(); });

        // Out-of-service branches 2-3 and 1-3 leave bus 3 an island without a reference bus
        auto island                 = data();
        island.branch[1].in_service = false;
        island.branch[2].in_service = false;

        success *= configurationRejected([&]()
                                         { return SystemT(island).verify(); });

        auto repeated          = data();
        repeated.bus[1].number = 1;

        success *= throws<std::invalid_argument>([&]()
                                                 { SystemT system(repeated); });

        return success.report(__func__);
      }

      /**
       * @brief MATPOWER rows become components on the system base
       *
       * Branch 2 and generator 2 are out of service, and branch 3 is a
       * phase-shifting transformer with line charging.
       */
      TestOutcome parseMatpowerCase()
      {
        TestStatus success = true;

        const auto data = parse(matpowerCase(MATPOWER_BRANCH, MATPOWER_GENCOST));

        success *= isEqual(data.va_base, 100.0e6);
        success *= data.bus.size() == 3;
        success *= data.bus[0].reference && !data.bus[1].reference && !data.bus[2].reference;
        success *= isEqual(data.bus[1].parameters.at(OPF::BusParameters::Vmin), 0.95);
        success *= isEqual(data.bus[1].parameters.at(OPF::BusParameters::Vmax), 1.05);

        success *= data.load.size() == 1 && data.load[0].id == "load_2";
        success *= isEqual(data.load[0].parameters.at(OPF::LoadParameters::P), 1.2);
        success *= isEqual(data.load[0].parameters.at(OPF::LoadParameters::Q), 0.3);
        success *= data.shunt.size() == 1 && data.shunt[0].id == "shunt_3";
        success *= isEqual(data.shunt[0].parameters.at(OPF::ShuntParameters::B), 0.1);

        success *= data.branch.size() == 3 && data.branch[2].id == "branch_3";
        success *= data.branch[0].in_service && !data.branch[1].in_service && data.branch[2].in_service;

        const auto& line = data.branch[0].parameters;

        success *= isEqual(line.at(OPF::BranchParameters::tap), ONE<RealT>);
        success *= isEqual(line.at(OPF::BranchParameters::Bmag), ZERO<RealT>);
        success *= isEqual(line.at(OPF::BranchParameters::Smax), 0.6);

        const auto& transformer = data.branch[2].parameters;

        success *= isEqual(transformer.at(OPF::BranchParameters::tap), 1.05);
        success *= isEqual(transformer.at(OPF::BranchParameters::phase), std::numbers::pi / 18.0);
        success *= isEqual(transformer.at(OPF::BranchParameters::Bmag), 0.01 / (1.05 * 1.05) - 0.01);
        success *= !transformer.contains(OPF::BranchParameters::Smax);

        success *= data.generator.size() == 3 && data.generator[2].id == "gen_3";
        success *= data.generator[0].in_service && !data.generator[1].in_service && data.generator[2].in_service;

        const auto& gen_1 = data.generator[0].parameters;

        success *= isEqual(gen_1.at(OPF::GeneratorParameters::Pmin), 0.1);
        success *= isEqual(gen_1.at(OPF::GeneratorParameters::Pmax), 1.5);
        success *= isEqual(gen_1.at(OPF::GeneratorParameters::Qmin), -0.5);
        success *= isEqual(gen_1.at(OPF::GeneratorParameters::Qmax), 0.5);
        success *= isEqual(gen_1.at(OPF::GeneratorParameters::c0), 5.0);
        success *= isEqual(gen_1.at(OPF::GeneratorParameters::c1), 1200.0);
        success *= isEqual(gen_1.at(OPF::GeneratorParameters::c2), 100.0);

        const auto& gen_3 = data.generator[2].parameters;

        success *= isEqual(gen_3.at(OPF::GeneratorParameters::c0), 7.0);
        success *= isEqual(gen_3.at(OPF::GeneratorParameters::c1), 1400.0);
        success *= isEqual(gen_3.at(OPF::GeneratorParameters::c2), 200.0);

        // An isolated bus, angle-difference limits, a piecewise-linear cost, a cubic cost, and a missing cost
        const std::string isolated     = MATPOWER_BUS + "\n4 4 0 0 0 0 1 1 0 230 1 1.04 0.96";
        const std::string angle_limits = MATPOWER_BRANCH + "\n1 3 0.02 0.2 0 0 0 0 0 0 1 -30 30";
        const std::string piecewise    = "1 0 0 2 0 0 100 1000;\n2 0 0 2 20 0;\n2 0 0 2 14 7";
        const std::string cubic        = "2 0 0 4 1 0.02 14 7;\n2 0 0 2 20 0;\n2 0 0 2 14 7";
        const std::string missing      = "2 0 0 3 0.01 12 5;\n2 0 0 2 20 0";

        success *= throws<std::invalid_argument>(parse, matpowerCase(MATPOWER_BRANCH, MATPOWER_GENCOST, isolated));
        success *= throws<std::invalid_argument>(parse, matpowerCase(angle_limits, MATPOWER_GENCOST));
        success *= throws<std::invalid_argument>(parse, matpowerCase(MATPOWER_BRANCH, piecewise));
        success *= throws<std::invalid_argument>(parse, matpowerCase(MATPOWER_BRANCH, cubic));
        success *= throws<std::invalid_argument>(parse, matpowerCase(MATPOWER_BRANCH, missing));

        return success.report(__func__);
      }

    private:
      static OPF::SystemModelData<RealT, IdxT> parse(const std::string& text)
      {
        std::istringstream input(text);
        return OPF::parseMatpowerCase(input);
      }

      /// MATPOWER case with the given branch, cost, and bus rows
      static std::string matpowerCase(const std::string& branch, const std::string& gencost, const std::string& bus = MATPOWER_BUS)
      {
        return "mpc.version = '2';\nmpc.baseMVA = 100;\nmpc.bus = [\n" + bus + "\n];\nmpc.gen = [\n" + MATPOWER_GEN
               + "\n];\nmpc.branch = [\n" + branch + "\n];\nmpc.gencost = [\n" + gencost + "\n];\n";
      }

      static inline const std::string MATPOWER_BUS = "1 3 0 0 0 0 1 1 0 230 1 1.06 0.94;\n"
                                                     "2 1 120 30 0 0 1 1 0 230 1 1.05 0.95;\n"
                                                     "3 2 0 0 0 10 1 1 0 230 1 1.04 0.96;";

      static inline const std::string MATPOWER_GEN = "1 0 0 50 -50 1 100 1 150 10\n"
                                                     "3 0 0 80 -40 1 100 0 90 0\n"
                                                     "3 0 0 80 -40 1 100 1 180 30";

      static inline const std::string MATPOWER_BRANCH = "1 2 0.01 0.1 0 60 0 0 0 0 1 -360 360\n"
                                                        "2 3 0.01 0.1 0 90 0 0 0 0 0 -360 360\n"
                                                        "1 3 0.02 0.2 0.02 0 0 0 1.05 10 1 0 0";

      static inline const std::string MATPOWER_GENCOST = "2 0.0 0.0 3 0.01 12 5;\n"
                                                         "2 0 0 2 20 0;\n"
                                                         "2 0 0 4 0 0.02 14 7";

      /// Values within `tol` relative to one plus the reference magnitude
      static bool isEqual(RealT value, RealT reference, RealT tol = std::numeric_limits<RealT>::epsilon())
      {
        if (Testing::isEqual(value, reference, tol))
        {
          return true;
        }
        std::cerr << std::setprecision(std::numeric_limits<RealT>::max_digits10) << value << " differs from "
                  << reference << " by " << std::abs(value - reference) / (ONE<RealT> + std::abs(reference)) / std::numeric_limits<RealT>::epsilon()
                  << " eps\n";
        return false;
      }

      /// Dependency maps with the same keys and values within `tol`
      static bool isEqual(const DependencyTracking::Variable::DependencyMap& actual,
                          const DependencyTracking::Variable::DependencyMap& expected,
                          RealT                                              tol)
      {
        bool equal = actual.size() == expected.size();
        for (const auto& [variable, value] : expected)
        {
          const auto entry = actual.find(variable);
          if (entry == actual.end())
          {
            std::cerr << "Derivative with respect to " << variable << " is missing\n";
            equal = false;
          }
          else
          {
            equal = isEqual(entry->second, value, tol) && equal;
          }
        }
        return equal;
      }

      static constexpr RealT SMAX         = 0.6;
      static constexpr RealT VMIN         = 0.95;
      static constexpr RealT VMAX         = 1.05;
      static constexpr RealT LOAD_P       = 1.2;
      static constexpr RealT LOAD_Q       = 0.3;
      static constexpr RealT SIGMA        = 0.7;
      static constexpr RealT LAMBDA_SCALE = 0.25;

      static constexpr RealT DERIVATIVE_TOL = 8 * std::numeric_limits<RealT>::epsilon();
      static constexpr RealT BALANCE_TOL    = 1.0e-12;

      /// Rows of the buses in order, then of the rated branch
      static inline const std::vector<IdxT> BALANCE_ROWS   = {0, 1, 3, 4, 6, 7};
      static inline const std::vector<IdxT> MAGNITUDE_ROWS = {2, 5, 8};
      static inline const std::vector<IdxT> FLOW_ROWS      = {9, 10};

      /// A point away from the flat start
      static inline const std::vector<RealT> POINT = {1.02, 0.0, 0.97, -0.07, 1.01, -0.03, 0.8, 0.2, 0.4, -0.1};

      static OPF::SystemModelData<RealT, IdxT> data()
      {
        OPF::SystemModelData<RealT, IdxT> data;

        for (IdxT number = 1; number <= 3; ++number)
        {
          auto& bus                                = data.bus.emplace_back();
          bus.number                               = number;
          bus.reference                            = number == 1;
          bus.parameters[OPF::BusParameters::Vmin] = VMIN;
          bus.parameters[OPF::BusParameters::Vmax] = VMAX;
        }

        addBranch(data, "branch_1_2", 1, 2, 0.01, 0.1, 0.02);
        addBranch(data, "branch_2_3", 2, 3, 0.01, 0.1, 0.0);
        addBranch(data, "branch_1_3", 1, 3, 0.02, 0.2, 0.0);
        data.branch[0].parameters[OPF::BranchParameters::Smax] = SMAX;

        addGenerator(data, "gen_1", 1, 10.0);
        addGenerator(data, "gen_3", 3, 20.0);

        auto& load                              = data.load.emplace_back();
        load.id                                 = "load_2";
        load.buses[OPF::LoadBuses::bus]         = 2;
        load.parameters[OPF::LoadParameters::P] = LOAD_P;
        load.parameters[OPF::LoadParameters::Q] = LOAD_Q;

        auto& shunt                               = data.shunt.emplace_back();
        shunt.id                                  = "shunt_3";
        shunt.buses[OPF::ShuntBuses::bus]         = 3;
        shunt.parameters[OPF::ShuntParameters::B] = 0.1;

        return data;
      }

      static void addBranch(OPF::SystemModelData<RealT, IdxT>& data,
                            const char*                        id,
                            IdxT                               bus1,
                            IdxT                               bus2,
                            RealT                              r,
                            RealT                              x,
                            RealT                              b)
      {
        auto& branch                                = data.branch.emplace_back();
        branch.id                                   = id;
        branch.buses[OPF::BranchBuses::bus1]        = bus1;
        branch.buses[OPF::BranchBuses::bus2]        = bus2;
        branch.parameters[OPF::BranchParameters::R] = r;
        branch.parameters[OPF::BranchParameters::X] = x;
        branch.parameters[OPF::BranchParameters::B] = b;
      }

      static void addGenerator(OPF::SystemModelData<RealT, IdxT>& data, const char* id, IdxT bus, RealT c1)
      {
        auto& generator                                      = data.generator.emplace_back();
        generator.id                                         = id;
        generator.buses[OPF::GeneratorBuses::bus]            = bus;
        generator.parameters[OPF::GeneratorParameters::Pmin] = 0.0;
        generator.parameters[OPF::GeneratorParameters::Pmax] = 2.0;
        generator.parameters[OPF::GeneratorParameters::Qmin] = -1.0;
        generator.parameters[OPF::GeneratorParameters::Qmax] = 1.0;
        generator.parameters[OPF::GeneratorParameters::c1]   = c1;
        generator.parameters[OPF::GeneratorParameters::c2]   = 0.01;
      }

      /// Move the variables to `POINT`
      template <typename ModelT>
      static void setVariables(ModelT& system)
      {
        auto* x = system.x().getData();
        for (IdxT j = 0; j < POINT.size(); ++j)
        {
          x[j] = POINT[j];
          if constexpr (std::is_same_v<typename ModelT::ScalarT, VariableT>)
          {
            x[j].setVariableNumber(j);
          }
        }
        system.x().setDataUpdated();
      }

      /// Hessian values at `sigma` and `lambda`, or none if the evaluation fails
      static std::vector<RealT> hessianValues(SystemT& system, RealT sigma, const std::vector<RealT>& lambda)
      {
        if (system.evaluateHessian(sigma, lambda.data()) != 0)
        {
          return {};
        }
        auto* hessian = system.getCsrHessian();
        return std::vector<RealT>(hessian->getValues(), hessian->getValues() + hessian->getNnz());
      }

      static bool optimize(SystemT& system)
      {
        Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
        app->Options()->SetIntegerValue("print_level", 0);
        app->Options()->SetNumericValue("bound_relax_factor", 0.0);
        app->Options()->SetNumericValue("tol", 1.0e-10);
        if (app->Initialize() != Ipopt::Solve_Succeeded)
        {
          return false;
        }

        Ipopt::SmartPtr<Ipopt::TNLP> problem = new AnalysisManager::IpoptInterface::SteadyState<RealT, IdxT>(&system);
        return app->OptimizeTNLP(problem) == Ipopt::Solve_Succeeded;
      }
    };
  } // namespace Testing
} // namespace GridKit
