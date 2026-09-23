#pragma once

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <IpIpoptApplication.hpp>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/OptimalPowerFlow/SystemModel.hpp>
#include <GridKit/Model/OptimalPowerFlow/SystemModelData.hpp>
#include <GridKit/Model/StateData.hpp>
#include <GridKit/Solver/Optimization/OptimizationProblem.hpp>
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
      /// Balance and magnitude rows at every bus, one reference row, and flow rows at the rated branch only
      TestOutcome allocation()
      {
        TestStatus success = true;

        SystemT system(data(), state());
        success *= system.allocate() == 0;

        // Three buses and two generators with two variables each
        success *= system.size() == 10;
        success *= system.sizeConstraints() == 12;

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
        success *= isEqual(g_lower[REFERENCE_ROW], ZERO<RealT>) && isEqual(g_upper[REFERENCE_ROW], ZERO<RealT>);
        for (const IdxT i : FLOW_ROWS)
        {
          success *= std::isinf(g_lower[i]) && isEqual(g_upper[i], SMAX * SMAX);
        }

        return success.report(__func__);
      }

      /// The lowest-numbered bus keeps its state angle, and the others move
      TestOutcome reference()
      {
        TestStatus success = true;

        SystemT system(data(), state());
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

        SystemT system(data(), state());
        system.allocate();
        setVariables(system);
        system.evaluateJacobian();
        const auto enzyme = MapFromCsr(system.getCsrJacobian());

        TrackedT tracked(data(), state());
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

        SystemT system(data(), state());
        system.allocate();
        setVariables(system);
        system.evaluateGradient();

        TrackedT tracked(data(), state());
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

        SystemT system(data(), state());
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

        SystemT system(data(), state());
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

      /// The solution state restarts the model at the solution
      TestOutcome solutionState()
      {
        TestStatus success = true;

        SystemT system(data(), state());
        system.allocate();
        optimize(system);

        SystemT restarted(data(), system.solutionState());
        success *= restarted.allocate() == 0;

        const RealT* x         = system.x().getData();
        const RealT* x_restart = restarted.x().getData();
        for (IdxT j = 0; j < system.size(); ++j)
        {
          success *= isEqual(x_restart[j], x[j], STATE_TOL);
        }

        return success.report(__func__);
      }

      /// Rows end at `;` or the line end, and scalars, cell arrays, and comments are skipped
      TestOutcome parseMatpowerData()
      {
        TestStatus success = true;

        const OPF::MatpowerData matpower = parse(R"(function mpc = test
          mpc.baseMVA = 100.0;
          mpc.bus_name = {
            'ONE';
          };
          mpc.gen = [ % 1 2 3;
            1 2.5	3;
            4 5 6
            7 8 9];
          mpc.gencost = [2 0 0 2 1 0];
        )");

        success *= matpower.matrices.size() == 2;
        success *= matpower.matrix("gen") == OPF::MatpowerMatrix{{1.0, 2.5, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
        success *= matpower.matrix("gencost") == OPF::MatpowerMatrix{{2.0, 0.0, 0.0, 2.0, 1.0, 0.0}};

        success *= throws<std::invalid_argument>([&]
                                                 { matpower.matrix("branch"); });
        success *= throws<std::invalid_argument>(parse, "mpc.bus = [\n 1 Inf;\n];");
        success *= throws<std::invalid_argument>(parse, "mpc.bus = [\n 1 2;\n");

        return success.report(__func__);
      }

      /**
       * @brief Buses match by number, in-service branches by their buses in
       * order, and in-service generators in order at buses with generators
       *
       * Branch 2-3 is listed as 3-2, and branch 1-3 has no rating. The
       * generator at bus 2 is a static injection, and the first one at bus 3
       * is out of service.
       */
      TestOutcome applyMatpowerData()
      {
        TestStatus success = true;

        OPF::SystemModelData<RealT, IdxT> network = data();
        network.branch[0].parameters.erase(OPF::BranchParameters::Smax);
        OPF::applyMatpowerData(network, parse(matpowerCase(MATPOWER_BRANCH, MATPOWER_GEN, MATPOWER_GENCOST)));

        success *= isEqual(network.bus[0].parameters.at(OPF::BusParameters::Vmin), 0.94);
        success *= isEqual(network.bus[0].parameters.at(OPF::BusParameters::Vmax), 1.06);
        success *= isEqual(network.bus[1].parameters.at(OPF::BusParameters::Vmin), 0.95);

        success *= isEqual(network.branch[0].parameters.at(OPF::BranchParameters::Smax), 0.6);
        success *= isEqual(network.branch[1].parameters.at(OPF::BranchParameters::Smax), 0.9);
        success *= !network.branch[2].parameters.contains(OPF::BranchParameters::Smax);

        const auto& gen_1 = network.generator[0].parameters;

        success *= isEqual(gen_1.at(OPF::GeneratorParameters::Pmin), 0.1);
        success *= isEqual(gen_1.at(OPF::GeneratorParameters::Pmax), 1.5);
        success *= isEqual(gen_1.at(OPF::GeneratorParameters::Qmin), -0.5);
        success *= isEqual(gen_1.at(OPF::GeneratorParameters::Qmax), 0.5);
        success *= isEqual(gen_1.at(OPF::GeneratorParameters::c0), 5.0);
        success *= isEqual(gen_1.at(OPF::GeneratorParameters::c1), 1200.0);
        success *= isEqual(gen_1.at(OPF::GeneratorParameters::c2), 100.0);

        const auto& gen_3 = network.generator[1].parameters;

        success *= isEqual(gen_3.at(OPF::GeneratorParameters::Pmin), 0.3);
        success *= isEqual(gen_3.at(OPF::GeneratorParameters::Pmax), 1.8);
        success *= isEqual(gen_3.at(OPF::GeneratorParameters::c0), 7.0);
        success *= isEqual(gen_3.at(OPF::GeneratorParameters::c1), 1400.0);
        success *= isEqual(gen_3.at(OPF::GeneratorParameters::c2), 200.0);

        // Extra and missing branches and generators, a piecewise-linear cost, and a cubic cost
        const std::string branch_extra = MATPOWER_BRANCH + "\n1 2 0.01 0.1 0 0 0 0 0 0 1";
        const std::string branch_short = "1 2 0.01 0.1 0.02 60 0 0 0 0 1";
        const std::string gen_extra    = "1 0 0 50 -50 1 100 1 150 10\n3 0 0 80 -40 1 100 1 90 0\n3 0 0 80 -40 1 100 1 180 30";
        const std::string gen_short    = "1 0 0 50 -50 1 100 1 150 10";
        const std::string costs        = "2 0 0 3 0.01 12 5;\n2 0 0 2 20 0;\n2 0 0 2 20 0";
        const std::string piecewise    = "1 0 0 2 0 0 100 10;\n2 0 0 2 20 0;\n2 0 0 2 20 0;\n2 0 0 2 14 7";
        const std::string cubic        = "2 0 0 3 0.01 12 5;\n2 0 0 2 20 0;\n2 0 0 2 20 0;\n2 0 0 4 1 0 14 7";

        success *= throws<std::invalid_argument>(apply, matpowerCase(branch_extra, MATPOWER_GEN, MATPOWER_GENCOST));
        success *= throws<std::invalid_argument>(apply, matpowerCase(branch_short, MATPOWER_GEN, MATPOWER_GENCOST));
        success *= throws<std::invalid_argument>(apply, matpowerCase(MATPOWER_BRANCH, gen_extra, costs));
        success *= throws<std::invalid_argument>(apply, matpowerCase(MATPOWER_BRANCH, gen_short, costs));
        success *= throws<std::invalid_argument>(apply, matpowerCase(MATPOWER_BRANCH, MATPOWER_GEN, piecewise));
        success *= throws<std::invalid_argument>(apply, matpowerCase(MATPOWER_BRANCH, MATPOWER_GEN, cubic));

        return success.report(__func__);
      }

    private:
      static OPF::MatpowerData parse(const std::string& text)
      {
        std::istringstream input(text);
        return OPF::parseMatpowerData(input);
      }

      static void apply(const std::string& text)
      {
        OPF::SystemModelData<RealT, IdxT> network = data();
        OPF::applyMatpowerData(network, parse(text));
      }

      /// MATPOWER case of the `data()` network with the given branch, generator, and cost rows
      static std::string matpowerCase(const std::string& branch, const std::string& gen, const std::string& gencost)
      {
        return "mpc.bus = [\n"
               "1 3 0 0 0 0 1 1 0 230 1 1.06 0.94;\n"
               "2 1 120 30 0 0 1 1 0 230 1 1.05 0.95;\n"
               "3 2 0 0 0 0 1 1 0 230 1 1.04 0.96;\n"
               "];\n"
               "mpc.branch = [\n"
               + branch + "\n];\nmpc.gen = [\n" + gen + "\n];\nmpc.gencost = [\n" + gencost + "\n];\n";
      }

      static inline const std::string MATPOWER_BRANCH = "1 2 0.01 0.1 0.02 60 0 0 0 0 1\n"
                                                        "3 2 0.01 0.1 0 90 0 0 0 0 1\n"
                                                        "1 3 0.02 0.2 0 0 0 0 0 0 1\n"
                                                        "1 3 0.02 0.2 0 50 0 0 0 0 0";

      static inline const std::string MATPOWER_GEN = "1 0 0 50 -50 1 100 1 150 10\n"
                                                     "2 20 0 0 0 1 100 1 20 20\n"
                                                     "3 0 0 80 -40 1 100 0 90 0\n"
                                                     "3 0 0 80 -40 1 100 1 180 30";

      static inline const std::string MATPOWER_GENCOST = "2 0 0 3 0.01 12 5;\n"
                                                         "2 0 0 2 20 0;\n"
                                                         "2 0 0 3 0 0 0;\n"
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

      // Measured worst errors: Jacobian 7.0 eps, gradient and Hessian exact,
      // solved balance 5.4e-13, restart 0.3 eps
      static constexpr RealT DERIVATIVE_TOL = 8 * std::numeric_limits<RealT>::epsilon();
      static constexpr RealT BALANCE_TOL    = 1.0e-12;
      static constexpr RealT STATE_TOL      = std::numeric_limits<RealT>::epsilon();

      /// Rows of the buses in order, then of the rated branch
      static inline const std::vector<IdxT> BALANCE_ROWS   = {0, 1, 4, 5, 7, 8};
      static inline const std::vector<IdxT> MAGNITUDE_ROWS = {2, 6, 9};
      static constexpr IdxT                 REFERENCE_ROW  = 3;
      static inline const std::vector<IdxT> FLOW_ROWS      = {10, 11};

      /// A point away from the flat start
      static inline const std::vector<RealT> POINT = {1.02, 0.0, 0.97, -0.07, 1.01, -0.03, 0.8, 0.2, 0.4, -0.1};

      static OPF::SystemModelData<RealT, IdxT> data()
      {
        OPF::SystemModelData<RealT, IdxT> data;

        for (IdxT number = 1; number <= 3; ++number)
        {
          auto& bus                                = data.bus.emplace_back();
          bus.number                               = number;
          bus.parameters[OPF::BusParameters::Vmin] = VMIN;
          bus.parameters[OPF::BusParameters::Vmax] = VMAX;
        }

        addBranch(data, "branch_1_2", 1, 2, 0.01, 0.1, 0.02);
        addBranch(data, "branch_2_3", 2, 3, 0.01, 0.1, 0.0);
        addBranch(data, "branch_1_3", 1, 3, 0.02, 0.2, 0.0);
        data.branch[0].parameters[OPF::BranchParameters::Smax] = SMAX;

        addGenerator(data, "gen_1", 1, 10.0);
        addGenerator(data, "gen_3", 3, 20.0);

        auto& load                      = data.load.emplace_back();
        load.id                         = "load_2";
        load.buses[OPF::LoadBuses::bus] = 2;

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

      /// Flat start with the load demand
      static Model::StateData state()
      {
        Model::StateData state;
        for (IdxT number = 1; number <= 3; ++number)
        {
          state.buses[Model::busKey(number)].values = {{"vr", ONE<RealT>}, {"vi", ZERO<RealT>}};
        }
        Model::setTerminalCurrent(state, "load_2", 2, 0, 1, -LOAD_P, -LOAD_Q);
        return state;
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

        Ipopt::SmartPtr<Ipopt::TNLP> problem = new AnalysisManager::IpoptInterface::OptimizationProblem<RealT, IdxT>(&system);
        return app->OptimizeTNLP(problem) == Ipopt::Solve_Succeeded;
      }
    };
  } // namespace Testing
} // namespace GridKit
