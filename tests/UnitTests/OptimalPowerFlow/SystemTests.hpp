#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/OptimalPowerFlow/SystemModel.hpp>
#include <GridKit/Model/OptimalPowerFlow/SystemModelData.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

#include "CoupledComponent.hpp"
#include "SystemTestData.hpp"
#include "TestHelpers.hpp"

namespace GridKit
{
  namespace Testing
  {
    namespace OPF = GridKit::OptimalPowerFlow;

    /**
     * @brief System assembly and derivatives on a three-bus
     * network
     *
     * Generator 1 is cheaper, but the rating of branch 1-2 forces generator
     * 3 to cover part of the load at bus 2.
     */
    class SystemTests
    {
      using RealT         = double;
      using IdxT          = size_t;
      using VariableT     = DependencyTracking::Variable;
      using DependencyMap = VariableT::DependencyMap;
      using SystemT       = OPF::SystemModel<RealT, IdxT>;
      using TrackedT      = OPF::SystemModel<VariableT, IdxT>;

    public:
      // The Hessian linearity sums reach 17 eps.
      static constexpr RealT kTol =
          static_cast<RealT>(32.0) * std::numeric_limits<RealT>::epsilon();

      /// Balance equalities and magnitude inequalities at every bus, and flow
      /// inequalities at every branch, unlimited at the unrated ones
      TestOutcome allocation()
      {
        TestStatus success = true;

        SystemT system(SystemTestData::data());
        success *= system.allocate() == 0 && system.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }

        // Three buses and two generators with two variables each
        success *= system.size() == 10;
        success *= system.sizeEqualities() == SystemTestData::BALANCE_ROWS.size();
        success *= system.sizeInequalities() == SystemTestData::MAGNITUDE_ROWS.size() + SystemTestData::FLOW_ROWS.size() + SystemTestData::UNRATED_ROWS.size();

        const RealT* d_lower = system.dLower().getData();
        const RealT* d_upper = system.dUpper().getData();
        for (const IdxT i : SystemTestData::MAGNITUDE_ROWS)
        {
          success *= isEqual(d_lower[i], SystemTestData::VMIN * SystemTestData::VMIN) && isEqual(d_upper[i], SystemTestData::VMAX * SystemTestData::VMAX);
        }
        for (const IdxT i : SystemTestData::FLOW_ROWS)
        {
          success *= std::isinf(d_lower[i]) && isEqual(d_upper[i], SystemTestData::SMAX * SystemTestData::SMAX);
        }
        for (const IdxT i : SystemTestData::UNRATED_ROWS)
        {
          success *= std::isinf(d_lower[i]) && std::isinf(d_upper[i]);
        }

        // The imaginary voltage of reference bus 1 is fixed at zero
        success *= isEqual(system.xLower().getData()[1], ZERO<RealT>) && isEqual(system.xUpper().getData()[1], ZERO<RealT>);

        return success.report(__func__);
      }

      /// Enzyme system Jacobians equal the dependency-tracking Jacobians
      TestOutcome jacobian()
      {
        TestStatus success = true;

        SystemT system(SystemTestData::data());
        success *= system.allocate() == 0 && system.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        setVariables(system, SystemTestData::POINT);
        success *= system.evaluateEqualityJacobian() == 0;
        success *= system.evaluateInequalityJacobian() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        const auto enzyme_c = MapFromCsr(system.getCsrEqualityJacobian());
        const auto enzyme_d = MapFromCsr(system.getCsrInequalityJacobian());

        TrackedT tracked(SystemTestData::data());
        success *= tracked.allocate() == 0 && tracked.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        success *= tracked.getCsrEqualityJacobian() != nullptr;
        success *= tracked.getCsrInequalityJacobian() != nullptr;
        setVariables(tracked, SystemTestData::POINT);
        success *= tracked.evaluateEqualityJacobian() == 0;
        success *= tracked.evaluateInequalityJacobian() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        const auto dependencies_c = MapFromCsr(tracked.getCsrEqualityJacobian());
        const auto dependencies_d = MapFromCsr(tracked.getCsrInequalityJacobian());

        success *= enzyme_c.size() == dependencies_c.size() && enzyme_d.size() == dependencies_d.size();
        if (!success)
        {
          return success.report(__func__);
        }
        for (IdxT i = 0; i < enzyme_c.size(); ++i)
        {
          success *= isEqual(enzyme_c[i], dependencies_c[i], kTol);
        }
        for (IdxT i = 0; i < enzyme_d.size(); ++i)
        {
          success *= isEqual(enzyme_d[i], dependencies_d[i], kTol);
        }

        return success.report(__func__);
      }

      /// Enzyme objective gradient equals the dependency-tracking gradient
      TestOutcome gradient()
      {
        TestStatus success = true;

        SystemT system(SystemTestData::data());
        success *= system.allocate() == 0 && system.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        setVariables(system, SystemTestData::POINT);
        success *= system.evaluateGradient() == 0;
        if (!success)
        {
          return success.report(__func__);
        }

        TrackedT tracked(SystemTestData::data());
        success *= tracked.allocate() == 0 && tracked.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        setVariables(tracked, SystemTestData::POINT);
        success *= tracked.evaluateGradient() == 0;
        if (!success)
        {
          return success.report(__func__);
        }

        const RealT* enzyme       = system.gradient().getData();
        const RealT* dependencies = tracked.gradient().getData();
        for (IdxT j = 0; j < system.size(); ++j)
        {
          success *= isEqual(enzyme[j], dependencies[j], kTol);
        }

        return success.report(__func__);
      }

      /// The Hessian keeps its pattern and is linear in the factors
      TestOutcome hessianLinearity()
      {
        TestStatus success = true;

        SystemT system(SystemTestData::data());
        success *= system.allocate() == 0 && system.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        setVariables(system, SystemTestData::POINT);

        const auto               expected_pattern = pattern(*system.getCsrHessian());
        const IdxT               nnz              = system.getCsrHessian()->getNnz();
        const std::vector<RealT> zeros_c(system.sizeEqualities(), ZERO<RealT>);
        const std::vector<RealT> zeros_d(system.sizeInequalities(), ZERO<RealT>);
        std::vector<RealT>       lambda_c(system.sizeEqualities());
        std::vector<RealT>       lambda_d(system.sizeInequalities());
        for (IdxT i = 0; i < lambda_c.size(); ++i)
        {
          lambda_c[i] = LAMBDA_SCALE * static_cast<RealT>(i + 1);
        }
        for (IdxT i = 0; i < lambda_d.size(); ++i)
        {
          lambda_d[i] = LAMBDA_SCALE * static_cast<RealT>(lambda_c.size() + i + 1);
        }

        const std::vector<RealT> objective    = hessianValues(system, SIGMA, zeros_c, zeros_d);
        success                              *= pattern(*system.getCsrHessian()) == expected_pattern;
        const std::vector<RealT> constraints  = hessianValues(system, ZERO<RealT>, lambda_c, lambda_d);
        success                              *= pattern(*system.getCsrHessian()) == expected_pattern;
        const std::vector<RealT> lagrangian   = hessianValues(system, SIGMA, lambda_c, lambda_d);
        success                              *= pattern(*system.getCsrHessian()) == expected_pattern;

        success *= system.getCsrHessian()->getNnz() == nnz;
        success *= objective.size() == nnz && constraints.size() == nnz && lagrangian.size() == nnz;
        if (!success)
        {
          return success.report(__func__);
        }
        for (IdxT k = 0; k < nnz; ++k)
        {
          success *= isEqual(lagrangian[k], objective[k] + constraints[k], kTol);
        }

        return success.report(__func__);
      }

      /// Inequality constraint functions without limits keep their rows, with
      /// infinite bounds
      TestOutcome unboundedInequalities()
      {
        TestStatus success = true;

        auto unbounded = SystemTestData::data();
        for (auto& bus : unbounded.bus)
        {
          bus.parameters.erase(OPF::BusParameters::Vmin);
          bus.parameters.erase(OPF::BusParameters::Vmax);
        }
        for (auto& branch : unbounded.branch)
        {
          branch.parameters.erase(OPF::BranchParameters::Smax);
        }

        SystemT system(unbounded);
        success *= system.allocate() == 0 && system.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        success *= system.sizeInequalities() == SystemTestData::MAGNITUDE_ROWS.size() + SystemTestData::FLOW_ROWS.size() + SystemTestData::UNRATED_ROWS.size();

        const RealT* d_lower = system.dLower().getData();
        const RealT* d_upper = system.dUpper().getData();
        for (IdxT i = 0; i < system.sizeInequalities(); ++i)
        {
          success *= std::isinf(d_lower[i]) && std::isinf(d_upper[i]);
        }
        success *= system.evaluateInequalities() == 0;
        success *= system.evaluateInequalityJacobian() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        success *= system.getCsrInequalityJacobian()->getNumRows() == system.sizeInequalities();

        TrackedT tracked(unbounded);
        success *= tracked.allocate() == 0 && tracked.initialize() == 0;
        success *= tracked.evaluateInequalityJacobian() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        success *= tracked.getCsrInequalityJacobian()->getNumRows() == tracked.sizeInequalities();

        return success.report(__func__);
      }

      /// Invalid bus numbers, connections, and islands are rejected
      TestOutcome configuration()
      {
        TestStatus success = true;

        auto no_reference                                             = SystemTestData::data();
        no_reference.bus[0].parameters[OPF::BusParameters::reference] = false;

        success *= configurationRejected([&]()
                                         { return SystemT(no_reference).verify(); });

        // Out-of-service branches 2-3 and 1-3 leave bus 3 an island without a reference bus
        auto island                                                    = SystemTestData::data();
        island.branch[1].parameters[OPF::BranchParameters::in_service] = false;
        island.branch[2].parameters[OPF::BranchParameters::in_service] = false;

        success *= configurationRejected([&]()
                                         { return SystemT(island).verify(); });

        auto repeated          = SystemTestData::data();
        repeated.bus[1].number = 1;

        success *= throws<std::invalid_argument>([&]()
                                                 { SystemT system(repeated); });

        auto invalid           = SystemTestData::data();
        invalid.bus[0].number  = INVALID_INDEX<IdxT>;
        success               *= configurationRejected([&]()
                                         { return SystemT(invalid).verify(); });

        auto missing = SystemTestData::data();
        missing.generator[0].ports.clear();
        success *= configurationRejected([&]()
                                         { return SystemT(missing).verify(); });

        auto duplicate          = SystemTestData::data();
        duplicate.branch[1].id  = duplicate.branch[0].id;
        success                *= configurationRejected([&]()
                                         { return SystemT(duplicate).verify(); });

        auto unknown                                = SystemTestData::data();
        unknown.load[0].ports[OPF::LoadPorts::bus]  = 99;
        success                                    *= configurationRejected([&]()
                                         { return SystemT(unknown).verify(); });

        return success.report(__func__);
      }

      /// Buses and devices from a JSON input file
      TestOutcome parseCase()
      {
        TestStatus success = true;

        const auto data = parse(jsonCase());

        success *= isEqual(data.va_base, 100.0e6);
        success *= data.bus.size() == 2 && std::get<bool>(data.bus[0].parameters.at(OPF::BusParameters::reference));
        success *= !data.bus[1].parameters.contains(OPF::BusParameters::reference);
        success *= isEqual(std::get<RealT>(data.bus[1].parameters.at(OPF::BusParameters::Vmin)), 0.95);
        success *= isEqual(std::get<RealT>(data.bus[1].parameters.at(OPF::BusParameters::B)), 0.1);

        success *= data.branch.size() == 1 && data.branch[0].ports.at(OPF::BranchPorts::bus2) == 2;
        success *= isEqual(std::get<RealT>(data.branch[0].parameters.at(OPF::BranchParameters::Smax)), 0.6);
        success *= !data.branch[0].parameters.contains(OPF::BranchParameters::tap);

        success *= data.generator.size() == 1 && !std::get<bool>(data.generator[0].parameters.at(OPF::GeneratorParameters::in_service));
        success *= data.load.size() == 1;

        // An unknown device class, parameter, and port
        success *= throws<std::invalid_argument>(parse, jsonCase(R"(, {"class": "Transformer", "id": "t", "ports": {}, "params": {}})"));
        success *= throws<std::invalid_argument>(parse, jsonCase(R"(, {"class": "Load", "id": "l", "ports": {"bus": 2}, "params": {"Pd": 1.0}})"));
        success *= throws<std::invalid_argument>(parse, jsonCase(R"(, {"class": "Load", "id": "l", "ports": {"node": 2}, "params": {}})"));

        // Every study starts flat and has no monitors, so init and mon are rejected
        success *= throws<std::invalid_argument>(parse, R"({"params": {"va_base": 1}, "buses": [{"class": "Bus", "number": 1, "params": {}, "init": {"Vr": 1.0}}], "devices": []})");
        success *= throws<std::invalid_argument>(parse, jsonCase(R"(, {"class": "Load", "id": "l", "ports": {"bus": 2}, "params": {}, "init": {}})"));
        success *= throws<std::invalid_argument>(parse, jsonCase(R"(, {"class": "Load", "id": "l", "ports": {"bus": 2}, "params": {}, "mon": []})"));

        // A boolean real and a numeric service status are read, and rejected by the model
        success *= throws<std::invalid_argument>([&]()
                                                 { SystemT system(parse(jsonCase(R"(, {"class": "Load", "id": "l", "ports": {"bus": 2}, "params": {"P": true}})"))); });
        success *= throws<std::invalid_argument>([&]()
                                                 { SystemT system(parse(jsonCase(R"(, {"class": "Load", "id": "l", "ports": {"bus": 2}, "params": {"in_service": 1}})"))); });

        const std::array<std::string, 4> invalid_numbers{
            "-1", "1.5", std::to_string(INVALID_INDEX<IdxT>), "1e30"};
        for (const auto& number : invalid_numbers)
        {
          const auto bus_case = R"({"params": {"va_base": 1}, "buses": [{"class": "Bus", "number": )"
                                + number + R"(, "params": {}}], "devices": []})";
          success *= throws<std::invalid_argument>(parse, bus_case);
          success *= throws<std::invalid_argument>(parse, jsonCase(R"(, {"class": "Load", "id": "l", "ports": {"bus": )" + number + R"(}, "params": {}})"));
        }

        for (const IdxT number : {IdxT{0}, INVALID_INDEX<IdxT> - 1})
        {
          const auto bus_case = R"({"params": {"va_base": 1}, "buses": [{"class": "Bus", "number": )"
                                + std::to_string(number) + R"(, "params": {}}], "devices": []})";
          success *= parse(bus_case).bus[0].number == number;
        }

        return success.report(__func__);
      }

      TestOutcome connections()
      {
        TestStatus success        = true;
        auto       owner_ptr      = std::make_unique<EqualityOwner>();
        auto       component_ptr  = std::make_unique<CoupledComponent>();
        auto&      owner          = *owner_ptr;
        auto&      component      = *component_ptr;
        success                  *= component.connect(0, owner, {1, 1, 0}) == 0;
        success                  *= component.connect(1, owner) == 0;
        success                  *= component.connect(2, owner) == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        success     *= component.connect(3, owner) != 0;
        success     *= component.connect(0, owner, {2, 0, 0}) != 0;
        success     *= component.connect(0, owner, {0, 2, 0}) != 0;
        success     *= component.connect(0, owner, {INVALID_INDEX<IdxT>, 0, 0}) != 0;
        success     *= component.ports()[0].offset.variables == 1;
        success     *= component.ports()[0].offset.equalities == 1;
        auto system  = makeSystem<RealT>(std::move(owner_ptr), std::move(component_ptr));
        success     *= system.allocate() == 0 && system.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        success *= component.connect(1, owner, {1, 0, 0}) != 0;
        return success.report(__func__);
      }

      TestOutcome objectiveAssembly()
      {
        TestStatus success  = true;
        auto       system   = coupledSystem();
        success            *= system.allocate() == 0 && system.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        setVariables(system, COUPLED_X);
        success *= system.evaluateObjective() == 0;
        success *= system.evaluateGradient() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        success               *= isEqual(system.f(), 18.5);
        const RealT* gradient  = system.gradient().getData();
        success               *= isEqual(gradient[0], 3.0);
        success               *= isEqual(gradient[1], 2.0);
        success               *= isEqual(gradient[2], 5.0);
        return success.report(__func__);
      }

      TestOutcome equalityAssembly()
      {
        TestStatus success  = true;
        auto       system   = coupledSystem();
        success            *= system.allocate() == 0 && system.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        setVariables(system, COUPLED_X);
        success *= system.evaluateEqualities() == 0;
        success *= system.evaluateEqualityJacobian() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        const RealT* c     = system.c().getData();
        success           *= isEqual(c[0], 20.0);
        success           *= isEqual(c[1], -4.0);
        success           *= isEqual(c[2], 15.0);
        const auto actual  = MapFromCsr(system.getCsrEqualityJacobian());
        success           *= isEqual(actual.at(0), DependencyMap{{0, 5.0}, {2, 2.0}}, kTol);
        success           *= isEqual(actual.at(1), DependencyMap{{0, 3.0}, {1, 2.0}}, kTol);
        success           *= isEqual(actual.at(2), DependencyMap{{1, 5.0}, {2, 3.0}}, kTol);
        return success.report(__func__);
      }

      TestOutcome inequalityAssembly()
      {
        TestStatus success  = true;
        auto       system   = coupledSystem();
        success            *= system.allocate() == 0 && system.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        setVariables(system, COUPLED_X);
        success *= system.evaluateInequalities() == 0;
        success *= system.evaluateInequalityJacobian() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        success           *= isEqual(system.d().getData()[0], 31.0);
        const auto actual  = MapFromCsr(system.getCsrInequalityJacobian());
        success           *= isEqual(actual.at(0), DependencyMap{{0, 3.0}, {1, 2.0}, {2, 10.0}}, kTol);
        return success.report(__func__);
      }

      TestOutcome hessianAssembly()
      {
        TestStatus success  = true;
        auto       system   = coupledSystem();
        success            *= system.allocate() == 0 && system.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        setVariables(system, COUPLED_X);
        success *= system.evaluateHessian(2.0, COUPLED_LAMBDA_C.data(), COUPLED_LAMBDA_D.data()) == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        const auto actual  = MapFromCsr(system.getCsrHessian());
        success           *= actual.at(0).empty();
        success           *= isEqual(actual.at(1), DependencyMap{{0, 18.0}}, kTol);
        success           *= isEqual(actual.at(2), DependencyMap{{0, 3.0}, {1, 7.0}, {2, 24.0}}, kTol);
        return success.report(__func__);
      }

      TestOutcome repeatedVariables()
      {
        TestStatus success  = true;
        auto       system   = coupledSystem(true);
        success            *= system.allocate() == 0 && system.initialize() == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        setVariables(system, COUPLED_X);
        success *= system.evaluateGradient() == 0;
        success *= system.evaluateHessian(2.0, COUPLED_LAMBDA_C.data(), COUPLED_LAMBDA_D.data()) == 0;
        if (!success)
        {
          return success.report(__func__);
        }
        const RealT* gradient  = system.gradient().getData();
        success               *= isEqual(gradient[0], 0.0);
        success               *= isEqual(gradient[1], 6.0);
        success               *= isEqual(gradient[2], 5.0);
        const auto actual      = MapFromCsr(system.getCsrHessian());
        success               *= actual.at(0).empty();
        success               *= isEqual(actual.at(1), DependencyMap{{1, 36.0}}, kTol);
        success               *= isEqual(actual.at(2), DependencyMap{{1, 10.0}, {2, 24.0}}, kTol);
        return success.report(__func__);
      }

      TestOutcome connectionOwnership()
      {
        TestStatus    success = true;
        EqualityOwner owner;
        auto          component = std::make_unique<CoupledComponent>();
        for (IdxT port = 0; port < CoupledComponent::NUM_PORT_VARIABLES.size(); ++port)
        {
          success *= component->connect(port, owner) == 0;
          if (!success)
          {
            return success.report(__func__);
          }
        }
        auto system  = makeSystem<RealT>(std::move(component));
        success     *= !system.verify().passed();
        success     *= system.allocate() != 0;
        return success.report(__func__);
      }

    private:
      static constexpr std::array<RealT, 3> COUPLED_X{2.0, 3.0, 5.0};
      static constexpr std::array<RealT, 3> COUPLED_LAMBDA_C{3.0, 5.0, 7.0};
      static constexpr std::array<RealT, 1> COUPLED_LAMBDA_D{11.0};

      static SystemT coupledSystem(bool repeated = false)
      {
        auto owner           = std::make_unique<EqualityOwner>();
        auto component       = std::make_unique<CoupledComponent>();
        IdxT variable_offset = 0;
        if (repeated)
        {
          variable_offset = 1;
        }
        if (component->connect(0, *owner, {1, 1, 0}) != 0
            || component->connect(1, *owner, {variable_offset, 0, 0}) != 0
            || component->connect(2, *owner) != 0)
        {
          throw std::invalid_argument("Invalid component test connections");
        }
        return makeSystem<RealT>(std::move(owner), std::move(component));
      }

      static constexpr RealT SIGMA        = 0.7;
      static constexpr RealT LAMBDA_SCALE = 0.25;

      static OPF::SystemModelData<RealT, IdxT> parse(const std::string& text)
      {
        std::istringstream input(text);
        return OPF::parseSystemModelData(input);
      }

      /// Two-bus case with `extra` appended to its devices
      static std::string jsonCase(const std::string& extra = "")
      {
        return R"({
          "params": {"freq_base": 60.0, "va_base": 100.0e6},
          "buses": [
            {"number": 1, "class": "Bus", "name": "1", "params": {"Vmin": 0.94, "Vmax": 1.06, "reference": true}},
            {"number": 2, "class": "Bus", "name": "2", "params": {"Vmin": 0.95, "Vmax": 1.05, "B": 0.1}}
          ],
          "devices": [
            {"class": "Branch", "id": "branch_1_2", "ports": {"bus1": 1, "bus2": 2}, "params": {"R": 0.01, "X": 0.1, "Smax": 0.6}},
            {"class": "Generator", "id": "gen_1", "ports": {"bus": 1}, "params": {"Pmax": 1.5, "c1": 12.0, "in_service": false}},
            {"class": "Load", "id": "load_2", "ports": {"bus": 2}, "params": {"P": 1.2, "Q": 0.3}})"
               + extra + "]}";
      }

      /// Hessian values at `sigma`, `lambda_c`, and `lambda_d`, or none if the evaluation fails
      static std::vector<RealT> hessianValues(SystemT&                  system,
                                              RealT                     sigma,
                                              const std::vector<RealT>& lambda_c,
                                              const std::vector<RealT>& lambda_d)
      {
        if (system.evaluateHessian(sigma, lambda_c.data(), lambda_d.data()) != 0)
        {
          return {};
        }
        auto* hessian = system.getCsrHessian();
        return std::vector<RealT>(hessian->getValues(), hessian->getValues() + hessian->getNnz());
      }
    };
  } // namespace Testing
} // namespace GridKit
