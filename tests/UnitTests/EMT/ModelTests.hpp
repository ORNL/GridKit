/**
 * @file ModelTests.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Residuals and Jacobians of EMT models on small cases.
 */

#pragma once

#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <memory>
#include <span>
#include <sstream>
#include <stdexcept>

#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/StateJSONParser.hpp>
#include <GridKit/Model/EMT/SystemModel.hpp>
#include <GridKit/Model/EMT/SystemModelData.hpp>
#include <GridKit/Testing/JacobianCheck.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    /// Residuals by variable name at a state, and the Jacobian against them
    class ModelTests
    {
      using RealT   = double;
      using IdxT    = size_t;
      using SystemT = EMT::SystemModel<RealT, IdxT>;
      using LoadZT  = EMT::LoadZ<RealT, IdxT>;

    public:
      /// Load of impedance 2 + 8/(s + 4) on a single-phase bus
      TestOutcome loadZ()
      {
        TestStatus success = true;

        constexpr const char* STATE = R"({
          "buses": {"1": {"v": [1.0]}},
          "devices": {"LD1": {"i": [-0.5], "z": {"w": [0.25]}}}
        })";

        const auto system = makeSystem(LOAD_CASE, STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        auto& load = system->device<LoadZT>("LD1");

        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(system->bus("1").getResidual(EMT::BusInternalVariables::v), {-0.5});
        success *= rowsEqual(load.getResidual(EMT::LoadZInternalVariables::i), {2.0});
        success *= rowsEqual(load.submodel(EMT::LoadZSubmodels::z).getResidual(EMT::VectorFitInternalVariables::w), {-1.5});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, kTol);

        return success.report(__func__);
      }

      /// Load of impedance 2 + (1 + 0.5j)/(s + 1 - 2j) + (1 - 0.5j)/(s + 1 + 2j) on a single-phase bus
      TestOutcome conjugatePoles()
      {
        TestStatus success = true;

        constexpr const char* CASE = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
          "devices": [{"class": "LoadZ", "id": "LD1", "ports": {"bus": 1},
                       "params": {"N": 1, "Z": {"D": [[2.0]], "poles": [[-1.0, 2.0], [-1.0, -2.0]],
                                                "residues": [[[[1.0, 0.5], [1.0, -0.5]]]]}}}]
        })";

        constexpr const char* STATE = R"({
          "buses": {"1": {"v": [1.0]}},
          "devices": {"LD1": {"i": [-0.5], "z": {"w": [0.25], "v": [0.125]}}}
        })";

        const auto system = makeSystem(CASE, STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        auto& load = system->device<LoadZT>("LD1");

        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(system->bus("1").getResidual(EMT::BusInternalVariables::v), {-0.5});
        success *= rowsEqual(load.getResidual(EMT::LoadZInternalVariables::i), {0.375});
        success *= rowsEqual(load.submodel(EMT::LoadZSubmodels::z).getResidual(EMT::VectorFitInternalVariables::w), {-1.0});
        success *= rowsEqual(load.submodel(EMT::LoadZSubmodels::z).getResidual(EMT::VectorFitInternalVariables::v), {0.375});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, kTol);

        return success.report(__func__);
      }

      /// A VectorFit with an unpaired complex pole, or a complex residue at a real pole, fails verification
      TestOutcome invalidVectorFit()
      {
        TestStatus success = true;

        constexpr const char* UNPAIRED = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
          "devices": [{"class": "LoadZ", "id": "LD1", "ports": {"bus": 1},
                       "params": {"N": 1, "Z": {"D": [[2.0]], "poles": [[-1.0, 2.0]], "residues": [[[[1.0, 0.5]]]]}}}]
        })";

        constexpr const char* COMPLEX_RESIDUE = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
          "devices": [{"class": "LoadZ", "id": "LD1", "ports": {"bus": 1},
                       "params": {"N": 1, "Z": {"D": [[2.0]], "poles": [[-4.0, 0.0]], "residues": [[[[8.0, 1.0]]]]}}}]
        })";

        success *= makeSystem(UNPAIRED, "{}") == nullptr;
        success *= makeSystem(COMPLEX_RESIDUE, "{}") == nullptr;

        return success.report(__func__);
      }

      /// Unknown parameters, ports, monitors and state entries, and a port without a bus, are rejected
      TestOutcome invalidInput()
      {
        TestStatus success = true;

        constexpr const char* UNKNOWN_PARAMETER = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
          "devices": [{"class": "LoadZ", "id": "LD1", "ports": {"bus": 1}, "params": {"N": 1, "R": 1.0}}]
        })";

        constexpr const char* UNKNOWN_PORT = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
          "devices": [{"class": "LoadZ", "id": "LD1", "ports": {"bus": 1, "bus2": 1}, "params": {"N": 1}}]
        })";

        constexpr const char* MISSING_PORT = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
          "devices": [{"class": "LoadZ", "id": "LD1", "ports": {}, "params": {"N": 1}}]
        })";

        constexpr const char* UNKNOWN_MONITOR = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
          "devices": [{"class": "LoadZ", "id": "LD1", "ports": {"bus": 1}, "params": {"N": 1}, "mon": ["v"]}]
        })";

        success *= rejected(UNKNOWN_PARAMETER, "{}");
        success *= rejected(UNKNOWN_PORT, "{}");
        success *= rejected(MISSING_PORT, "{}");
        success *= rejected(UNKNOWN_MONITOR, "{}");
        success *= rejected(LOAD_CASE, R"({"devices": {"LD1": {"v": [1.0]}}})");

        return success.report(__func__);
      }

      /// Bus voltage at its base voltage, load current at the power base, and memory state at the current over its pole
      TestOutcome absoluteTolerances()
      {
        TestStatus success = true;

        const auto system = makeSystem(LOAD_CASE, "{}");

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        constexpr RealT REL_TOL = 1.0e-3;

        auto& load = system->device<LoadZT>("LD1");

        success *= system->setAbsoluteTolerance(REL_TOL) == 0;
        success *= rowsEqual(system->bus("1").absoluteTolerance(EMT::BusInternalVariables::v), {1.0});
        success *= rowsEqual(load.absoluteTolerance(EMT::LoadZInternalVariables::i), {1.0});
        success *= rowsEqual(load.submodel(EMT::LoadZSubmodels::z).absoluteTolerance(EMT::VectorFitInternalVariables::w), {0.25});

        return success.report(__func__);
      }

    private:
      /// Load of impedance 2 + 8/(s + 4) on a single-phase bus of base voltage 1 kV, with a power base of 1 MVA
      static constexpr const char* LOAD_CASE = R"({
        "params": {"freq_base": 60.0, "va_base": 1.0e6},
        "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
        "devices": [{"class": "LoadZ", "id": "LD1", "ports": {"bus": 1},
                     "params": {"N": 1, "Z": {"D": [[2.0]], "poles": [[-4.0, 0.0]], "residues": [[[[8.0, 0.0]]]]}}}]
      })";

      static constexpr RealT ALPHA = 0.5;
      static constexpr RealT kTol  = 1.0e-11;

      /// System of the case at the state, or none if it cannot be allocated
      static std::unique_ptr<SystemT> makeSystem(const char* case_json, const char* state_json)
      {
        std::istringstream input(case_json);
        auto               system = std::make_unique<SystemT>(EMT::parseSystemModelData(input));
        if (system->allocate() != 0 || system->initialize() != 0)
        {
          return nullptr;
        }
        EMT::applyState(nlohmann::json::parse(state_json), *system);
        return system;
      }

      /// Building the case at the state throws `std::invalid_argument`
      static bool rejected(const char* case_json, const char* state_json)
      {
        try
        {
          makeSystem(case_json, state_json);
        }
        catch (const std::invalid_argument&)
        {
          return true;
        }
        return false;
      }

      static bool rowsEqual(std::span<const RealT> rows, std::initializer_list<RealT> expected)
      {
        return std::ranges::equal(rows, expected, [](RealT row, RealT value)
                                  { return isEqual(row, value); });
      }
    };
  } // namespace Testing
} // namespace GridKit
