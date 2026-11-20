/**
 * @file ModelTests.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Residuals and Jacobians of EMT models on small cases.
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <format>
#include <initializer_list>
#include <limits>
#include <memory>
#include <numbers>
#include <span>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <nlohmann/json.hpp>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>
#include <GridKit/Model/EMT/Operators/Modulation/PWM/PwmImpl.hpp>
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
      using RealT      = double;
      using IdxT       = size_t;
      using SystemT    = EMT::SystemModel<RealT, IdxT>;
      using LoadZT     = EMT::LoadZ<RealT, IdxT>;
      using VectorFitT = EMT::VectorFit<RealT, IdxT>;

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

      /// All eight combinations of R, L and C on one bus
      TestOutcome loadRLC()
      {
        TestStatus success = true;

        constexpr const char* CASE = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 8, "kv": 1.0}}],
          "devices": [{"class": "LoadRLC", "id": "LD1", "terminals": {"bus": 1},
                       "params": {"N": 8, "R": [2, 4, null, null, 5, null, 2, null],
                                          "L": [3, null, null, 2, 4, null, null, 6],
                                          "C": [null, null, null, null, 2, 3, 4, 5]}}]
        })";

        constexpr const char* STATE = R"({
          "buses": {"1": {"v": [2, 4, 6, 8, 10, 12, 14, 16]}},
          "devices": {"LD1": {"i_l": [1, -4, -1, 2]}}
        })";

        const auto system = makeSystem(CASE, STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        constexpr RealT REL_TOL = 1.0e-3;

        // Voltages up to 16 raise the round-off of the difference, which passes at 1e-10 and not at 1e-11
        constexpr RealT LOAD_TOL = 1.0e-10;

        auto& load = system->device<EMT::LoadRLC<RealT, IdxT>>("LD1");

        // One inductor current for each of the four phases with an inductor
        success *= system->size() == 12;
        success *= load.yp().setToConst(0.5) == 0;
        success *= system->bus("1").yp().setToConst(0.25) == 0;
        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(load.getResidual(EMT::LoadRLCInternalVariables::i_l), {3.5, 9.0, 12.0, 19.0});
        success *= rowsEqual(system->bus("1").getResidual(EMT::BusInternalVariables::v), {0.0, -1.0, 0.0, -4.0, -3.5, -0.75, -8.0, 0.75});
        success *= load.evaluateOutputs() == 0;
        success *= rowsEqual(load.getOutput(EMT::LoadRLCOutputs::i), {0.0, -1.0, 0.0, -4.0, -3.5, -0.75, -8.0, 0.75});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, LOAD_TOL);

        success *= system->setAbsoluteTolerance(REL_TOL) == 0;
        success *= rowsEqual(load.absoluteTolerance(EMT::LoadRLCInternalVariables::i_l), {1.0, 1.0, 1.0, 1.0});

        return success.report(__func__);
      }

      /// A load whose arrays do not have N entries, or hold an entry that is not positive, fails verification
      TestOutcome invalidLoadRLC()
      {
        TestStatus success = true;

        static constexpr const char* CASE = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
          "devices": [{"class": "LoadRLC", "id": "LD1", "terminals": {"bus": 1}, "params": {"N": 1, "R": [2.0]}}]
        })";

        const auto invalid = [](auto change)
        {
          auto load = nlohmann::json::parse(CASE);
          change(load["devices"][0]["params"]);
          return makeSystem(load.dump().c_str(), "{}") == nullptr;
        };

        success *= makeSystem(CASE, "{}") != nullptr;
        success *= invalid([](auto& params)
                           { params["R"] = {2.0, 2.0}; });
        success *= invalid([](auto& params)
                           { params["L"] = nlohmann::json::array(); });
        success *= invalid([](auto& params)
                           { params["R"] = {0.0}; });
        success *= invalid([](auto& params)
                           { params["L"] = {-1.0}; });
        success *= invalid([](auto& params)
                           { params["C"] = {0.0}; });

        return success.report(__func__);
      }

      /// Load of impedance 2 + (1 + 0.5j)/(s + 1 - 2j) + (1 - 0.5j)/(s + 1 + 2j) on a single-phase bus
      TestOutcome conjugatePoles()
      {
        TestStatus success = true;

        constexpr const char* CASE = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
          "devices": [{"class": "LoadZ", "id": "LD1", "terminals": {"bus": 1},
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

      /// A VectorFit with an unpaired complex pole, a complex residue at a real pole, or dimensions other than its parent's fails verification
      TestOutcome invalidVectorFit()
      {
        TestStatus success = true;

        constexpr const char* UNPAIRED = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
          "devices": [{"class": "LoadZ", "id": "LD1", "terminals": {"bus": 1},
                       "params": {"N": 1, "Z": {"D": [[2.0]], "poles": [[-1.0, 2.0]], "residues": [[[[1.0, 0.5]]]]}}}]
        })";

        constexpr const char* COMPLEX_RESIDUE = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
          "devices": [{"class": "LoadZ", "id": "LD1", "terminals": {"bus": 1},
                       "params": {"N": 1, "Z": {"D": [[2.0]], "poles": [[-4.0, 0.0]], "residues": [[[[8.0, 1.0]]]]}}}]
        })";

        // The data of an operator give its dimensions, which its parent checks
        constexpr const char* MISSHAPED = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
          "devices": [{"class": "LoadZ", "id": "LD1", "terminals": {"bus": 1},
                       "params": {"N": 1, "Z": {"D": [[2.0, 0.0], [0.0, 2.0]]}}}]
        })";

        constexpr const char* OMITTED = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
          "devices": [{"class": "LoadZ", "id": "LD1", "terminals": {"bus": 1}, "params": {"N": 1}}]
        })";

        success *= makeSystem(UNPAIRED, "{}") == nullptr;
        success *= makeSystem(COMPLEX_RESIDUE, "{}") == nullptr;
        success *= makeSystem(MISSHAPED, "{}") == nullptr;
        success *= makeSystem(OMITTED, "{}") == nullptr;

        return success.report(__func__);
      }

      /// Unknown parameters, terminals, monitors and state entries, and a terminal without a bus, are rejected
      TestOutcome invalidInput()
      {
        TestStatus success = true;

        constexpr const char* UNKNOWN_PARAMETER = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
          "devices": [{"class": "LoadZ", "id": "LD1", "terminals": {"bus": 1}, "params": {"N": 1, "R": 1.0}}]
        })";

        constexpr const char* UNKNOWN_TERMINAL = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
          "devices": [{"class": "LoadZ", "id": "LD1", "terminals": {"bus": 1, "bus2": 1}, "params": {"N": 1}}]
        })";

        constexpr const char* MISSING_TERMINAL = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
          "devices": [{"class": "LoadZ", "id": "LD1", "terminals": {}, "params": {"N": 1}}]
        })";

        constexpr const char* UNKNOWN_MONITOR = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
          "devices": [{"class": "LoadZ", "id": "LD1", "terminals": {"bus": 1}, "params": {"N": 1}, "mon": ["v"]}]
        })";

        success *= rejected(UNKNOWN_PARAMETER, "{}");
        success *= rejected(UNKNOWN_TERMINAL, "{}");
        success *= rejected(MISSING_TERMINAL, "{}");
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

      /// A monitor is an internal variable or an output, of a device or of one of its submodels
      TestOutcome monitors()
      {
        TestStatus success = true;

        constexpr const char* STATE = R"({
          "buses": {"1": {"v": [1.0]}},
          "devices": {"LD1": {"i": [-0.5], "z": {"w": [0.25]}}}
        })";

        auto monitored                 = nlohmann::json::parse(LOAD_CASE);
        monitored["devices"][0]["mon"] = {"i", "z.w", "z.y"};
        const auto system              = makeSystem(monitored.dump().c_str(), STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        // The output of the impedance is 2 i + 8 w
        auto& z  = system->device<LoadZT>("LD1").submodel(EMT::LoadZSubmodels::z);
        success *= z.outputSize(EMT::VectorFitOutputs::y) == 1;
        success *= z.evaluateOutputs() == 0;
        success *= rowsEqual(z.getOutput(EMT::VectorFitOutputs::y), {1.0});
        success *= system->monitoring();

        // A name that the model or its submodel does not have is rejected
        const auto rejects = [](const char* section, const char* name)
        {
          auto unknown = nlohmann::json::parse(LOAD_CASE);
          if (std::string_view{section} == "mon")
          {
            unknown["devices"][0][section] = {name};
          }
          else
          {
            unknown["devices"][0][section][name] = 1;
          }
          return rejected(unknown.dump().c_str(), "{}");
        };
        success *= rejects("mon", "y");
        success *= rejects("mon", "z.i");
        success *= rejects("mon", "q.y");
        success *= rejects("inputs", "u");
        success *= rejects("outputs", "y");

        return success.report(__func__);
      }

      /// An operator in a case reads the bus voltage and adds its output to the current balance, here an admittance of 0.5
      TestOutcome operatorDevice()
      {
        TestStatus success = true;

        const auto system = makeSystem(ADMITTANCE_CASE, R"({"buses": {"1": {"v": [2.0]}}})");

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        auto& admittance = system->device<VectorFitT>("Y1");

        // An operator without poles has no variable
        success *= system->size() == 1;
        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(system->bus("1").getResidual(EMT::BusInternalVariables::v), {1.0});
        success *= admittance.evaluateOutputs() == 0;
        success *= rowsEqual(admittance.getOutput(EMT::VectorFitOutputs::y), {1.0});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, kTol);

        return success.report(__func__);
      }

      /// An input that the case leaves out is held: it reads a signal that the state gives and the system holds
      TestOutcome heldInput()
      {
        TestStatus success = true;

        auto held = nlohmann::json::parse(ADMITTANCE_CASE);
        held["devices"][0].erase("inputs");
        const auto system = makeSystem(held.dump().c_str(), R"({"devices": {"Y1": {"u": [0.25]}}})");

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        // The bus voltage and the held input
        constexpr IdxT HELD = 1;

        auto& admittance  = system->device<VectorFitT>("Y1");
        success          *= system->size() == 2;
        success          *= rowsEqual(system->held(admittance, EMT::VectorFitInputs::u), {0.25});

        // The residual of a held signal is minus its derivative
        system->yp().getData()[HELD]  = 0.5;
        success                      *= system->evaluateResidual() == 0;
        success                      *= rowsEqual(system->bus("1").getResidual(EMT::BusInternalVariables::v), {0.125});
        success                      *= isEqual(system->getResidual().getData()[HELD], -0.5);
        success                      *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, kTol);

        // The state cannot set an input that reads the bus
        success *= rejected(ADMITTANCE_CASE, R"({"devices": {"Y1": {"u": [1.0]}}})");

        return success.report(__func__);
      }

      /// An input reads a variable of a device that the case lists after it, in place
      TestOutcome variableReference()
      {
        TestStatus success = true;

        static constexpr const char* READER = R"({"class": "VectorFit", "id": "Y1", "inputs": {"u": "LD1.i"}, "outputs": {"y": 1},
                                                  "params": {"D": [[0.5]]}})";

        auto read = nlohmann::json::parse(LOAD_CASE);
        read["devices"].insert(read["devices"].begin(), nlohmann::json::parse(READER));

        const auto system = makeSystem(read.dump().c_str(), R"({
          "buses": {"1": {"v": [1.0]}},
          "devices": {"LD1": {"i": [-0.5], "z": {"w": [0.25]}}}
        })");

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        // The bus voltage, the load current and its memory state: the reader adds no variable.
        // The bus collects the load current and half of it again from the reader.
        success *= system->size() == 3;
        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(system->bus("1").getResidual(EMT::BusInternalVariables::v), {-0.75});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, kTol);

        // An input reads a device or a bus that the case has, and a variable or an output that it has
        const auto rejects = [&read](const char* reference)
        {
          auto unknown                         = read;
          unknown["devices"][0]["inputs"]["u"] = reference;
          return rejected(unknown.dump().c_str(), "{}");
        };
        success *= rejects("LD2.i");
        success *= rejects("LD1.q");
        success *= rejects("2");

        // Only an input reads an output
        auto terminal                               = read;
        terminal["devices"][1]["terminals"]["bus"]  = "Y1.y";
        success                                    *= rejected(terminal.dump().c_str(), "{}");

        return success.report(__func__);
      }

      /// Two inputs read the output of an operator through the one signal that follows it
      TestOutcome outputReference()
      {
        TestStatus success = true;

        static constexpr const char* CASE = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
          "devices": [{"class": "VectorFit", "id": "Y2", "inputs": {"u": "Y1.y"}, "outputs": {"y": 1}, "params": {"D": [[4.0]]}},
                      {"class": "VectorFit", "id": "Y3", "inputs": {"u": "Y1.y"}, "outputs": {"y": 1}, "params": {"D": [[1.0]]}},
                      {"class": "VectorFit", "id": "Y1", "inputs": {"u": 1}, "params": {"D": [[0.5]]}}]
        })";

        const auto system = makeSystem(CASE, R"({"buses": {"1": {"v": [2.0]}}})");

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        constexpr RealT REL_TOL = 1.0e-3;

        // The bus voltage and the signal of `Y1.y`, which both readers share
        constexpr IdxT SIGNAL  = 1;
        success               *= system->size() == 2;

        // The signal is 0.5 and the output it follows is 1, and the readers add 4 + 1 times the signal to the bus
        system->y().getData()[SIGNAL]  = 0.5;
        success                       *= system->evaluateResidual() == 0;
        success                       *= isEqual(system->getResidual().getData()[SIGNAL], 0.5);
        success                       *= rowsEqual(system->bus("1").getResidual(EMT::BusInternalVariables::v), {2.5});
        success                       *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, kTol);

        // A signal takes the tolerance of the output it follows, here the relative tolerance
        success *= system->setAbsoluteTolerance(REL_TOL) == 0;
        success *= isEqual(system->absoluteTolerance().getData()[SIGNAL], REL_TOL);

        // An input fits the output it reads
        auto misfit                          = nlohmann::json::parse(CASE);
        misfit["devices"][0]["params"]["D"]  = {{4.0, 4.0}};
        success                             *= rejected(misfit.dump().c_str(), "{}");

        // An output is added to one residual, whichever device the case lists first
        auto twice                           = nlohmann::json::parse(CASE);
        twice["devices"][2]["outputs"]["y"]  = 1;
        success                             *= rejected(twice.dump().c_str(), "{}");
        std::swap(twice["devices"][0], twice["devices"][2]);
        success *= rejected(twice.dump().c_str(), "{}");

        return success.report(__func__);
      }

      /// Pmu on a three-phase bus at (0, 1, -1) with currents (3, 0, -3): squared voltages 2, active power 3 and
      /// reactive power 3 sqrt(3)
      TestOutcome pmu()
      {
        TestStatus success = true;

        const auto system = makeSystem(PMU_CASE, PMU_STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        constexpr RealT REL_TOL    = 1.0e-3;
        constexpr RealT VQ         = TWO<RealT> / std::numbers::sqrt3_v<RealT>;
        constexpr RealT OMEGA_BASE = TWO<RealT> * std::numbers::pi_v<RealT> * 60.0;

        // The frequency base raises the round-off of the difference, which passes at 1e-8 and not at 1e-9
        constexpr RealT PMU_TOL = 1.0e-8;

        auto& pmu = system->device<EMT::Pmu<RealT, IdxT>>("PMU1");

        // Three bus voltages, three Pmu and two Pll variables, and the three held currents: the Pmu connects its own
        // input, and nothing here reads its outputs
        constexpr IdxT VMS  = 3;
        success            *= system->size() == 11;

        // Each mean moves toward its instantaneous value over half a cycle of the frequency base
        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(pmu.getResidual(EMT::PmuInternalVariables::vms), {-240.0});
        success *= rowsEqual(pmu.getResidual(EMT::PmuInternalVariables::p), {360.0});
        success *= rowsEqual(pmu.getResidual(EMT::PmuInternalVariables::q), {360.0 * std::numbers::sqrt3_v<RealT>});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, PMU_TOL);

        success *= pmu.evaluateOutputs() == 0;
        success *= rowsEqual(pmu.getOutput(EMT::PmuOutputs::vrms), {2.0});
        success *= rowsEqual(pmu.getOutput(EMT::PmuOutputs::theta), {0.25});
        success *= rowsEqual(pmu.getOutput(EMT::PmuOutputs::omega), {OMEGA_BASE + 0.5 * VQ * std::cos(0.25) + 0.5});

        // The tolerances follow the base voltage of the bus, 2 kV, and the system power base, 1 MVA
        success *= system->setAbsoluteTolerance(REL_TOL) == 0;
        success *= rowsEqual(std::span<const RealT>{system->absoluteTolerance().getData() + VMS, 3}, {4.0, 1000.0, 1000.0});

        return success.report(__func__);
      }

      /// Operators read the outputs of a Pmu, which read its Pll, through the signals that follow them
      TestOutcome pmuReference()
      {
        TestStatus success = true;

        auto read = nlohmann::json::parse(PMU_CASE);
        read["buses"].push_back(nlohmann::json::parse(R"({"number": 2, "class": "Bus", "params": {"N": 1, "kv": 1.0}})"));
        for (const char* output : {"omega", "theta", "vrms"})
        {
          read["devices"].push_back(nlohmann::json::parse(std::format(
              R"({{"class": "VectorFit", "id": "{0}", "inputs": {{"u": "PMU1.{0}"}}, "outputs": {{"y": 2}}, "params": {{"D": [[0.5]]}}}})",
              output)));
        }
        const auto system = makeSystem(read.dump().c_str(), PMU_STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        constexpr RealT REL_TOL    = 1.0e-3;
        constexpr RealT VQ         = TWO<RealT> / std::numbers::sqrt3_v<RealT>;
        constexpr RealT OMEGA_BASE = TWO<RealT> * std::numbers::pi_v<RealT> * 60.0;

        // The frequency base raises the round-off of the difference, which passes at 1e-8 and not at 1e-9
        constexpr RealT PMU_TOL = 1.0e-8;

        // Four bus voltages, three Pmu and two Pll variables, the three held currents, and the signals of `omega`, `theta`
        // and `vrms`
        constexpr IdxT OMEGA  = 12;
        constexpr IdxT THETA  = 13;
        constexpr IdxT VRMS   = 14;
        success              *= system->size() == 15;

        // The signals start at zero, so their residuals are the outputs
        success *= system->evaluateResidual() == 0;
        success *= isEqual(system->getResidual().getData()[OMEGA], OMEGA_BASE + 0.5 * VQ * std::cos(0.25) + 0.5);
        success *= isEqual(system->getResidual().getData()[THETA], 0.25);
        success *= isEqual(system->getResidual().getData()[VRMS], 2.0);
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, PMU_TOL);

        // The signals take the tolerances of the outputs
        success *= system->setAbsoluteTolerance(REL_TOL) == 0;
        success *= isEqual(system->absoluteTolerance().getData()[OMEGA], REL_TOL * OMEGA_BASE);
        success *= isEqual(system->absoluteTolerance().getData()[THETA], REL_TOL);
        success *= isEqual(system->absoluteTolerance().getData()[VRMS], 2.0);

        // A dead bus has no RMS voltage
        const auto dead  = makeSystem(read.dump().c_str(), "{}");
        success         *= dead != nullptr && dead->evaluateResidual() == 0 && isEqual(dead->getResidual().getData()[VRMS], 0.0);

        return success.report(__func__);
      }

      /// A layout is set once, and an absent submodel is an error
      TestOutcome layout()
      {
        TestStatus success = true;

        struct Probe final : EMT::Component<RealT, IdxT>
        {
          Model::ConfigurationChecks verify() const override
          {
            return Model::ConfigurationChecks("Probe");
          }

          void setLayoutTwice()
          {
            this->setLayout({1}, {}, {}, {});
            this->setLayout({1}, {}, {}, {});
          }

          void absentSubmodel()
          {
            this->submodel<Probe>(EMT::LoadZSubmodels::z);
          }
        };

        success *= throws<std::logic_error>([]
                                            { Probe{}.setLayoutTwice(); });
        success *= throws<std::invalid_argument>([]
                                                 { Probe{}.absentSubmodel(); });

        return success.report(__func__);
      }

      /// Unit source behind a unit conductance on a single-phase bus, at time zero
      TestOutcome voltageSource()
      {
        TestStatus success = true;

        constexpr const char* CASE = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
          "devices": [{"class": "VoltageSource", "id": "VS1", "terminals": {"bus": 1},
                       "params": {"N": 1, "E": [1.0], "phi": [0.0], "omega": 1.0, "Y": {"D": [[1.0]]}}}]
        })";

        constexpr const char* STATE = R"({
          "buses": {"1": {"v": [1.0]}},
          "devices": {"VS1": {"u": [0.5], "i": [0.125]}}
        })";

        const auto system = makeSystem(CASE, STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        auto& source = system->device("VS1");

        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(system->bus("1").getResidual(EMT::BusInternalVariables::v), {0.125});
        success *= rowsEqual(source.getResidual(EMT::VoltageSourceInternalVariables::u), {std::numbers::sqrt2_v<RealT> - 1.5});
        success *= rowsEqual(source.getResidual(EMT::VoltageSourceInternalVariables::i), {0.375});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, kTol);

        // The admittance has the dimensions of the source
        auto misshaped                               = nlohmann::json::parse(CASE);
        misshaped["devices"][0]["params"]["Y"]["D"]  = {{1.0, 0.0}, {0.0, 1.0}};
        success                                     *= makeSystem(misshaped.dump().c_str(), "{}") == nullptr;

        return success.report(__func__);
      }

      /// Series resistance 1 and shunt conductance 0.5 between single-phase buses
      TestOutcome lineLumped()
      {
        TestStatus success = true;

        constexpr const char* CASE = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}},
                    {"number": 2, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
          "devices": [{"class": "LineLumped", "id": "L1", "terminals": {"bus1": 1, "bus2": 2},
                       "params": {"N": 1, "K": 1, "conductors": [1], "dx": 2.0, "Zp": {"D": [[0.5]]}, "Yp": {"D": [[0.25]]}}}]
        })";

        constexpr const char* STATE = R"({
          "buses": {"1": {"v": [2.0]}, "2": {"v": [1.0]}},
          "devices": {"L1": {"i12": [0.5]}}
        })";

        const auto system = makeSystem(CASE, STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        auto& line = system->device<EMT::LineLumped<RealT, IdxT>>("L1");

        // Half of the shunt conductance, 0.25, draws from each bus
        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(system->bus("1").getResidual(EMT::BusInternalVariables::v), {-1.0});
        success *= rowsEqual(system->bus("2").getResidual(EMT::BusInternalVariables::v), {0.25});
        success *= rowsEqual(line.getResidual(EMT::LineLumpedInternalVariables::i12), {-0.5});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, kTol);

        // The shunt currents are the outputs of the shunt operators
        auto& y1  = line.submodel(EMT::LineLumpedSubmodels::y1);
        auto& y2  = line.submodel(EMT::LineLumpedSubmodels::y2);
        success  *= y1.evaluateOutputs() == 0 && y2.evaluateOutputs() == 0;
        success  *= rowsEqual(y1.getOutput(EMT::VectorFitOutputs::y), {-0.5});
        success  *= rowsEqual(y2.getOutput(EMT::VectorFitOutputs::y), {-0.25});

        // A line without a shunt admittance has no shunt operators
        auto series = nlohmann::json::parse(CASE);
        series["devices"][0]["params"].erase("Yp");
        const auto series_system = makeSystem(series.dump().c_str(), STATE);

        success *= series_system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        success *= series_system->size() == 3;
        success *= series_system->evaluateResidual() == 0;
        success *= rowsEqual(series_system->bus("1").getResidual(EMT::BusInternalVariables::v), {-0.5});
        success *= rowsEqual(series_system->bus("2").getResidual(EMT::BusInternalVariables::v), {0.5});
        success *= jacobianMatchesResidual(*series_system, ZERO<RealT>, ALPHA, kTol);
        success *= rejected(series.dump().c_str(), R"({"devices": {"L1": {"y1": {"w": []}}}})");

        return success.report(__func__);
      }

      /// Two conductors on one phase, each of series resistance 1 and shunt conductance 0.5
      TestOutcome bundledLineLumped()
      {
        TestStatus success = true;

        constexpr const char* CASE = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}},
                    {"number": 2, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
          "devices": [{"class": "LineLumped", "id": "L1", "terminals": {"bus1": 1, "bus2": 2},
                       "params": {"N": 1, "K": 2, "conductors": [1, 1], "dx": 2.0,
                                  "Zp": {"D": [[0.5, 0.0], [0.0, 0.5]]}, "Yp": {"D": [[0.25, 0.0], [0.0, 0.25]]}}}]
        })";

        constexpr const char* STATE = R"({
          "buses": {"1": {"v": [2.0]}, "2": {"v": [1.0]}},
          "devices": {"L1": {"i12": [0.5, 0.25]}}
        })";

        const auto system = makeSystem(CASE, STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        auto& line = system->device("L1");

        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(system->bus("1").getResidual(EMT::BusInternalVariables::v), {-1.75});
        success *= rowsEqual(system->bus("2").getResidual(EMT::BusInternalVariables::v), {0.25});
        success *= rowsEqual(line.getResidual(EMT::LineLumpedInternalVariables::i12), {-0.5, -0.75});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, kTol);

        return success.report(__func__);
      }

      /// A line whose conductors or operators do not match its phase and conductor counts fails verification
      TestOutcome invalidLineLumped()
      {
        TestStatus success = true;

        static constexpr const char* CASE = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}},
                    {"number": 2, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
          "devices": [{"class": "LineLumped", "id": "L1", "terminals": {"bus1": 1, "bus2": 2},
                       "params": {"N": 1, "K": 1, "conductors": [1], "dx": 2.0, "Zp": {"D": [[0.5]]}, "Yp": {"D": [[0.25]]}}}]
        })";

        const auto invalid = [](auto change)
        {
          auto line = nlohmann::json::parse(CASE);
          change(line["devices"][0]["params"]);
          return makeSystem(line.dump().c_str(), "{}") == nullptr;
        };

        success *= makeSystem(CASE, "{}") != nullptr;
        success *= invalid([](auto& params)
                           { params["conductors"] = {2}; });
        success *= invalid([](auto& params)
                           { params["conductors"] = {1, 1}; });
        success *= invalid([](auto& params)
                           { params["Zp"]["D"] = {{0.5, 0.0}, {0.0, 0.5}}; });
        success *= invalid([](auto& params)
                           { params["Yp"]["D"] = {{0.25, 0.0}, {0.0, 0.25}}; });

        return success.report(__func__);
      }

      /// A closed switch joins the voltages of its single-phase buses
      TestOutcome closedSwitch()
      {
        TestStatus success = true;

        const auto system = makeSystem(SWITCH_CASE, R"({
          "buses": {"1": {"v": [2.0]}, "2": {"v": [1.0]}},
          "devices": {"S1": {"i12": [0.5]}}
        })");

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(system->bus("1").getResidual(EMT::BusInternalVariables::v), {-0.5});
        success *= rowsEqual(system->bus("2").getResidual(EMT::BusInternalVariables::v), {0.5});
        success *= rowsEqual(system->device("S1").getResidual(EMT::SwitchInternalVariables::i12), {-1.0});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, kTol);

        return success.report(__func__);
      }

      /// An open switch blocks the current between its single-phase buses
      TestOutcome openSwitch()
      {
        TestStatus success = true;

        const auto system = makeSystem(SWITCH_CASE, R"({
          "buses": {"1": {"v": [2.0]}, "2": {"v": [1.0]}},
          "devices": {"S1": {"open": [1.0]}}
        })");

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        system->device("S1").y(EMT::SwitchInternalVariables::i12)[0] = 0.5;

        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(system->bus("1").getResidual(EMT::BusInternalVariables::v), {-0.5});
        success *= rowsEqual(system->bus("2").getResidual(EMT::BusInternalVariables::v), {0.5});
        success *= rowsEqual(system->device("S1").getResidual(EMT::SwitchInternalVariables::i12), {0.5});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, kTol);

        return success.report(__func__);
      }

      /// A value written to `open` waits while the contact carries current, and is set when the current crosses zero
      TestOutcome switchWaitsForZero()
      {
        TestStatus success = true;

        const auto system = makeSystem(SWITCH_CASE, R"({
          "devices": {"S1": {"i12": [0.5], "open": [1.0]}}
        })");

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        auto&      contact = system->device<EMT::Switch<RealT, IdxT>>("S1");
        const auto open    = system->held(contact, EMT::SwitchInputs::open);

        RealT current  = ZERO<RealT>;
        success       *= open[0] == ZERO<RealT>;
        success       *= system->sizeRoots() == 1;
        success       *= system->evaluateRoots(system->y().getData(), &current) == 0;
        success       *= current == 0.5;

        system->applyRoots({true});
        success *= open[0] == ONE<RealT>;
        success *= contact.y(EMT::SwitchInternalVariables::i12)[0] == ZERO<RealT>;
        success *= system->sizeRoots() == 0;

        return success.report(__func__);
      }

      /// Two switches read one command, a held signal that the case declares after them
      TestOutcome sharedCommand()
      {
        TestStatus success = true;

        constexpr const char* CASE = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}},
                    {"number": 2, "class": "Bus", "params": {"N": 1, "kv": 1.0}},
                    {"number": 3, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
          "devices": [{"class": "Switch", "id": "S1", "terminals": {"bus1": 1, "bus2": 2}, "inputs": {"open": "K1.s"}, "params": {"N": 1}},
                      {"class": "Switch", "id": "S2", "terminals": {"bus1": 2, "bus2": 3}, "inputs": {"open": "K1.s"}, "params": {"N": 1}},
                      {"class": "Signal", "id": "K1", "params": {"N": 1, "held": true}}]
        })";

        const auto system = makeSystem(CASE, R"({
          "buses": {"1": {"v": [3.0]}, "2": {"v": [2.0]}, "3": {"v": [1.0]}},
          "devices": {"S1": {"i12": [0.5]}, "S2": {"i12": [0.25]}, "K1": {"s": [1.0]}}
        })");

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        // Three bus voltages, two switch currents and the command, which both switches read in place.
        // The command is 1, so both switches are open.
        success *= system->size() == 6;
        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(system->device("S1").getResidual(EMT::SwitchInternalVariables::i12), {0.5});
        success *= rowsEqual(system->device("S2").getResidual(EMT::SwitchInternalVariables::i12), {0.25});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, kTol);

        return success.report(__func__);
      }

      /// Wye to delta transformer of unit ratings with winding 1 saturated on phase a
      TestOutcome transformer()
      {
        TestStatus success = true;

        constexpr const char* CASE = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 3, "kv": 0.001}},
                    {"number": 2, "class": "Bus", "params": {"N": 3, "kv": 0.001}}],
          "devices": [{"class": "Transformer", "id": "T1", "terminals": {"bus1": 1, "bus2": 2},
                       "params": {"S": 1.0, "V1": 1.0, "V2": 1.0, "f": 60.0, "R": 0.5, "X": 0.5, "I0": 0.5, "P0": 0.1, "tap": 1.1, "split": 0.35,
                                  "P2": [[1.0, 0.0, -1.0], [-1.0, 1.0, 0.0], [0.0, -1.0, 1.0]]}}]
        })";

        constexpr const char* STATE = R"({
          "buses": {"1": {"v": [1.0, 0.0, 0.0]}, "2": {"v": [0.0, 0.0, 0.0]}},
          "devices": {"T1": {"i12": [0.5, 0.0, 0.0], "psi1": [2.0, 0.0, 0.0]}}
        })";

        const auto system = makeSystem(CASE, STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        using Variables = EMT::TransformerInternalVariables;
        using Outputs   = EMT::TransformerOutputs;

        auto& transformer = system->device<EMT::Transformer<RealT, IdxT>>("T1");

        // The curvature of the saturation ramp raises the error of the difference, which passes at 1e-10 and not at 1e-11
        constexpr RealT TRANSFORMER_TOL = 1.0e-10;

        const RealT PEAK   = std::sqrt(TWO<RealT> / THREE<RealT>);
        const RealT DELTA  = std::numbers::sqrt2_v<RealT> / THREE<RealT>;
        success           *= transformer.evaluateOutputs() == 0;
        const RealT e1     = transformer.getOutput(Outputs::e1)[0];
        const RealT e2     = transformer.getOutput(Outputs::e2)[0];
        const RealT im1    = transformer.getOutput(Outputs::im1)[0];
        const RealT im2    = transformer.getOutput(Outputs::im2)[0];
        const RealT iw1    = transformer.getOutput(Outputs::iw1)[0];
        const RealT iw2    = transformer.getOutput(Outputs::iw2)[0];

        // Winding equations of phase a, with half of R in each winding and the core split 0.35 to 0.65.
        // Winding 1 is past its knee at 1.2, where the slope changes from 1/Lm to 1/Lsat = 1.
        const RealT INV_LM  = std::sqrt(0.5 * 0.5 - 0.1 * 0.1);
        const RealT g1      = 2.0 * INV_LM + (ONE<RealT> - INV_LM) * 0.8;
        success            *= system->size() == 15;
        success            *= isEqual(e1, ONE<RealT> / PEAK - 0.25 * iw1);
        success            *= isEqual(e2, -1.1 * 0.25 * iw2);
        success            *= isEqual(im2, 0.65 * 0.1 * e2);
        success            *= isEqual(iw1, im1 + 0.5);
        success            *= isEqual(iw2, 1.1 * (im2 - 0.5));
        success            *= isEqual(im1, 0.35 * (g1 + 0.1 * e1));

        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(system->bus("1").getResidual(EMT::BusInternalVariables::v), {-PEAK * iw1, 0.0, 0.0});
        success *= rowsEqual(system->bus("2").getResidual(EMT::BusInternalVariables::v), {-DELTA * iw2, DELTA * iw2, 0.0});
        success *= rowsEqual(transformer.getResidual(Variables::i12), {e2 - e1, 0.0, 0.0});
        success *= rowsEqual(transformer.getResidual(Variables::psi1), {-e1, 0.0, 0.0});
        success *= rowsEqual(transformer.getResidual(Variables::psi2), {-e2, 0.0, 0.0});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, TRANSFORMER_TOL);

        return success.report(__func__);
      }

      /// Machine of unit ratings at angle zero and synchronous speed on a dead bus, saturated on both axes
      TestOutcome machine()
      {
        TestStatus success = true;

        const auto system = makeSystem(MACHINE_CASE, MACHINE_STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        const RealT OMEGA_BASE = TWO<RealT> * std::numbers::pi_v<RealT> * 60.0;
        const RealT PEAK       = std::sqrt(TWO<RealT> / THREE<RealT>);
        const RealT QUADRATURE = 0.2 * std::numbers::sqrt3_v<RealT>;

        using Variables = EMT::MachineInternalVariables;
        using Outputs   = EMT::MachineOutputs;

        auto& machine = system->device<EMT::Machine<RealT, IdxT>>("M1");

        // The saturation function of an axis is (|psia| - 0.8)^2, so the magnetizing inductances are divided by
        // 1.04 at air-gap flux 1.0 in d and by 1.25 at -1.3 in q.
        // The air-gap fluxes are then 2 (0.9 + 1.4 + 0.96) / 6.52 = 1.0 and -2 (1.5 + 1.5 + 1.7125) / 7.25 = -1.3.
        success *= system->size() == 16;
        success *= machine.evaluateOutputs() == 0;
        success *= rowsEqual(machine.getOutput(Outputs::id), {0.2});
        success *= rowsEqual(machine.getOutput(Outputs::iq), {0.4});
        success *= rowsEqual(machine.getOutput(Outputs::i0), {0.0});
        success *= rowsEqual(machine.getOutput(Outputs::ifd), {0.8});
        success *= rowsEqual(machine.getOutput(Outputs::i1d), {-0.08});
        success *= rowsEqual(machine.getOutput(Outputs::i1q), {-0.4});
        success *= rowsEqual(machine.getOutput(Outputs::i2q), {-0.825});
        success *= rowsEqual(machine.getOutput(Outputs::te), {0.66});
        success *= rowsEqual(machine.getOutput(Outputs::is), {0.2, -0.1 + QUADRATURE, -0.1 - QUADRATURE});

        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(system->bus("1").getResidual(EMT::BusInternalVariables::v),
                             {0.2 * PEAK, (-0.1 + QUADRATURE) * PEAK, (-0.1 - QUADRATURE) * PEAK});
        success *= rowsEqual(machine.getResidual(Variables::theta), {-OMEGA_BASE});
        success *= rowsEqual(machine.getResidual(Variables::omega), {-0.34});
        success *= rowsEqual(machine.getResidual(Variables::psid), {1.5});
        success *= rowsEqual(machine.getResidual(Variables::psiq), {0.9});
        success *= rowsEqual(machine.getResidual(Variables::psi0), {0.0});
        success *= rowsEqual(machine.getResidual(Variables::psifd), {-0.1});
        success *= rowsEqual(machine.getResidual(Variables::psi1d), {-0.08});
        success *= rowsEqual(machine.getResidual(Variables::psi1q), {-0.4});
        success *= rowsEqual(machine.getResidual(Variables::psi2q), {-0.825});
        success *= rowsEqual(machine.getResidual(Variables::psiad), {0.0});
        success *= rowsEqual(machine.getResidual(Variables::psiaq), {0.0});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, MACHINE_TOL);

        return success.report(__func__);
      }

      /// The iron-loss resistance of a machine draws the bus voltage over it from the stator currents
      TestOutcome machineIronLoss()
      {
        TestStatus success = true;

        auto data                           = nlohmann::json::parse(MACHINE_CASE);
        auto state                          = nlohmann::json::parse(MACHINE_STATE);
        data["devices"][0]["params"]["Rfe"] = 2.0;
        state["buses"]["1"]["v"]            = {1.0, -0.5, -0.5};
        const auto system                   = makeSystem(data.dump().c_str(), state.dump().c_str());

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        const RealT PEAK       = std::sqrt(TWO<RealT> / THREE<RealT>);
        const RealT QUADRATURE = 0.2 * std::numbers::sqrt3_v<RealT>;

        // The conductance of a phase is S / (Rfe V^2) = 0.5 S
        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(system->bus("1").getResidual(EMT::BusInternalVariables::v),
                             {0.2 * PEAK - 0.5, (-0.1 + QUADRATURE) * PEAK + 0.25, (-0.1 - QUADRATURE) * PEAK + 0.25});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, MACHINE_TOL);

        return success.report(__func__);
      }

      /// A Pll reads the stator currents of a machine, an output, through the signal that follows it
      TestOutcome machineCurrentReference()
      {
        TestStatus success = true;

        auto read = nlohmann::json::parse(MACHINE_CASE);
        read["devices"].push_back(nlohmann::json::parse(
            R"({"class": "Pll", "id": "PLL1", "inputs": {"v": "M1.is"}, "params": {"Kp": 0.5, "Ki": 0.25}})"));
        const auto system = makeSystem(read.dump().c_str(), MACHINE_STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        constexpr RealT REL_TOL = 1.0e-3;

        const RealT QUADRATURE = 0.2 * std::numbers::sqrt3_v<RealT>;

        // Three bus voltages, eleven machine and two Pll variables, the two held inputs, and the three signals of `is`.
        // The other outputs of the machine add no variable.
        constexpr IdxT SIGNAL  = 18;
        success               *= system->size() == 21;

        // The signals start at zero, so their residuals are the stator currents
        success *= system->evaluateResidual() == 0;
        const std::span<const RealT> signal{system->getResidual().getData() + SIGNAL, 3};
        success *= rowsEqual(signal, {0.2, -0.1 + QUADRATURE, -0.1 - QUADRATURE});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, MACHINE_TOL);

        // A signal takes the tolerance of the output it follows, here the relative tolerance
        success *= system->setAbsoluteTolerance(REL_TOL) == 0;
        const std::span<const RealT> tol{system->absoluteTolerance().getData() + SIGNAL, 3};
        success *= rowsEqual(tol, {REL_TOL, REL_TOL, REL_TOL});

        // An input fits the output it reads
        read["devices"][1]["inputs"]["v"]  = "M1.te";
        success                           *= rejected(read.dump().c_str(), "{}");

        return success.report(__func__);
      }

      /// Pll at angle zero on a three-phase bus at (0, 1, -1), whose q-axis voltage is 2/sqrt(3)
      TestOutcome pll()
      {
        TestStatus success = true;

        const auto system = makeSystem(PLL_CASE, PLL_STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        constexpr RealT VQ         = TWO<RealT> / std::numbers::sqrt3_v<RealT>;
        constexpr RealT OMEGA_BASE = TWO<RealT> * std::numbers::pi_v<RealT> * 60.0;

        // The frequency base raises the round-off of the difference, which passes at 1e-8 and not at 1e-9
        constexpr RealT PLL_TOL = 1.0e-8;

        auto& pll = system->device<EMT::Pll<RealT, IdxT>>("PLL1");

        // The frequency is an output, which nothing here reads
        success *= system->size() == 5;
        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(system->bus("1").getResidual(EMT::BusInternalVariables::v), {0.0, 0.0, 0.0});
        success *= rowsEqual(pll.getResidual(EMT::PllInternalVariables::theta), {OMEGA_BASE + 0.5 * VQ + 0.5});
        success *= rowsEqual(pll.getResidual(EMT::PllInternalVariables::x), {0.25 * VQ});
        success *= pll.evaluateOutputs() == 0;
        success *= rowsEqual(pll.getOutput(EMT::PllOutputs::omega), {OMEGA_BASE + 0.5 * VQ + 0.5});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, PLL_TOL);

        return success.report(__func__);
      }

      /// An operator reads the frequency of a Pll, an output with Enzyme derivatives, through the signal that follows it
      TestOutcome pllFrequencyReference()
      {
        TestStatus success = true;

        auto read = nlohmann::json::parse(PLL_CASE);
        read["buses"].push_back(nlohmann::json::parse(R"({"number": 2, "class": "Bus", "params": {"N": 1, "kv": 1.0}})"));
        read["devices"].push_back(nlohmann::json::parse(
            R"({"class": "VectorFit", "id": "Y1", "inputs": {"u": "PLL1.omega"}, "outputs": {"y": 2}, "params": {"D": [[0.5]]}})"));
        const auto system = makeSystem(read.dump().c_str(), PLL_STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        constexpr RealT REL_TOL    = 1.0e-3;
        constexpr RealT VQ         = TWO<RealT> / std::numbers::sqrt3_v<RealT>;
        constexpr RealT OMEGA_BASE = TWO<RealT> * std::numbers::pi_v<RealT> * 60.0;

        // The frequency base raises the round-off of the difference, which passes at 1e-8 and not at 1e-9
        constexpr RealT PLL_TOL = 1.0e-8;

        // Four bus voltages, the two Pll variables and the signal of `omega`
        constexpr IdxT SIGNAL  = 6;
        success               *= system->size() == 7;

        // The signal starts at zero, so its residual is the frequency
        success *= system->evaluateResidual() == 0;
        success *= isEqual(system->getResidual().getData()[SIGNAL], OMEGA_BASE + 0.5 * VQ + 0.5);
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, PLL_TOL);

        // The signal takes the tolerance of the frequency
        success *= system->setAbsoluteTolerance(REL_TOL) == 0;
        success *= isEqual(system->absoluteTolerance().getData()[SIGNAL], REL_TOL * OMEGA_BASE);

        return success.report(__func__);
      }

      /// The mean of each switching function over a carrier period is its duty ratio, for smooth and sharp edges
      TestOutcome pwmMean()
      {
        TestStatus success = true;

        constexpr IdxT SAMPLES = 400;

        const std::array<RealT, 3> m{0.5, -0.25, 0.0};
        for (const RealT fc : {1000.0, 10.0})
        {
          EMT::PwmData<RealT, IdxT> data;
          data.N         = 3;
          data.fc        = fc;
          data.injection = false;
          const EMT::Pwm<RealT, IdxT> pwm(data);

          success *= pwm.verify().passed();

          std::array<RealT, 3> mean{};
          for (IdxT k = 0; k < SAMPLES; ++k)
          {
            std::array<RealT, 3> s{};
            pwm.switching(m.data(), static_cast<RealT>(k) / (static_cast<RealT>(SAMPLES) * fc), s.data());
            for (IdxT n = 0; n < 3; ++n)
            {
              mean[n] += s[n] / static_cast<RealT>(SAMPLES);
            }
          }

          // Rounding accumulates over the samples
          const RealT                tol = static_cast<RealT>(SAMPLES) * std::numeric_limits<RealT>::epsilon();
          const std::array<RealT, 3> duty{0.75, 0.375, 0.5};
          for (IdxT n = 0; n < 3; ++n)
          {
            success *= isEqual(mean[n], duty[n], tol);
          }
        }

        return success.report(__func__);
      }

      /// The smoothing scale times the carrier period selects how the switching functions are evaluated
      TestOutcome pwmSeries()
      {
        TestStatus success = true;

        const auto series = [](RealT mu)
        {
          EMT::PwmData<RealT, IdxT> data;
          data.N  = 3;
          data.fc = 1000.0;
          data.mu = mu;
          return EMT::Pwm<RealT, IdxT>(data).series();
        };

        success *= series(240.0) == EMT::PwmSeries::duty;
        success *= series(1000.0) == EMT::PwmSeries::harmonics;
        success *= series(50000.0) == EMT::PwmSeries::pulses;

        return success.report(__func__);
      }

      /// Bridge at angle zero and modulation (0.5, 0) between DC terminals at 4 and 0, with duty ratios (0.75, 0.375, 0.375)
      TestOutcome bridge()
      {
        TestStatus success = true;

        const auto system = makeSystem(BRIDGE_CASE, BRIDGE_STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        constexpr RealT TWO_THIRDS = TWO<RealT> / THREE<RealT>;

        // The DC voltage of 4 raises the round-off of the difference, which passes at 1e-10 and not at 1e-11
        constexpr RealT BRIDGE_TOL = 1.0e-10;

        auto& converter = system->device("BR1");

        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(system->bus("1").getResidual(EMT::BusInternalVariables::v), {0.5, 0.0, 0.0});
        success *= rowsEqual(system->bus("2").getResidual(EMT::BusInternalVariables::v), {-0.75});
        success *= rowsEqual(system->bus("3").getResidual(EMT::BusInternalVariables::v), {-0.25});
        success *= rowsEqual(converter.getResidual(EMT::BridgeInternalVariables::i), {0.0, -1.5, -1.5});
        success *= rowsEqual(converter.getResidual(EMT::BridgeInternalVariables::vo), {-0.5, 0.0, 0.0});
        success *= rowsEqual(converter.getResidual(EMT::BridgeInternalVariables::ig), {-0.5, 0.0, 0.0});
        success *= rowsEqual(converter.getResidual(EMT::BridgeInternalVariables::idq), {TWO_THIRDS, 0.0});
        success *= rowsEqual(converter.getResidual(EMT::BridgeInternalVariables::vodq), {TWO<RealT> * TWO_THIRDS, 0.0});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, BRIDGE_TOL);

        return success.report(__func__);
      }

      /// The Jacobian of a bridge matches its residual for carrier harmonics and for resolved pulses, away from a carrier edge
      TestOutcome bridgeSwitching()
      {
        TestStatus success = true;

        // A quarter of a carrier period after its start
        constexpr RealT TIME = 2.5e-4;

        // The resolved edges raise the error of the difference, which passes at 1e-9 and not at 1e-10
        constexpr RealT SWITCHING_TOL = 1.0e-9;

        for (const RealT mu : {1000.0, 50000.0})
        {
          auto switching                                 = nlohmann::json::parse(BRIDGE_CASE);
          switching["devices"][0]["params"]["pwm"]["mu"] = mu;
          const auto system                              = makeSystem(switching.dump().c_str(), BRIDGE_STATE);

          success *= system != nullptr && jacobianMatchesResidual(*system, TIME, ALPHA, SWITCHING_TOL);
        }

        return success.report(__func__);
      }

      /// A 1 A peak rating and axis (0.6, 0.8) project the orders (0.6, 0.8) to (1, 0) A
      TestOutcome currentControl()
      {
        TestStatus success = true;

        constexpr const char* CASE = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 0.01}},
                    {"number": 2, "class": "Bus", "params": {"N": 1, "kv": 0.01}}],
          "devices": [{"class": "CurrentControl", "id": "CC1", "inputs": {"vpos": 1, "vneg": 2},
                       "params": {"V": 816.496580927726, "mva": 0.001, "Kp": 2.0, "Ki": 4.0, "L": 0.5}}]
        })";

        // The voltages, the frequency, the direction and the orders are held inputs
        constexpr const char* STATE = R"({
          "buses": {"1": {"v": [10.0]}, "2": {"v": [0.0]}},
          "devices": {"CC1": {"xdq": [0.0, 0.0], "mdq": [0.5, 0.0], "vdq": [1.0, 0.0], "omega": [2.0],
                              "axis": [0.6, 0.8], "ip": [0.6], "iq": [0.8]}}
        })";

        const auto system = makeSystem(CASE, STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        // The DC voltage of 10 raises the round-off of the difference, which passes at 1e-10 and not at 1e-11
        constexpr RealT CONTROL_TOL = 1.0e-10;

        auto& control = system->device("CC1");

        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(control.getResidual(EMT::CurrentControlInternalVariables::xdq), {4.0, 0.0});
        success *= rowsEqual(control.getResidual(EMT::CurrentControlInternalVariables::mdq), {0.1, 0.0});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, CONTROL_TOL);

        return success.report(__func__);
      }

      /// A reactive-current order of 0.3 is reduced by twice the voltage above 1.1, and one below the bound passes through
      TestOutcome hvrcm()
      {
        TestStatus success = true;

        constexpr const char* CASE = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [],
          "devices": [{"class": "HVRCM", "id": "HV1", "params": {"Vhvmax": 1.1, "Khv": 2.0, "Qmin": -0.5}}]
        })";

        const auto system = makeSystem(CASE, R"({"devices": {"HV1": {"vt": [1.2], "iqcmd": [0.3]}}})");

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        using HvrcmT = EMT::HVRCM<RealT, IdxT>;

        auto& management = system->device<HvrcmT>("HV1");

        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(management.getResidual(EMT::HVRCMInternalVariables::iq), {0.3 - 2.0 * Math::ramp(1.2 - 1.1)});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, kTol);

        system->held(management, EMT::HVRCMInputs::iqcmd)[0]  = -0.8;
        success                                              *= system->evaluateResidual() == 0;
        success                                              *= rowsEqual(management.getResidual(EMT::HVRCMInternalVariables::iq), {-0.8});

        return success.report(__func__);
      }

      /// At voltage 0.5 a ceiling of 0.4 recovers toward its target and limits an order of 0.9 to the available current
      TestOutcome lvrt()
      {
        TestStatus success = true;

        constexpr const char* CASE = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [],
          "devices": [{"class": "LVRT", "id": "LV1",
                       "params": {"TM": 0.5, "VL0": 0.2, "VL1": 0.8, "IL1": 1.2, "Rpmax": 1.0, "VA0": 0.1, "VA1": 0.6}}]
        })";

        const auto system = makeSystem(CASE, R"({"devices": {"LV1": {"il": [0.4], "vt": [0.5], "ipcmd": [0.9]}}})");

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        auto& limit = system->device("LV1");

        const RealT target    = Math::linseg(0.5, 0.2, 0.8, 1.2);
        const RealT available = 0.4 * Math::linseg(0.5, 0.1, 0.6, 1.0);

        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(limit.getResidual(EMT::LVRTInternalVariables::il), {(target - 0.4) / 0.5});
        success *= rowsEqual(limit.getResidual(EMT::LVRTInternalVariables::ip), {available});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, kTol);

        return success.report(__func__);
      }

      /// Regca at a terminal voltage of 0.75 and angle zero, with orders of 0.5 that its ceiling and its rate limit hold back
      TestOutcome regca()
      {
        TestStatus success = true;

        constexpr const char* CASE = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 3, "kv": 1.0}}],
          "devices": [{"class": "Regca", "id": "GEN", "terminals": {"bus": 1},
                       "params": {"V": 1.0, "mva": 1.0e-6, "Tg": 0.5, "TM": 0.5, "Rqmax": 0.25, "Rqmin": -0.25, "Rpmax": 10.0,
                                  "IL1": 1.0, "VL0": 0.5, "VL1": 1.0, "VA0": 0.25, "VA1": 1.25, "Vhvmax": 1.25, "Khv": 2.0,
                                  "Qmin": -1.0}}]
        })";

        constexpr const char* STATE = R"({"devices": {"GEN": {
          "il": [0.75], "ip": [0.25], "iq": [0.125], "i": [0.5, -0.25, -0.25],
          "vt": [0.75], "ipcmd": [0.5], "iqcmd": [0.5]}}})";

        const auto system = makeSystem(CASE, STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        using RegcaT    = EMT::Regca<RealT, IdxT>;
        using Variables = EMT::RegcaInternalVariables;

        constexpr RealT REL_TOL = 1.0e-3;

        // Peak phase current at the rating, and the reactive current on phases b and c at angle zero
        const RealT peak       = std::sqrt(TWO<RealT> / THREE<RealT>);
        const RealT quadrature = HALF<RealT> * std::numbers::sqrt3_v<RealT> * 0.125;

        auto& converter = system->device<RegcaT>("GEN");

        // The ceiling falls toward 0.5 and lets 0.375 of the active order through. The reactive current would rise at
        // 0.75 and is held to its rate limit.
        success *= system->size() == 13;
        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(converter.getResidual(Variables::il), {-0.5});
        success *= rowsEqual(converter.getResidual(Variables::ip), {0.25});
        success *= rowsEqual(converter.getResidual(Variables::iq), {0.25});
        success *= rowsEqual(converter.getResidual(Variables::i),
                             {-0.5 + peak * 0.25, 0.25 + peak * (-0.125 - quadrature), 0.25 + peak * (-0.125 + quadrature)});
        success *= rowsEqual(system->bus("1").getResidual(EMT::BusInternalVariables::v), {0.5, -0.25, -0.25});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, kTol);

        // The currents into the grid take the tolerance of the system power base at the base voltage of the bus
        constexpr IdxT I  = 6;
        success          *= system->setAbsoluteTolerance(REL_TOL) == 0;
        success          *= rowsEqual(std::span<const RealT>{system->absoluteTolerance().getData() + I, 3}, {1.0, 1.0, 1.0});

        // A rate limit of zero is none
        auto free                              = nlohmann::json::parse(CASE);
        free["devices"][0]["params"]["Rqmax"]  = 0.0;
        const auto unlimited                   = makeSystem(free.dump().c_str(), STATE);
        success                               *= unlimited != nullptr && unlimited->evaluateResidual() == 0;
        if (success)
        {
          success *= rowsEqual(unlimited->device<RegcaT>("GEN").getResidual(Variables::iq), {0.75});
        }

        auto invalid                              = nlohmann::json::parse(CASE);
        invalid["devices"][0]["params"]["Rqmin"]  = 1.0;
        success                                  *= makeSystem(invalid.dump().c_str(), "{}") == nullptr;

        return success.report(__func__);
      }

      /// A bus at (0, 1, -1) has q-axis voltage 2/sqrt(3) at angle zero, and coordinates (3, 4) on a base of 5 have direction (0.6, 0.8)
      TestOutcome park()
      {
        TestStatus success = true;

        constexpr const char* CASE = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 3, "kv": 1.0}}],
          "devices": [{"class": "Park", "id": "PK1", "inputs": {"abc": 1}, "params": {"base": 5.0}}]
        })";

        const auto system = makeSystem(CASE, R"({"buses": {"1": {"v": [0.0, 1.0, -1.0]}}})");

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        constexpr RealT VQ = TWO<RealT> / std::numbers::sqrt3_v<RealT>;

        // Coordinates up to 4 raise the round-off of the difference, which passes at 1e-10 and not at 1e-11
        constexpr RealT PARK_TOL = 1.0e-10;

        using ParkT = EMT::Park<RealT, IdxT>;

        auto& transform = system->device<ParkT>("PK1");

        // With no coordinates yet, the direction is the d axis
        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(transform.getResidual(EMT::ParkInternalVariables::dq), {0.0, VQ});
        success *= rowsEqual(transform.getResidual(EMT::ParkInternalVariables::axis), {1.0, 0.0});

        const auto dq  = transform.y(EMT::ParkInternalVariables::dq);
        dq[0]          = 3.0;
        dq[1]          = 4.0;
        success       *= system->evaluateResidual() == 0;
        success       *= rowsEqual(transform.getResidual(EMT::ParkInternalVariables::dq), {-3.0, VQ - 4.0});
        success       *= rowsEqual(transform.getResidual(EMT::ParkInternalVariables::axis), {0.6, 0.8});
        success       *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, PARK_TOL);

        return success.report(__func__);
      }

      /// A converter whose control and bridge read each other and a Pll: variables in place, the frequency through a signal
      TestOutcome converter()
      {
        TestStatus success = true;

        constexpr const char* CASE = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e4},
          "buses": [{"number": 1, "class": "Bus", "params": {"N": 3, "kv": 0.4}},
                    {"number": 2, "class": "Bus", "params": {"N": 1, "kv": 0.4}},
                    {"number": 3, "class": "Bus", "params": {"N": 1, "kv": 0.4}}],
          "devices": [{"class": "CurrentControl", "id": "CC1",
                       "inputs": {"vpos": 2, "vneg": 3, "vdq": "BR1.vodq", "idq": "BR1.idq", "omega": "PLL1.omega"},
                       "params": {"V": 400.0, "mva": 0.01, "Kp": 4.0, "Ki": 100.0, "L": 0.002}},
                      {"class": "Bridge", "id": "BR1", "terminals": {"bus": 1, "vpos": 2, "vneg": 3},
                       "inputs": {"theta": "PLL1.theta", "mdq": "CC1.mdq"},
                       "params": {"N": 3, "Rs": 0.05, "Ls": 0.002, "C": 2.0e-5, "Rg": 0.02, "Lg": 0.0005, "pwm": {"N": 3, "fc": 10000.0}}},
                      {"class": "Pll", "id": "PLL1", "inputs": {"v": 1}, "params": {"Kp": 0.5, "Ki": 20.0}}]
        })";

        constexpr const char* STATE = R"({
          "buses": {"1": {"v": [320.0, -160.0, -160.0]}, "2": {"v": [400.0]}, "3": {"v": [-400.0]}},
          "devices": {"CC1": {"mdq": [0.8, 0.0], "axis": [1.0, 0.0], "ip": [1.0]},
                      "BR1": {"i": [8.0, -4.0, -4.0], "vo": [320.0, -160.0, -160.0], "ig": [8.0, -4.0, -4.0],
                              "idq": [8.0, 0.0], "vodq": [320.0, 0.0]}}
        })";

        const auto system = makeSystem(CASE, STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        // Voltages of 400 and the frequency base raise the round-off of the difference, which passes at 1e-8 and not at 1e-9
        constexpr RealT CONVERTER_TOL = 1.0e-8;

        // Five bus voltages, four control, thirteen bridge and two Pll variables,
        // the four held inputs `axis`, `ip` and `iq`, and the signal of the frequency
        success *= system->size() == 29;
        success *= system->evaluateResidual() == 0;
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, CONVERTER_TOL);

        return success.report(__func__);
      }

      /// A valve at 0.5 opens toward a reference of 0.375 less a speed deviation of 0.125 over a droop of 0.25
      TestOutcome tgov1()
      {
        TestStatus success = true;

        constexpr const char* CASE = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [],
          "devices": [{"class": "Tgov1", "id": "GOV",
                       "params": {"Trate": 2.0, "S": 1.0e6, "R": 0.25, "T1": 0.5, "T2": 1.0, "T3": 2.0, "Dt": 2.0}}]
        })";

        const auto system = makeSystem(
            CASE, R"({"devices": {"GOV": {"PTX": [0.25], "PV": [0.5], "PM": [0.75], "omega": [1.125], "pref": [0.375]}}})");

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        using Variables = EMT::Tgov1InternalVariables;

        auto& governor = system->device("GOV");

        // The machine rating is half the turbine rating
        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(governor.getResidual(Variables::PTX), {0.125});
        success *= rowsEqual(governor.getResidual(Variables::PV), {Math::antiwindup(0.5, 0.5, 0.0, 1.0) / 0.5});
        success *= rowsEqual(governor.getResidual(Variables::PM), {-0.375});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, kTol);

        auto invalid                          = nlohmann::json::parse(CASE);
        invalid["devices"][0]["params"]["R"]  = 0.0;
        success                              *= makeSystem(invalid.dump().c_str(), "{}") == nullptr;

        return success.report(__func__);
      }

      /// A fuel demand of 0.375 is the lower of the load and temperature demands, with a reference of 1.5 on twice the turbine rating
      TestOutcome gastPti()
      {
        TestStatus success = true;

        constexpr const char* CASE = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [],
          "devices": [{"class": "GastPti", "id": "GOV",
                       "params": {"Trate": 2.0, "S": 1.0e6, "R": 0.25, "T1": 0.5, "T2": 0.25, "T3": 2.0,
                                  "At": 1.0, "Kt": 2.0, "Dturb": 2.0}}]
        })";

        constexpr const char* STATE = R"({"devices": {"GOV": {
          "XVALVE": [0.25], "XFLOW": [0.5], "XTEMP": [0.75], "VLOAD": [0.375], "VTEMP": [2.0], "VLV": [0.5], "PMECH": [0.75],
          "omega": [1.125], "pref": [1.5]}}})";

        const auto system = makeSystem(CASE, STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        using Variables = EMT::GastPtiInternalVariables;

        auto& governor = system->device("GOV");

        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(governor.getResidual(Variables::XVALVE), {Math::antiwindup(0.25, 0.25, 0.0, 1.0) / 0.5});
        success *= rowsEqual(governor.getResidual(Variables::XFLOW), {-1.0});
        success *= rowsEqual(governor.getResidual(Variables::XTEMP), {-0.125});
        success *= rowsEqual(governor.getResidual(Variables::VLOAD), {-0.03125});
        success *= rowsEqual(governor.getResidual(Variables::VTEMP), {-0.5});
        success *= rowsEqual(governor.getResidual(Variables::VLV), {-0.5 + Math::min(0.375, 2.0)});
        success *= rowsEqual(governor.getResidual(Variables::PMECH), {-0.125});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, kTol);

        auto invalid                          = nlohmann::json::parse(CASE);
        invalid["devices"][0]["params"]["R"]  = 0.0;
        success                              *= makeSystem(invalid.dump().c_str(), "{}") == nullptr;

        return success.report(__func__);
      }

      /// A head of 0.75 under a dam head of 1 accelerates the flow, and a flow of 0.5 at a gate power of 0.5 lowers the head
      TestOutcome hygov()
      {
        TestStatus success = true;

        constexpr const char* CASE = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [],
          "devices": [{"class": "Hygov", "id": "GOV",
                       "params": {"Trate": 2.0, "S": 1.0e6, "Rperm": 0.5, "Rtemp": 0.5, "Tr": 2.0, "Tf": 0.5, "Tg": 0.5,
                                  "Velm": 0.25, "Tw": 2.0, "At": 2.0, "Dturb": 0.5, "Qnl": 0.125, "Tn": 0.5, "Tnp": 1.0}}]
        })";

        constexpr const char* STATE = R"({"devices": {"GOV": {
          "XN": [0.25], "XF": [0.25], "C": [0.5], "G": [0.25], "Q": [0.5], "OMEGADB": [0.125], "EF": [0.5], "FC": [0.125],
          "RC": [0.125], "PGV": [0.5], "H": [0.75], "PMECH": [0.5], "omega": [1.125], "pref": [1.0], "paux": [0.5]}}})";

        const auto system = makeSystem(CASE, STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        using Variables = EMT::HygovInternalVariables;

        auto& governor = system->device("GOV");

        // The references add to 0.75 on the turbine rating, and the speed lead-lag output is 0.1875
        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(governor.getResidual(Variables::XN), {-0.125});
        success *= rowsEqual(governor.getResidual(Variables::XF), {0.5});
        success *= rowsEqual(governor.getResidual(Variables::C), {Math::antiwindup(0.5, 0.125, 0.0, 1.0)});
        success *= rowsEqual(governor.getResidual(Variables::G), {0.5});
        success *= rowsEqual(governor.getResidual(Variables::Q), {0.125});
        success *= rowsEqual(governor.getResidual(Variables::EF), {-0.1875});
        success *= rowsEqual(governor.getResidual(Variables::FC), {0.5625});
        success *= rowsEqual(governor.getResidual(Variables::H), {-0.0625});
        success *= rowsEqual(governor.getResidual(Variables::PMECH), {0.296875});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, kTol);

        auto invalid                            = nlohmann::json::parse(CASE);
        invalid["devices"][0]["params"]["db2"]  = 0.125;
        success                                *= makeSystem(invalid.dump().c_str(), "{}") == nullptr;

        return success.report(__func__);
      }

      /// A terminal voltage of 1.25 against a reference of 1.5 and three added signals leaves a voltage error residual of 0.203125
      TestOutcome sexsPti()
      {
        TestStatus success = true;

        constexpr const char* CASE = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [],
          "devices": [{"class": "SexsPti", "id": "EX",
                       "params": {"V": 2.0, "K": 8.0, "Ta": 1.0, "Tb": 2.0, "Te": 0.5, "Efdmin": -4.0, "Efdmax": 4.0}}]
        })";

        constexpr const char* STATE = R"({"devices": {"EX": {
          "VR": [0.25], "EFD": [1.0], "VTR": [0.125],
          "vt": [2.5], "vref": [1.5], "vs": [0.0625], "vuel": [0.03125], "voel": [-0.015625]}}})";

        const auto system = makeSystem(CASE, STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        using Variables = EMT::SexsPtiInternalVariables;

        auto& exciter = system->device("EX");

        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(exciter.getResidual(Variables::VR), {-0.1875});
        success *= rowsEqual(exciter.getResidual(Variables::EFD), {Math::antiwindup(1.0, -3.0, -4.0, 4.0)});
        success *= rowsEqual(exciter.getResidual(Variables::VTR), {0.203125});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, kTol);

        auto invalid = nlohmann::json::parse(CASE);
        invalid["devices"][0]["params"].erase("Tb");
        success *= makeSystem(invalid.dump().c_str(), "{}") == nullptr;

        return success.report(__func__);
      }

      /// Saturation points (1, 0) and (2, 0.5) begin the saturation at 1 with unit gain, and speed 1.125 raises a field voltage of 1.5
      TestOutcome ieeet1()
      {
        TestStatus success = true;

        constexpr const char* CASE = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [],
          "devices": [{"class": "Ieeet1", "id": "EX",
                       "params": {"V": 2.0, "Tr": 0.5, "Ka": 8.0, "Ta": 0.5, "Ke": 0.5, "Te": 0.5, "Kf": 0.125, "Tf": 2.0,
                                  "Vrmin": -4.0, "Vrmax": 4.0, "E1": 1.0, "Se1": 0.0, "E2": 2.0, "Se2": 0.5,
                                  "Ispdlim": 1.0}}]
        })";

        constexpr const char* STATE = R"({"devices": {"EX": {
          "VTS": [1.0], "VR": [1.0], "EFDP": [1.5], "VFX": [0.25], "VTR": [0.25], "VF": [0.0625], "VE": [0.125],
          "EFD": [1.25], "KSAT": [0.5],
          "vt": [2.5], "omega": [1.125], "vref": [1.5], "vs": [0.0625], "vuel": [0.03125], "voel": [-0.015625]}}})";

        const auto system = makeSystem(CASE, STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        using Variables = EMT::Ieeet1InternalVariables;

        auto& exciter = system->device("EX");

        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(exciter.getResidual(Variables::VTS), {0.5});
        success *= rowsEqual(exciter.getResidual(Variables::VR), {Math::antiwindup(1.0, 2.0, -4.0, 4.0)});
        success *= rowsEqual(exciter.getResidual(Variables::EFDP), {0.25});
        success *= rowsEqual(exciter.getResidual(Variables::VFX), {0.03125});
        success *= rowsEqual(exciter.getResidual(Variables::VTR), {0.265625});
        success *= rowsEqual(exciter.getResidual(Variables::VF), {-0.4375});
        success *= rowsEqual(exciter.getResidual(Variables::VE), {0.375});
        success *= rowsEqual(exciter.getResidual(Variables::EFD), {0.4375});
        success *= rowsEqual(exciter.getResidual(Variables::KSAT), {-0.5 + Math::qramp(0.5)});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, kTol);

        auto invalid                           = nlohmann::json::parse(CASE);
        invalid["devices"][0]["params"]["Ke"]  = 0.0;
        success                               *= makeSystem(invalid.dump().c_str(), "{}") == nullptr;

        return success.report(__func__);
      }

      /// With `UEL` 2 the limiter signal adds to the voltage error and the lead-lag output passes the gate
      TestOutcome esdc1a()
      {
        TestStatus success = true;

        constexpr const char* CASE = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [],
          "devices": [{"class": "Esdc1a", "id": "EX",
                       "params": {"V": 2.0, "Tr": 0.5, "Ka": 8.0, "Ta": 0.5, "Tb": 2.0, "Tc": 1.0, "Vrmin": -4.0,
                                  "Vrmax": 4.0, "Ke": 0.5, "Te": 0.5, "Kf": 0.25, "Tf1": 2.0, "Spdmlt": true, "E1": 1.0,
                                  "Se1": 0.0, "E2": 2.0, "Se2": 0.5, "UEL": 2, "exclim": false}}]
        })";

        constexpr const char* STATE = R"({"devices": {"EX": {
          "EFDP": [1.5], "VC": [1.0], "VR": [1.0], "VF": [0.0625], "XLL": [0.125], "EV": [0.25], "VLL": [0.125],
          "VHV": [0.25], "SE": [0.5], "VFE": [0.75], "EFD": [1.25],
          "vt": [2.5], "omega": [1.125], "vref": [1.5], "vs": [0.0625], "vuel": [0.03125]}}})";

        const auto system = makeSystem(CASE, STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        using Variables = EMT::Esdc1aInternalVariables;

        auto& exciter = system->device("EX");

        // The field-voltage state rises at (1 - 0.75) / 0.5 without its lower limit
        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(exciter.getResidual(Variables::EFDP), {0.5});
        success *= rowsEqual(exciter.getResidual(Variables::VC), {0.5});
        success *= rowsEqual(exciter.getResidual(Variables::VR), {Math::antiwindup(1.0, 1.0, -4.0, 4.0) / 0.5});
        success *= rowsEqual(exciter.getResidual(Variables::VF), {0.03125});
        success *= rowsEqual(exciter.getResidual(Variables::XLL), {0.0625});
        success *= rowsEqual(exciter.getResidual(Variables::EV), {0.28125});
        success *= rowsEqual(exciter.getResidual(Variables::VLL), {0.0625});
        success *= rowsEqual(exciter.getResidual(Variables::VHV), {-0.125});
        success *= rowsEqual(exciter.getResidual(Variables::SE), {-0.5 + Math::qramp(0.5)});
        success *= rowsEqual(exciter.getResidual(Variables::VFE), {0.5});
        success *= rowsEqual(exciter.getResidual(Variables::EFD), {0.4375});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, kTol);

        auto invalid                            = nlohmann::json::parse(CASE);
        invalid["devices"][0]["params"]["UEL"]  = 4;
        success                                *= makeSystem(invalid.dump().c_str(), "{}") == nullptr;

        return success.report(__func__);
      }

      /// A second-order notch filter 1 + 0.5 s + 0.25 s^2 with input 0.5, two lead-lags, a washout and limits of 0.125
      TestOutcome ieeest()
      {
        TestStatus success = true;

        constexpr const char* CASE = R"({
          "params": {"freq_base": 60.0, "va_base": 1.0e6},
          "buses": [],
          "devices": [{"class": "Ieeest", "id": "PSS",
                       "params": {"A1": 0.5, "A2": 0.25, "A5": 1.0, "A6": 0.5, "T1": 1.0, "T2": 2.0, "T3": 1.0, "T4": 2.0,
                                  "T5": 2.0, "T6": 4.0, "Ks": 2.0, "Lsmin": -0.125, "Lsmax": 0.125}}]
        })";

        constexpr const char* STATE = R"({"devices": {"PSS": {
          "X1": [0.25], "X2": [0.125], "X5": [0.25], "X6": [0.125], "X7": [0.0625], "V4": [0.5], "V5": [0.3125],
          "V6": [0.25], "V7": [0.0625], "VSS": [0.03125], "u": [0.5]}}})";

        const auto system = makeSystem(CASE, STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        using Variables = EMT::IeeestInternalVariables;

        auto& stabilizer = system->device("PSS");

        // The second derivative of the filtered signal is (0.5 - 0.25 - 0.5 * 0.125) / 0.25
        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(stabilizer.getResidual(Variables::X1), {0.125});
        success *= rowsEqual(stabilizer.getResidual(Variables::X2), {0.75});
        success *= rowsEqual(stabilizer.getResidual(Variables::X3), {0.0});
        success *= rowsEqual(stabilizer.getResidual(Variables::X4), {0.0});
        success *= rowsEqual(stabilizer.getResidual(Variables::X5), {0.25});
        success *= rowsEqual(stabilizer.getResidual(Variables::X6), {0.1875});
        success *= rowsEqual(stabilizer.getResidual(Variables::X7), {0.1875});
        success *= rowsEqual(stabilizer.getResidual(Variables::V4), {0.25});
        success *= rowsEqual(stabilizer.getResidual(Variables::V5), {0.125});
        success *= rowsEqual(stabilizer.getResidual(Variables::V6), {-0.0625});
        success *= rowsEqual(stabilizer.getResidual(Variables::V7), {0.5});
        success *= rowsEqual(stabilizer.getResidual(Variables::VSS), {-0.03125 + Math::clamp(0.0625, -0.125, 0.125)});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, kTol);

        auto invalid                               = nlohmann::json::parse(CASE);
        invalid["devices"][0]["params"]["Tdelay"]  = 0.125;
        success                                   *= makeSystem(invalid.dump().c_str(), "{}") == nullptr;

        return success.report(__func__);
      }

      /// A terminal voltage of 1, within the voltage band, with an active power of 0.5625 and a reactive power of -0.28125
      TestOutcome reecb()
      {
        TestStatus success = true;

        constexpr const char* CASE = R"({
          "params": {"freq_base": 60.0, "va_base": 2.0},
          "buses": [],
          "devices": [{"class": "Reecb", "id": "REEC",
                       "params": {"V": 1.0, "mva": 1.0e-6, "Vref0": 1.0, "Trv": 0.5, "Tp": 0.5, "Tiq": 0.5, "Tpord": 0.5,
                                  "Imax": 1.0}}]
        })";

        constexpr const char* STATE = R"({"devices": {"REEC": {
          "VMEAS": [0.75], "PMEAS": [0.25], "QV": [0.25], "PORD": [0.5], "VSAFE": [0.75], "SDIP": [1.0],
          "IQV": [0.125], "QREF": [0.125], "EQ": [0.125], "VPIQ": [1.0], "RPORD": [0.25], "ILCAP": [0.75],
          "IQMAX": [0.875], "IPMAX": [0.5], "IQRAW": [0.25], "IQCMD": [0.5], "IPCMD": [0.5],
          "qext": [0.125], "pref": [0.375]}}})";

        const auto system = makeSystem(CASE, STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        using ReecbT    = EMT::Reecb<RealT, IdxT>;
        using Variables = EMT::ReecbInternalVariables;

        auto& control = system->device<ReecbT>("REEC");

        system->held(control, EMT::ReecbInputs::vt)[0] = ONE<RealT>;
        system->held(control, EMT::ReecbInputs::p)[0]  = 0.5625;
        system->held(control, EMT::ReecbInputs::q)[0]  = -0.28125;

        const auto near = [&control](Variables variable, RealT value)
        { return isEqual(control.getResidual(variable)[0], value, kTol); };

        // The system power base is twice the device rating, and the reactive current has priority
        success *= system->evaluateResidual() == 0;
        success *= near(Variables::VMEAS, 0.5);
        success *= near(Variables::PMEAS, 0.625);
        success *= near(Variables::QREF, 0.125);
        success *= near(Variables::EQ, 0.28125);
        success *= near(Variables::IQMAX, 0.125);
        success *= near(Variables::IPMAX, 0.25);
        success *= near(Variables::IQRAW, 0.125);
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, PLANT_TOL);

        auto invalid = nlohmann::json::parse(CASE);
        invalid["devices"][0]["params"].erase("Vref0");
        success *= makeSystem(invalid.dump().c_str(), "{}") == nullptr;

        return success.report(__func__);
      }

      /// A line drop of 1.28 + j 0.96 turns a unit voltage with powers (0.5, -0.25) into a compensated voltage of magnitude 1
      TestOutcome repca()
      {
        TestStatus success = true;

        constexpr const char* CASE = R"({
          "params": {"freq_base": 60.0, "va_base": 2.0},
          "buses": [],
          "devices": [{"class": "Repca", "id": "PLANT",
                       "params": {"V": 1.0, "mva": 1.0e-6, "VcompFlag": false, "Freqflag": true, "Tfltr": 0.5, "Tfv": 2.0,
                                  "Tft": 1.0, "Tp": 0.5, "Tlag": 2.0, "Kc": 0.5, "Rc": 1.28, "Xc": 0.96}}]
        })";

        constexpr const char* STATE = R"({"devices": {"PLANT": {
          "VMEAS": [0.75], "QMEAS": [-0.5], "XQPI": [0.125], "XQLAG": [0.25], "PMEAS": [0.25], "XPPI": [0.25],
          "PREF": [0.5], "VLDC": [0.75], "VDROOP": [0.5], "VCTRL": [0.375], "SFRZ": [1.0], "ERQ": [0.125],
          "ERQDB": [0.125], "ERQLIM": [0.125], "QPI": [0.75], "QEXT": [0.1875], "EP": [0.125], "EPLIM": [0.125],
          "PPI": [1.0], "PEXT": [0.125],
          "vref": [1.0], "pref": [0.25], "freqref": [1.0]}}})";

        const auto system = makeSystem(CASE, STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        using RepcaT    = EMT::Repca<RealT, IdxT>;
        using Variables = EMT::RepcaInternalVariables;

        auto& plant = system->device<RepcaT>("PLANT");

        system->held(plant, EMT::RepcaInputs::v)[0]     = ONE<RealT>;
        system->held(plant, EMT::RepcaInputs::p)[0]     = 0.5;
        system->held(plant, EMT::RepcaInputs::q)[0]     = -0.25;
        system->held(plant, EMT::RepcaInputs::omega)[0] = TWO<RealT> * std::numbers::pi_v<RealT> * 60.0;

        const auto near = [&plant](Variables variable, RealT value)
        { return isEqual(plant.getResidual(variable)[0], value, kTol); };

        // The powers are 0.5 and -0.25, the frequency is at its reference, and the commands are on twice the device rating
        success *= system->evaluateResidual() == 0;
        success *= near(Variables::VMEAS, -0.75);
        success *= near(Variables::QMEAS, 0.5);
        success *= near(Variables::XQLAG, 0.25);
        success *= near(Variables::PMEAS, 0.5);
        success *= near(Variables::PREF, 0.25);
        success *= near(Variables::VLDC, 0.25);
        success *= near(Variables::VDROOP, 0.375);
        success *= near(Variables::VCTRL, 0.125);
        success *= near(Variables::ERQ, 0.125);
        success *= near(Variables::QEXT, 0.25);
        success *= near(Variables::EF, 0.0);
        success *= near(Variables::EP, 0.125);
        success *= near(Variables::PEXT, 0.25);

        // The quotient in the line drop raises the truncation of the difference, which passes at 1e-9 and not at 1e-10
        constexpr RealT REPCA_TOL  = 1.0e-9;
        success                   *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, REPCA_TOL);

        auto invalid                             = nlohmann::json::parse(CASE);
        invalid["devices"][0]["params"]["Qmin"]  = 2.0;
        success                                 *= makeSystem(invalid.dump().c_str(), "{}") == nullptr;

        return success.report(__func__);
      }

      /// A machine reads its field voltage and mechanical power from an exciter and a governor, and the exciter the voltage of a Pmu
      TestOutcome machineControls()
      {
        TestStatus success = true;

        auto read                    = nlohmann::json::parse(MACHINE_CASE);
        read["devices"][0]["inputs"] = {{"efd", "EX.EFD"}, {"pm", "GOV.PM"}};
        read["devices"].push_back(nlohmann::json::parse(
            R"({"class": "Tgov1", "id": "GOV", "inputs": {"omega": "M1.omega"}, "params": {"Trate": 1.0e-6, "S": 1.0}})"));
        read["devices"].push_back(nlohmann::json::parse(
            R"({"class": "SexsPti", "id": "EX", "inputs": {"vt": "PMU.vrms"},
                "params": {"V": 1.0, "K": 8.0, "Ta": 1.0, "Tb": 2.0, "Te": 0.5, "Efdmin": -4.0, "Efdmax": 4.0}})"));
        read["devices"].push_back(nlohmann::json::parse(
            R"({"class": "Pmu", "id": "PMU", "inputs": {"v": 1}, "params": {"pll": {"Kp": 0.5, "Ki": 0.25}}})"));

        constexpr const char* STATE = R"({
          "buses": {"1": {"v": [0.0, 1.0, -1.0]}},
          "devices": {"M1": {"omega": [1.0], "psid": [0.9], "psiq": [-1.5], "psifd": [1.4], "psi1d": [0.96],
                             "psi1q": [-1.5], "psi2q": [-1.7125], "psiad": [1.0], "psiaq": [-1.3]},
                      "GOV": {"PTX": [0.5], "PV": [0.5], "PM": [1.0]},
                      "EX": {"EFD": [2.0], "vref": [1.5]},
                      "PMU": {"vms": [2.0]}}
        })";

        const auto system = makeSystem(read.dump().c_str(), STATE);

        success *= system != nullptr;
        if (!success)
        {
          return success.report(__func__);
        }

        // Three bus voltages, eleven machine, three governor, three exciter, three Pmu and two Pll variables, eight held
        // inputs, and the signal of `vrms`
        constexpr IdxT SIGNAL  = 33;
        success               *= system->size() == 34;

        // The machine accelerates as with a held mechanical power of 1, and the signal starts at zero
        success *= system->evaluateResidual() == 0;
        success *= rowsEqual(system->device("M1").getResidual(EMT::MachineInternalVariables::omega), {-0.34});
        const std::span<const RealT> signal{system->getResidual().getData() + SIGNAL, 1};
        success *= rowsEqual(signal, {std::numbers::sqrt2_v<RealT>});
        success *= jacobianMatchesResidual(*system, ZERO<RealT>, ALPHA, MACHINE_TOL);

        return success.report(__func__);
      }

    private:
      /// Bridge of unit filter parameters between a three-phase bus and DC terminals, with its angle and modulation held
      static constexpr const char* BRIDGE_CASE = R"({
        "params": {"freq_base": 60.0, "va_base": 1.0e6},
        "buses": [{"number": 1, "class": "Bus", "params": {"N": 3, "kv": 0.001}},
                  {"number": 2, "class": "Bus", "params": {"N": 1, "kv": 0.004}},
                  {"number": 3, "class": "Bus", "params": {"N": 1, "kv": 0.004}}],
        "devices": [{"class": "Bridge", "id": "BR1", "terminals": {"bus": 1, "vpos": 2, "vneg": 3},
                     "params": {"N": 3, "Rs": 1.0, "Ls": 1.0, "C": 1.0, "Rg": 1.0, "Lg": 1.0,
                                "pwm": {"N": 3, "fc": 1000.0, "injection": false}}}]
      })";

      /// Angle zero and modulation (0.5, 0) between DC terminals at 4 and 0
      static constexpr const char* BRIDGE_STATE = R"({
        "buses": {"1": {"v": [1.0, 0.0, 0.0]}, "2": {"v": [4.0]}, "3": {"v": [0.0]}},
        "devices": {"BR1": {"i": [1.0, 0.0, 0.0], "vo": [2.0, 0.0, 0.0], "ig": [0.5, 0.0, 0.0], "mdq": [0.5, 0.0]}}
      })";

      /// Machine of unit ratings on a three-phase bus, with its field voltage and mechanical power held
      static constexpr const char* MACHINE_CASE = R"({
        "params": {"freq_base": 60.0, "va_base": 1.0e6},
        "buses": [{"number": 1, "class": "Bus", "params": {"N": 3, "kv": 0.001}}],
        "devices": [{"class": "Machine", "id": "M1", "terminals": {"bus": 1},
                     "params": {"S": 1.0, "V": 1.0, "f": 60.0, "H": 2.0, "Ll": 0.5, "Lmd": 2.0, "Lmq": 1.0,
                                "Rfd": 0.5, "Llfd": 0.5, "R1d": 1.0, "Ll1d": 0.5, "R1q": 1.0, "Ll1q": 0.5,
                                "R2q": 1.0, "Ll2q": 0.5, "S10": 0.04, "S12": 0.16}}]
      })";

      /// Angle zero and synchronous speed on a dead bus, with stator currents 0.2 and 0.4 in d and q
      static constexpr const char* MACHINE_STATE = R"({
        "buses": {"1": {"v": [0.0, 0.0, 0.0]}},
        "devices": {"M1": {"omega": [1.0], "psid": [0.9], "psiq": [-1.5], "psifd": [1.4], "psi1d": [0.96],
                           "psi1q": [-1.5], "psi2q": [-1.7125], "psiad": [1.0], "psiaq": [-1.3], "efd": [2.0], "pm": [1.0]}}
      })";

      /// Switch between two single-phase buses, with its command left open
      static constexpr const char* SWITCH_CASE = R"({
        "params": {"freq_base": 60.0, "va_base": 1.0e6},
        "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}},
                  {"number": 2, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
        "devices": [{"class": "Switch", "id": "S1", "terminals": {"bus1": 1, "bus2": 2}, "params": {"N": 1}}]
      })";

      /// Pll on a three-phase bus
      static constexpr const char* PLL_CASE = R"({
        "params": {"freq_base": 60.0, "va_base": 1.0e6},
        "buses": [{"number": 1, "class": "Bus", "params": {"N": 3, "kv": 1.0}}],
        "devices": [{"class": "Pll", "id": "PLL1", "inputs": {"v": 1}, "params": {"Kp": 0.5, "Ki": 0.25}}]
      })";

      static constexpr const char* PMU_CASE = R"({
        "params": {"freq_base": 60.0, "va_base": 1.0e6},
        "buses": [{"number": 1, "class": "Bus", "params": {"N": 3, "kv": 2.0}}],
        "devices": [{"class": "Pmu", "id": "PMU1", "inputs": {"v": 1}, "params": {"pll": {"Kp": 0.5, "Ki": 0.25}}}]
      })";

      /// Mean-square voltage 4 on a bus at (0, 1, -1) with held currents (3, 0, -3), and the Pll at angle 0.25
      static constexpr const char* PMU_STATE = R"({
        "buses": {"1": {"v": [0.0, 1.0, -1.0]}},
        "devices": {"PMU1": {"vms": [4.0], "i": [3.0, 0.0, -3.0], "pll": {"theta": [0.25], "x": [0.5]}}}
      })";

      /// Angle zero on a bus at (0, 1, -1), whose q-axis voltage is 2/sqrt(3)
      static constexpr const char* PLL_STATE = R"({
        "buses": {"1": {"v": [0.0, 1.0, -1.0]}},
        "devices": {"PLL1": {"x": [0.5]}}
      })";

      /// Operator of gain 0.5 that reads the voltage of a single-phase bus and adds its output to the current balance
      static constexpr const char* ADMITTANCE_CASE = R"({
        "params": {"freq_base": 60.0, "va_base": 1.0e6},
        "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
        "devices": [{"class": "VectorFit", "id": "Y1", "inputs": {"u": 1}, "outputs": {"y": 1}, "params": {"D": [[0.5]]}}]
      })";

      /// Load of impedance 2 + 8/(s + 4) on a single-phase bus of base voltage 1 kV, with a power base of 1 MVA
      static constexpr const char* LOAD_CASE = R"({
        "params": {"freq_base": 60.0, "va_base": 1.0e6},
        "buses": [{"number": 1, "class": "Bus", "params": {"N": 1, "kv": 1.0}}],
        "devices": [{"class": "LoadZ", "id": "LD1", "terminals": {"bus": 1},
                     "params": {"N": 1, "Z": {"D": [[2.0]], "poles": [[-4.0, 0.0]], "residues": [[[[8.0, 0.0]]]]}}}]
      })";

      static constexpr RealT ALPHA = 0.5;
      static constexpr RealT kTol  = 1.0e-11;

      /// The frequency base of a machine raises the round-off of the difference, which passes at 1e-8 and not at 1e-9
      static constexpr RealT MACHINE_TOL = 1.0e-8;

      /// The per-unit bases of Reecb raise the round-off of the difference, which passes at 1e-10 and not at 1e-11
      static constexpr RealT PLANT_TOL = 1.0e-10;

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
