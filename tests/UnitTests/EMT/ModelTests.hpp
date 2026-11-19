/**
 * @file ModelTests.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Residuals and Jacobians of EMT models on small cases.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <format>
#include <initializer_list>
#include <memory>
#include <numbers>
#include <span>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/ComponentImpl.hpp>
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

    private:
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
