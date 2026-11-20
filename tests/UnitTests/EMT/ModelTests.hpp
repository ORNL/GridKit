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

    private:
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
