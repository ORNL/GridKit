/**
 * @file SystemTests.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief EMT system model tests: JSON parsing, assembly, and time integration.
 *
 */
#pragma once

#include <cmath>
#include <complex>
#include <numbers>
#include <sstream>
#include <string>

#include <GridKit/Definitions.hpp>
#include <GridKit/Model/EMT/Component/Switch/Switch.hpp>
#include <GridKit/Model/EMT/SystemModel.hpp>
#include <GridKit/Model/EMT/SystemModelData.hpp>
#include <GridKit/Testing/Testing.hpp>

#ifdef GRIDKIT_ENABLE_SUNDIALS
#include <GridKit/Solver/Dynamic/Ida.hpp>
#endif

namespace GridKit
{
  namespace Testing
  {
    /**
     * @brief Tests for the EMT system model on a one-bus source and load case.
     */
    template <typename ScalarT, typename IdxT>
    class SystemTests
    {
    public:
      using RealT = ScalarT;

      SystemTests()  = default;
      ~SystemTests() = default;

    private:
      /// One bus, one sinusoidal source behind a series R-L, one R-L load.
      /// The per-phase matrices are diagonal so each phase solves the same
      /// scalar phasor circuit.
      static std::string caseJson()
      {
        return R"({
          "header": {
            "case_name": "EMT two-device smoke case",
            "case_description": "One bus with a voltage source and an impedance load",
            "case_comments": "Used by EMT SystemTests"
          },
          "devices": [
            { "class": "Bus", "id": "bus_1" },
            {
              "class": "VoltageSource",
              "id": "source_1",
              "params": {
                "E": [100.0, 100.0, 100.0],
                "phi": [0.0, -2.0943951023931953, 2.0943951023931953],
                "omega": 376.99111843077515,
                "Rs": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                "Ls": [[0.01, 0.0, 0.0], [0.0, 0.01, 0.0], [0.0, 0.0, 0.01]]
              },
              "inputs": { "bus": "bus_1" }
            },
            {
              "class": "LoadZ",
              "id": "load_1",
              "params": {
                "R": [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
                "L": [[0.04, 0.0, 0.0], [0.0, 0.04, 0.0], [0.0, 0.0, 0.04]]
              },
              "inputs": { "bus": "bus_1" }
            }
          ]
        })";
      }

      /// Source at bus 1, lumped line from bus 1 to bus 2, load at bus 2.
      /// The per-phase matrices are diagonal so each phase solves the same
      /// scalar phasor ladder circuit.
      static std::string threeBusCaseJson()
      {
        return R"({
          "header": {
            "case_name": "EMT source-line-load case",
            "case_description": "Voltage source, lumped line, and impedance load",
            "case_comments": "Used by EMT SystemTests"
          },
          "devices": [
            { "class": "Bus", "id": "bus_1" },
            { "class": "Bus", "id": "bus_2" },
            {
              "class": "VoltageSource",
              "id": "source_1",
              "params": {
                "E": [100.0, 100.0, 100.0],
                "phi": [0.0, -2.0943951023931953, 2.0943951023931953],
                "omega": 376.99111843077515,
                "Rs": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                "Ls": [[0.01, 0.0, 0.0], [0.0, 0.01, 0.0], [0.0, 0.0, 0.01]]
              },
              "inputs": { "bus": "bus_1" }
            },
            {
              "class": "LineLumped",
              "id": "line_1_2",
              "params": {
                "N": 3,
                "K": 3,
                "conductors": [1, 2, 3],
                "dx": 1.0,
                "Rp": [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]],
                "Lp": [[0.02, 0.0, 0.0], [0.0, 0.02, 0.0], [0.0, 0.0, 0.02]],
                "Gp": [[1.0e-4, 0.0, 0.0], [0.0, 1.0e-4, 0.0], [0.0, 0.0, 1.0e-4]],
                "Cp": [[1.0e-5, 0.0, 0.0], [0.0, 1.0e-5, 0.0], [0.0, 0.0, 1.0e-5]]
              },
              "inputs": { "bus1": "bus_1", "bus2": "bus_2" }
            },
            {
              "class": "LoadZ",
              "id": "load_2",
              "params": {
                "R": [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
                "L": [[0.04, 0.0, 0.0], [0.0, 0.04, 0.0], [0.0, 0.0, 0.04]]
              },
              "inputs": { "bus": "bus_2" }
            }
          ]
        })";
      }

      /// Source at bus 1, initially open switch from bus 1 to bus 2, load at
      /// bus 2.
      static std::string switchCaseJson()
      {
        return R"({
          "header": {
            "case_name": "EMT switch energization case",
            "case_description": "Voltage source energizing an impedance load through a switch",
            "case_comments": "Used by EMT SystemTests"
          },
          "devices": [
            { "class": "Bus", "id": "bus_1" },
            { "class": "Bus", "id": "bus_2" },
            {
              "class": "VoltageSource",
              "id": "source_1",
              "params": {
                "E": [100.0, 100.0, 100.0],
                "phi": [0.0, -2.0943951023931953, 2.0943951023931953],
                "omega": 376.99111843077515,
                "Rs": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                "Ls": [[0.01, 0.0, 0.0], [0.0, 0.01, 0.0], [0.0, 0.0, 0.01]]
              },
              "inputs": { "bus": "bus_1" }
            },
            {
              "class": "Switch",
              "id": "switch_1_2",
              "params": { "open": true },
              "inputs": { "bus1": "bus_1", "bus2": "bus_2" }
            },
            {
              "class": "LoadZ",
              "id": "load_1",
              "params": {
                "R": [[100.0, 0.0, 0.0], [0.0, 100.0, 0.0], [0.0, 0.0, 100.0]]
              },
              "inputs": { "bus": "bus_1" }
            },
            {
              "class": "LoadZ",
              "id": "load_2",
              "params": {
                "R": [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
              },
              "inputs": { "bus": "bus_2" }
            }
          ]
        })";
      }

      /// Synchronous machine at bus 1 serving a purely resistive load. The
      /// bus initial voltage and machine power schedule are consistent with
      /// the load, so the initialized trajectory is an exact steady state.
      static std::string machineCaseJson()
      {
        return R"({
          "header": {
            "case_name": "EMT machine flat-start case",
            "case_description": "Synchronous machine serving a resistive load",
            "case_comments": "Used by EMT SystemTests"
          },
          "devices": [
            {
              "class": "Bus",
              "id": "bus_1",
              "init": {
                "va": 11267.65281680262,
                "vb": -5633.82640840131,
                "vc": -5633.82640840131
              }
            },
            {
              "class": "Machine",
              "id": "machine_1",
              "params": {
                "N": 3,
                "S": 100.0e6,
                "V": 13800.0,
                "f": 60.0,
                "H": 3.7,
                "F": 0.0,
                "Rs": 0.003,
                "Ll": 0.15,
                "Lmd": 1.66,
                "Lmq": 1.61,
                "L0": 0.15,
                "Rfd": 0.0006,
                "Llfd": 0.165,
                "R1d": 0.0284,
                "Ll1d": 0.1713,
                "R1q": 0.0062,
                "Ll1q": 0.7252,
                "R2q": 0.0237,
                "Ll2q": 0.125,
                "S10": 0.1,
                "S12": 0.5,
                "p0": 50.0e6,
                "q0": 0.0
              },
              "inputs": { "bus": "bus_1" }
            },
            {
              "class": "LoadZ",
              "id": "load_1",
              "params": {
                "R": [[3.8088, 0.0, 0.0], [0.0, 3.8088, 0.0], [0.0, 0.0, 3.8088]]
              },
              "inputs": { "bus": "bus_1" }
            }
          ]
        })";
      }

      /// The machine flat-start case with a TGOV1 governor wired through
      /// signals: machine speed to the governor, governor mechanical
      /// power back to the machine.
      static std::string machineGovernorCaseJson()
      {
        return R"({
          "header": {
            "case_name": "EMT governed machine flat-start case",
            "case_description": "Synchronous machine with a TGOV1 governor serving a resistive load",
            "case_comments": "Used by EMT SystemTests"
          },
          "signals": [
            { "id": "speed_1" },
            { "id": "pmech_1" }
          ],
          "devices": [
            {
              "class": "Bus",
              "id": "bus_1",
              "init": {
                "va": 11267.65281680262,
                "vb": -5633.82640840131,
                "vc": -5633.82640840131
              }
            },
            {
              "class": "Machine",
              "id": "machine_1",
              "params": {
                "N": 3,
                "S": 100.0e6,
                "V": 13800.0,
                "f": 60.0,
                "H": 3.7,
                "F": 0.0,
                "Rs": 0.003,
                "Ll": 0.15,
                "Lmd": 1.66,
                "Lmq": 1.61,
                "L0": 0.15,
                "Rfd": 0.0006,
                "Llfd": 0.165,
                "R1d": 0.0284,
                "Ll1d": 0.1713,
                "R1q": 0.0062,
                "Ll1q": 0.7252,
                "R2q": 0.0237,
                "Ll2q": 0.125,
                "S10": 0.1,
                "S12": 0.5,
                "p0": 50.0e6,
                "q0": 0.0
              },
              "inputs": { "bus": "bus_1", "pm": "pmech_1" },
              "outputs": { "speed": "speed_1" }
            },
            {
              "class": "Tgov1",
              "id": "governor_1",
              "params": {
                "Trate": 100.0,
                "R": 0.05,
                "T1": 0.5,
                "T2": 2.5,
                "T3": 7.5,
                "Pvmax": 1.0,
                "Pvmin": 0.0,
                "Dt": 0.0
              },
              "inputs": { "speed": "speed_1" },
              "outputs": { "pmech": "pmech_1" }
            },
            {
              "class": "LoadZ",
              "id": "load_1",
              "params": {
                "R": [[3.8088, 0.0, 0.0], [0.0, 3.8088, 0.0], [0.0, 0.0, 3.8088]]
              },
              "inputs": { "bus": "bus_1" }
            }
          ]
        })";
      }

      /// The twin-circuit pair: a parallel R-L branch is exactly a
      /// one-real-pole rational impedance, Z(s) = R || sL
      /// = R - (R^2/L)/(s + R/L). With R = 10 and L = 0.04: D = 10 I,
      /// p = -250, A = -2500 I.
      static std::string twinPrimitiveCaseJson()
      {
        return R"({
          "header": {
            "case_name": "EMT twin circuit, primitive form",
            "case_description": "Parallel R and L loads built from primitives",
            "case_comments": "Used by EMT SystemTests"
          },
          "devices": [
            { "class": "Bus", "id": "bus_1" },
            {
              "class": "VoltageSource",
              "id": "source_1",
              "params": {
                "E": [100.0, 100.0, 100.0],
                "phi": [0.0, -2.0943951023931953, 2.0943951023931953],
                "omega": 376.99111843077515,
                "Rs": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                "Ls": [[0.01, 0.0, 0.0], [0.0, 0.01, 0.0], [0.0, 0.0, 0.01]]
              },
              "inputs": { "bus": "bus_1" }
            },
            {
              "class": "LoadZ",
              "id": "load_r",
              "params": {
                "R": [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
              },
              "inputs": { "bus": "bus_1" }
            },
            {
              "class": "LoadZ",
              "id": "load_l",
              "params": {
                "L": [[0.04, 0.0, 0.0], [0.0, 0.04, 0.0], [0.0, 0.0, 0.04]]
              },
              "inputs": { "bus": "bus_1" }
            }
          ]
        })";
      }

      static std::string twinRationalCaseJson()
      {
        return R"({
          "header": {
            "case_name": "EMT twin circuit, rational form",
            "case_description": "The parallel R and L branch as a rational impedance",
            "case_comments": "Used by EMT SystemTests"
          },
          "devices": [
            { "class": "Bus", "id": "bus_1" },
            {
              "class": "VoltageSource",
              "id": "source_1",
              "params": {
                "E": [100.0, 100.0, 100.0],
                "phi": [0.0, -2.0943951023931953, 2.0943951023931953],
                "omega": 376.99111843077515,
                "Rs": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                "Ls": [[0.01, 0.0, 0.0], [0.0, 0.01, 0.0], [0.0, 0.0, 0.01]]
              },
              "inputs": { "bus": "bus_1" }
            },
            {
              "class": "LoadZ",
              "id": "load_rl",
              "params": {},
              "inputs": { "bus": "bus_1" },
              "submodels": {
                "Z": {
                  "D": [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
                  "poles": [[-250.0, 0.0]],
                  "residues": [
                    [[[-2500.0, 0.0], [0.0, 0.0], [0.0, 0.0]],
                     [[0.0, 0.0], [-2500.0, 0.0], [0.0, 0.0]],
                     [[0.0, 0.0], [0.0, 0.0], [-2500.0, 0.0]]]
                  ]
                }
              }
            }
          ]
        })";
      }

      static std::string twinLineMatrixCaseJson()
      {
        return R"({
          "header": {
            "case_name": "EMT twin line, matrix form",
            "case_description": "Matrix-parameter line with explicit series-RL shunt loads at each terminal",
            "case_comments": "Used by EMT SystemTests"
          },
          "devices": [
            { "class": "Bus", "id": "bus_1" },
            { "class": "Bus", "id": "bus_2" },
            {
              "class": "VoltageSource",
              "id": "source_1",
              "params": {
                "E": [100.0, 100.0, 100.0],
                "phi": [0.0, -2.0943951023931953, 2.0943951023931953],
                "omega": 376.99111843077515,
                "Rs": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                "Ls": [[0.01, 0.0, 0.0], [0.0, 0.01, 0.0], [0.0, 0.0, 0.01]]
              },
              "inputs": { "bus": "bus_1" }
            },
            {
              "class": "LineLumped",
              "id": "line_1_2",
              "params": {
                "N": 3,
                "K": 3,
                "conductors": [1, 2, 3],
                "dx": 1.0,
                "Rp": [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]],
                "Lp": [[0.02, 0.0, 0.0], [0.0, 0.02, 0.0], [0.0, 0.0, 0.02]],
                "Gp": [[1.0e-4, 0.0, 0.0], [0.0, 1.0e-4, 0.0], [0.0, 0.0, 1.0e-4]],
                "Cp": [[1.0e-5, 0.0, 0.0], [0.0, 1.0e-5, 0.0], [0.0, 0.0, 1.0e-5]]
              },
              "inputs": { "bus1": "bus_1", "bus2": "bus_2" }
            },
            {
              "class": "LoadZ",
              "id": "load_2",
              "params": {
                "R": [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
                "L": [[0.04, 0.0, 0.0], [0.0, 0.04, 0.0], [0.0, 0.0, 0.04]]
              },
              "inputs": { "bus": "bus_2" }
            },
            {
              "class": "LoadZ",
              "id": "load_sh1",
              "params": {
                "R": [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
                "L": [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.1]]
              },
              "inputs": { "bus": "bus_1" }
            },
            {
              "class": "LoadZ",
              "id": "load_sh2",
              "params": {
                "R": [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
                "L": [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.1]]
              },
              "inputs": { "bus": "bus_2" }
            }
          ]
        })";
      }

      static std::string twinLineRationalCaseJson()
      {
        return R"({
          "header": {
            "case_name": "EMT twin line, rational form",
            "case_description": "Line with rational series and shunt fits absorbing the terminal RL shunts",
            "case_comments": "Used by EMT SystemTests"
          },
          "devices": [
            { "class": "Bus", "id": "bus_1" },
            { "class": "Bus", "id": "bus_2" },
            {
              "class": "VoltageSource",
              "id": "source_1",
              "params": {
                "E": [100.0, 100.0, 100.0],
                "phi": [0.0, -2.0943951023931953, 2.0943951023931953],
                "omega": 376.99111843077515,
                "Rs": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                "Ls": [[0.01, 0.0, 0.0], [0.0, 0.01, 0.0], [0.0, 0.0, 0.01]]
              },
              "inputs": { "bus": "bus_1" }
            },
            {
              "class": "LineLumped",
              "id": "line_1_2",
              "params": {
                "N": 3,
                "K": 3,
                "conductors": [1, 2, 3],
                "dx": 2.0
              },
              "inputs": { "bus1": "bus_1", "bus2": "bus_2" },
              "submodels": {
                "Zp": {
                  "D": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                  "E": [[0.01, 0.0, 0.0], [0.0, 0.01, 0.0], [0.0, 0.0, 0.01]]
                },
                "Yp": {
                  "D": [[5.0e-5, 0.0, 0.0], [0.0, 5.0e-5, 0.0], [0.0, 0.0, 5.0e-5]],
                  "E": [[5.0e-6, 0.0, 0.0], [0.0, 5.0e-6, 0.0], [0.0, 0.0, 5.0e-6]],
                  "poles": [[-100.0, 0.0]],
                  "residues": [
                    [[[10.0, 0.0], [0.0, 0.0], [0.0, 0.0]],
                     [[0.0, 0.0], [10.0, 0.0], [0.0, 0.0]],
                     [[0.0, 0.0], [0.0, 0.0], [10.0, 0.0]]]
                  ]
                }
              }
            },
            {
              "class": "LoadZ",
              "id": "load_2",
              "params": {
                "R": [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
                "L": [[0.04, 0.0, 0.0], [0.0, 0.04, 0.0], [0.0, 0.0, 0.04]]
              },
              "inputs": { "bus": "bus_2" }
            }
          ]
        })";
      }

      static std::string twinSourceCaseJson()
      {
        return R"({
          "header": {
            "case_name": "EMT twin source, rational form",
            "case_description": "The series R and L source branch as a rational admittance",
            "case_comments": "Used by EMT SystemTests"
          },
          "devices": [
            { "class": "Bus", "id": "bus_1" },
            {
              "class": "VoltageSource",
              "id": "source_1",
              "params": {
                "E": [100.0, 100.0, 100.0],
                "phi": [0.0, -2.0943951023931953, 2.0943951023931953],
                "omega": 376.99111843077515
              },
              "inputs": { "bus": "bus_1" },
              "submodels": {
                "Y": {
                  "poles": [[-100.0, 0.0]],
                  "residues": [
                    [[[100.0, 0.0], [0.0, 0.0], [0.0, 0.0]],
                     [[0.0, 0.0], [100.0, 0.0], [0.0, 0.0]],
                     [[0.0, 0.0], [0.0, 0.0], [100.0, 0.0]]]
                  ]
                }
              }
            },
            {
              "class": "LoadZ",
              "id": "load_r",
              "params": {
                "R": [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
              },
              "inputs": { "bus": "bus_1" }
            },
            {
              "class": "LoadZ",
              "id": "load_l",
              "params": {
                "L": [[0.04, 0.0, 0.0], [0.0, 0.04, 0.0], [0.0, 0.0, 0.04]]
              },
              "inputs": { "bus": "bus_1" }
            }
          ]
        })";
      }

    public:
      /**
       * @brief Case JSON parses into the expected data containers.
       */
      TestOutcome parse()
      {
        TestStatus success = true;

        std::istringstream stream(caseJson());
        const auto         data = GridKit::EMT::parseSystemModelData(stream);

        success *= (data.case_name == "EMT two-device smoke case");
        success *= (data.bus.size() == 1);
        success *= (data.voltage_source.size() == 1);
        success *= (data.loadz.size() == 1);
        success *= (data.bus[0].id == "bus_1");
        success *= (data.voltage_source[0].id == "source_1");
        success *= (data.voltage_source[0].inputs.at(
                        GridKit::EMT::VoltageSourceInputs::bus)
                    == "bus_1");
        success *= (data.loadz[0].id == "load_1");
        success *= (data.loadz[0].inputs.at(GridKit::EMT::LoadZInputs::bus)
                    == "bus_1");

        std::istringstream governed_stream(machineGovernorCaseJson());
        const auto         governed = GridKit::EMT::parseSystemModelData(governed_stream);

        success *= (governed.signal.size() == 2);
        success *= (governed.signal[0].id == "speed_1");
        success *= (governed.signal[1].id == "pmech_1");
        success *= (governed.machine.size() == 1);
        success *= (governed.gov.size() == 1);
        success *= (governed.machine[0].inputs.at(GridKit::EMT::MachineInputs::bus)
                    == "bus_1");
        success *= (governed.machine[0].inputs.at(GridKit::EMT::MachineInputs::pm)
                    == "pmech_1");
        success *= (governed.machine[0].outputs.at(GridKit::EMT::MachineOutputs::speed)
                    == "speed_1");
        success *= (governed.gov[0].inputs.at(
                        GridKit::EMT::Controller::Tgov1Inputs::speed)
                    == "speed_1");
        success *= (governed.gov[0].outputs.at(
                        GridKit::EMT::Controller::Tgov1Outputs::pmech)
                    == "pmech_1");

        try
        {
          std::istringstream duplicate_stream(R"({
            "header": {
              "case_name": "duplicate IDs",
              "case_description": "Parser rejection coverage",
              "case_comments": ""
            },
            "devices": [
              { "class": "Bus", "id": "duplicate" },
              { "class": "Bus", "id": "duplicate" }
            ]
          })");
          (void) GridKit::EMT::parseSystemModelData(duplicate_stream);
          success *= false;
        }
        catch (const std::runtime_error&)
        {
        }

        return success.report(__func__);
      }

      /**
       * @brief Data-driven system assembles, verifies, and sizes correctly.
       */
      TestOutcome system()
      {
        TestStatus success = true;

        std::istringstream stream(caseJson());
        const auto         data = GridKit::EMT::parseSystemModelData(stream);

        GridKit::EMT::SystemModel<ScalarT, IdxT> sys(data);
        sys.allocate();

        success *= (sys.size() == 12);
        success *= (sys.verify() == 0);

        sys.tagDifferentiable();
        const auto& tag  = sys.tag();
        // Bus voltages and source voltages are algebraic; branch currents are
        // differential.
        success         *= (tag[0] == false);
        success         *= (tag[1] == false);
        success         *= (tag[2] == false);
        success         *= (tag[3] == false);
        success         *= (tag[4] == false);
        success         *= (tag[5] == false);
        success         *= (tag[6] == true);
        success         *= (tag[7] == true);
        success         *= (tag[8] == true);
        success         *= (tag[9] == true);
        success         *= (tag[10] == true);
        success         *= (tag[11] == true);

#ifdef GRIDKIT_ENABLE_ENZYME
        success *= sys.hasJacobian();
        success *= (sys.getCsrJacobian() != nullptr);
        success *= (sys.nnz() > 0);
#endif

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_SUNDIALS
      /**
       * @brief Integrate to sinusoidal steady state and compare against the
       * analytic per-phase phasor solution.
       */
      TestOutcome steadyState()
      {
        TestStatus success = true;

        std::istringstream stream(caseJson());
        const auto         data = GridKit::EMT::parseSystemModelData(stream);

        GridKit::EMT::SystemModel<ScalarT, IdxT> sys(data);
        sys.allocate();

        AnalysisManager::Sundials::Ida<ScalarT, IdxT> ida(&sys);
        ida.setMaxSteps(100000);
        ida.setTolerance(1.0e-9, 1.0e-9);
        ida.configureSimulation();
        ida.initializeSimulation(0.0, true);

        const RealT t_final = 0.3;
        ida.runSimulation(t_final);

        // Analytic per-phase phasor solution of the series circuit
        const RealT               omega = 376.99111843077515;
        const RealT               E_rms = 100.0;
        const RealT               sqrt2 = std::numbers::sqrt2_v<RealT>;
        const std::complex<RealT> Z_total{1.0 + 10.0, omega * (0.01 + 0.04)};
        const std::complex<RealT> Zs{1.0, omega * 0.01};
        const std::complex<RealT> rotation = std::exp(std::complex<RealT>{0.0, omega * t_final});

        const auto* y = sys.y().getData();

        const std::array<RealT, 3> phi{0.0, -2.0943951023931953, 2.0943951023931953};
        for (size_t n = 0; n < 3; ++n)
        {
          const std::complex<RealT> e_pk = sqrt2 * E_rms * std::exp(std::complex<RealT>{0.0, phi[n]});
          const std::complex<RealT> i_pk = e_pk / Z_total;
          const std::complex<RealT> v_pk = e_pk - Zs * i_pk;

          const RealT v_expected      = (v_pk * rotation).real();
          const RealT i_src_expected  = (i_pk * rotation).real();
          const RealT i_load_expected = -i_src_expected;

          // Tolerance is the measured error floor at the 1.0e-9 solver
          // tolerance, phase-c current 8.97e-10, with headroom.
          success *= isEqual(y[n], v_expected, 2.0e-9);
          success *= isEqual(y[6 + n], i_src_expected, 2.0e-9);
          success *= isEqual(y[9 + n], i_load_expected, 2.0e-9);
        }

        return success.report(__func__);
      }

      /**
       * @brief Integrate the source-line-load case to sinusoidal steady state
       * and compare against the analytic per-phase phasor ladder solution.
       *
       * The line's shunt capacitance makes both bus voltages differential.
       */
      TestOutcome threeBusSteadyState()
      {
        TestStatus success = true;

        std::istringstream stream(threeBusCaseJson());
        const auto         data = GridKit::EMT::parseSystemModelData(stream);

        GridKit::EMT::SystemModel<ScalarT, IdxT> sys(data);
        sys.allocate();

        success *= (sys.size() == 24);

        sys.tagDifferentiable();
        const auto& tag = sys.tag();
        // Both bus voltages are differential through the line's shunt
        // capacitance.
        for (size_t j = 0; j < 6; ++j)
        {
          success *= (tag[j] == true);
        }

        AnalysisManager::Sundials::Ida<ScalarT, IdxT> ida(&sys);
        ida.setMaxSteps(1000000);
        ida.setTolerance(1.0e-9, 1.0e-9);
        ida.configureSimulation();
        ida.initializeSimulation(0.0, true);

        const RealT t_final = 0.3;
        ida.runSimulation(t_final);

        // Analytic per-phase phasor ladder solution
        const RealT               omega = 376.99111843077515;
        const RealT               E_rms = 100.0;
        const RealT               sqrt2 = std::numbers::sqrt2_v<RealT>;
        const std::complex<RealT> Zs{1.0, omega * 0.01};
        const std::complex<RealT> Zl{2.0, omega * 0.02};
        const std::complex<RealT> Yh = std::complex<RealT>{1.0e-4, omega * 1.0e-5} / RealT{2.0};
        const std::complex<RealT> Zload{10.0, omega * 0.04};
        const std::complex<RealT> rotation = std::exp(std::complex<RealT>{0.0, omega * t_final});

        // Nodal equations:
        //   (1/Zs + Yh + 1/Zl) v1 - (1/Zl) v2 = e/Zs
        //   -(1/Zl) v1 + (1/Zl + Yh + 1/Zload) v2 = 0
        const std::complex<RealT> a11 = RealT{1.0} / Zs + Yh + RealT{1.0} / Zl;
        const std::complex<RealT> a12 = -RealT{1.0} / Zl;
        const std::complex<RealT> a21 = -RealT{1.0} / Zl;
        const std::complex<RealT> a22 = RealT{1.0} / Zl + Yh + RealT{1.0} / Zload;
        const std::complex<RealT> det = a11 * a22 - a12 * a21;

        const auto* y = sys.y().getData();

        const std::array<RealT, 3> phi{0.0, -2.0943951023931953, 2.0943951023931953};
        for (size_t n = 0; n < 3; ++n)
        {
          const std::complex<RealT> e_pk = sqrt2 * E_rms * std::exp(std::complex<RealT>{0.0, phi[n]});
          const std::complex<RealT> b1   = e_pk / Zs;
          const std::complex<RealT> v1   = (b1 * a22) / det;
          const std::complex<RealT> v2   = (-a21 * b1) / det;
          const std::complex<RealT> i12  = (v1 - v2) / Zl;
          const std::complex<RealT> isrc = (e_pk - v1) / Zs;
          const std::complex<RealT> ish1 = -Yh * v1;
          const std::complex<RealT> ish2 = -Yh * v2;
          const std::complex<RealT> ild  = -v2 / Zload;

          // Layout: bus1 v [0,3), bus2 v [3,6), source e [6,9) i [9,12),
          // line i12 [12,15) i_sh1 [15,18) i_sh2 [18,21), load i [21,24)
          // Tolerance is the measured error floor at the 1.0e-9 solver
          // tolerance with headroom, phase-a bus 1 voltage 1.5e-8, dominated
          // by integration error over the stiff shunt dynamics.
          success *= isEqual(y[n], (v1 * rotation).real(), 3.0e-8);
          success *= isEqual(y[3 + n], (v2 * rotation).real(), 3.0e-8);
          success *= isEqual(y[9 + n], (isrc * rotation).real(), 3.0e-8);
          success *= isEqual(y[12 + n], (i12 * rotation).real(), 3.0e-8);
          success *= isEqual(y[15 + n], (ish1 * rotation).real(), 3.0e-8);
          success *= isEqual(y[18 + n], (ish2 * rotation).real(), 3.0e-8);
          success *= isEqual(y[21 + n], (ild * rotation).real(), 3.0e-8);
        }

        return success.report(__func__);
      }

      /**
       * @brief Energize the load by closing the switch mid-simulation.
       *
       * Closing the switch changes the Jacobian sparsity pattern, so the
       * system Jacobian structure is invalidated and the linear solver is
       * reconfigured before the integrator restarts.
       */
      TestOutcome switchEnergization()
      {
        TestStatus success = true;

        std::istringstream stream(switchCaseJson());
        const auto         data = GridKit::EMT::parseSystemModelData(stream);

        GridKit::EMT::SystemModel<ScalarT, IdxT> sys(data);
        sys.allocate();

        // Component order: bus 1, bus 2, source, load 1, load 2, switch
        auto* sw  = dynamic_cast<GridKit::EMT::Switch<ScalarT, IdxT>*>(sys.getComponent(5));
        success  *= (sw != nullptr);
        success  *= (sw == sys.getSwitch("switch_1_2"));
        success  *= sw->isOpen();

        AnalysisManager::Sundials::Ida<ScalarT, IdxT> ida(&sys);
        ida.setMaxSteps(1000000);
        ida.setTolerance(1.0e-9, 1.0e-9);
        ida.configureSimulation();
        ida.initializeSimulation(0.0, true);

        const RealT t_close = 0.05;
        ida.runSimulation(t_close);

        // Before the close, the load side is de-energized.
        // Layout: bus1 v [0,3), bus2 v [3,6), source e [6,9) i [9,12),
        // load 1 i [12,15), load 2 i [15,18), switch i12 [18,21)
        const auto* y = sys.y().getData();
        for (size_t j = 0; j < 3; ++j)
        {
          success *= isEqual(y[3 + j], RealT{0.0}, 1.0e-12);  // bus 2 voltage
          success *= isEqual(y[15 + j], RealT{0.0}, 1.0e-12); // load 2 current
          success *= isEqual(y[18 + j], RealT{0.0}, 1.0e-12); // switch current
        }

        // Close the switch: rediscover the Jacobian structure, rebuild the
        // linear solver, and reinitialize the integrator.
        sw->setOpen(false);
        sys.resetJacobianStructure();
        sys.evaluateResidual();
        sys.evaluateJacobian();
        ida.configureLinearSolver();
        ida.initializeSimulation(t_close, true);

        const RealT t_final = 0.4;
        ida.runSimulation(t_final);

        // Analytic per-phase phasor solution of the energized circuit with
        // the two resistive loads in parallel behind the source impedance
        const RealT               omega = 376.99111843077515;
        const RealT               E_rms = 100.0;
        const RealT               sqrt2 = std::numbers::sqrt2_v<RealT>;
        const RealT               R_par = 100.0 * 10.0 / 110.0;
        const std::complex<RealT> Z_total{1.0 + R_par, omega * 0.01};
        const std::complex<RealT> rotation = std::exp(std::complex<RealT>{0.0, omega * t_final});

        const std::array<RealT, 3> phi{0.0, -2.0943951023931953, 2.0943951023931953};
        for (size_t n = 0; n < 3; ++n)
        {
          const std::complex<RealT> e_pk  = sqrt2 * E_rms * std::exp(std::complex<RealT>{0.0, phi[n]});
          const std::complex<RealT> i_pk  = e_pk / Z_total;
          const std::complex<RealT> v1_pk = R_par * i_pk;

          const RealT v1_expected      = (v1_pk * rotation).real();
          const RealT i_src_expected   = (i_pk * rotation).real();
          const RealT i_load2_expected = -v1_expected / 10.0;

          // Tolerance is the measured error floor at the 1.0e-9 solver
          // tolerance with headroom.
          success *= isEqual(y[n], v1_expected, 3.0e-8);
          success *= isEqual(y[3 + n], v1_expected, 3.0e-8);
          success *= isEqual(y[9 + n], i_src_expected, 3.0e-8);
          success *= isEqual(y[15 + n], i_load2_expected, 3.0e-8);
          success *= isEqual(y[18 + n], -i_load2_expected, 3.0e-8);
        }

        return success.report(__func__);
      }

      /**
       * @brief Machine flat start: an exactly initialized steady state stays
       * flat over one second of simulation.
       */
      TestOutcome machineFlatStart()
      {
        TestStatus success = true;

        std::istringstream stream(machineCaseJson());
        const auto         data = GridKit::EMT::parseSystemModelData(stream);

        GridKit::EMT::SystemModel<ScalarT, IdxT> sys(data);
        sys.allocate();

        // Layout: bus v [0,3), machine [3,27), load i [27,30)
        success *= (sys.size() == 30);

        // The initialized state satisfies the assembled residual exactly.
        sys.initialize();
        sys.evaluateResidual();
        const auto* f             = sys.getResidual().getData();
        RealT       residual_norm = 0.0;
        for (IdxT j = 0; j < sys.size(); ++j)
        {
          residual_norm += f[j] * f[j];
        }
        residual_norm  = std::sqrt(residual_norm);
        // The residual mixes per-unit machine rows with SI current-balance
        // rows on a 50 MW operating point, so the exact-initialization
        // floor is set by cancellation in the ampere-scale rows.
        success       *= (residual_norm < 1.0e-6);

        const auto* y   = sys.y().getData();
        const RealT te0 = y[3 + 20];

        AnalysisManager::Sundials::Ida<ScalarT, IdxT> ida(&sys);
        ida.setMaxSteps(1000000);
        ida.setTolerance(1.0e-9, 1.0e-9);
        ida.configureSimulation();
        ida.initializeSimulation(0.0, true);

        const RealT t_final = 1.0;
        ida.runSimulation(t_final);

        // The trajectory stays flat: synchronous speed, constant torque, and
        // the terminal voltage amplitude holds.
        //
        // Tolerances are the measured drift floors at the 1.0e-9 solver
        // tolerance with headroom.
        success *= isEqual(y[3 + 1], RealT{1.0}, 1.0e-7);
        success *= isEqual(y[3 + 20], te0, 1.0e-6);

        const RealT v_peak       = 11267.65281680262;
        const RealT va           = y[0];
        const RealT vb           = y[1];
        const RealT vc           = y[2];
        const RealT v_amplitude  = std::sqrt((TWO<RealT> / THREE<RealT>) *(va * va + vb * vb + vc * vc));
        success                 *= isEqual(v_amplitude, v_peak, 1.0e-6);

        return success.report(__func__);
      }

      /**
       * @brief Governed machine flat start: the TGOV1 governor initializes
       * from the machine dispatch and the trajectory stays flat.
       */
      TestOutcome machineGovernorFlatStart()
      {
        TestStatus success = true;

        std::istringstream stream(machineGovernorCaseJson());
        const auto         data = GridKit::EMT::parseSystemModelData(stream);

        GridKit::EMT::SystemModel<ScalarT, IdxT> sys(data);
        sys.allocate();

        // Layout: bus v [0,3), machine [3,27), load i [27,30), governor [30,33)
        success *= (sys.size() == 33);

        sys.initialize();
        sys.evaluateResidual();
        const auto* f             = sys.getResidual().getData();
        RealT       residual_norm = 0.0;
        for (IdxT j = 0; j < sys.size(); ++j)
        {
          residual_norm += f[j] * f[j];
        }
        residual_norm  = std::sqrt(residual_norm);
        success       *= (residual_norm < 1.0e-6);

        const auto* y   = sys.y().getData();
        const RealT te0 = y[3 + 20];
        const RealT pm0 = y[30 + 2];

        // The governor mechanical power matches the machine dispatch.
        success *= isEqual(pm0, te0, 1.0e-12);

        AnalysisManager::Sundials::Ida<ScalarT, IdxT> ida(&sys);
        ida.setMaxSteps(1000000);
        ida.setTolerance(1.0e-9, 1.0e-9);
        ida.configureSimulation();
        ida.initializeSimulation(0.0, true);

        const RealT t_final = 1.0;
        ida.runSimulation(t_final);

        // Tolerances are the measured drift floors at the 1.0e-9 solver
        // tolerance with headroom.
        success *= isEqual(y[3 + 1], RealT{1.0}, 1.0e-7);
        success *= isEqual(y[3 + 20], te0, 1.0e-6);
        success *= isEqual(y[30 + 2], pm0, 1.0e-6);

        return success.report(__func__);
      }

      /**
       * @brief Twin-circuit trajectory test for the rational operator.
       *
       * The same source energizes a parallel R-L branch twice: once built
       * from primitive loads and once as a one-real-pole rational impedance
       * with hand-derived coefficients. The trajectories of the shared
       * variables must agree at the measured solver floor mid-transient and
       * at steady state, and the rational branch current must equal the sum
       * of the primitive branch currents.
       */
      TestOutcome twinCircuit()
      {
        TestStatus success = true;

        std::istringstream primitive_stream(twinPrimitiveCaseJson());
        const auto         primitive_data = GridKit::EMT::parseSystemModelData(primitive_stream);
        std::istringstream rational_stream(twinRationalCaseJson());
        const auto         rational_data = GridKit::EMT::parseSystemModelData(rational_stream);

        GridKit::EMT::SystemModel<ScalarT, IdxT> primitive(primitive_data);
        primitive.allocate();
        GridKit::EMT::SystemModel<ScalarT, IdxT> rational(rational_data);
        rational.allocate();

        // Layouts: bus v [0,3), source e [3,6) i [6,9); then either the two
        // primitive load currents [9,15) or the rational current [9,12) and
        // one memory-state triple [12,15)
        success *= (primitive.size() == 15);
        success *= (rational.size() == 15);
        success *= (rational.verify() == 0);

        AnalysisManager::Sundials::Ida<ScalarT, IdxT> ida_primitive(&primitive);
        ida_primitive.setMaxSteps(1000000);
        ida_primitive.setTolerance(1.0e-9, 1.0e-9);
        ida_primitive.configureSimulation();
        ida_primitive.initializeSimulation(0.0, true);

        AnalysisManager::Sundials::Ida<ScalarT, IdxT> ida_rational(&rational);
        ida_rational.setMaxSteps(1000000);
        ida_rational.setTolerance(1.0e-9, 1.0e-9);
        ida_rational.configureSimulation();
        ida_rational.initializeSimulation(0.0, true);

        const auto* y_primitive = primitive.y().getData();
        const auto* y_rational  = rational.y().getData();

        // Tolerance is the measured trajectory-agreement floor at the 1.0e-9
        // solver tolerance with headroom.
        const RealT agreement_tol = 1.0e-6;

        for (const RealT t_check : {0.01, 0.3})
        {
          ida_primitive.runSimulation(t_check);
          ida_rational.runSimulation(t_check);

          for (size_t j = 0; j < 9; ++j)
          {
            success *= isEqual(y_rational[j], y_primitive[j], agreement_tol);
          }
          for (size_t n = 0; n < 3; ++n)
          {
            success *= isEqual(y_rational[9 + n],
                               y_primitive[9 + n] + y_primitive[12 + n],
                               agreement_tol);
          }
        }

        return success.report(__func__);
      }

      /**
       * @brief Rational line agrees with an equivalent matrix-parameter line.
       *
       * The rational line carries a feedthrough-only series fit and a shunt
       * fit with one real pole, scaled by a segment length of two. The twin
       * builds the same electrical circuit from the matrix-parameter line and
       * explicit series-RL loads at each terminal absorbing the shunt pole
       * term. The shared trajectories must agree at the measured solver
       * floor, and each rational shunt current must equal the matrix shunt
       * current plus the terminal load current.
       */
      TestOutcome twinLine()
      {
        TestStatus success = true;

        std::istringstream matrix_stream(twinLineMatrixCaseJson());
        const auto         matrix_data = GridKit::EMT::parseSystemModelData(matrix_stream);
        std::istringstream rational_stream(twinLineRationalCaseJson());
        const auto         rational_data = GridKit::EMT::parseSystemModelData(rational_stream);

        GridKit::EMT::SystemModel<ScalarT, IdxT> matrix(matrix_data);
        matrix.allocate();
        GridKit::EMT::SystemModel<ScalarT, IdxT> rational(rational_data);
        rational.allocate();

        // Layouts: bus voltages [0,6), source e [6,9) i [9,12), line i12
        // [12,15) ish1 [15,18) ish2 [18,21); then either load_2 [21,24) and
        // the terminal loads [24,30), or the two shunt memory-state triples
        // [21,27) and load_2 [27,30)
        success *= (matrix.size() == 30);
        success *= (rational.size() == 30);
        success *= (rational.verify() == 0);

        AnalysisManager::Sundials::Ida<ScalarT, IdxT> ida_matrix(&matrix);
        ida_matrix.setMaxSteps(1000000);
        ida_matrix.setTolerance(1.0e-9, 1.0e-9);
        ida_matrix.configureSimulation();
        ida_matrix.initializeSimulation(0.0, true);

        AnalysisManager::Sundials::Ida<ScalarT, IdxT> ida_rational(&rational);
        ida_rational.setMaxSteps(1000000);
        ida_rational.setTolerance(1.0e-9, 1.0e-9);
        ida_rational.configureSimulation();
        ida_rational.initializeSimulation(0.0, true);

        const auto* y_matrix   = matrix.y().getData();
        const auto* y_rational = rational.y().getData();

        // Tolerance is the measured trajectory-agreement floor at the 1.0e-9
        // solver tolerance with headroom.
        const RealT agreement_tol = 1.0e-6;

        for (const RealT t_check : {0.01, 0.3})
        {
          ida_matrix.runSimulation(t_check);
          ida_rational.runSimulation(t_check);

          for (size_t j = 0; j < 15; ++j)
          {
            success *= isEqual(y_rational[j], y_matrix[j], agreement_tol);
          }
          for (size_t n = 0; n < 3; ++n)
          {
            success *= isEqual(y_rational[15 + n],
                               y_matrix[15 + n] + y_matrix[24 + n],
                               agreement_tol);
            success *= isEqual(y_rational[18 + n],
                               y_matrix[18 + n] + y_matrix[27 + n],
                               agreement_tol);
            success *= isEqual(y_rational[27 + n], y_matrix[21 + n], agreement_tol);
          }
        }

        return success.report(__func__);
      }

      /**
       * @brief Rational-admittance source agrees with the series-matrix form.
       *
       * The series R and L source branch is exactly the one-real-pole
       * admittance \f$ (1/L_s)/(s + R_s/L_s) \f$ with no feedthrough. The
       * twin drives the primitive-load circuit of the twin-circuit test from
       * the same source written both ways. The shared trajectories must
       * agree at the measured solver floor, the branch voltage must equal
       * the EMF minus the terminal voltage, and the memory state scaled by
       * the residue must recover the primitive branch current.
       */
      TestOutcome twinSource()
      {
        TestStatus success = true;

        std::istringstream primitive_stream(twinPrimitiveCaseJson());
        const auto         primitive_data = GridKit::EMT::parseSystemModelData(primitive_stream);
        std::istringstream rational_stream(twinSourceCaseJson());
        const auto         rational_data = GridKit::EMT::parseSystemModelData(rational_stream);

        GridKit::EMT::SystemModel<ScalarT, IdxT> primitive(primitive_data);
        primitive.allocate();
        GridKit::EMT::SystemModel<ScalarT, IdxT> rational(rational_data);
        rational.allocate();

        // Layouts: bus v [0,3), source e [3,6); then either the primitive
        // branch current [6,9) and the loads [9,15), or the branch voltage
        // [6,9), one memory-state triple [9,12), and the loads [12,18)
        success *= (primitive.size() == 15);
        success *= (rational.size() == 18);
        success *= (rational.verify() == 0);

        AnalysisManager::Sundials::Ida<ScalarT, IdxT> ida_primitive(&primitive);
        ida_primitive.setMaxSteps(1000000);
        ida_primitive.setTolerance(1.0e-9, 1.0e-9);
        ida_primitive.configureSimulation();
        ida_primitive.initializeSimulation(0.0, true);

        AnalysisManager::Sundials::Ida<ScalarT, IdxT> ida_rational(&rational);
        ida_rational.setMaxSteps(1000000);
        ida_rational.setTolerance(1.0e-9, 1.0e-9);
        ida_rational.configureSimulation();
        ida_rational.initializeSimulation(0.0, true);

        const auto* y_primitive = primitive.y().getData();
        const auto* y_rational  = rational.y().getData();

        // Tolerance is the measured trajectory-agreement floor at the 1.0e-9
        // solver tolerance with headroom.
        const RealT agreement_tol = 1.0e-6;

        // The residue folded into the memory-state output
        const RealT residue = 100.0;

        for (const RealT t_check : {0.01, 0.3})
        {
          ida_primitive.runSimulation(t_check);
          ida_rational.runSimulation(t_check);

          for (size_t j = 0; j < 6; ++j)
          {
            success *= isEqual(y_rational[j], y_primitive[j], agreement_tol);
          }
          for (size_t n = 0; n < 3; ++n)
          {
            success *= isEqual(y_rational[6 + n],
                               y_rational[3 + n] - y_rational[n],
                               agreement_tol);
            success *= isEqual(residue * y_rational[9 + n],
                               y_primitive[6 + n],
                               agreement_tol);
            success *= isEqual(y_rational[12 + n], y_primitive[9 + n], agreement_tol);
            success *= isEqual(y_rational[15 + n], y_primitive[12 + n], agreement_tol);
          }
        }

        return success.report(__func__);
      }
#endif
    };

  } // namespace Testing
} // namespace GridKit
