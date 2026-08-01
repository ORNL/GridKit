#pragma once

#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string>
#include <variant>

#include <GridKit/Model/OPF/SystemData.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    class SystemDataTests
    {
    public:
      TestOutcome parseCompleteSystem()
      {
        const std::string document = R"(
          {
            "header": {
              "format_version": 0,
              "format_revision": 1,
              "case_name": "Complete OPF parser test",
              "case_date_time": "2026-07-31T12:00:00-05:00",
              "case_description": null,
              "case_comments": "Exercises every supported device"
            },
            "params": {
              "freq_base": 60,
              "va_base": 100000000.0
            },
            "buses": [
              {
                "class": "Bus",
                "id": "Bus_1",
                "params": {
                  "number": 1,
                  "kv": 230.0,
                  "vmin": 0.9,
                  "vmax": 1.1
                }
              },
              {
                "class": "Bus",
                "id": "Bus_2",
                "params": {
                  "number": 2,
                  "kv": 230
                }
              }
            ],
            "devices": [
              {
                "class": "Branch",
                "buses": {"from": 1, "to": 2},
                "id": "BR1",
                "params": {
                  "R": 0.01,
                  "X": 0.1,
                  "G": 0.0,
                  "B": 0.02,
                  "smax": 2.5
                }
              },
              {
                "class": "Generator",
                "buses": {"bus": 1},
                "id": "G1",
                "params": {
                  "pmin": null,
                  "qmin": -0.4,
                  "qmax": 0.4,
                  "mva": 100.0,
                  "c1": 12.0
                }
              },
              {
                "class": "Load",
                "buses": {"bus": 2},
                "id": "L1",
                "params": {
                  "pmin": -2.0,
                  "pmax": 0.0,
                  "qmax": null
                }
              },
              {
                "class": "Shunt",
                "buses": {"bus": 2},
                "id": "SH1",
                "params": {"G": 0.01, "B": -0.05}
              }
            ]
          })";

        std::istringstream stream(document);
        const auto         system_data = OPF::parseSystemData(stream);

        TestStatus success  = true;
        success            *= system_data.header.format_version == 0U;
        success            *= system_data.header.format_revision == 1U;
        success            *= system_data.header.case_name == "Complete OPF parser test";
        success            *= system_data.header.case_date_time
                   == "2026-07-31T12:00:00-05:00";
        success *= !system_data.header.case_description.has_value();
        success *= system_data.header.case_comments
                   == "Exercises every supported device";
        success *= system_data.params.freq_base == 60.0;
        success *= system_data.params.va_base == 100000000.0;

        success *= system_data.buses.size() == 2;
        if (system_data.buses.size() == 2)
        {
          success *= system_data.buses[0].number == 1;
          success *= system_data.buses[0].vmin == 0.9;
          success *= system_data.buses[0].vmax == 1.1;
          success *= system_data.buses[1].number == 2;
          success *= !system_data.buses[1].vmin.has_value();
          success *= !system_data.buses[1].vmax.has_value();
        }

        success *= system_data.devices.size() == 4;
        if (system_data.devices.size() == 4)
        {
          const auto* branch = std::get_if<OPF::SystemData<>::BranchDataT>(
              &system_data.devices[0]);
          const auto* generator = std::get_if<OPF::SystemData<>::GeneratorDataT>(
              &system_data.devices[1]);
          const auto* load = std::get_if<OPF::SystemData<>::LoadDataT>(
              &system_data.devices[2]);
          const auto* shunt = std::get_if<OPF::SystemData<>::ShuntDataT>(
              &system_data.devices[3]);

          success *= branch != nullptr;
          success *= generator != nullptr;
          success *= load != nullptr;
          success *= shunt != nullptr;

          if (branch != nullptr)
          {
            success *= branch->id == "BR1";
            success *= branch->from == 1;
            success *= branch->to == 2;
            success *= branch->smax == 2.5;
          }
          if (generator != nullptr)
          {
            success *= generator->id == "G1";
            success *= !generator->pmin.has_value();
            success *= !generator->pmax.has_value();
            success *= generator->qmin == -0.4;
            success *= generator->qmax == 0.4;
            success *= generator->c0 == 0.0;
            success *= generator->c1 == 12.0;
            success *= generator->c2 == 0.0;
          }
          if (load != nullptr)
          {
            success *= load->id == "L1";
            success *= load->pmin == -2.0;
            success *= load->pmax == 0.0;
            success *= !load->qmin.has_value();
            success *= !load->qmax.has_value();
          }
          if (shunt != nullptr)
          {
            success *= shunt->id == "SH1";
            success *= shunt->bus == 2;
            success *= shunt->G == 0.01;
            success *= shunt->B == -0.05;
          }
        }

        return success.report(__func__);
      }

      TestOutcome parseRepositoryCases()
      {
        TestStatus success = true;
        try
        {
          const auto two_bus = OPF::parseSystemData(
              std::filesystem::path{"TwoBusBasic.opf.json"});
          const auto two_bus_gensal = OPF::parseSystemData(
              std::filesystem::path{"TwoBusGensal.opf.json"});
          const auto three_bus = OPF::parseSystemData(
              std::filesystem::path{"ThreeBusBasic.opf.json"});
          const auto texas = OPF::parseSystemData(
              std::filesystem::path{"texas.opf.json"});

          success *= two_bus.buses.size() == 2;
          success *= two_bus.devices.size() == 3;
          success *= two_bus_gensal.buses.size() == 2;
          success *= !two_bus_gensal.devices.empty();
          success *= three_bus.buses.size() == 3;
          success *= !three_bus.devices.empty();
          success *= texas.buses.size() == 2000;
          success *= texas.devices.size() == 5145;
        }
        catch (const std::exception&)
        {
          success = false;
        }
        return success.report(__func__);
      }

      TestOutcome rejectInvalidDocuments()
      {
        TestStatus success  = true;
        success            *= parserRejects("{]");
        success            *= parserRejects(R"({
          "header": {"format_version": 0, "format_revision": 1, "case_name": "x"},
          "params": {"freq_base": 60, "va_base": 100},
          "buses": [], "devices": [], "extra": true
        })");
        success            *= parserRejects(R"({
          "header": {"format_version": 0.0, "format_revision": 1, "case_name": "x"},
          "params": {"freq_base": 60, "va_base": 100},
          "buses": [], "devices": []
        })");
        success            *= parserRejects(R"({
          "header": {"format_version": 1, "format_revision": 0, "case_name": "x"},
          "params": {"freq_base": 60, "va_base": 100},
          "buses": [], "devices": []
        })");
        success            *= parserRejects(R"({
          "header": {"format_version": 0, "format_revision": 1, "case_name": "x"},
          "params": {"freq_base": 60, "va_base": 100},
          "buses": [{"class": "Reference", "id": "B1",
                      "params": {"number": 1, "kv": 115}}],
          "devices": []
        })");
        success            *= parserRejects(R"({
          "header": {"format_version": 0, "format_revision": 1, "case_name": "x"},
          "params": {"freq_base": 60, "va_base": 100},
          "buses": [{"class": "Bus", "id": "B1",
                      "params": {"number": -1, "kv": 115}}],
          "devices": []
        })");
        success            *= parserRejects(R"({
          "header": {"format_version": 0, "format_revision": 1, "case_name": "x"},
          "params": {"freq_base": 60, "va_base": 100},
          "buses": [],
          "devices": [{"class": "Transformer", "id": "T1",
                        "buses": {"from": 1, "to": 2}, "params": {}}]
        })");
        success            *= parserRejects(R"({
          "header": {"format_version": 0, "format_revision": 1, "case_name": "x"},
          "params": {"freq_base": 60, "va_base": 100},
          "buses": [],
          "devices": [{"class": "Branch", "id": "BR1",
                        "buses": {"from": 1, "to": 2},
                        "params": {"R": 0, "X": 0.1, "G": 0}}]
        })");
        success            *= parserRejectsWithMessage(
            R"({
              "header": {"format_version": 0, "format_revision": 1, "case_name": "x"},
              "params": {"freq_base": 60, "va_base": 100},
              "buses": [],
              "devices": [{"class": "Generator", "id": "G1",
                            "buses": {"bus": 1},
                            "params": {"mva": 100, "pmaxx": 2}}]
            })",
            "devices[0] (Generator \"G1\").params.pmaxx");
        success *= parserRejects(R"({
          "header": {"format_version": 0, "format_revision": 1, "case_name": "x"},
          "params": {"freq_base": 60, "va_base": 100},
          "buses": [],
          "devices": [{"class": "Generator", "id": "G1",
                        "buses": {"bus": 1},
                        "params": {"mva": 100, "c0": null}}]
        })");

        return success.report(__func__);
      }

      TestOutcome rejectMissingFile()
      {
        TestStatus success = false;
        try
        {
          static_cast<void>(OPF::parseSystemData(
              std::filesystem::path{"missing-opf-system-data-file.json"}));
        }
        catch (const std::runtime_error&)
        {
          success = true;
        }
        return success.report(__func__);
      }

    private:
      static bool parserRejects(const std::string& document)
      {
        std::istringstream stream(document);
        try
        {
          static_cast<void>(OPF::parseSystemData(stream));
        }
        catch (const std::invalid_argument&)
        {
          return true;
        }
        return false;
      }

      static bool parserRejectsWithMessage(const std::string& document,
                                           const std::string& expected_message)
      {
        std::istringstream stream(document);
        try
        {
          static_cast<void>(OPF::parseSystemData(stream));
        }
        catch (const std::invalid_argument& error)
        {
          return std::string(error.what()).find(expected_message) != std::string::npos;
        }
        return false;
      }
    };

  } // namespace Testing
} // namespace GridKit
