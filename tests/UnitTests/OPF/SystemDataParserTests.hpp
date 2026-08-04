#pragma once

#include <array>
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
    class SystemDataParserTests
    {
    public:
      TestOutcome parseCompleteSystem()
      {
        std::istringstream stream(completeDocument());
        const auto         data = OPF::parseSystemData(stream);

        TestStatus success  = true;
        success            *= data.header.format_version == 0U;
        success            *= data.header.format_revision == 1U;
        success            *= data.header.case_name == "Complete parser test";
        success            *= data.header.case_date_time == "2026-08-03T12:00:00-05:00";
        success            *= data.header.case_description == "Every supported field";
        success            *= data.header.case_comments == "Parser contract";
        success            *= data.params.freq_base == 60.0;
        success            *= data.params.va_base == 100000000.0;

        success *= data.buses.size() == 2;
        if (data.buses.size() == 2)
        {
          success *= data.buses[0].id == "Bus_1";
          success *= data.buses[0].number == 1;
          success *= data.buses[0].kv == 230.0;
          success *= data.buses[0].vmin == 0.9;
          success *= data.buses[0].vmax == 1.1;
          success *= data.buses[1].id == "Bus_2";
          success *= data.buses[1].number == 2;
        }

        success *= data.devices.size() == 4;
        if (data.devices.size() == 4)
        {
          const auto* branch = std::get_if<OPF::SystemData<>::BranchDataT>(
              &data.devices[0]);
          const auto* generator = std::get_if<OPF::SystemData<>::GeneratorDataT>(
              &data.devices[1]);
          const auto* load = std::get_if<OPF::SystemData<>::LoadDataT>(
              &data.devices[2]);
          const auto* shunt = std::get_if<OPF::SystemData<>::ShuntDataT>(
              &data.devices[3]);

          success *= branch != nullptr;
          success *= generator != nullptr;
          success *= load != nullptr;
          success *= shunt != nullptr;

          if (branch != nullptr)
          {
            success *= branch->id == "BR1";
            success *= branch->from == 1;
            success *= branch->to == 2;
            success *= branch->R == 0.01;
            success *= branch->X == 0.1;
            success *= branch->G == 0.001;
            success *= branch->B == 0.02;
            success *= branch->smax == 2.5;
          }
          if (generator != nullptr)
          {
            success *= generator->id == "G1";
            success *= generator->bus == 1;
            success *= generator->pmin == 0.0;
            success *= generator->pmax == 1.2;
            success *= generator->qmin == -0.4;
            success *= generator->qmax == 0.4;
            success *= generator->mva == 120.0;
            success *= generator->c0 == 1.0;
            success *= generator->c1 == 12.0;
            success *= generator->c2 == 0.02;
          }
          if (load != nullptr)
          {
            success *= load->id == "L1";
            success *= load->bus == 2;
            success *= load->pmin == -2.0;
            success *= load->pmax == 0.0;
            success *= load->qmin == -1.0;
            success *= load->qmax == 0.0;
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

      TestOutcome parseOptionalFieldsAndDefaults()
      {
        const std::string document = R"({
          "header": {
            "format_version": 0,
            "format_revision": 1,
            "case_name": "Optional fields",
            "case_date_time": null,
            "case_description": null
          },
          "params": {"freq_base": 50, "va_base": 1},
          "buses": [
            {"class": "Bus", "id": "B1",
             "params": {"number": 0, "kv": 115, "vmin": null}}
          ],
          "devices": [
            {"class": "Branch", "buses": {"from": 0, "to": 1}, "id": "BR",
             "params": {"R": 0, "X": 0.1, "G": 0, "B": 0, "smax": null}},
            {"class": "Generator", "buses": {"bus": 0}, "id": "G",
             "params": {"pmin": null, "qmax": null, "mva": 100}},
            {"class": "Load", "buses": {"bus": 0}, "id": "L",
             "params": {"pmax": null}},
            {"class": "Shunt", "buses": {"bus": 0}, "id": "S",
             "params": {"G": 0, "B": 0}}
          ]
        })";

        std::istringstream stream(document);
        const auto         data = OPF::parseSystemData(stream);

        TestStatus success  = true;
        success            *= !data.header.case_date_time.has_value();
        success            *= !data.header.case_description.has_value();
        success            *= !data.header.case_comments.has_value();
        success            *= !data.buses[0].vmin.has_value();
        success            *= !data.buses[0].vmax.has_value();

        const auto& branch = std::get<OPF::SystemData<>::BranchDataT>(
            data.devices[0]);
        const auto& generator = std::get<OPF::SystemData<>::GeneratorDataT>(
            data.devices[1]);
        const auto& load = std::get<OPF::SystemData<>::LoadDataT>(data.devices[2]);

        success *= !branch.smax.has_value();
        success *= !generator.pmin.has_value();
        success *= !generator.pmax.has_value();
        success *= !generator.qmin.has_value();
        success *= !generator.qmax.has_value();
        success *= generator.c0 == 0.0;
        success *= generator.c1 == 0.0;
        success *= generator.c2 == 0.0;
        success *= !load.pmin.has_value();
        success *= !load.pmax.has_value();
        success *= !load.qmin.has_value();
        success *= !load.qmax.has_value();

        return success.report(__func__);
      }

      TestOutcome rejectUnknownFields()
      {
        TestStatus success  = true;
        success            *= parserRejects(replaceOnce(
            minimalDocument(), "\"devices\": []", "\"devices\": [], \"extra\": true"));
        success            *= parserRejects(replaceOnce(
            minimalDocument(), "\"case_name\": \"Minimal\"", "\"case_name\": \"Minimal\", \"extra\": true"));
        success            *= parserRejects(replaceOnce(
            minimalDocument(), "\"va_base\": 100", "\"va_base\": 100, \"extra\": 0"));
        success            *= parserRejects(replaceOnce(
            completeDocument(), "\"id\": \"Bus_1\"", "\"id\": \"Bus_1\", \"extra\": 0"));
        success            *= parserRejects(replaceOnce(
            completeDocument(), "\"vmax\": 1.1", "\"vmax\": 1.1, \"extra\": 0"));
        success            *= parserRejects(replaceOnce(
            completeDocument(), "\"id\": \"BR1\"", "\"id\": \"BR1\", \"extra\": 0"));
        success            *= parserRejects(replaceOnce(
            completeDocument(), "\"from\": 1, \"to\": 2", "\"from\": 1, \"to\": 2, \"extra\": 0"));
        success            *= parserRejectsWithMessage(
            replaceOnce(completeDocument(), "\"c2\": 0.02", "\"c2\": 0.02, \"pmaxx\": 2"),
            "devices[1] (Generator \"G1\").params.pmaxx");

        return success.report(__func__);
      }

      TestOutcome rejectInvalidTypesAndVersions()
      {
        TestStatus success  = true;
        success            *= parserRejects("[]");
        success            *= parserRejects("{]");
        success            *= parserRejects(R"({
          "header": null,
          "params": {"freq_base": 60, "va_base": 100},
          "buses": [],
          "devices": []
        })");
        success            *= parserRejects(replaceOnce(
            minimalDocument(), "\"format_version\": 0", "\"format_version\": 0.0"));
        success            *= parserRejects(replaceOnce(
            minimalDocument(), "\"format_version\": 0", "\"format_version\": 1"));
        success            *= parserRejects(replaceOnce(
            minimalDocument(), "\"format_revision\": 1", "\"format_revision\": 2"));
        success            *= parserRejects(replaceOnce(
            minimalDocument(), "\"freq_base\": 60", "\"freq_base\": \"60\""));
        success            *= parserRejects(replaceOnce(
            minimalDocument(), "\"va_base\": 100", "\"va_base\": 1e999"));
        success            *= parserRejects(replaceOnce(
            completeDocument(),
            "\"case_date_time\": \"2026-08-03T12:00:00-05:00\"",
            "\"case_date_time\": true"));
        success            *= parserRejects(replaceOnce(
            minimalDocument(), "\"buses\": []", "\"buses\": {}"));
        success            *= parserRejects(replaceOnce(
            completeDocument(), "\"number\": 1", "\"number\": 1.0"));
        success            *= parserRejects(replaceOnce(
            completeDocument(), "\"number\": 1", "\"number\": -1"));
        success            *= parserRejects(replaceOnce(
            completeDocument(), "\"kv\": 230.0", "\"kv\": \"230\""));
        success            *= parserRejects(replaceOnce(
            completeDocument(), "\"vmin\": 0.9", "\"vmin\": false"));
        success            *= parserRejects(replaceOnce(
            completeDocument(), "\"class\": \"Bus\"", "\"class\": \"Reference\""));
        success            *= parserRejects(replaceOnce(
            completeDocument(), "\"class\": \"Branch\"", "\"class\": \"Transformer\""));
        success            *= parserRejects(replaceOnce(
            completeDocument(), "\"from\": 1", "\"from\": -1"));
        success            *= parserRejects(replaceOnce(
            completeDocument(), "\"c0\": 1.0", "\"c0\": null"));

        return success.report(__func__);
      }

      TestOutcome rejectMissingRequiredFields()
      {
        TestStatus success  = true;
        success            *= parserRejects(replaceOnce(
            minimalDocument(), "\"params\": {\"freq_base\": 60, \"va_base\": 100},", ""));
        success            *= parserRejects(replaceOnce(
            minimalDocument(), "\"case_name\": \"Minimal\"", "\"case_name_missing\": \"Minimal\""));
        success            *= parserRejects(replaceOnce(
            completeDocument(), "\"B\": 0.02,", ""));
        success            *= parserRejects(replaceOnce(
            completeDocument(), "\"mva\": 120.0,", ""));
        success            *= parserRejects(replaceOnce(
            completeDocument(), "\"G\": 0.01,", ""));

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

      TestOutcome parseCuratedFixtures()
      {
        const std::array<std::string, 8> fixture_names{{"TwoBusBasic.opf.json",
                                                        "TwoBusGensal.opf.json",
                                                        "TwoBusIeeet1.opf.json",
                                                        "TwoBusTgov1.opf.json",
                                                        "ThreeBusBasic.opf.json",
                                                        "ThreeBusClassical.opf.json",
                                                        "ThreeBusConstantSource.opf.json",
                                                        "ThreeBusZipLoad.opf.json"}};

        TestStatus success = true;
        for (const auto& fixture_name : fixture_names)
        {
          const auto path = std::filesystem::path{GRIDKIT_OPF_CASE_DIRECTORY}
                            / fixture_name;
          const auto data  = OPF::parseSystemData(path);
          success         *= !data.header.case_name.empty();
          success         *= !data.buses.empty();
          success         *= !data.devices.empty();
        }
        return success.report(__func__);
      }

    private:
      static std::string minimalDocument()
      {
        return R"({
          "header": {
            "format_version": 0,
            "format_revision": 1,
            "case_name": "Minimal"
          },
          "params": {"freq_base": 60, "va_base": 100},
          "buses": [],
          "devices": []
        })";
      }

      static std::string completeDocument()
      {
        return R"({
          "header": {
            "format_version": 0,
            "format_revision": 1,
            "case_name": "Complete parser test",
            "case_date_time": "2026-08-03T12:00:00-05:00",
            "case_description": "Every supported field",
            "case_comments": "Parser contract"
          },
          "params": {"freq_base": 60, "va_base": 100000000.0},
          "buses": [
            {"class": "Bus", "id": "Bus_1",
             "params": {"number": 1, "kv": 230.0, "vmin": 0.9, "vmax": 1.1}},
            {"class": "Bus", "id": "Bus_2",
             "params": {"number": 2, "kv": 230}}
          ],
          "devices": [
            {"class": "Branch", "buses": {"from": 1, "to": 2}, "id": "BR1",
             "params": {"R": 0.01, "X": 0.1, "G": 0.001, "B": 0.02,
                        "smax": 2.5}},
            {"class": "Generator", "buses": {"bus": 1}, "id": "G1",
             "params": {"pmin": 0.0, "pmax": 1.2, "qmin": -0.4,
                        "qmax": 0.4, "mva": 120.0, "c0": 1.0,
                        "c1": 12.0, "c2": 0.02}},
            {"class": "Load", "buses": {"bus": 2}, "id": "L1",
             "params": {"pmin": -2.0, "pmax": 0.0, "qmin": -1.0,
                        "qmax": 0.0}},
            {"class": "Shunt", "buses": {"bus": 2}, "id": "SH1",
             "params": {"G": 0.01, "B": -0.05}}
          ]
        })";
      }

      static std::string replaceOnce(std::string        document,
                                     const std::string& original,
                                     const std::string& replacement)
      {
        const auto position = document.find(original);
        if (position == std::string::npos)
        {
          throw std::logic_error("Parser test replacement target was not found");
        }
        document.replace(position, original.size(), replacement);
        return document;
      }

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
