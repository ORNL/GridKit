/**
 * @file PartitionTests.hpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 * @brief Tests for applying partition files to flat cases.
 */

#pragma once

#include <exception>
#include <filesystem>
#include <limits>
#include <sstream>
#include <string>

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/Partition/PartitionData.hpp>
#include <GridKit/Model/PhasorDynamics/Partition/PartitionPlacement.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class PartitionTests
    {
    private:
      using ComponentT  = PhasorDynamics::Component<ScalarT, IdxT>;
      using RealT       = typename ComponentT::RealT;
      using DataT       = PhasorDynamics::SystemModelData<RealT, IdxT>;
      using SignalDataT = typename DataT::SignalDataT;
      using PlacementT  = PhasorDynamics::PartitionPlacement<RealT, IdxT>;
      using SystemT     = PhasorDynamics::SystemModel<ScalarT, IdxT>;
      using BusT        = PhasorDynamics::BusBase<ScalarT, IdxT>;
      using Log         = ::GridKit::Utilities::Logger;

    public:
      PartitionTests()  = default;
      ~PartitionTests() = default;

      /**
       * @brief Split ThreeBusBasic in two and compare with the flat system.
       *
       * Devices go to the first partition listing all their buses.
       */
      TestOutcome twoPartitions()
      {
        TestStatus success = true;

        const DataT flat_data =
            PhasorDynamics::parseSystemModelData(std::string("ThreeBusBasic.case.json"));
        const DataT top_data = apply(flat_data, R"({
          "name": "root",
          "partitions": [
            { "name": "A", "external": [0, 1] },
            { "name": "B", "internal": [2], "external": [0, 1] }
          ]
        })");

        success *= top_data.bus.size() == 2;
        success *= top_data.branch.empty();
        success *= top_data.genrou.empty();
        success *= top_data.loadz.empty();
        success *= top_data.case_name == flat_data.case_name;
        success *= top_data.monitor_sink.size() == flat_data.monitor_sink.size();
        success *= top_data.subsystems.size() == 2;
        if (top_data.subsystems.size() == 2)
        {
          const DataT& a = top_data.subsystems[0];
          const DataT& b = top_data.subsystems[1];

          success *= a.id == "0";
          success *= b.id == "1";
          success *= a.level == 1;
          success *= b.level == 1;
          success *= a.bus.empty();
          success *= a.branch.size() == 1;
          success *= a.genrou.size() == 1;
          success *= a.loadz.empty();
          success *= b.bus.size() == 1;
          success *= b.bus[0].bus_id == 2;
          success *= b.branch.size() == 2;
          success *= b.genrou.size() == 1;
          success *= b.loadz.size() == 1;
          success *= b.freq_base == flat_data.freq_base;
          success *= b.va_base == flat_data.va_base;

          SystemT flat(flat_data);
          SystemT top(top_data);
          compareSystems(flat_data, flat, top, success);
        }

        return success.report(__func__);
      }

      /**
       * @brief A group partition nests subsystems and keeps external buses at the top.
       */
      TestOutcome multiLevel()
      {
        TestStatus success = true;

        const DataT flat_data =
            PhasorDynamics::parseSystemModelData(std::string("ThreeBusBasic.case.json"));
        const DataT top_data = apply(flat_data, R"({
          "name": "root",
          "partitions": [
            { "name": "A", "external": [0, 1] },
            {
              "name": "group",
              "partitions": [
                { "name": "B", "internal": [2], "external": [0, 1] }
              ]
            }
          ]
        })");

        success *= top_data.bus.size() == 2;
        success *= top_data.subsystems.size() == 2;
        if (top_data.subsystems.size() == 2)
        {
          const DataT& group = top_data.subsystems[1];

          success *= group.bus.empty();
          success *= group.branch.empty();
          success *= group.subsystems.size() == 1;
          if (group.subsystems.size() == 1)
          {
            const DataT& b = group.subsystems[0];

            success *= b.id == "1.0";
            success *= b.level == 2;
            success *= b.bus.size() == 1;
            success *= b.branch.size() == 2;

            SystemT flat(flat_data);
            SystemT top(top_data);
            compareSystems(flat_data, flat, top, success);
          }
        }

        return success.report(__func__);
      }

      /**
       * @brief Bus lists and child partitions are optional.
       */
      TestOutcome optionalKeys()
      {
        TestStatus success = true;

        const DataT flat_data =
            PhasorDynamics::parseSystemModelData(std::string("ThreeBusBasic.case.json"));

        // No external list
        const DataT no_external  = apply(flat_data, R"({
          "name": "root",
          "partitions": [
            { "name": "A", "internal": [0, 1, 2] }
          ]
        })");
        success                 *= no_external.bus.empty();
        success                 *= no_external.subsystems.size() == 1;
        if (no_external.subsystems.size() == 1)
        {
          success *= no_external.subsystems[0].bus.size() == 3;
          success *= no_external.subsystems[0].branch.size() == 3;
        }

        // A root partition without children keeps the case flat.
        const DataT no_children  = apply(flat_data, R"({
          "name": "root",
          "internal": [0, 1, 2],
          "partitions": []
        })");
        success                 *= no_children.subsystems.empty();
        success                 *= no_children.bus.size() == 3;
        success                 *= no_children.branch.size() == 3;
        success                 *= no_children.genrou.size() == 2;
        success                 *= no_children.loadz.size() == 1;

        // A partition with no lists holds nothing.
        const DataT empty_leaf  = apply(flat_data, R"({
          "name": "root",
          "partitions": [
            { "name": "A", "internal": [0, 1, 2] },
            { "name": "B" }
          ]
        })");
        success                *= empty_leaf.subsystems.size() == 2;
        if (empty_leaf.subsystems.size() == 2)
        {
          success *= empty_leaf.subsystems[1].bus.empty();
          success *= empty_leaf.subsystems[1].branch.empty();
        }

        return success.report(__func__);
      }

      /**
       * @brief Devices connected by a signal go to the first partition listing all their buses.
       *
       * The generator on bus 1 fits partition A on its own, but its signal
       * peer on bus 2 fits only partition B.
       */
      TestOutcome signalGroupPlacement()
      {
        TestStatus success = true;

        DataT flat_data =
            PhasorDynamics::parseSystemModelData(std::string("ThreeBusBasic.case.json"));

        SignalDataT signal;
        signal.signal_id = 10;
        flat_data.signal.push_back(signal);
        flat_data.genrou[0].signal_outputs[PhasorDynamics::GenrouSignalOutputs::speed] = 10;
        flat_data.genrou[1].signal_inputs[PhasorDynamics::GenrouSignalInputs::pmech]   = 10;

        const DataT top_data = apply(flat_data, R"({
          "name": "root",
          "partitions": [
            { "name": "A", "external": [0, 1] },
            { "name": "B", "internal": [2], "external": [0, 1] }
          ]
        })");

        success *= top_data.signal.empty();
        success *= top_data.subsystems.size() == 2;
        if (top_data.subsystems.size() == 2)
        {
          success *= top_data.subsystems[0].genrou.empty();
          success *= top_data.subsystems[0].signal.empty();
          success *= top_data.subsystems[1].genrou.size() == 2;
          success *= top_data.subsystems[1].signal.size() == 1;
        }

        return success.report(__func__);
      }

      /**
       * @brief Badly formed partition files are rejected.
       */
      TestOutcome invalidFiles()
      {
        TestStatus success = true;

        const DataT flat_data =
            PhasorDynamics::parseSystemModelData(std::string("ThreeBusBasic.case.json"));

        const auto previous_verbosity = Log::verbosity();
        // Suppress expected errors from the rejected files below.
        // Use EVERYTHING to inspect those diagnostics.
        Log::setVerbosity(Log::Verbosity::NONE);

        // File structure
        success *= applyThrows(flat_data, R"({ "internal": [0, 1, 2] })");
        success *= applyThrows(flat_data, R"({ "name": "", "internal": [0, 1, 2] })");
        success *= applyThrows(flat_data, R"({
          "name": "root",
          "internal": [0, 1, 2],
          "description": "extra"
        })");
        success *= applyThrows(flat_data, R"({ "name": "root", "partitions": 1 })");
        success *= applyThrows(flat_data, R"({ "name": "root", "internal": 0 })");
        success *= applyThrows(flat_data, R"({ "name": "root", "internal": [-1, 0, 1, 2] })");
        success *= applyThrows(flat_data, R"({ "name": "root", "internal": ["0", 1, 2] })");

        // Bus lists on a partition with children are not supported.
        success *= applyThrows(flat_data, R"({
          "name": "root",
          "external": [0],
          "partitions": [{ "name": "A", "internal": [1, 2] }]
        })");

        // Every bus listed, and only case buses
        success *= applyThrows(flat_data, R"({ "name": "root", "internal": [0, 1] })");
        success *= applyThrows(flat_data, R"({ "name": "root", "internal": [0, 1, 2, 7] })");

        // Internal and external lists
        success *= applyThrows(flat_data, R"({ "name": "root", "internal": [0, 0, 1, 2] })");
        success *= applyThrows(flat_data, R"({ "name": "root", "external": [0, 0, 1, 2] })");
        success *= applyThrows(flat_data, R"({
          "name": "root",
          "internal": [0, 1, 2],
          "external": [0]
        })");
        success *= applyThrows(flat_data, R"({
          "name": "root",
          "partitions": [
            { "name": "A", "internal": [0, 1, 2] },
            { "name": "B", "internal": [2] }
          ]
        })");
        success *= applyThrows(flat_data, R"({
          "name": "root",
          "partitions": [
            { "name": "A", "internal": [0, 1, 2] },
            { "name": "B", "external": [2] }
          ]
        })");
        success *= applyThrows(flat_data, R"({
          "name": "root",
          "partitions": [
            { "name": "A", "external": [2] },
            { "name": "B", "internal": [0, 1, 2] }
          ]
        })");

        // A branch between buses of two partitions
        success *= applyThrows(flat_data, R"({
          "name": "root",
          "partitions": [
            { "name": "A", "internal": [0] },
            { "name": "B", "internal": [1, 2] }
          ]
        })");

        // Devices connected by a signal with no partition listing all their buses
        DataT no_common_partition = flat_data;
        no_common_partition.branch.pop_back();
        no_common_partition.genrou[0].signal_outputs[PhasorDynamics::GenrouSignalOutputs::speed] = 10;
        no_common_partition.genrou[1].signal_inputs[PhasorDynamics::GenrouSignalInputs::pmech]   = 10;

        success *= applyThrows(no_common_partition, R"({
          "name": "root",
          "partitions": [
            { "name": "A", "internal": [1], "external": [0] },
            { "name": "B", "internal": [2], "external": [0] }
          ]
        })");

        // A device using a signal of a bus
        DataT bus_signal                                                              = flat_data;
        bus_signal.bus[2].signal_outputs[PhasorDynamics::BusSignalOutputs::vr]        = 20;
        bus_signal.genrou[1].signal_inputs[PhasorDynamics::GenrouSignalInputs::pmech] = 20;

        success *= applyThrows(bus_signal, R"({ "name": "root", "internal": [0, 1, 2] })");

        // A device with no bus, directly or by signals
        DataT no_bus = flat_data;
        no_bus.constant_source.emplace_back();
        success *= applyThrows(no_bus, R"({ "name": "root", "internal": [0, 1, 2] })");

        // A device on a bus outside the case
        DataT unknown_bus                          = flat_data;
        unknown_bus.loadz[0].buses.begin()->second = 99;

        success *= applyThrows(unknown_bus, R"({ "name": "root", "internal": [0, 1, 2] })");

        // A case that already has subsystems
        DataT hierarchical = flat_data;
        hierarchical.subsystems.emplace_back();
        success *= applyThrows(hierarchical, R"({ "name": "root", "internal": [0, 1, 2] })");

        Log::setVerbosity(previous_verbosity);

        return success.report(__func__);
      }

      /**
       * @brief Apply the shipped WECC240 partition file and compare with the flat system.
       */
      TestOutcome wecc240()
      {
        TestStatus success = true;

        const DataT flat_data =
            PhasorDynamics::parseSystemModelData(std::string("WECC240.case.json"));
        const PhasorDynamics::PartitionData<IdxT> partition =
            PhasorDynamics::parsePartitionData(std::filesystem::path("WECC240.partition.json"));
        const DataT top_data = PlacementT::apply(flat_data, partition);

        success *= top_data.bus.size() == 12;
        success *= top_data.subsystems.size() == 4;
        if (top_data.subsystems.size() == 4)
        {
          success *= top_data.subsystems[0].bus.size() == 58;
          success *= top_data.subsystems[1].bus.size() == 83;
          success *= top_data.subsystems[2].bus.size() == 57;
          success *= top_data.subsystems[3].bus.size() == 33;

          size_t num_branches = 0;
          size_t num_genrous  = 0;
          size_t num_signals  = top_data.signal.size();
          for (const auto& subsystem : top_data.subsystems)
          {
            num_branches += subsystem.branch.size();
            num_genrous  += subsystem.genrou.size();
            num_signals  += subsystem.signal.size();
          }
          success *= num_branches == flat_data.branch.size();
          success *= num_genrous == flat_data.genrou.size();
          success *= num_signals == flat_data.signal.size();

          SystemT flat(flat_data);
          SystemT top(top_data);
          compareSystems(flat_data, flat, top, success);
        }

        return success.report(__func__);
      }

    private:
      /// Apply the partition file contents in `text` to `data`
      static DataT apply(const DataT& data, const std::string& text)
      {
        std::stringstream stream(text);
        return PlacementT::apply(data, PhasorDynamics::parsePartitionData(stream));
      }

      /// Whether applying the partition file contents in `text` to `data` throws
      static bool applyThrows(const DataT& data, const std::string& text)
      {
        bool throws = false;
        // Wrong JSON types throw exceptions of the JSON library.
        try
        {
          apply(data, text);
        }
        catch (const std::exception&)
        {
          throws = true;
        }
        return throws;
      }

      /// Compare variables and residuals by bus ID and component label
      static void compareSystems(const DataT& flat_data,
                                 SystemT&     flat,
                                 SystemT&     top,
                                 TestStatus&  success)
      {
        // External buses may sum device contributions in a different order.
        const RealT tol = 100 * std::numeric_limits<RealT>::epsilon();

        for (auto* system : {&flat, &top})
        {
          success *= system->allocate() == 0;
          success *= system->initialize() == 0;
          success *= system->evaluateResidual() == 0;
        }

        success *= top.size() == flat.size();
        if (top.size() == flat.size())
        {
          const auto* y_flat = flat.y().getData();
          const auto* y      = top.y().getData();
          const auto* f_flat = flat.getResidual().getData();
          const auto* f      = top.getResidual().getData();
          for (const auto& bus_data : flat_data.bus)
          {
            BusT* flat_bus = flat.getBus(bus_data.bus_id);
            BusT* bus      = top.getBus(bus_data.bus_id);
            for (IdxT variable = 0; variable < bus->size(); ++variable)
            {
              const IdxT y_flat_index = flat_bus->getVariableIndex(variable);
              const IdxT y_index      = bus->getVariableIndex(variable);
              const IdxT f_flat_index = flat_bus->getResidualIndex(variable);
              const IdxT f_index      = bus->getResidualIndex(variable);

              success *= isEqual(y[y_index], y_flat[y_flat_index], tol);
              success *= isEqual(f[f_index], f_flat[f_flat_index], tol);
            }
          }

          for (IdxT component_index = 0; component_index < flat.numComponents(); ++component_index)
          {
            ComponentT* flat_component = flat.getComponent(component_index);
            ComponentT* component      = top.getComponent(flat.getComponentLabel(component_index));
            for (IdxT variable = 0; variable < flat_component->size(); ++variable)
            {
              const IdxT y_flat_index = flat_component->getVariableIndex(variable);
              const IdxT y_index      = component->getVariableIndex(variable);
              const IdxT f_flat_index = flat_component->getResidualIndex(variable);
              const IdxT f_index      = component->getResidualIndex(variable);

              success *= isEqual(y[y_index], y_flat[y_flat_index]);
              success *= isEqual(f[f_index], f_flat[f_flat_index]);
            }
          }
        }
      }
    };
  } // namespace Testing
} // namespace GridKit
