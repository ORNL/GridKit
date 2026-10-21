#pragma once

#include <cstddef>
#include <set>
#include <string_view>
#include <type_traits>

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentLibrary.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit
{
  namespace Testing
  {
    enum class ExampleParameters : size_t
    {
      gain,
      time_constant,
    };

    enum class ExampleBuses : size_t
    {
      bus,
    };

    enum class ExampleSignalInputs : size_t
    {
    };

    enum class ExampleSignalOutputs : size_t
    {
      output,
    };

    enum class ExampleMonitorableVariables : size_t
    {
      state,
    };

    template <typename RealT, typename IdxT>
    struct ExampleData : PhasorDynamics::ComponentData<RealT,
                                                       IdxT,
                                                       ExampleParameters,
                                                       ExampleBuses,
                                                       ExampleSignalInputs,
                                                       ExampleSignalOutputs,
                                                       ExampleMonitorableVariables>
    {
      static constexpr typename ExampleData::Manifest MANIFEST{
          .class_name  = "Example",
          .description = "Example component",
          .parameters  = {{
              {ExampleParameters::gain, "p.u.", "K", "Gain"},
              {ExampleParameters::time_constant, "s", "T", "Time constant"},
          }},
          .buses       = {{
              {ExampleBuses::bus, "Terminal bus"},
          }},
          .inputs      = {},
          .outputs     = {{
              {ExampleSignalOutputs::output, "Output signal"},
          }},
          .monitors    = {{
              {ExampleMonitorableVariables::state, "p.u.", "State"},
          }},
      };

      static_assert(MANIFEST.complete());
    };

    using ExampleDataT = ExampleData<double, size_t>;

    /// Single-entry bus table left out, so its entry is value-initialized and undescribed
    inline constexpr ExampleDataT::Manifest INCOMPLETE_MANIFEST{
        .class_name  = "Example",
        .description = "Example component",
        .parameters  = ExampleDataT::MANIFEST.parameters,
        .outputs     = ExampleDataT::MANIFEST.outputs,
        .monitors    = ExampleDataT::MANIFEST.monitors,
    };
    static_assert(!INCOMPLETE_MANIFEST.complete());

    /// Checks `ComponentManifest` entries and the `ComponentLibrary` model list.
    class ComponentManifestTests
    {
    public:
      TestOutcome manifest()
      {
        TestStatus success = true;

        const auto& manifest = ExampleDataT::MANIFEST;

        success *= (manifest.parameters.size() == Utilities::enum_size<ExampleParameters>());
        success *= (manifest.parameters[1].name() == "time_constant");
        success *= (manifest.buses[0].name() == "bus");
        success *= (manifest.outputs[0].name() == "output");
        success *= (manifest.monitors[0].name() == "state");

        // Compiles the empty-enumeration path of name(), the loop body never runs
        success *= manifest.inputs.empty();
        for (const auto& entry : manifest.inputs)
        {
          success *= !entry.name().empty();
        }

        return success.report(__func__);
      }

      TestOutcome library()
      {
        TestStatus success = true;

        std::size_t                visited = 0;
        std::set<std::string_view> class_names;
        PhasorDynamics::ComponentLibrary<>::forEach(
            [&]<typename ModelT, typename DataT>(std::type_identity<ModelT>, std::type_identity<DataT>)
            {
              class_names.insert(DataT::MANIFEST.class_name);
              ++visited;
            });

        // Every model has its own class, so the parser can dispatch on it
        success *= (class_names.size() == visited);

        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
