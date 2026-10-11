#pragma once

#include <cstddef>
#include <tuple>
#include <type_traits>

#include <GridKit/Model/PhasorDynamics/ComponentDefinition.hpp>
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

    using ExampleDefinition = PhasorDynamics::ComponentDefinition<ExampleParameters,
                                                                  ExampleBuses,
                                                                  ExampleSignalInputs,
                                                                  ExampleSignalOutputs,
                                                                  ExampleMonitorableVariables>;

    inline constexpr ExampleDefinition EXAMPLE_DEFINITION{
        .device_class          = "Example",
        .description           = "Example component",
        .parameters            = {{
            {ExampleParameters::gain, "p.u.", "K", "Gain"},
            {ExampleParameters::time_constant, "s", "T", "Time constant"},
        }},
        .buses                 = {{
            {ExampleBuses::bus, "Terminal bus"},
        }},
        .signal_inputs         = {},
        .signal_outputs        = {{
            {ExampleSignalOutputs::output, "Output signal"},
        }},
        .monitorable_variables = {{
            {ExampleMonitorableVariables::state, "p.u.", "State"},
        }},
    };
    static_assert(EXAMPLE_DEFINITION.complete());

    /// Single-entry bus section left out, so its row is value-initialized and undescribed
    inline constexpr ExampleDefinition INCOMPLETE_DEFINITION{
        .device_class          = "Example",
        .description           = "Example component",
        .parameters            = EXAMPLE_DEFINITION.parameters,
        .signal_outputs        = EXAMPLE_DEFINITION.signal_outputs,
        .monitorable_variables = EXAMPLE_DEFINITION.monitorable_variables,
    };
    static_assert(!INCOMPLETE_DEFINITION.complete());

    /// Checks `ComponentDefinition` rows and the component definition registry.
    class ComponentDefinitionTests
    {
    public:
      TestOutcome definition()
      {
        TestStatus success = true;

        success *= (EXAMPLE_DEFINITION.parameters.size() == Utilities::enum_size<ExampleParameters>());
        success *= (EXAMPLE_DEFINITION.parameters[1].name() == "time_constant");
        success *= (EXAMPLE_DEFINITION.buses[0].name() == "bus");
        success *= (EXAMPLE_DEFINITION.signal_outputs[0].name() == "output");
        success *= (EXAMPLE_DEFINITION.monitorable_variables[0].name() == "state");

        // Compiles the empty-enumeration path of name(), the loop body never runs
        success *= EXAMPLE_DEFINITION.signal_inputs.empty();
        for (const auto& row : EXAMPLE_DEFINITION.signal_inputs)
        {
          success *= !row.name().empty();
        }

        return success.report(__func__);
      }

      TestOutcome registry()
      {
        TestStatus success = true;

        std::size_t visited = 0;
        PhasorDynamics::forEachComponentDefinition([&](const auto& definition)
                                                   {
                                                     success *= !definition.device_class.empty();
                                                     ++visited; });

        using RegistryT  = std::remove_cvref_t<decltype(PhasorDynamics::COMPONENT_DEFINITIONS)>;
        success         *= (visited == std::tuple_size_v<RegistryT>);

        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
