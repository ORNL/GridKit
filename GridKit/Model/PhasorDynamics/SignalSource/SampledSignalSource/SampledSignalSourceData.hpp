#pragma once

#include <filesystem>
#include <string>
#include <utility>
#include <vector>

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace SignalSource
    {
      enum class SampledSignalSourceParameters
      {
        scale,
        offset,
      };

      enum class SampledSignalSourcePorts
      {
        output,
      };

      enum class SampledSignalSourceMonitorableVariables
      {
        out,
      };

      template <typename real_type, typename index_type>
      struct SampledSignalSourceData : public ComponentData<real_type,
                                                            index_type,
                                                            SampledSignalSourceParameters,
                                                            SampledSignalSourcePorts,
                                                            SampledSignalSourceMonitorableVariables>
      {
        using RealT                = real_type;
        using IdxT                 = index_type;
        using Parameters           = SampledSignalSourceParameters;
        using Ports                = SampledSignalSourcePorts;
        using MonitorableVariables = SampledSignalSourceMonitorableVariables;

        SampledSignalSourceData() = default;

        std::string                          source_type{"samples"};
        std::filesystem::path                file;
        std::string                          time_column{"t"};
        std::string                          value_column{"u"};
        std::vector<std::pair<RealT, RealT>> samples;
      };
    } // namespace SignalSource
  } // namespace PhasorDynamics
} // namespace GridKit
