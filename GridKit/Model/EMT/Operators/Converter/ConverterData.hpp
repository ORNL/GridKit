#pragma once

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    enum class ConverterParameters
    {
      NONE,
    };

    enum class ConverterInputs : size_t
    {
      sa,
      sb,
      sc,
      vdc,
      SIZE,
    };

    enum class ConverterOutputs : size_t
    {
      voa,
      vob,
      voc,
      SIZE,
    };

    enum class ConverterMonitorableVariables
    {
      voa,
      vob,
      voc,
    };

    template <typename real_type, typename index_type>
    struct ConverterData : public ComponentData<real_type,
                                                index_type,
                                                ConverterParameters,
                                                ConverterInputs,
                                                ConverterOutputs,
                                                ConverterMonitorableVariables>
    {
      ConverterData() = default;

      using Parameters           = ConverterParameters;
      using Inputs               = ConverterInputs;
      using Outputs              = ConverterOutputs;
      using MonitorableVariables = ConverterMonitorableVariables;
    };
  } // namespace EMT
} // namespace GridKit
