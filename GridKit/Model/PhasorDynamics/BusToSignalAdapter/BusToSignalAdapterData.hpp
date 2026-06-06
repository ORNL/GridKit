/**
 * @file BusToSignalAdapterData.hpp
 * @author Philip Fackler (facklerpw@ornl.gov)
 *
 * @brief Data structure for BusToSignalAdapter Data
 *
 */
#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Parameter keys for BusToSignalAdapter component
     *
     * These enum values serve as keys for the parameters map in ComponentData.
     */
    enum class BusToSignalAdapterParameters
    {
      NONE,
    };

    /**
     * @brief BusToSignalAdapter terminals
     */
    enum class BusToSignalAdapterTerminals : size_t
    {
      bus,
      SIZE
    };

    /**
     * @brief BusToSignalAdapter input ports
     */
    enum class BusToSignalAdapterInputPorts : size_t
    {
      ir,
      ii,
      SIZE
    };

    /**
     * @brief BusToSignalAdapter output ports
     */
    enum class BusToSignalAdapterOutputPorts : size_t
    {
      vr,
      vi,
      SIZE
    };

    /**
     * @brief Placeholder enum for BusToSignalAdapter monitorable variables
     */
    enum class BusToSignalAdapterMonitorableVariables
    {
      NONE,
    };

    /**
     * @brief Modeling data for BusToSignalAdapter using ComponentData base
     *
     * @tparam RealT Real number type (e.g., double)
     * @tparam IdxT  Index type (e.g., size_t)
     */
    template <typename real_type, typename index_type>
    struct BusToSignalAdapterData
      : public ComponentData<real_type,
                             index_type,
                             BusToSignalAdapterParameters,
                             BusToSignalAdapterTerminals,
                             BusToSignalAdapterInputPorts,
                             BusToSignalAdapterOutputPorts,
                             BusToSignalAdapterMonitorableVariables>
    {
      BusToSignalAdapterData() = default;

      using Parameters           = BusToSignalAdapterParameters;
      using Terminals            = BusToSignalAdapterTerminals;
      using InputPorts           = BusToSignalAdapterInputPorts;
      using OutputPorts          = BusToSignalAdapterOutputPorts;
      using MonitorableVariables = BusToSignalAdapterMonitorableVariables;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
