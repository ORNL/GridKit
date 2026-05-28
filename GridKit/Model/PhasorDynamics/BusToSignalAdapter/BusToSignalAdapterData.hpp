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
     * @brief BusToSignalAdapter ports
     */
    enum class BusToSignalAdapterPorts
    {
      bus,
      vr,
      vi,
      ir,
      ii,
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
    template <typename RealT, typename IdxT>
    struct BusToSignalAdapterData
      : public ComponentData<RealT,
                             IdxT,
                             BusToSignalAdapterParameters,
                             BusToSignalAdapterPorts,
                             BusToSignalAdapterMonitorableVariables>
    {
      BusToSignalAdapterData() = default;

      using Parameters           = BusToSignalAdapterParameters;
      using Ports                = BusToSignalAdapterPorts;
      using MonitorableVariables = BusToSignalAdapterMonitorableVariables;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
