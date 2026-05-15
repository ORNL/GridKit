/**
 * @file CoSimData.hpp
 * @author Philip Fackler (facklerpw@ornl.gov)
 *
 * @brief Data structure for CoSim Data
 *
 */
#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Connector
    {
      /**
       * @brief Parameter keys for CoSim component
       *
       * These enum values serve as keys for the parameters map in ComponentData.
       */
      enum class CoSimParameters
      {
        NONE,
      };

      /**
       * @brief CoSim ports
       */
      enum class CoSimPorts
      {
        bus,
        vr,
        vi,
        ir,
        ii,
      };

      /**
       * @brief Placeholder enum for CoSim monitorable variables
       */
      enum class CoSimMonitorableVariables
      {
        NONE,
      };

      /**
       * @brief Modeling data for CoSim Connector using ComponentData base
       *
       * @tparam RealT Real number type (e.g., double)
       * @tparam IdxT  Index type (e.g., size_t)
       */
      template <typename RealT, typename IdxT>
      struct CoSimData : public ComponentData<RealT, IdxT, CoSimParameters, CoSimPorts, CoSimMonitorableVariables>
      {
        CoSimData() = default;

        using Parameters           = CoSimParameters;
        using Ports                = CoSimPorts;
        using MonitorableVariables = CoSimMonitorableVariables;
      };

    } // namespace Connector
  } // namespace PhasorDynamics
} // namespace GridKit
