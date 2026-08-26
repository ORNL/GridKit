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
    enum class BusToSignalAdapterParameters : size_t
    {
    };

    /**
     * @brief BusToSignalAdapter buses
     */
    enum class BusToSignalAdapterBuses : size_t
    {
      bus,
    };

    /**
     * @brief BusToSignalAdapter signal inputs
     */
    enum class BusToSignalAdapterSignalInputs : size_t
    {
      vr_in,
      vi_in,
      ir_in,
      ii_in
    };

    /**
     * @brief BusToSignalAdapter signal outputs
     */
    enum class BusToSignalAdapterSignalOutputs : size_t
    {
      vr_out,
      vi_out,
      ir_out,
      ii_out
    };

    /**
     * @brief Placeholder enum for BusToSignalAdapter monitorable variables
     */
    enum class BusToSignalAdapterMonitorableVariables : size_t
    {
    };

    /**
     * @brief Modeling data for BusToSignalAdapter using ComponentData base
     *
     * @tparam RealT Real number type (e.g., double)
     * @tparam IdxT  Index type (e.g., size_t)
     */
    template <typename real_type, typename index_type>
    using BusToSignalAdapterData =
        ComponentData<real_type,
                      index_type,
                      BusToSignalAdapterParameters,
                      BusToSignalAdapterBuses,
                      BusToSignalAdapterSignalInputs,
                      BusToSignalAdapterSignalOutputs,
                      BusToSignalAdapterMonitorableVariables>;
  } // namespace PhasorDynamics
} // namespace GridKit
