/**
 * @file GasBusData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for gas buses.
 */
#pragma once

#include <string>

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Gas
    {
      /// Variables able to be monitored for a gas bus
      enum class GasBusMonitorableVariables : size_t
      {
        p, ///< \f$p\f$ Bus pressure [Pa]
      };

      /**
       * @brief Contains modeling data for a GasBus or GasBusInfinite
       *
       * The `device_class` field holds the bus class from the case file.
       *
       * @tparam real_type  Real parameter data type
       * @tparam index_type Integer parameter data type
       */
      template <typename real_type, typename index_type>
      struct GasBusData : public ComponentData<real_type,
                                               index_type,
                                               Utilities::EmptyEnum,
                                               Utilities::EmptyEnum,
                                               Utilities::EmptyEnum,
                                               Utilities::EmptyEnum,
                                               GasBusMonitorableVariables>
      {
        using RealT = real_type;
        using IdxT  = index_type;

        std::string name;      ///< A name given to this gas bus
        RealT       p0{0.0};   ///< Initial bus pressure [Pa]
        IdxT        bus_id{0}; ///< The unique ID of the gas bus
      };
    } // namespace Gas
  } // namespace PhasorDynamics
} // namespace GridKit
