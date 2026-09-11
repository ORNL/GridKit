/**
 * @file BusData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data for buses (nodes)
 *
 */
#pragma once

#include <map>
#include <set>
#include <string>
#include <type_traits>
#include <variant>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Indices of the variables able to be monitored on this component
    enum class BusMonitorableVariables : size_t
    {
      Vr,
      Vi,
      Vm,
      Va
    };

    /// Parameters for a bus
    enum class BusParameters : size_t
    {
      kv, ///< Voltage base [kV]
    };

    /**
     * @brief Contains modeling data for a Bus
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     *
     * @todo Decide on naming scheme for model parameters.
     */
    template <typename real_type, typename index_type>
    struct BusData
    {
      using RealT = real_type;
      using IdxT  = index_type;

      using Parameters     = BusParameters;
      using ParameterValue = std::variant<bool, RealT, IdxT>;

      std::string name; ///< A name given to this bus

      RealT Vr0{1.0}; ///< Initial value for the real bus voltage
      RealT Vi0{0.0}; ///< Initial value for the imaginary bus voltage

      IdxT bus_id{0}; ///< The unique ID of the bus

      /// Enumeration over the kinds of bus this data structure can be for
      enum class BusType
      {
        INVALID,
        DEFAULT,
        SLACK,
      };

      BusType bus_type{BusType::INVALID}; ///< The kind of bus this data is for

      std::map<Parameters, ParameterValue> parameters; ///< Mapping of parameters to parameter values

      /// Alias
      using MonitorableVariables = BusMonitorableVariables;

      /// Set indicating the variables being monitored
      std::set<MonitorableVariables> monitored_variables;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
