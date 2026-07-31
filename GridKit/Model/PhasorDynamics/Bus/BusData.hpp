/**
 * @file BusData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data for buses (nodes)
 *
 */
#pragma once

#include <map>
#include <optional>
#include <set>
#include <string>
#include <type_traits>
#include <variant>

#include <GridKit/Model/StateData.hpp>

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
    enum class BusParameters
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

      // TODO: Bus-level freq_base and va_base are parsed but not applied.
      // Prefer removing them as bus parameters in a future cleanup.
      std::optional<RealT> freq_base; ///< Override for the system-wide base frequency
      std::optional<RealT> va_base;   ///< Override for the system-wide power base

      /// Initial operating state supplied separately from model parameters
      std::optional<Model::BusState> initial_state;

      /// Alias
      using MonitorableVariables = BusMonitorableVariables;

      /// Set indicating the variables being monitored
      std::set<MonitorableVariables> monitored_variables;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
