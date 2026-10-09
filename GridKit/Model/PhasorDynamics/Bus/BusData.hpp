/**
 * @file BusData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data for buses (nodes)
 *
 */
#pragma once

#include <cstddef>

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Parameters for a bus
    enum class BusParameters : size_t
    {
      kv, ///< Voltage base [kV]
    };

    /// Signal inlets that a bus could have.
    enum class BusSignalInputs : size_t
    {
      // TODO: this is a bit of a hack but there's no better solution without
      //       reworking how buses are parsed
      vr,
      vi,
      ir,
      ii,
    };

    enum class BusSignalOutputs : size_t
    {
      // TODO: likewise
      vr,
      vi,
      ir,
      ii,
    };

    /// Indices of the variables able to be monitored on this component
    enum class BusMonitorableVariables : size_t
    {
      Vr,
      Vi,
      Vm,
      Va
    };

    /**
     * @brief Contains modeling data for a Bus
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     */
    template <typename real_type, typename index_type>
    struct BusData : public ComponentData<real_type,
                                          index_type,
                                          BusParameters,
                                          Utilities::EmptyEnum,
                                          BusSignalInputs,
                                          BusSignalOutputs,
                                          BusMonitorableVariables>
    {
      using RealT = real_type;
      using IdxT  = index_type;

      std::string name; ///< A name given to this bus

      // TODO: these should become parameters
      RealT Vr0{1.0}; ///< Initial value for the real bus voltage
      RealT Vi0{0.0}; ///< Initial value for the imaginary bus voltage

      IdxT bus_id{0}; ///< The unique ID of the bus

      /// Enumeration over the kinds of bus this data structure can be for
      enum class BusType
      {
        INVALID,
        DEFAULT,
        SLACK,
        SIGNAL_VOLTAGE_OUT, ///< Bus with voltage signal outputs and current signal inputs
        SIGNAL_VOLTAGE_IN,  ///< Bus with voltage signal inputs and current signal outputs
      };

      BusType bus_type{BusType::INVALID}; ///< The kind of bus this data is for
    };

  } // namespace PhasorDynamics
} // namespace GridKit
