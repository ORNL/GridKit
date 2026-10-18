/**
 * @file BusData.hpp
 * @brief Modeling data for optimal power flow buses.
 */

#pragma once

#include <cstddef>
#include <map>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /// Parameters for a bus
    enum class BusParameters : size_t
    {
      Vmin, ///< \f$V^{\min}\f$ Voltage magnitude lower limit [p.u.]
      Vmax, ///< \f$V^{\max}\f$ Voltage magnitude upper limit [p.u.]
      G,    ///< \f$G\f$ Shunt conductance to ground [p.u.]
      B,    ///< \f$B\f$ Shunt susceptance to ground [p.u.]
    };

    /**
     * @brief Contains modeling data for a Bus
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer data type
     */
    template <typename real_type, typename index_type>
    struct BusData
    {
      using RealT      = real_type;
      using IdxT       = index_type;
      using Parameters = BusParameters;

      /// Bus number
      IdxT number{0};

      /// Angle reference, with zero imaginary voltage
      bool reference{false};

      /// Mapping of parameters to parameter values
      std::map<Parameters, RealT> parameters;
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
