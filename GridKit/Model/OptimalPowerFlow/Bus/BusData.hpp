/**
 * @file BusData.hpp
 * @brief Modeling data for optimal power flow buses.
 */

#pragma once

#include <cstddef>
#include <format>
#include <string>

#include <GridKit/Model/OptimalPowerFlow/ComponentData.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /// Parameters for a bus
    enum class BusParameters : size_t
    {
      Vmin,      ///< \f$V^{\min}\f$ Voltage magnitude lower limit [p.u.]
      Vmax,      ///< \f$V^{\max}\f$ Voltage magnitude upper limit [p.u.]
      G,         ///< \f$G\f$ Shunt conductance to ground [p.u.]
      B,         ///< \f$B\f$ Shunt susceptance to ground [p.u.]
      reference, ///< Angle reference, with zero imaginary voltage
    };

    /// Data of a `Bus`, identified by its number
    template <typename real_type, typename index_type>
    struct BusData
    {
      using RealT      = real_type;
      using IdxT       = index_type;
      using Parameters = BusParameters;

      /// Bus number
      IdxT number{0};

      ParameterMap<Parameters, RealT, IdxT> parameters;

      /// Bus `id`, from its number
      std::string id() const
      {
        return std::format("bus {}", number);
      }
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
