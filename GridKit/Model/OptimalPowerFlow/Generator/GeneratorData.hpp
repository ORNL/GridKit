/**
 * @file GeneratorData.hpp
 * @brief Modeling data for optimal power flow generators.
 */

#pragma once

#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/ComponentData.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /// Parameters for a generator
    enum class GeneratorParameters : size_t
    {
      Pmin, ///< \f$P^{\min}\f$ Active power lower limit [p.u.]
      Pmax, ///< \f$P^{\max}\f$ Active power upper limit [p.u.]
      Qmin, ///< \f$Q^{\min}\f$ Reactive power lower limit [p.u.]
      Qmax, ///< \f$Q^{\max}\f$ Reactive power upper limit [p.u.]
      c0,   ///< \f$c_0\f$ Constant cost coefficient
      c1,   ///< \f$c_1\f$ Linear cost coefficient
      c2,   ///< \f$c_2\f$ Quadratic cost coefficient
    };

    /// Ports of a generator
    enum class GeneratorPorts : size_t
    {
      bus, ///< Port bus
    };

    /**
     * @brief Contains modeling data for a Generator
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer data type
     */
    template <typename real_type, typename index_type>
    using GeneratorData = ComponentData<real_type, index_type, GeneratorParameters, GeneratorPorts>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
