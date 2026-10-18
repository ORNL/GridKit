/**
 * @file LoadData.hpp
 * @brief Modeling data for optimal power flow loads.
 */

#pragma once

#include <cstddef>

#include <GridKit/Model/OptimalPowerFlow/ComponentData.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /// Parameters for a load, each a demand at 1 p.u. voltage
    enum class LoadParameters : size_t
    {
      P,          ///< \f$P\f$ Constant-power active demand [p.u.]
      Q,          ///< \f$Q\f$ Constant-power reactive demand [p.u.]
      Pi,         ///< \f$P_I\f$ Constant-current active demand [p.u.]
      Qi,         ///< \f$Q_I\f$ Constant-current reactive demand [p.u.]
      Pz,         ///< \f$P_Z\f$ Constant-impedance active demand [p.u.]
      Qz,         ///< \f$Q_Z\f$ Constant-impedance reactive demand [p.u.]
      in_service, ///< \f$u\f$ Service status
    };

    /// Ports of a load
    enum class LoadPorts : size_t
    {
      bus, ///< Port bus
    };

    /**
     * @brief Contains modeling data for a Load
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer data type
     */
    template <typename real_type, typename index_type>
    using LoadData = ComponentData<real_type, index_type, LoadParameters, LoadPorts>;
  } // namespace OptimalPowerFlow
} // namespace GridKit
