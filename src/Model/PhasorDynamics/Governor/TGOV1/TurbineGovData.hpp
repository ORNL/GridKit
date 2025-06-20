/**
 * @file GenrouData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for TGOV1 (transmission lines)
 *
 */
#pragma once

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Contains modeling data for a Genrou generator model.
     *
     * @tparam RealT Real parameter data type
     * @tparam IdxT  Integer parameter data type
     *
     * Integer parameters are of the same type as matrix and vector indices.
     *
     * @todo Decide on naming scheme for model parameters.
     */
    template <typename RealT, typename IdxT>
    struct GenrouData
    {
      RealT R{0.0};     ///< R
      RealT Pvmin{0.0}; ///< Min Valve Power
      RealT Pvmax{0.0}; ///< Max Valve Power
      RealT T1{0.0};    ///< 
      RealT T2{0.0};    ///< 
      RealT T3{0.0};    ///< 
      RealT Dt{0.0};    ///< 
    };
  } // namespace PhasorDynamics
} // namespace GridKit
