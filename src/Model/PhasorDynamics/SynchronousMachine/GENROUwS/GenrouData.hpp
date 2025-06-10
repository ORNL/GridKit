/**
 * @file GenrouData.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Modeling data for branches (transmission lines)
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
      IdxT unit_id{0}; ///< Unique unit ID

      RealT p0{0.0};    ///< Initial active power
      RealT q0{0.0};    ///< Initial reactive power
      RealT H{0.0};     ///< Rotor inertia
      RealT D{0.0};     ///< Damping coefficient
      RealT Ra{0.0};    ///< Winding resistance
      RealT Tdop{0.0};  ///< Open circuit direct axis transient time
      RealT Tdopp{0.0}; ///< Open circuit direct axis sub-transient time
      RealT Tqop{0.0};  ///< Open circuit quadrature axis transient
      RealT Tqopp{0.0}; ///< Open circuit quadrature axis sub-transient time
      RealT Xd{0.0};    ///< Direct axis synchronous reactance
      RealT Xdp{0.0};   ///< Direct axis transient reactance
      RealT Xdpp{0.0};  ///< Direct axis sub-transient reactance
      RealT Xq{0.0};    ///< Quadrature axis synchronous reactance
      RealT Xqp{0.0};   ///< Quadrature axis transient reactance
      RealT Xqpp{0.0};  ///< Quadrature axis sub-transient reactance
      RealT Xl{0.0};    ///< Stator leakage reactance
      RealT S10{0.0};   ///< Saturation factor at 1.0 pu flux
      RealT S12{0.0};   ///< Saturation factor at 1.2 pu flux

      IdxT bus_id{0}; ///< Unique ID of the connecting bus
    };
  } // namespace PhasorDynamics
} // namespace GridKit
