/**
 * @file Tgov1Data.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for TGOV1
 *
 */
#pragma once

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      /**
       * @brief Contains modeling data for a TGOV1 Governor model.
       *
       * @tparam RealT Real parameter data type
       * @tparam IdxT  Integer parameter data type
       *
       * Integer parameters are of the same type as matrix and vector indices.
       *
       * @todo Decide on naming scheme for model parameters.
       */
      template <typename RealT, typename IdxT>
      struct Tgov1Data
      {
        RealT R{0.05};    ///< Droop Constant
        RealT T1{0.5};    ///< Valve Time Delay
        RealT T2{2.5};    ///< Turbine Numerator Time Constant
        RealT T3{7.5};    ///< Turbine Delay
        RealT Pvmax{1.0}; ///< Max Valve Power
        RealT Pvmin{0.0}; ///< Min Valve Power
        RealT Dt{0.0};    ///<
      };
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
