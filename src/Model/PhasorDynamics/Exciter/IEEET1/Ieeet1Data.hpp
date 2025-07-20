/**
 * @file Ieeet1Data.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * 
 * @brief Data structure for IEEET1 Data
 *
 */
#pragma once

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /**
       * @brief Contains modeling data for a IEEET1 Exciter model.
       *
       * @tparam RealT Real parameter data type
       * @tparam IdxT  Integer parameter data type
       *
       * Integer parameters are of the same type as matrix and vector indices.
       *
       * @todo Decide on naming scheme for model parameters.
       */
      template <typename RealT, typename IdxT>
      struct Ieeet1Data
      {
        RealT Tr{0};        ///< Time constant for voltage sensing
        RealT Ka{50};       ///< Coefficient for voltage regulation
        RealT Ta{0.04};     ///< Time constant for voltage regulation
        RealT Ke{-0.06};    ///< Coefficient for excitation system
        RealT Te{0.6};      ///< Time constant for excitation system
        RealT Kf{0.09};     ///< Coefficient for feedback
        RealT Tf{1.46};     ///< Time constant for feedback
        RealT Vrmin{-1};    ///< LL to voltage regulation
        RealT Vrmax{1};     ///< HH to voltage regulation
        RealT E1{2.8};      ///< Saturation parameter
        RealT E2{3.373};    ///< Saturation parameter
        RealT Se1{0.04};    ///< Saturation parameter
        RealT Se2{0.33};    ///< Saturation parameter
        RealT Ispdlim{0.0}; ///< Speed limit flag indicator

        // NOTE signal data format TBD
        // IdxT signal_efd{0};
        // IdxT signal_vref{0};
        // IdxT signal_speed{0};
      };
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit