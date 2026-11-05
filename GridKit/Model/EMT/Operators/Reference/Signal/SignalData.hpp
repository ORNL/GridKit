/**
 * @file SignalData.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Modeling data for EMT signals.
 */

#pragma once

#include <cstddef>
#include <tuple>

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Parameters of a `Signal`
    enum class SignalParameters : size_t
    {
      N,    ///< \f$N\f$ Number of signals
      held, ///< \f$\mathrm{held}\f$ Whether the signals are held
    };

    /// Data of a `Signal`
    template <typename RealT, typename IdxT>
    struct SignalData
    {
      using Parameters = SignalParameters;

      IdxT N{0};
      bool held{false};

      static constexpr std::tuple PARAMETERS{Parameter{Parameters::N, &SignalData::N}, Parameter{Parameters::held, &SignalData::held}};
    };
  } // namespace EMT
} // namespace GridKit
