/**
 * @file BusImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT bus.
 */

#pragma once

#include <algorithm>

#include <GridKit/Model/EMT/Bus/Bus.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    Bus<scalar_type, index_type>::Bus(const ModelDataT& data)
      : N_(data.N),
        kv_(data.kv)
    {
      this->setLayout({N_}, {});
    }

    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Bus<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Bus");
      checks.check(N_ > 0, "N must be positive");
      checks.check(kv_ > ZERO<RealT>, "kv must be positive");
      return checks;
    }

    /// Voltages at the bus base voltage
    template <typename scalar_type, typename index_type>
    int Bus<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      std::ranges::fill(this->absoluteTolerance(BusInternalVariables::v), rel_tol * KILO<RealT> * kv_);
      return 0;
    }
  } // namespace EMT
} // namespace GridKit
