/**
 * @file HygovImpl.hpp
 * @brief Definition of an optimal power flow HYGOV governor.
 */

#pragma once

#include <cmath>
#include <type_traits>

#include <GridKit/Model/OptimalPowerFlow/ComponentImpl.hpp>
#include <GridKit/Model/OptimalPowerFlow/Governor/HYGOV/Hygov.hpp>
#include <GridKit/Model/ParameterReader.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Governor
    {
      template <typename scalar_type, typename index_type>
      Hygov<scalar_type, index_type>::Hygov(const ModelDataT& data)
        : ComponentT(std::type_identity<Hygov>{}, data.id)
      {
        initializeParameters(data);
      }

      template <typename scalar_type, typename index_type>
      Model::ConfigurationChecks Hygov<scalar_type, index_type>::verify() const
      {
        Model::ConfigurationChecks checks("Hygov");
        checks.check(Trate_ > ZERO<RealT>, "Trate must be positive");
        checks.check(Gmin_ < Gmax_, "Gmin must be less than Gmax");
        checks.check(At_ > ZERO<RealT>, "At must be positive");
        checks.check(Hdam_ > ZERO<RealT>, "Hdam must be positive");
        return checks;
      }

      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::tagInequality()
      {
        const auto GATE = static_cast<size_t>(HygovInternalConstraints::GATE);

        tag_[GATE] = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Hygov<scalar_type, index_type>::initialize()
      {
        const auto GATE = static_cast<size_t>(HygovInternalConstraints::GATE);

        RealT* g_lower = g_lower_.getData();
        RealT* g_upper = g_upper_.getData();

        g_lower[GATE] = Gmin_;
        g_upper[GATE] = Gmax_;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      [[gnu::always_inline]] void
      Hygov<scalar_type, index_type>::internalConstraints([[maybe_unused]] const ScalarT* x,
                                                          const ScalarT*                  external,
                                                          ScalarT*                        g) const
      {
        const auto PM   = static_cast<size_t>(HygovExternalVariables::PM);
        const auto GATE = static_cast<size_t>(HygovInternalConstraints::GATE);

        g[GATE] = (this->toComponentBase(external[PM]) / (At_ * Hdam_) + Qnl_) / std::sqrt(Hdam_);
      }

      template <typename scalar_type, typename index_type>
      void Hygov<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Params = HygovParameters;

        const Model::ParameterReader reader(data, "Hygov");

        Trate_ = reader.readReal(Params::Trate).value_or(Trate_);
        Gmax_  = reader.readReal(Params::Gmax).value_or(Gmax_);
        Gmin_  = reader.readReal(Params::Gmin).value_or(Gmin_);
        At_    = reader.readReal(Params::At).value_or(At_);
        Qnl_   = reader.readReal(Params::Qnl).value_or(Qnl_);
        Hdam_  = reader.readReal(Params::Hdam).value_or(Hdam_);
        setDerivedParameters();
      }

      template <typename scalar_type, typename index_type>
      void Hygov<scalar_type, index_type>::setDerivedParameters()
      {
        this->setComponentBase(MEGA<RealT> * Trate_);
      }
    } // namespace Governor
  } // namespace OptimalPowerFlow
} // namespace GridKit
