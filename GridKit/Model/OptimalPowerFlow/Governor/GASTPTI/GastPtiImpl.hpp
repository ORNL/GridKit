/**
 * @file GastPtiImpl.hpp
 * @brief Definition of an optimal power flow GASTPTI governor.
 */

#pragma once

#include <algorithm>
#include <type_traits>

#include <GridKit/Model/OptimalPowerFlow/ComponentImpl.hpp>
#include <GridKit/Model/OptimalPowerFlow/Governor/GASTPTI/GastPti.hpp>
#include <GridKit/Model/ParameterReader.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Governor
    {
      template <typename scalar_type, typename index_type>
      GastPti<scalar_type, index_type>::GastPti(const ModelDataT& data)
        : ComponentT(std::type_identity<GastPti>{}, data.id)
      {
        initializeParameters(data);
      }

      template <typename scalar_type, typename index_type>
      Model::ConfigurationChecks GastPti<scalar_type, index_type>::verify() const
      {
        Model::ConfigurationChecks checks("GastPti");
        checks.check(Trate_ > ZERO<RealT>, "Trate must be positive");
        checks.check(Vmin_ <= Vmax_, "Vmin must not exceed Vmax");
        return checks;
      }

      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::tagInequality()
      {
        const auto XF = static_cast<size_t>(GastPtiInternalConstraints::XF);

        tag_[XF] = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int GastPti<scalar_type, index_type>::initialize()
      {
        const auto XF = static_cast<size_t>(GastPtiInternalConstraints::XF);

        RealT* g_lower = g_lower_.getData();
        RealT* g_upper = g_upper_.getData();

        g_lower[XF] = Vmin_;
        g_upper[XF] = std::min(Vmax_, At_);
        return 0;
      }

      template <typename scalar_type, typename index_type>
      [[gnu::always_inline]] void
      GastPti<scalar_type, index_type>::internalConstraints([[maybe_unused]] const ScalarT* x,
                                                            const ScalarT*                  external,
                                                            ScalarT*                        g) const
      {
        const auto PM = static_cast<size_t>(GastPtiExternalVariables::PM);
        const auto XF = static_cast<size_t>(GastPtiInternalConstraints::XF);

        g[XF] = this->toComponentBase(external[PM]);
      }

      template <typename scalar_type, typename index_type>
      void GastPti<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Params = GastPtiParameters;

        const Model::ParameterReader reader(data, "GastPti");

        Trate_ = reader.readReal(Params::Trate).value_or(Trate_);
        Vmax_  = reader.readReal(Params::Vmax).value_or(Vmax_);
        Vmin_  = reader.readReal(Params::Vmin).value_or(Vmin_);
        At_    = reader.readReal(Params::At).value_or(At_);
        setDerivedParameters();
      }

      template <typename scalar_type, typename index_type>
      void GastPti<scalar_type, index_type>::setDerivedParameters()
      {
        this->setComponentBase(MEGA<RealT> * Trate_);
      }
    } // namespace Governor
  } // namespace OptimalPowerFlow
} // namespace GridKit
