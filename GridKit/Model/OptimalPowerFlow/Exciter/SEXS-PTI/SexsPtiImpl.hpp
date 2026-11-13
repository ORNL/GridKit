/**
 * @file SexsPtiImpl.hpp
 * @brief Definition of an optimal power flow SEXS-PTI exciter.
 */

#pragma once

#include <type_traits>

#include <GridKit/Model/OptimalPowerFlow/ComponentImpl.hpp>
#include <GridKit/Model/OptimalPowerFlow/Exciter/SEXS-PTI/SexsPti.hpp>
#include <GridKit/Model/ParameterReader.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Exciter
    {
      template <typename scalar_type, typename index_type>
      SexsPti<scalar_type, index_type>::SexsPti(const ModelDataT& data)
        : ComponentT(std::type_identity<SexsPti>{}, data.id)
      {
        initializeParameters(data);
      }

      template <typename scalar_type, typename index_type>
      Model::ConfigurationChecks SexsPti<scalar_type, index_type>::verify() const
      {
        Model::ConfigurationChecks checks("SexsPti");
        checks.check(Efdmin_ < Efdmax_, "Efdmin must be less than Efdmax");
        return checks;
      }

      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::tagInequality()
      {
        const auto EFD = static_cast<size_t>(SexsPtiInternalConstraints::EFD);

        tag_[EFD] = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int SexsPti<scalar_type, index_type>::initialize()
      {
        const auto EFD = static_cast<size_t>(SexsPtiInternalConstraints::EFD);

        RealT* g_lower = g_lower_.getData();
        RealT* g_upper = g_upper_.getData();

        g_lower[EFD] = Efdmin_;
        g_upper[EFD] = Efdmax_;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      [[gnu::always_inline]] void
      SexsPti<scalar_type, index_type>::internalConstraints([[maybe_unused]] const ScalarT* x,
                                                            const ScalarT*                  external,
                                                            ScalarT*                        g) const
      {
        g[static_cast<size_t>(SexsPtiInternalConstraints::EFD)] = external[static_cast<size_t>(SexsPtiExternalVariables::EFD)];
      }

      template <typename scalar_type, typename index_type>
      void SexsPti<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Params = SexsPtiParameters;

        const Model::ParameterReader reader(data, "SexsPti");

        Efdmin_ = reader.readReal(Params::Efdmin).value_or(Efdmin_);
        Efdmax_ = reader.readReal(Params::Efdmax).value_or(Efdmax_);
      }
    } // namespace Exciter
  } // namespace OptimalPowerFlow
} // namespace GridKit
