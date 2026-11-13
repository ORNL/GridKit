/**
 * @file Tgov1Impl.hpp
 * @brief Definition of an optimal power flow TGOV1 governor.
 */

#pragma once

#include <type_traits>

#include <GridKit/Model/OptimalPowerFlow/ComponentImpl.hpp>
#include <GridKit/Model/OptimalPowerFlow/Governor/Tgov1/Tgov1.hpp>
#include <GridKit/Model/ParameterReader.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Governor
    {
      template <typename scalar_type, typename index_type>
      Tgov1<scalar_type, index_type>::Tgov1(const ModelDataT& data)
        : ComponentT(std::type_identity<Tgov1>{}, data.id)
      {
        initializeParameters(data);
      }

      template <typename scalar_type, typename index_type>
      Model::ConfigurationChecks Tgov1<scalar_type, index_type>::verify() const
      {
        Model::ConfigurationChecks checks("Tgov1");
        checks.check(Trate_ > ZERO<RealT>, "Trate must be positive");
        checks.check(Pvmin_ <= Pvmax_, "Pvmin must not exceed Pvmax");
        return checks;
      }

      template <typename scalar_type, typename index_type>
      int Tgov1<scalar_type, index_type>::tagInequality()
      {
        const auto PV = static_cast<size_t>(Tgov1InternalConstraints::PV);

        tag_[PV] = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Tgov1<scalar_type, index_type>::initialize()
      {
        const auto PV = static_cast<size_t>(Tgov1InternalConstraints::PV);

        RealT* g_lower = g_lower_.getData();
        RealT* g_upper = g_upper_.getData();

        g_lower[PV] = Pvmin_;
        g_upper[PV] = Pvmax_;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      [[gnu::always_inline]] void
      Tgov1<scalar_type, index_type>::internalConstraints([[maybe_unused]] const ScalarT* x,
                                                          const ScalarT*                  external,
                                                          ScalarT*                        g) const
      {
        const auto PM = static_cast<size_t>(Tgov1ExternalVariables::PM);
        const auto PV = static_cast<size_t>(Tgov1InternalConstraints::PV);

        g[PV] = this->toComponentBase(external[PM]);
      }

      template <typename scalar_type, typename index_type>
      void Tgov1<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Params = Tgov1Parameters;

        const Model::ParameterReader reader(data, "Tgov1");

        Trate_ = reader.readReal(Params::Trate).value_or(Trate_);
        Pvmax_ = reader.readReal(Params::Pvmax).value_or(Pvmax_);
        Pvmin_ = reader.readReal(Params::Pvmin).value_or(Pvmin_);
        setDerivedParameters();
      }

      template <typename scalar_type, typename index_type>
      void Tgov1<scalar_type, index_type>::setDerivedParameters()
      {
        this->setComponentBase(MEGA<RealT> * Trate_);
      }
    } // namespace Governor
  } // namespace OptimalPowerFlow
} // namespace GridKit
