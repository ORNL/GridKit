/**
 * @file Ieeet1Impl.hpp
 * @brief Definition of an optimal power flow IEEET1 exciter.
 */

#pragma once

#include <cmath>
#include <type_traits>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentImpl.hpp>
#include <GridKit/Model/OptimalPowerFlow/Exciter/IEEET1/Ieeet1.hpp>
#include <GridKit/Model/ParameterReader.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Exciter
    {
      template <typename scalar_type, typename index_type>
      Ieeet1<scalar_type, index_type>::Ieeet1(const ModelDataT& data)
        : ComponentT(std::type_identity<Ieeet1>{}, data.id)
      {
        initializeParameters(data);
      }

      template <typename scalar_type, typename index_type>
      Model::ConfigurationChecks Ieeet1<scalar_type, index_type>::verify() const
      {
        const bool unsaturated = Se1_ == ZERO<RealT> && Se2_ == ZERO<RealT>;
        const bool saturated   = E1_ > ZERO<RealT> && E2_ > ZERO<RealT>
                               && Se1_ >= ZERO<RealT> && Se2_ >= ZERO<RealT>
                               && (E2_ - E1_) * (Se2_ - Se1_) > ZERO<RealT>;

        Model::ConfigurationChecks checks("Ieeet1");
        checks.check(Vrmin_ <= Vrmax_, "Vrmin must not exceed Vrmax");
        checks.check(unsaturated || saturated, "E1, E2, Se1, and Se2 must give a saturation curve");
        return checks;
      }

      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::tagInequality()
      {
        const auto VR = static_cast<size_t>(Ieeet1InternalConstraints::VR);

        tag_[VR] = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::initialize()
      {
        const auto VR = static_cast<size_t>(Ieeet1InternalConstraints::VR);

        RealT* g_lower = g_lower_.getData();
        RealT* g_upper = g_upper_.getData();

        g_lower[VR] = Vrmin_;
        g_upper[VR] = Vrmax_;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      [[gnu::always_inline]] void
      Ieeet1<scalar_type, index_type>::internalConstraints([[maybe_unused]] const ScalarT* x,
                                                           const ScalarT*                  external,
                                                           ScalarT*                        g) const
      {
        const auto EFD = static_cast<size_t>(Ieeet1ExternalVariables::EFD);
        const auto VR  = static_cast<size_t>(Ieeet1InternalConstraints::VR);

        const ScalarT efd  = external[EFD];
        const ScalarT ksat = SB_ * Math::qramp(efd - SA_);

        g[VR] = ke_on_ * (Ke_ * efd + ksat) + (ONE<RealT> - ke_on_) * Vrmax_ / 10.0;
      }

      template <typename scalar_type, typename index_type>
      void Ieeet1<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Params = Ieeet1Parameters;

        const Model::ParameterReader reader(data, "Ieeet1");

        Ke_    = reader.readReal(Params::Ke).value_or(Ke_);
        Vrmin_ = reader.readReal(Params::Vrmin).value_or(Vrmin_);
        Vrmax_ = reader.readReal(Params::Vrmax).value_or(Vrmax_);
        E1_    = reader.readReal(Params::E1).value_or(E1_);
        E2_    = reader.readReal(Params::E2).value_or(E2_);
        Se1_   = reader.readReal(Params::Se1).value_or(Se1_);
        Se2_   = reader.readReal(Params::Se2).value_or(Se2_);
        setDerivedParameters();
      }

      /// Masks and the saturation coefficients of the PhasorDynamics model
      template <typename scalar_type, typename index_type>
      void Ieeet1<scalar_type, index_type>::setDerivedParameters()
      {
        ke_on_ = ZERO<RealT>;
        if (Ke_ != ZERO<RealT>)
        {
          ke_on_ = ONE<RealT>;
        }

        SA_ = ZERO<RealT>;
        SB_ = ZERO<RealT>;
        if (Se1_ == ZERO<RealT> && Se2_ == ZERO<RealT>)
        {
          return;
        }

        if (Se1_ == ZERO<RealT>)
        {
          SA_ = E1_;
          SB_ = Se2_ * E2_ / ((E2_ - E1_) * (E2_ - E1_));
          return;
        }

        if (Se2_ == ZERO<RealT>)
        {
          SA_ = E2_;
          SB_ = Se1_ * E1_ / ((E1_ - E2_) * (E1_ - E2_));
          return;
        }

        const RealT C = std::sqrt(Se2_ * E2_ / (Se1_ * E1_));

        SA_ = (C * E1_ - E2_) / (C - ONE<RealT>);
        SB_ = Se1_ * E1_ / ((E1_ - SA_) * (E1_ - SA_));
      }
    } // namespace Exciter
  } // namespace OptimalPowerFlow
} // namespace GridKit
