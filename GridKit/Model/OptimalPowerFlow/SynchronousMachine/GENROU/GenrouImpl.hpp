/**
 * @file GenrouImpl.hpp
 * @brief Definition of an optimal power flow round-rotor machine.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentImpl.hpp>
#include <GridKit/Model/OptimalPowerFlow/SynchronousMachine/GENROU/Genrou.hpp>
#include <GridKit/Model/ParameterReader.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    Genrou<scalar_type, index_type>::Genrou(const ModelDataT& data)
      : ComponentT(std::type_identity<Genrou>{}, data.id)
    {
      initializeParameters(data);
    }

    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Genrou<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Genrou");
      checks.check(Xd_ > Xl_, "Xd must exceed Xl");
      return checks;
    }

    template <typename scalar_type, typename index_type>
    int Genrou<scalar_type, index_type>::tagInequality()
    {
      const auto MECHANICAL = static_cast<size_t>(GenrouInternalConstraints::MECHANICAL);
      const auto FIELD      = static_cast<size_t>(GenrouInternalConstraints::FIELD);

      tag_[MECHANICAL] = false;
      tag_[FIELD]      = false;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Genrou<scalar_type, index_type>::initialize()
    {
      const auto IR = static_cast<size_t>(GenrouInternalVariables::IR);
      const auto II = static_cast<size_t>(GenrouInternalVariables::II);

      RealT* x_lower = x_lower_.getData();
      RealT* x_upper = x_upper_.getData();

      x_lower[IR] = -ONE<RealT>;
      x_upper[IR] = ONE<RealT>;
      x_lower[II] = -ONE<RealT>;
      x_upper[II] = ONE<RealT>;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    [[gnu::always_inline]] void
    Genrou<scalar_type, index_type>::internalConstraints(const ScalarT* x,
                                                         const ScalarT* external,
                                                         ScalarT*       g) const
    {
      const auto IR         = static_cast<size_t>(GenrouInternalVariables::IR);
      const auto II         = static_cast<size_t>(GenrouInternalVariables::II);
      const auto VR         = static_cast<size_t>(GenrouExternalVariables::VR);
      const auto VI         = static_cast<size_t>(GenrouExternalVariables::VI);
      const auto PM         = static_cast<size_t>(GenrouExternalVariables::PM);
      const auto EFD        = static_cast<size_t>(GenrouExternalVariables::EFD);
      const auto MECHANICAL = static_cast<size_t>(GenrouInternalConstraints::MECHANICAL);
      const auto FIELD      = static_cast<size_t>(GenrouInternalConstraints::FIELD);

      const ScalarT ir = x[IR];
      const ScalarT ii = x[II];
      const ScalarT vr = external[VR];
      const ScalarT vi = external[VI];

      const ScalarT psir   = vr + Ra_ * ir - Xdpp_ * ii;
      const ScalarT psii   = vi + Ra_ * ii + Xdpp_ * ir;
      const ScalarT ksat   = SB_ * Math::qramp(std::sqrt(psir * psir + psii * psii) - SA_);
      const ScalarT kq     = ONE<RealT> + Xqd_ * ksat;
      const ScalarT er     = kq * psir - (Xq_ - Xdpp_) * ii;
      const ScalarT ei     = kq * psii + (Xq_ - Xdpp_) * ir;
      const ScalarT e      = std::sqrt(er * er + ei * ei);
      const ScalarT id     = (ir * ei - ii * er) / e;
      const ScalarT psidpp = (psir * er + psii * ei) / e;

      g[MECHANICAL] = -external[PM] + this->toSystemBase(vr * ir + vi * ii + Ra_ * (ir * ir + ii * ii));
      g[FIELD]      = -external[EFD] + (ONE<RealT> + ksat) * psidpp + (Xd_ - Xdpp_) * id;
    }

    template <typename scalar_type, typename index_type>
    [[gnu::always_inline]] void
    Genrou<scalar_type, index_type>::externalConstraints(const ScalarT* x,
                                                         const ScalarT* external,
                                                         ScalarT*       g) const
    {
      const auto IR = static_cast<size_t>(GenrouInternalVariables::IR);
      const auto II = static_cast<size_t>(GenrouInternalVariables::II);
      const auto VR = static_cast<size_t>(GenrouExternalVariables::VR);
      const auto VI = static_cast<size_t>(GenrouExternalVariables::VI);
      const auto DP = static_cast<size_t>(GenrouExternalConstraints::DP);
      const auto DQ = static_cast<size_t>(GenrouExternalConstraints::DQ);

      g[DP] = this->toSystemBase(external[VR] * x[IR] + external[VI] * x[II]);
      g[DQ] = this->toSystemBase(external[VI] * x[IR] - external[VR] * x[II]);
    }

    template <typename scalar_type, typename index_type>
    void Genrou<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Params = GenrouParameters;

      const Model::ParameterReader reader(data, "Genrou");

      Ra_   = reader.readReal(Params::Ra).value_or(Ra_);
      Xd_   = reader.readReal(Params::Xd).value_or(Xd_);
      Xq_   = reader.readReal(Params::Xq).value_or(Xq_);
      Xdpp_ = reader.readReal(Params::Xdpp).value_or(Xdpp_);
      Xl_   = reader.readReal(Params::Xl).value_or(Xl_);
      S10_  = reader.readReal(Params::S10).value_or(S10_);
      S12_  = reader.readReal(Params::S12).value_or(S12_);
      mva_  = reader.readReal(Params::mva).value_or(mva_);
      setDerivedParameters();
    }

    /// Machine base and the saturation coefficients of the PhasorDynamics model
    template <typename scalar_type, typename index_type>
    void Genrou<scalar_type, index_type>::setDerivedParameters()
    {
      this->setComponentBase(MEGA<RealT> * mva_);

      Xqd_ = (Xq_ - Xl_) / (Xd_ - Xl_);

      SA_ = ZERO<RealT>;
      SB_ = ZERO<RealT>;
      if (S12_ != ZERO<RealT>)
      {
        const RealT r = std::sqrt(S10_ / S12_);

        SA_ = std::min((1.2 * r + ONE<RealT>) / (r + ONE<RealT>), (1.2 * r - ONE<RealT>) / (r - ONE<RealT>) );
        SB_ = S12_ / ((SA_ - 1.2) * (SA_ - 1.2));
      }
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
