/**
 * @file GensalImpl.hpp
 * @brief Definition of an optimal power flow salient-pole machine.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/OptimalPowerFlow/ComponentImpl.hpp>
#include <GridKit/Model/OptimalPowerFlow/SynchronousMachine/GENSAL/Gensal.hpp>
#include <GridKit/Model/ParameterReader.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    Gensal<scalar_type, index_type>::Gensal(const ModelDataT& data)
      : ComponentT(std::type_identity<Gensal>{}, data.id)
    {
      initializeParameters(data);
    }

    /// Every parameter set is valid
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Gensal<scalar_type, index_type>::verify() const
    {
      return Model::ConfigurationChecks("Gensal");
    }

    template <typename scalar_type, typename index_type>
    int Gensal<scalar_type, index_type>::tagInequality()
    {
      const auto MECHANICAL = static_cast<size_t>(GensalInternalConstraints::MECHANICAL);
      const auto FIELD      = static_cast<size_t>(GensalInternalConstraints::FIELD);

      tag_[MECHANICAL] = false;
      tag_[FIELD]      = false;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Gensal<scalar_type, index_type>::initialize()
    {
      const auto IR = static_cast<size_t>(GensalInternalVariables::IR);
      const auto II = static_cast<size_t>(GensalInternalVariables::II);

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
    Gensal<scalar_type, index_type>::internalConstraints(const ScalarT* x,
                                                         const ScalarT* external,
                                                         ScalarT*       g) const
    {
      const auto IR         = static_cast<size_t>(GensalInternalVariables::IR);
      const auto II         = static_cast<size_t>(GensalInternalVariables::II);
      const auto VR         = static_cast<size_t>(GensalExternalVariables::VR);
      const auto VI         = static_cast<size_t>(GensalExternalVariables::VI);
      const auto PM         = static_cast<size_t>(GensalExternalVariables::PM);
      const auto EFD        = static_cast<size_t>(GensalExternalVariables::EFD);
      const auto MECHANICAL = static_cast<size_t>(GensalInternalConstraints::MECHANICAL);
      const auto FIELD      = static_cast<size_t>(GensalInternalConstraints::FIELD);

      const ScalarT ir = x[IR];
      const ScalarT ii = x[II];
      const ScalarT vr = external[VR];
      const ScalarT vi = external[VI];

      const ScalarT er   = vr + Ra_ * ir - Xq_ * ii;
      const ScalarT ei   = vi + Ra_ * ii + Xq_ * ir;
      const ScalarT e    = std::sqrt(er * er + ei * ei);
      const ScalarT id   = (ir * ei - ii * er) / e;
      const ScalarT eqp  = e - (Xq_ - Xdp_) * id;
      const ScalarT ksat = SB_ * Math::qramp(eqp - SA_);

      g[MECHANICAL] = -external[PM] + this->toSystemBase(vr * ir + vi * ii + Ra_ * (ir * ir + ii * ii));
      g[FIELD]      = -external[EFD] + (ONE<RealT> + ksat) * eqp + (Xd_ - Xdp_) * id;
    }

    template <typename scalar_type, typename index_type>
    [[gnu::always_inline]] void
    Gensal<scalar_type, index_type>::externalConstraints(const ScalarT* x,
                                                         const ScalarT* external,
                                                         ScalarT*       g) const
    {
      const auto IR = static_cast<size_t>(GensalInternalVariables::IR);
      const auto II = static_cast<size_t>(GensalInternalVariables::II);
      const auto VR = static_cast<size_t>(GensalExternalVariables::VR);
      const auto VI = static_cast<size_t>(GensalExternalVariables::VI);
      const auto DP = static_cast<size_t>(GensalExternalConstraints::DP);
      const auto DQ = static_cast<size_t>(GensalExternalConstraints::DQ);

      g[DP] = this->toSystemBase(external[VR] * x[IR] + external[VI] * x[II]);
      g[DQ] = this->toSystemBase(external[VI] * x[IR] - external[VR] * x[II]);
    }

    template <typename scalar_type, typename index_type>
    void Gensal<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Params = GensalParameters;

      const Model::ParameterReader reader(data, "Gensal");

      Ra_  = reader.readReal(Params::Ra).value_or(Ra_);
      Xd_  = reader.readReal(Params::Xd).value_or(Xd_);
      Xdp_ = reader.readReal(Params::Xdp).value_or(Xdp_);
      Xq_  = reader.readReal(Params::Xq).value_or(Xq_);
      S10_ = reader.readReal(Params::S10).value_or(S10_);
      S12_ = reader.readReal(Params::S12).value_or(S12_);
      mva_ = reader.readReal(Params::mva).value_or(mva_);
      setDerivedParameters();
    }

    /// Machine base and the saturation coefficients of the PhasorDynamics model
    template <typename scalar_type, typename index_type>
    void Gensal<scalar_type, index_type>::setDerivedParameters()
    {
      this->setComponentBase(MEGA<RealT> * mva_);

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
