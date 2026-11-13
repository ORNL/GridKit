/**
 * @file GenClassicalImpl.hpp
 * @brief Definition of an optimal power flow classical machine.
 */

#pragma once

#include <cmath>
#include <type_traits>

#include <GridKit/Model/OptimalPowerFlow/ComponentImpl.hpp>
#include <GridKit/Model/OptimalPowerFlow/SynchronousMachine/GenClassical/GenClassical.hpp>
#include <GridKit/Model/ParameterReader.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    GenClassical<scalar_type, index_type>::GenClassical(const ModelDataT& data)
      : ComponentT(std::type_identity<GenClassical>{}, data.id)
    {
      initializeParameters(data);
    }

    /// Every parameter set is valid
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks GenClassical<scalar_type, index_type>::verify() const
    {
      return Model::ConfigurationChecks("GenClassical");
    }

    template <typename scalar_type, typename index_type>
    int GenClassical<scalar_type, index_type>::tagInequality()
    {
      const auto MECHANICAL = static_cast<size_t>(GenClassicalInternalConstraints::MECHANICAL);
      const auto FIELD      = static_cast<size_t>(GenClassicalInternalConstraints::FIELD);

      tag_[MECHANICAL] = false;
      tag_[FIELD]      = false;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int GenClassical<scalar_type, index_type>::initialize()
    {
      const auto IR = static_cast<size_t>(GenClassicalInternalVariables::IR);
      const auto II = static_cast<size_t>(GenClassicalInternalVariables::II);

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
    GenClassical<scalar_type, index_type>::internalConstraints(const ScalarT* x,
                                                               const ScalarT* external,
                                                               ScalarT*       g) const
    {
      const auto IR         = static_cast<size_t>(GenClassicalInternalVariables::IR);
      const auto II         = static_cast<size_t>(GenClassicalInternalVariables::II);
      const auto VR         = static_cast<size_t>(GenClassicalExternalVariables::VR);
      const auto VI         = static_cast<size_t>(GenClassicalExternalVariables::VI);
      const auto PM         = static_cast<size_t>(GenClassicalExternalVariables::PM);
      const auto EFD        = static_cast<size_t>(GenClassicalExternalVariables::EFD);
      const auto MECHANICAL = static_cast<size_t>(GenClassicalInternalConstraints::MECHANICAL);
      const auto FIELD      = static_cast<size_t>(GenClassicalInternalConstraints::FIELD);

      const ScalarT ir = x[IR];
      const ScalarT ii = x[II];
      const ScalarT vr = external[VR];
      const ScalarT vi = external[VI];

      const ScalarT er = vr + Ra_ * ir - Xdp_ * ii;
      const ScalarT ei = vi + Ra_ * ii + Xdp_ * ir;

      g[MECHANICAL] = -external[PM] + this->toSystemBase(vr * ir + vi * ii + Ra_ * (ir * ir + ii * ii));
      g[FIELD]      = -external[EFD] + std::sqrt(er * er + ei * ei);
    }

    template <typename scalar_type, typename index_type>
    [[gnu::always_inline]] void
    GenClassical<scalar_type, index_type>::externalConstraints(const ScalarT* x,
                                                               const ScalarT* external,
                                                               ScalarT*       g) const
    {
      const auto IR = static_cast<size_t>(GenClassicalInternalVariables::IR);
      const auto II = static_cast<size_t>(GenClassicalInternalVariables::II);
      const auto VR = static_cast<size_t>(GenClassicalExternalVariables::VR);
      const auto VI = static_cast<size_t>(GenClassicalExternalVariables::VI);
      const auto DP = static_cast<size_t>(GenClassicalExternalConstraints::DP);
      const auto DQ = static_cast<size_t>(GenClassicalExternalConstraints::DQ);

      g[DP] = this->toSystemBase(external[VR] * x[IR] + external[VI] * x[II]);
      g[DQ] = this->toSystemBase(external[VI] * x[IR] - external[VR] * x[II]);
    }

    template <typename scalar_type, typename index_type>
    void GenClassical<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Params = GenClassicalParameters;

      const Model::ParameterReader reader(data, "GenClassical");

      Ra_  = reader.readReal(Params::Ra).value_or(Ra_);
      Xdp_ = reader.readReal(Params::Xdp).value_or(Xdp_);
      mva_ = reader.readReal(Params::mva).value_or(mva_);
      setDerivedParameters();
    }

    /// Machine base
    template <typename scalar_type, typename index_type>
    void GenClassical<scalar_type, index_type>::setDerivedParameters()
    {
      this->setComponentBase(MEGA<RealT> * mva_);
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
