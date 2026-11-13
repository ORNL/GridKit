/**
 * @file RepcaImpl.hpp
 * @brief Definition of an optimal power flow REPCA controller.
 */

#pragma once

#include <type_traits>

#include <GridKit/Model/OptimalPowerFlow/ComponentImpl.hpp>
#include <GridKit/Model/OptimalPowerFlow/Controller/REPCA/Repca.hpp>
#include <GridKit/Model/ParameterReader.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Controller
    {
      template <typename scalar_type, typename index_type>
      Repca<scalar_type, index_type>::Repca(const ModelDataT& data)
        : ComponentT(std::type_identity<Repca>{}, data.id)
      {
        initializeParameters(data);
      }

      template <typename scalar_type, typename index_type>
      Model::ConfigurationChecks Repca<scalar_type, index_type>::verify() const
      {
        Model::ConfigurationChecks checks("Repca");
        checks.check(mva_ > ZERO<RealT>, "mva must be positive");
        checks.check(Qmin_ <= Qmax_, "Qmin must not exceed Qmax");
        checks.check(Pmin_ <= Pmax_, "Pmin must not exceed Pmax");
        return checks;
      }

      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::tagInequality()
      {
        const auto Q = static_cast<size_t>(RepcaInternalConstraints::Q);
        const auto P = static_cast<size_t>(RepcaInternalConstraints::P);

        tag_[Q] = true;
        tag_[P] = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Repca<scalar_type, index_type>::initialize()
      {
        const auto Q = static_cast<size_t>(RepcaInternalConstraints::Q);
        const auto P = static_cast<size_t>(RepcaInternalConstraints::P);

        RealT* g_lower = g_lower_.getData();
        RealT* g_upper = g_upper_.getData();

        g_lower[Q] = Qmin_;
        g_upper[Q] = Qmax_;
        g_lower[P] = Pmin_;
        g_upper[P] = Pmax_;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      [[gnu::always_inline]] void
      Repca<scalar_type, index_type>::internalConstraints([[maybe_unused]] const ScalarT* x,
                                                          const ScalarT*                  external,
                                                          ScalarT*                        g) const
      {
        const auto PBR = static_cast<size_t>(RepcaExternalVariables::PBR);
        const auto QBR = static_cast<size_t>(RepcaExternalVariables::QBR);
        const auto Q   = static_cast<size_t>(RepcaInternalConstraints::Q);
        const auto P   = static_cast<size_t>(RepcaInternalConstraints::P);

        g[Q] = this->toComponentBase(external[QBR]);
        g[P] = freq_on_ * this->toComponentBase(external[PBR]);
      }

      template <typename scalar_type, typename index_type>
      void Repca<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Params = RepcaParameters;

        const Model::ParameterReader reader(data, "Repca");

        mva_      = reader.readReal(Params::mva).value_or(mva_);
        Freqflag_ = reader.readBool(Params::Freqflag).value_or(Freqflag_);
        Qmax_     = reader.readReal(Params::Qmax).value_or(Qmax_);
        Qmin_     = reader.readReal(Params::Qmin).value_or(Qmin_);
        Pmax_     = reader.readReal(Params::Pmax).value_or(Pmax_);
        Pmin_     = reader.readReal(Params::Pmin).value_or(Pmin_);
        setDerivedParameters();
      }

      template <typename scalar_type, typename index_type>
      void Repca<scalar_type, index_type>::setDerivedParameters()
      {
        this->setComponentBase(MEGA<RealT> * mva_);

        freq_on_ = ZERO<RealT>;
        if (Freqflag_)
        {
          freq_on_ = ONE<RealT>;
        }
      }
    } // namespace Controller
  } // namespace OptimalPowerFlow
} // namespace GridKit
