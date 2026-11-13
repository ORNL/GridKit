/**
 * @file ReecbImpl.hpp
 * @brief Definition of an optimal power flow REECB controller.
 */

#pragma once

#include <type_traits>

#include <GridKit/Model/OptimalPowerFlow/ComponentImpl.hpp>
#include <GridKit/Model/OptimalPowerFlow/Controller/REECB/Reecb.hpp>
#include <GridKit/Model/ParameterReader.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Controller
    {
      template <typename scalar_type, typename index_type>
      Reecb<scalar_type, index_type>::Reecb(const ModelDataT& data)
        : ComponentT(std::type_identity<Reecb>{}, data.id)
      {
        initializeParameters(data);
      }

      template <typename scalar_type, typename index_type>
      Model::ConfigurationChecks Reecb<scalar_type, index_type>::verify() const
      {
        Model::ConfigurationChecks checks("Reecb");
        checks.check(mva_ > ZERO<RealT>, "mva must be positive");
        checks.check(Pmin_ <= Pmax_, "Pmin must not exceed Pmax");
        checks.check(Qmin_ <= Qmax_, "Qmin must not exceed Qmax");
        checks.check(Imax_ > ZERO<RealT>, "Imax must be positive");
        return checks;
      }

      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::tagInequality()
      {
        const auto P  = static_cast<size_t>(ReecbInternalConstraints::P);
        const auto Q  = static_cast<size_t>(ReecbInternalConstraints::Q);
        const auto I2 = static_cast<size_t>(ReecbInternalConstraints::I2);

        tag_[P]  = true;
        tag_[Q]  = true;
        tag_[I2] = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Reecb<scalar_type, index_type>::initialize()
      {
        const auto P  = static_cast<size_t>(ReecbInternalConstraints::P);
        const auto Q  = static_cast<size_t>(ReecbInternalConstraints::Q);
        const auto I2 = static_cast<size_t>(ReecbInternalConstraints::I2);

        RealT* g_lower = g_lower_.getData();
        RealT* g_upper = g_upper_.getData();

        g_lower[P]  = Pmin_;
        g_upper[P]  = Pmax_;
        g_lower[Q]  = Qmin_;
        g_upper[Q]  = Qmax_;
        g_upper[I2] = Imax_ * Imax_;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      [[gnu::always_inline]] void
      Reecb<scalar_type, index_type>::internalConstraints([[maybe_unused]] const ScalarT* x,
                                                          const ScalarT*                  external,
                                                          ScalarT*                        g) const
      {
        const auto IPCMD = static_cast<size_t>(ReecbExternalVariables::IPCMD);
        const auto IQCMD = static_cast<size_t>(ReecbExternalVariables::IQCMD);
        const auto PE    = static_cast<size_t>(ReecbExternalVariables::PE);
        const auto QGEN  = static_cast<size_t>(ReecbExternalVariables::QGEN);
        const auto P     = static_cast<size_t>(ReecbInternalConstraints::P);
        const auto Q     = static_cast<size_t>(ReecbInternalConstraints::Q);
        const auto I2    = static_cast<size_t>(ReecbInternalConstraints::I2);

        const ScalarT ip = this->toComponentBase(external[IPCMD]);
        const ScalarT iq = this->toComponentBase(external[IQCMD]);

        g[P]  = this->toComponentBase(external[PE]);
        g[Q]  = this->toComponentBase(external[QGEN]);
        g[I2] = ip * ip + iq * iq;
      }

      template <typename scalar_type, typename index_type>
      void Reecb<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Params = ReecbParameters;

        const Model::ParameterReader reader(data, "Reecb");

        mva_  = reader.readReal(Params::mva).value_or(mva_);
        Pmax_ = reader.readReal(Params::Pmax).value_or(Pmax_);
        Pmin_ = reader.readReal(Params::Pmin).value_or(Pmin_);
        Qmax_ = reader.readReal(Params::Qmax).value_or(Qmax_);
        Qmin_ = reader.readReal(Params::Qmin).value_or(Qmin_);
        Imax_ = reader.readReal(Params::Imax).value_or(Imax_);
        setDerivedParameters();
      }

      template <typename scalar_type, typename index_type>
      void Reecb<scalar_type, index_type>::setDerivedParameters()
      {
        this->setComponentBase(MEGA<RealT> * mva_);
      }
    } // namespace Controller
  } // namespace OptimalPowerFlow
} // namespace GridKit
