/**
 * @file RegcaImpl.hpp
 * @brief Definition of an optimal power flow REGCA converter.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>

#include <GridKit/Model/OptimalPowerFlow/ComponentImpl.hpp>
#include <GridKit/Model/OptimalPowerFlow/Converter/REGCA/Regca.hpp>
#include <GridKit/Model/ParameterReader.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace Converter
    {
      template <typename scalar_type, typename index_type>
      Regca<scalar_type, index_type>::Regca(const ModelDataT& data)
        : ComponentT(std::type_identity<Regca>{}, data.id)
      {
        initializeParameters(data);
      }

      template <typename scalar_type, typename index_type>
      Model::ConfigurationChecks Regca<scalar_type, index_type>::verify() const
      {
        Model::ConfigurationChecks checks("Regca");
        checks.check(mva_ > ZERO<RealT>, "mva must be positive");
        checks.check(IL1_ >= ZERO<RealT>, "IL1 must be nonnegative");
        checks.check(VA1_ < Vhvmax_, "VA1 must be less than Vhvmax");
        return checks;
      }

      template <typename scalar_type, typename index_type>
      int Regca<scalar_type, index_type>::tagInequality()
      {
        const auto ACTIVE_COMMAND   = static_cast<size_t>(RegcaInternalConstraints::ACTIVE_COMMAND);
        const auto REACTIVE_COMMAND = static_cast<size_t>(RegcaInternalConstraints::REACTIVE_COMMAND);
        const auto ACTIVE_POWER     = static_cast<size_t>(RegcaInternalConstraints::ACTIVE_POWER);
        const auto REACTIVE_POWER   = static_cast<size_t>(RegcaInternalConstraints::REACTIVE_POWER);
        const auto VT2              = static_cast<size_t>(RegcaInternalConstraints::VT2);

        tag_[ACTIVE_COMMAND]   = false;
        tag_[REACTIVE_COMMAND] = false;
        tag_[ACTIVE_POWER]     = false;
        tag_[REACTIVE_POWER]   = false;
        tag_[VT2]              = true;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Regca<scalar_type, index_type>::initialize()
      {
        const auto IP  = static_cast<size_t>(RegcaInternalVariables::IP);
        const auto VT2 = static_cast<size_t>(RegcaInternalConstraints::VT2);

        RealT* x_upper = x_upper_.getData();
        RealT* g_lower = g_lower_.getData();
        RealT* g_upper = g_upper_.getData();

        RealT Vmin = VA1_;
        if (sL_)
        {
          x_upper[IP] = IL1_;
          Vmin        = std::max(VA1_, VL1_);
        }

        g_lower[VT2] = Vmin * Vmin;
        g_upper[VT2] = Vhvmax_ * Vhvmax_;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      [[gnu::always_inline]] void
      Regca<scalar_type, index_type>::internalConstraints(const ScalarT* x,
                                                          const ScalarT* external,
                                                          ScalarT*       g) const
      {
        const auto IP               = static_cast<size_t>(RegcaInternalVariables::IP);
        const auto IQ               = static_cast<size_t>(RegcaInternalVariables::IQ);
        const auto VR               = static_cast<size_t>(RegcaExternalVariables::VR);
        const auto VI               = static_cast<size_t>(RegcaExternalVariables::VI);
        const auto IPCMD            = static_cast<size_t>(RegcaExternalVariables::IPCMD);
        const auto IQCMD            = static_cast<size_t>(RegcaExternalVariables::IQCMD);
        const auto PBR              = static_cast<size_t>(RegcaExternalVariables::PBR);
        const auto QBR              = static_cast<size_t>(RegcaExternalVariables::QBR);
        const auto ACTIVE_COMMAND   = static_cast<size_t>(RegcaInternalConstraints::ACTIVE_COMMAND);
        const auto REACTIVE_COMMAND = static_cast<size_t>(RegcaInternalConstraints::REACTIVE_COMMAND);
        const auto ACTIVE_POWER     = static_cast<size_t>(RegcaInternalConstraints::ACTIVE_POWER);
        const auto REACTIVE_POWER   = static_cast<size_t>(RegcaInternalConstraints::REACTIVE_POWER);
        const auto VT2              = static_cast<size_t>(RegcaInternalConstraints::VT2);

        const ScalarT ip  = x[IP];
        const ScalarT iq  = x[IQ];
        const ScalarT vt2 = external[VR] * external[VR] + external[VI] * external[VI];
        const ScalarT vt  = std::sqrt(vt2);

        g[ACTIVE_COMMAND]   = -external[IPCMD] + this->toSystemBase(ip);
        g[REACTIVE_COMMAND] = -external[IQCMD] + this->toSystemBase(iq);
        g[ACTIVE_POWER]     = -external[PBR] + this->toSystemBase(vt * ip);
        g[REACTIVE_POWER]   = -external[QBR] + this->toSystemBase(vt * iq);
        g[VT2]              = vt2;
      }

      template <typename scalar_type, typename index_type>
      [[gnu::always_inline]] void
      Regca<scalar_type, index_type>::externalConstraints(const ScalarT* x,
                                                          const ScalarT* external,
                                                          ScalarT*       g) const
      {
        const auto IP = static_cast<size_t>(RegcaInternalVariables::IP);
        const auto IQ = static_cast<size_t>(RegcaInternalVariables::IQ);
        const auto VR = static_cast<size_t>(RegcaExternalVariables::VR);
        const auto VI = static_cast<size_t>(RegcaExternalVariables::VI);
        const auto DP = static_cast<size_t>(RegcaExternalConstraints::DP);
        const auto DQ = static_cast<size_t>(RegcaExternalConstraints::DQ);

        const ScalarT vt = std::sqrt(external[VR] * external[VR] + external[VI] * external[VI]);

        g[DP] = this->toSystemBase(vt * x[IP]);
        g[DQ] = this->toSystemBase(vt * x[IQ]);
      }

      template <typename scalar_type, typename index_type>
      void Regca<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Params = RegcaParameters;

        const Model::ParameterReader reader(data, "Regca");

        mva_    = reader.readReal(Params::mva).value_or(mva_);
        sL_     = reader.readBool(Params::sL).value_or(sL_);
        IL1_    = reader.readReal(Params::IL1).value_or(IL1_);
        VL1_    = reader.readReal(Params::VL1).value_or(VL1_);
        VA1_    = reader.readReal(Params::VA1).value_or(VA1_);
        Vhvmax_ = reader.readReal(Params::Vhvmax).value_or(Vhvmax_);
        setDerivedParameters();
      }

      template <typename scalar_type, typename index_type>
      void Regca<scalar_type, index_type>::setDerivedParameters()
      {
        this->setComponentBase(MEGA<RealT> * mva_);
      }
    } // namespace Converter
  } // namespace OptimalPowerFlow
} // namespace GridKit
