/**
 * @file RegcaImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT renewable energy generator/converter model.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/Converter/Regca/Regca.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct a converter, its rate-limit selections and its peak current.
     *
     * Layout: internal variables `il`, `ip` and `iq` of size 1 and `i` of
     * size 3, terminal `bus` of size 3, and inputs `vt`, `theta`, `ipcmd` and
     * `iqcmd` of size 1. `verify` reports invalid data.
     */
    template <typename scalar_type, typename index_type>
    Regca<scalar_type, index_type>::Regca(const ModelDataT& data)
      : ComponentT(std::type_identity<Regca>{}),
        V_(data.V),
        mva_(data.mva),
        Tg_(data.Tg),
        TM_(data.TM),
        Rqmax_(data.Rqmax),
        Rqmin_(data.Rqmin),
        Rpmax_(data.Rpmax),
        IL1_(data.IL1),
        VL0_(data.VL0),
        VL1_(data.VL1),
        VA0_(data.VA0),
        VA1_(data.VA1),
        Vhvmax_(data.Vhvmax),
        Khv_(data.Khv),
        Qmin_(data.Qmin)
    {
      if (Rqmax_ > ZERO<RealT>)
      {
        rq_max_on_ = ONE<RealT>;
      }
      if (Rqmin_ < ZERO<RealT>)
      {
        rq_min_on_ = ONE<RealT>;
      }
      if (V_ > ZERO<RealT>)
      {
        current_peak_ = std::sqrt(TWO<RealT> / THREE<RealT>) * MEGA<RealT> * mva_ / V_;
      }
    }

    /**
     * @brief Check the ratings, the time constants, the voltage breakpoints,
     *        the ceiling target, the rate limits and the voltage response.
     *
     * @return The failed checks.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Regca<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Regca");
      checks.check(V_ > ZERO<RealT> && mva_ > ZERO<RealT> && std::isfinite(current_peak_), "V and mva must be positive");
      checks.check(Tg_ > ZERO<RealT> && TM_ > ZERO<RealT>, "Tg and TM must be positive");
      checks.check(Rqmax_ >= ZERO<RealT> && Rqmin_ <= ZERO<RealT>, "Rqmin <= 0 <= Rqmax is required");
      checks.check(VL0_ >= ZERO<RealT> && VL1_ > VL0_, "VL0 must be nonnegative and below VL1");
      checks.check(VA0_ >= ZERO<RealT> && VA1_ > VA0_, "VA0 must be nonnegative and below VA1");
      checks.check(IL1_ > ZERO<RealT> && Rpmax_ > ZERO<RealT>, "IL1 and Rpmax must be positive");
      checks.check(Vhvmax_ > ZERO<RealT>, "Vhvmax must be positive");
      checks.check(Khv_ >= ZERO<RealT>, "Khv must be nonnegative");
      checks.check(Qmin_ <= ZERO<RealT>, "Qmin must be nonpositive");
      return checks;
    }

    /**
     * @brief Set the tolerance of the currents into the grid from the base
     *        voltage of the bus and the system power base.
     *
     * @param[in] rel_tol Relative tolerance of the integrator.
     * @return 0.
     *
     * @pre The system has bound this model and set the tolerance of its bus.
     */
    template <typename scalar_type, typename index_type>
    int Regca<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      const RealT v_base = this->externalTolerance(this->terminalOffset(RegcaTerminals::bus)) / rel_tol;

      this->absoluteTolerance(RegcaInternalVariables::il)[0] = rel_tol;
      this->absoluteTolerance(RegcaInternalVariables::ip)[0] = rel_tol;
      this->absoluteTolerance(RegcaInternalVariables::iq)[0] = rel_tol;
      std::ranges::fill(this->absoluteTolerance(RegcaInternalVariables::i), rel_tol * this->va_system_base_ / v_base);
      return 0;
    }

    /**
     * @brief Write the residuals of the current ceiling, of the active and
     *        reactive currents and of the currents into the grid.
     *
     * The active-current order is limited at low terminal voltage and the
     * reactive-current order reduced at high terminal voltage. Each current
     * follows its order with the time constant `Tg`, the reactive current
     * within its rate limits.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    Regca<scalar_type, index_type>::internalResidual(const T* y, const T* yp, const T* w, const T*, T* f) const
    {
      const T vt = w[VT] / V_;

      const T target    = Math::linseg(vt, VL0_, VL1_, IL1_);
      const T available = y[IL] * Math::linseg(vt, VA0_, VA1_, ONE<RealT>);
      const T ip_order  = w[IP_CMD] - std::fdim(w[IP_CMD], available);

      const T absorbed = Khv_ * Math::ramp(vt - Vhvmax_);
      const T floor    = std::min(w[IQ_CMD], T{Qmin_});
      const T iq_order = floor + std::fdim(w[IQ_CMD] - floor, absorbed);
      const T iq_rate  = (iq_order - y[IQ]) / Tg_;

      f[IL] = -yp[IL] + std::min((target - y[IL]) / TM_, T{Rpmax_});
      f[IP] = -yp[IP] + (ip_order - y[IP]) / Tg_;
      f[IQ] = -yp[IQ] + iq_rate + rq_max_on_ * (std::min(iq_rate, T{Rqmax_}) - iq_rate)
              + rq_min_on_ * (std::max(iq_rate, T{Rqmin_}) - iq_rate);

#pragma clang loop unroll(full)
      for (size_t n = 0; n < PHASES; ++n)
      {
        const T angle = w[THETA] + GAMMA[n];

        f[I + n] = -y[I + n] + current_peak_ * (y[IP] * std::cos(angle) + y[IQ] * std::sin(angle));
      }
    }

    /**
     * @brief Add the currents into the grid to the current balance of the bus.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    Regca<scalar_type, index_type>::externalResidual(const T* y, const T*, const T*, const T*, T* f) const
    {
#pragma clang loop unroll(full)
      for (size_t n = 0; n < PHASES; ++n)
      {
        f[BUS + n] = y[I + n];
      }
    }
  } // namespace EMT
} // namespace GridKit
