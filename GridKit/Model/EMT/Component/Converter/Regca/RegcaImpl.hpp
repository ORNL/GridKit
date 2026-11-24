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
#include <GridKit/Model/RegcaLimits.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct a converter, its rate-limit selections and its peak current.
     *
     * Layout: internal variables `vmeas`, `ip` and `iq` of size 1 and `i` of
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
        Qmin_(data.Qmin),
        lvpl_(data.sL && data.IL1 >= ZERO<RealT>)
    {
      if (V_ > ZERO<RealT>)
      {
        current_peak_ = std::sqrt(TWO<RealT> / THREE<RealT>) * MEGA<RealT> * mva_ / V_;
      }
    }

    template <typename scalar_type, typename index_type>
    int Regca<scalar_type, index_type>::initialize()
    {
      if (this->time() == ZERO<RealT>)
      {
        const auto* y   = this->y_.getData();
        const RealT iq0 = std::max(Qmin_, y[IQ] - Khv_ * std::max(ZERO<RealT>, y[VM] - Vhvmax_));
        rq_max_on_      = iq0 > ZERO<RealT> && Rqmax_ > ZERO<RealT>;
        rq_min_on_      = iq0 < ZERO<RealT> && Rqmin_ < ZERO<RealT>;
      }
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Regca<scalar_type, index_type>::evaluateRootsAt(RealT, const scalar_type* y, RealT* roots)
    {
      roots[0] = y[VM] - VL1_;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    bool Regca<scalar_type, index_type>::applyRootEvents(RealT, std::span<const int> crossed)
    {
      auto& ip = this->y_.getData()[IP];
      if (crossed[0] < 0 && ip > IL1_)
      {
        ip = IL1_;
        return true;
      }
      return false;
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
      checks.check(Tg_ > ZERO<RealT> && TM_ >= ZERO<RealT>, "Tg must be positive and TM nonnegative");
      checks.check(!lvpl_ || (VL0_ >= ZERO<RealT> && VL1_ > VL0_), "VL0 must be nonnegative and below VL1");
      checks.check(VA0_ >= ZERO<RealT> && VA1_ > VA0_, "VA0 must be nonnegative and below VA1");
      checks.check(Rpmax_ >= ZERO<RealT>, "Rpmax must be nonnegative");
      checks.check(Vhvmax_ > ZERO<RealT>, "Vhvmax must be positive");
      checks.check(Khv_ >= ZERO<RealT> && Qmin_ <= ZERO<RealT>, "Khv must be nonnegative and Qmin nonpositive");
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

      this->absoluteTolerance(RegcaInternalVariables::vmeas)[0] = rel_tol;
      this->absoluteTolerance(RegcaInternalVariables::ip)[0]    = rel_tol;
      this->absoluteTolerance(RegcaInternalVariables::iq)[0]    = rel_tol;
      std::ranges::fill(this->absoluteTolerance(RegcaInternalVariables::i), rel_tol * this->va_system_base_ / v_base);
      return 0;
    }

    /**
     * @brief Converter current states and their algebraic network limits.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    Regca<scalar_type, index_type>::internalResidual(const T* y, const T* yp, const T* w, const T*, T* f) const
    {
      const T vt = w[VT] / V_;

      const T fp   = (w[IP_CMD] - y[IP]) / Tg_;
      const T fq   = (w[IQ_CMD] - y[IQ]) / Tg_;
      const T rate = Model::rrpwr(y[IP], fp, Rpmax_);
      const T ip   = y[IP] * Math::linseg(vt, VA0_, VA1_, ONE<RealT>);
      const T iq   = std::max(T{Qmin_}, y[IQ] - Khv_ * std::fdim(vt, Vhvmax_));

      f[VM]                    = vt - y[VM] - TM_ * yp[VM];
      const RealT slope        = lvpl_ ? IL1_ / (VL1_ - VL0_) : ZERO<RealT>;
      const T     ceiling      = slope * Math::ramp(y[VM] - VL0_);
      const T     ceiling_rate = slope * Math::sigmoid(y[VM] - VL0_) * yp[VM];
      const T     limited      = Model::awmax(y[IP], rate, ceiling, ceiling_rate);
      f[IP]                    = -yp[IP] + rate + lvpl_ * (y[VM] <= VL1_) * (limited - rate);
      f[IQ]                    = fq - yp[IQ] + rq_max_on_ * (std::min(fq, T{Rqmax_}) - fq)
              + rq_min_on_ * (std::max(fq, T{Rqmin_}) - fq);

#pragma clang loop unroll(full)
      for (size_t n = 0; n < PHASES; ++n)
      {
        const T angle = w[THETA] + GAMMA[n];

        f[I + n] = -y[I + n] + current_peak_ * (ip * std::cos(angle) + iq * std::sin(angle));
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
