/**
 * @file BridgeImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT converter bridge.
 *
 * The model equations are in the `README.md` of this directory, and the
 * kernel arguments are defined in `ComponentImpl.hpp`.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <numbers>

#include <GridKit/Model/EMT/Component/Converter/Bridge/Bridge.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>
#include <GridKit/Model/EMT/Operators/Modulation/PWM/PwmImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct a bridge with three legs.
     *
     * Layout: internal variables `i`, `vo` and `ig` of size 3 and `idq` and
     * `vodq` of size 2, terminal `bus` of size 3 and terminals `vpos` and
     * `vneg` of size 1, input `theta` of size 1 and input `mdq` of size 2.
     * `verify` reports invalid data.
     */
    template <typename scalar_type, typename index_type>
    Bridge<scalar_type, index_type>::Bridge(const ModelDataT& data)
      : ComponentT(std::type_identity<Bridge>{}),
        N_(data.N),
        Rs_(data.Rs),
        Ls_(data.Ls),
        C_(data.C),
        Rg_(data.Rg),
        Lg_(data.Lg),
        pwm_(data.pwm)
    {
      for (size_t n = 0; n < LEGS; ++n)
      {
        gamma_[n] = -TWO<RealT> * std::numbers::pi_v<RealT> * static_cast<RealT>(n) / THREE<RealT>;
      }
    }

    /**
     * @brief Check the leg count, the filter parameters and the modulation.
     *
     * @return The failed checks.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Bridge<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Bridge");
      checks.check(N_ == LEGS, "N must be 3");
      checks.check(Rs_ >= ZERO<RealT> && Rg_ >= ZERO<RealT>, "Rs and Rg must be nonnegative");
      checks.check(Ls_ > ZERO<RealT> && Lg_ > ZERO<RealT>, "Ls and Lg must be positive");
      checks.check(C_ > ZERO<RealT>, "C must be positive");
      const Model::ConfigurationChecks pwm_checks = pwm_.verify();
      for (const auto& error : pwm_checks.errors())
      {
        checks.fail(error);
      }
      return checks;
    }

    /**
     * @brief Set the voltage tolerances to the bus tolerance, and the current
     *        tolerances to `rel_tol` times the system power base over the
     *        bus base voltage.
     *
     * @param[in] rel_tol Relative tolerance of the integrator.
     * @return 0.
     *
     * @pre The bus at terminal `bus` has set its tolerances with the same `rel_tol`.
     */
    template <typename scalar_type, typename index_type>
    int Bridge<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      const RealT v_tol  = this->externalTolerance(this->terminalOffset(BridgeTerminals::bus));
      const RealT v_base = v_tol / rel_tol;
      const RealT i_tol  = rel_tol * this->va_system_base_ / v_base;

      std::ranges::fill(this->absoluteTolerance(BridgeInternalVariables::i), i_tol);
      std::ranges::fill(this->absoluteTolerance(BridgeInternalVariables::vo), v_tol);
      std::ranges::fill(this->absoluteTolerance(BridgeInternalVariables::ig), i_tol);
      std::ranges::fill(this->absoluteTolerance(BridgeInternalVariables::idq), i_tol);
      std::ranges::fill(this->absoluteTolerance(BridgeInternalVariables::vodq), v_tol);
      return 0;
    }

    /**
     * @brief Call `evaluate` with the kernels of the series that the
     *        modulation selects.
     *
     * The series is chosen here, outside the kernels, since Enzyme does not
     * differentiate a choice among them.
     *
     * @tparam EvaluateT Callable taking the kernels and returning a status.
     * @return The status of `evaluate`.
     */
    template <typename scalar_type, typename index_type>
    template <typename EvaluateT>
    int Bridge<scalar_type, index_type>::evaluate(EvaluateT evaluate) const
    {
      switch (pwm_.series())
      {
      case PwmSeries::pulses:
        return evaluate(Kernels<PwmSeries::pulses>{*this});
      case PwmSeries::harmonics:
        return evaluate(Kernels<PwmSeries::harmonics>{*this});
      case PwmSeries::duty:
        break;
      }
      return evaluate(Kernels<PwmSeries::duty>{*this});
    }

    /**
     * @brief Leg angles and switching functions at the current time.
     *
     * The modulation signals are rotated to the legs at angle `theta`.
     *
     * @tparam series Series that evaluates the switching functions.
     * @tparam T Scalar type of the kernel.
     * @param[in] w External values of the bridge.
     */
    template <typename scalar_type, typename index_type>
    template <PwmSeries series, typename T>
    [[gnu::always_inline]] auto Bridge<scalar_type, index_type>::legs(const T* w) const -> Legs<T>
    {
      Legs<T>             legs;
      std::array<T, LEGS> m;
#pragma clang loop unroll(full)
      for (size_t n = 0; n < LEGS; ++n)
      {
        legs.cosine[n] = std::cos(w[THETA] + gamma_[n]);
        legs.sine[n]   = std::sin(w[THETA] + gamma_[n]);
        m[n]           = w[MD] * legs.cosine[n] - w[MQ] * legs.sine[n];
      }
      pwm_.template switching<series>(m.data(), this->time(), legs.s.data());
      return legs;
    }

    /**
     * @brief Write the filter equations of each leg and the dq projections.
     *
     * The indices are constants, since Enzyme needs the layout at compile time.
     */
    template <typename scalar_type, typename index_type>
    template <PwmSeries series, typename T>
    [[gnu::always_inline]] void
    Bridge<scalar_type, index_type>::internalResidual(const T* y, const T* yp, const T* w, const T*, T* f) const
    {
      constexpr RealT TWO_THIRDS = TWO<RealT> / THREE<RealT>;

      const auto [cosine, sine, s] = legs<series>(w);

#pragma clang loop unroll(full)
      for (size_t n = 0; n < LEGS; ++n)
      {
        const T e = (ONE<RealT> - s[n]) * w[VNEG] + s[n] * w[VPOS];

        f[I + n]  = Ls_ * yp[I + n] + Rs_ * y[I + n] + y[VO + n] - e;
        f[VO + n] = C_ * yp[VO + n] + y[IG + n] - y[I + n];
        f[IG + n] = Lg_ * yp[IG + n] + Rg_ * y[IG + n] + w[BUS + n] - y[VO + n];
      }

      f[ID]  = -y[ID] + TWO_THIRDS * (y[I] * cosine[0] + y[I + 1] * cosine[1] + y[I + 2] * cosine[2]);
      f[IQ]  = -y[IQ] - TWO_THIRDS * (y[I] * sine[0] + y[I + 1] * sine[1] + y[I + 2] * sine[2]);
      f[VOD] = -y[VOD] + TWO_THIRDS * (y[VO] * cosine[0] + y[VO + 1] * cosine[1] + y[VO + 2] * cosine[2]);
      f[VOQ] = -y[VOQ] - TWO_THIRDS * (y[VO] * sine[0] + y[VO + 1] * sine[1] + y[VO + 2] * sine[2]);
    }

    /**
     * @brief Write the grid-side currents to the residual of the bus voltage,
     *        and the DC currents to the residuals of the DC terminal voltages.
     */
    template <typename scalar_type, typename index_type>
    template <PwmSeries series, typename T>
    [[gnu::always_inline]] void
    Bridge<scalar_type, index_type>::externalResidual(const T* y, const T*, const T* w, const T*, T* f) const
    {
      const auto s = legs<series>(w).s;

#pragma clang loop unroll(full)
      for (size_t n = 0; n < LEGS; ++n)
      {
        f[BUS + n] = y[IG + n];
      }
      f[VPOS] = -(s[0] * y[I] + s[1] * y[I + 1] + s[2] * y[I + 2]);
      f[VNEG] = -((ONE<RealT> - s[0]) * y[I] + (ONE<RealT> - s[1]) * y[I + 1] + (ONE<RealT> - s[2]) * y[I + 2]);
    }
  } // namespace EMT
} // namespace GridKit
