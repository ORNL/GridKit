/**
 * @file MachineImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT synchronous machine.
 *
 * The model equations are in the `README.md` of this directory, and the
 * kernel arguments are defined in `ComponentImpl.hpp`.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <numbers>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/Source/Machine/Machine.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>
#include <GridKit/Model/EMT/Operators/Reference/ThreePhase.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct a machine and its per-unit bases.
     *
     * Layout: eleven internal variables of size 1, terminal `bus` of size 3,
     * inputs `efd` and `pm` of size 1, and eight outputs of size 1 before
     * output `is` of size 3. `verify` reports invalid data.
     */
    template <typename scalar_type, typename index_type>
    Machine<scalar_type, index_type>::Machine(const ModelDataT& data)
      : ComponentT(std::type_identity<Machine>{}),
        S_(data.S),
        V_(data.V),
        f_b_(data.f),
        H_(data.H),
        F_(data.F),
        Rs_(data.Rs),
        Ll_(data.Ll),
        Lmd_(data.Lmd),
        Lmq_(data.Lmq),
        L0_(data.L0),
        Rfd_(data.Rfd),
        Llfd_(data.Llfd),
        R1d_(data.R1d),
        Ll1d_(data.Ll1d),
        R1q_(data.R1q),
        Ll1q_(data.Ll1q),
        R2q_(data.R2q),
        Ll2q_(data.Ll2q),
        S10_(data.S10),
        S12_(data.S12),
        Rfe_(data.Rfe)
    {
      if (L0_ == ZERO<RealT>)
      {
        L0_ = Ll_;
      }

      const auto axis = [this](RealT lm, RealT l1, RealT l2) -> Axis
      {
        if (Ll_ <= ZERO<RealT> || lm <= ZERO<RealT> || l1 <= ZERO<RealT> || l2 <= ZERO<RealT>)
        {
          return {};
        }
        const RealT stator = ONE<RealT> / Ll_;
        const RealT rotor1 = ONE<RealT> / l1;
        const RealT rotor2 = ONE<RealT> / l2;
        return {ONE<RealT> / lm, stator, rotor1, rotor2, stator + rotor1 + rotor2};
      };
      d_ = axis(Lmd_, Llfd_, Ll1d_);
      q_ = axis(Lmq_, Ll1q_, Ll2q_);
      if (L0_ > ZERO<RealT>)
      {
        inv_L0_ = ONE<RealT> / L0_;
      }

      omega_b_ = TWO<RealT> * std::numbers::pi_v<RealT> * f_b_;
      if (S_ > ZERO<RealT> && V_ > ZERO<RealT> && Lmd_ > ZERO<RealT>)
      {
        Vpk_ = std::sqrt(TWO<RealT> / THREE<RealT>) * V_;
        Ipk_ = std::numbers::sqrt2_v<RealT> * S_ / (std::numbers::sqrt3_v<RealT> * V_);
        kfd_ = Rfd_ / Lmd_;
        if (Rfe_ > ZERO<RealT>)
        {
          Gfe_ = S_ / (Rfe_ * V_ * V_);
        }
      }

      if (S12_ > ZERO<RealT> && S10_ >= ZERO<RealT>)
      {
        // Per-unit flux of the saturation factor `S12`
        constexpr RealT FLUX = 1.2;

        const RealT s = std::sqrt(S10_ / S12_);
        SA_           = std::min((FLUX * s + ONE<RealT>) / (s + ONE<RealT>), (FLUX * s - ONE<RealT>) / (s - ONE<RealT>) );
        SB_           = S12_ / ((SA_ - FLUX) * (SA_ - FLUX));
      }
    }

    /**
     * @brief Check the ratings and the per-unit parameters.
     *
     * @return The failed checks.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Machine<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("Machine");
      checks.check(S_ > ZERO<RealT> && V_ > ZERO<RealT> && f_b_ > ZERO<RealT>, "S, V and f must be positive");
      checks.check(H_ > ZERO<RealT>, "H must be positive");
      checks.check(F_ >= ZERO<RealT> && Rs_ >= ZERO<RealT>, "F and Rs must be nonnegative");
      checks.check(Ll_ > ZERO<RealT> && Lmd_ > ZERO<RealT> && Lmq_ > ZERO<RealT> && L0_ > ZERO<RealT>, "Ll, Lmd, Lmq and L0 must be positive");
      checks.check(Llfd_ > ZERO<RealT> && Ll1d_ > ZERO<RealT> && Ll1q_ > ZERO<RealT> && Ll2q_ > ZERO<RealT>,
                   "Llfd, Ll1d, Ll1q and Ll2q must be positive");
      checks.check(Rfd_ > ZERO<RealT> && R1d_ > ZERO<RealT> && R1q_ > ZERO<RealT> && R2q_ > ZERO<RealT>, "Rfd, R1d, R1q and R2q must be positive");
      checks.check(S10_ >= ZERO<RealT> && S12_ >= ZERO<RealT>, "S10 and S12 must be nonnegative");
      checks.check(Rfe_ >= ZERO<RealT>, "Rfe must be nonnegative");
      return checks;
    }

    /**
     * @brief Phase currents a, b and c of the currents `id`, `iq` and `i0`
     *        at angle `theta`.
     *
     * @tparam T Scalar type of the kernel.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] std::array<T, Machine<scalar_type, index_type>::PHASES>
    Machine<scalar_type, index_type>::phaseCurrents(T theta, T id, T iq, T i0)
    {
      const T c     = std::cos(theta);
      const T s     = std::sin(theta);
      const T alpha = id * c - iq * s;
      const T beta  = id * s + iq * c;
      return {i0 + alpha,
              i0 - (alpha - std::numbers::sqrt3_v<RealT> * beta) / TWO<RealT>,
              i0 - (alpha + std::numbers::sqrt3_v<RealT> * beta) / TWO<RealT>};
    }

    /**
     * @brief Reciprocal saturation factor of an axis, even in its air-gap
     *        flux linkage `psi`.
     *
     * @tparam T Scalar type of the kernel.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] T Machine<scalar_type, index_type>::inverseSaturation(T psi) const
    {
      return ONE<RealT> + SB_ * (Math::qramp(psi - SA_) + Math::qramp(-psi - SA_));
    }

    /**
     * @brief Winding quantities, in MachineOutputs order, from the flux
     *        linkages.
     *
     * @tparam T Scalar type of the kernel.
     * @param[in] y Internal values of the machine.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] auto Machine<scalar_type, index_type>::windingValues(const T* y) const -> std::array<T, NUM_WINDING_VALUES>
    {
      const T    id = d_.inv_Ll * (y[PSIAD] - y[PSID]);
      const T    iq = q_.inv_Ll * (y[PSIAQ] - y[PSIQ]);
      const T    i0 = -inv_L0_ * y[PSI0];
      const auto is = phaseCurrents(y[THETA], id, iq, i0);
      return {id,
              iq,
              i0,
              d_.inv_L1 * (y[PSIFD] - y[PSIAD]),
              d_.inv_L2 * (y[PSI1D] - y[PSIAD]),
              q_.inv_L1 * (y[PSI1Q] - y[PSIAQ]),
              q_.inv_L2 * (y[PSI2Q] - y[PSIAQ]),
              y[PSID] * iq - y[PSIQ] * id,
              is[0],
              is[1],
              is[2]};
    }

    /**
     * @brief Write the rotor, the flux linkage and the air-gap flux equations.
     *
     * The indices are constants, since Enzyme needs the layout at compile time.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    Machine<scalar_type, index_type>::internalResidual(const T* y, const T* yp, const T* w, const T*, T* f) const
    {
      const auto a        = windingValues(y);
      const auto [vd, vq] = park(w + BUS, y[THETA]);
      const T v0          = (w[BUS] + w[BUS + 1] + w[BUS + 2]) / (THREE<RealT> * Vpk_);

      f[THETA] = yp[THETA] - omega_b_ * y[OMEGA];
      f[OMEGA] = TWO<RealT> * H_ * yp[OMEGA] - w[PM] / y[OMEGA] + a[TE] + F_ * y[OMEGA];
      f[PSID]  = yp[PSID] / omega_b_ - y[OMEGA] * y[PSIQ] - Rs_ * a[ID] - vd / Vpk_;
      f[PSIQ]  = yp[PSIQ] / omega_b_ + y[OMEGA] * y[PSID] - Rs_ * a[IQ] - vq / Vpk_;
      f[PSI0]  = yp[PSI0] / omega_b_ - Rs_ * a[I0] - v0;
      f[PSIFD] = yp[PSIFD] / omega_b_ + Rfd_ * a[IFD] - kfd_ * w[EFD];
      f[PSI1D] = yp[PSI1D] / omega_b_ + R1d_ * a[I1D];
      f[PSI1Q] = yp[PSI1Q] / omega_b_ + R1q_ * a[I1Q];
      f[PSI2Q] = yp[PSI2Q] / omega_b_ + R2q_ * a[I2Q];

      f[PSIAD] = y[PSIAD] - d_.airGapFlux(y[PSID], y[PSIFD], y[PSI1D], inverseSaturation(y[PSIAD]));
      f[PSIAQ] = y[PSIAQ] - q_.airGapFlux(y[PSIQ], y[PSI1Q], y[PSI2Q], inverseSaturation(y[PSIAQ]));
    }

    /**
     * @brief Write the stator currents, less the iron-loss currents, to the
     *        residual of the bus voltage.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    Machine<scalar_type, index_type>::externalResidual(const T* y, const T*, const T* w, const T*, T* f) const
    {
      const auto a = windingValues(y);
#pragma clang loop unroll(full)
      for (size_t n = 0; n < PHASES; ++n)
      {
        f[BUS + n] = Ipk_ * a[IS + n] - Gfe_ * w[BUS + n];
      }
    }

    /**
     * @brief Write winding quantity `output`.
     *
     * The winding values are in MachineOutputs order, and only `is`, the
     * last, has more than one entry.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    Machine<scalar_type, index_type>::output(OutputsT output, const T* y, const T*, const T*, const T*, T* f) const
    {
      const auto a = windingValues(y);
      const auto k = static_cast<size_t>(output);

      f[0] = a[k];
      if (output == MachineOutputs::is)
      {
        f[1] = a[k + 1];
        f[2] = a[k + 2];
      }
    }
  } // namespace EMT
} // namespace GridKit
