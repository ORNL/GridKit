/**
 * @file TransformerImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT transformer.
 *
 * The model equations are in the `README.md` of this directory, and the
 * kernel arguments are defined in `ComponentImpl.hpp`.
 */

#pragma once

#include <array>
#include <cmath>
#include <numbers>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/Transformer/Transformer.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct a transformer and its per-unit bases.
     *
     * Layout: internal variables `i12`, `psi1` and `psi2` of size 3,
     * terminals `bus1` and `bus2` of size 3, and six outputs of size 3.
     * `verify` reports invalid data.
     */
    template <typename scalar_type, typename index_type>
    Transformer<scalar_type, index_type>::Transformer(const ModelDataT& data)
      : ComponentT(std::type_identity<Transformer>{}),
        N_(data.N),
        S_(data.S),
        V1_(data.V1),
        V2_(data.V2),
        f_b_(data.f),
        P1_(data.P1),
        P2_(data.P2),
        tap_(data.tap),
        phase_(data.phase),
        R_(data.R),
        X_(data.X),
        I0_(data.I0),
        P0_(data.P0),
        knee_(data.knee),
        Lsat_(data.Lsat),
        split_(data.split)
    {
      if (Lsat_ == ZERO<RealT>)
      {
        Lsat_ = TWO<RealT> * X_;
      }

      // An invalid map is replaced by the identity, which `verify` reports
      connectionMap(P1_, w1_.P);
      connectionMap(P2_, w2_.P);

      // A map without a phase shift keeps its exact entries
      if (phase_ != ZERO<RealT>)
      {
        shiftPhase(w1_.P, phase_);
      }

      omega_b_ = TWO<RealT> * std::numbers::pi_v<RealT> * f_b_;
      Rw_      = HALF<RealT> * R_;
      if (S_ > ZERO<RealT>)
      {
        Gc_ = P0_ / S_;
      }

      const auto initializeWinding = [this](Winding& winding, RealT voltage, RealT turns, RealT beta)
      {
        const RealT Vw = windingVoltage(winding.P, voltage);
        if (Vw > ZERO<RealT>)
        {
          winding.Vpk = std::numbers::sqrt2_v<RealT> * Vw;
          winding.Ipk = std::numbers::sqrt2_v<RealT> * S_ / (THREE<RealT> * Vw);
        }
        winding.turns = turns;
        winding.beta  = beta;
        winding.G     = beta * Gc_;
        winding.gain  = turns / (ONE<RealT> + turns * turns * Rw_ * winding.G);
      };
      initializeWinding(w1_, V1_, ONE<RealT>, split_);
      initializeWinding(w2_, V2_, tap_, ONE<RealT> - split_);

      const RealT Im2 = I0_ * I0_ - Gc_ * Gc_;
      if (Im2 > ZERO<RealT> && Lsat_ > ZERO<RealT>)
      {
        inv_Lm_ = std::sqrt(Im2);
        slope_  = ONE<RealT> / Lsat_ - inv_Lm_;
      }
    }

    /**
     * @brief Check the ratings, the per-unit parameters and the connection maps.
     *
     * @return The failed checks.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Transformer<scalar_type, index_type>::verify() const
    {
      MapT map;

      Model::ConfigurationChecks checks("Transformer");
      checks.check(N_ == PHASES, "N must be 3");
      checks.check(S_ > ZERO<RealT> && V1_ > ZERO<RealT> && V2_ > ZERO<RealT> && f_b_ > ZERO<RealT>, "S, V1, V2 and f must be positive");
      checks.check(tap_ > ZERO<RealT>, "tap must be positive");
      checks.check(R_ >= ZERO<RealT>, "R must be nonnegative");
      checks.check(X_ > ZERO<RealT>, "X must be positive");
      checks.check(P0_ >= ZERO<RealT>, "P0 must be nonnegative");
      checks.check(I0_ > Gc_, "I0 must exceed P0 / S");
      checks.check(knee_ > ZERO<RealT> && Lsat_ > ZERO<RealT>, "knee and Lsat must be positive");
      checks.check(split_ >= ZERO<RealT> && split_ <= ONE<RealT>, "split must be between 0 and 1");
      checks.check(connectionMap(P1_, map) && connectionMap(P2_, map),
                   "each column of P1 and P2 must have one 1, at most one -1 and zeros elsewhere");
      checks.check(slope_ > ZERO<RealT>, "the magnetizing inductance must exceed Lsat");
      return checks;
    }

    /**
     * @brief Copy a connection map to row-major storage.
     *
     * @param[in] matrix Map of 3 rows and 3 columns, empty if omitted.
     * @param[out] map Row-major map, the identity if `matrix` is omitted or invalid.
     * @return True if `matrix` is omitted, or each of its columns has one 1,
     *         at most one -1 and zeros elsewhere.
     */
    template <typename scalar_type, typename index_type>
    bool Transformer<scalar_type, index_type>::connectionMap(const Matrix<RealT>& matrix, MapT& map)
    {
      map.fill(ZERO<RealT>);
      for (size_t n = 0; n < PHASES; ++n)
      {
        map[PHASES * n + n] = ONE<RealT>;
      }
      if (matrix.empty())
      {
        return true;
      }

      bool valid = matrix.size() == PHASES;
      for (const auto& row : matrix)
      {
        valid = valid && row.size() == PHASES;
      }
      if (!valid)
      {
        return false;
      }

      for (size_t n = 0; n < PHASES; ++n)
      {
        size_t plus  = 0;
        size_t minus = 0;
        size_t zero  = 0;
        for (size_t m = 0; m < PHASES; ++m)
        {
          if (matrix[m][n] == ONE<RealT>)
          {
            ++plus;
          }
          if (matrix[m][n] == MINUS_ONE<RealT>)
          {
            ++minus;
          }
          if (matrix[m][n] == ZERO<RealT>)
          {
            ++zero;
          }
        }
        valid = valid && plus == 1 && minus <= 1 && plus + minus + zero == PHASES;
      }
      if (!valid)
      {
        return false;
      }

      for (size_t m = 0; m < PHASES; ++m)
      {
        for (size_t n = 0; n < PHASES; ++n)
        {
          map[PHASES * m + n] = matrix[m][n];
        }
      }
      return true;
    }

    /**
     * @brief Rotate the bus side of a connection map about the zero-sequence axis.
     *
     * @param[in,out] map Row-major connection map of the terminal.
     * @param[in] angle Angle by which the positive-sequence voltage of the bus
     *                  leads that of the winding.
     */
    template <typename scalar_type, typename index_type>
    void Transformer<scalar_type, index_type>::shiftPhase(MapT& map, RealT angle)
    {
      constexpr RealT                     SHIFT = TWO<RealT> * std::numbers::pi_v<RealT> / THREE<RealT>;
      constexpr std::array<RealT, PHASES> GAMMA{ZERO<RealT>, -SHIFT, SHIFT};

      const MapT winding = map;
      for (size_t m = 0; m < PHASES; ++m)
      {
        for (size_t n = 0; n < PHASES; ++n)
        {
          RealT entry = ZERO<RealT>;
          for (size_t k = 0; k < PHASES; ++k)
          {
            const RealT rotation  = (ONE<RealT> + TWO<RealT> * std::cos(angle + GAMMA[m] - GAMMA[k])) / THREE<RealT>;
            entry                += rotation * winding[PHASES * k + n];
          }
          map[PHASES * m + n] = entry;
        }
      }
    }

    /**
     * @brief Rated RMS voltage across a winding.
     *
     * @param[in] map Row-major connection map of the terminal.
     * @param[in] V Rated line-to-line RMS voltage of the terminal.
     * @return The phase voltage scaled by the first column of `map` applied
     *         to the positive-sequence unit phasors.
     */
    template <typename scalar_type, typename index_type>
    typename Transformer<scalar_type, index_type>::RealT
    Transformer<scalar_type, index_type>::windingVoltage(const MapT& map, RealT V)
    {
      constexpr RealT SHIFT = TWO<RealT> * std::numbers::pi_v<RealT> / THREE<RealT>;

      const RealT real = map[0] + (map[PHASES] + map[2 * PHASES]) * std::cos(SHIFT);
      const RealT imag = (map[2 * PHASES] - map[PHASES]) * std::sin(SHIFT);
      return V / std::numbers::sqrt3_v<RealT> * std::hypot(real, imag);
    }

    /**
     * @brief Magnetizing current of flux linkage `psi`, a two-slope curve
     *        with its knees at plus and minus `knee_`.
     *
     * @tparam T Scalar type of the kernel.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] T Transformer<scalar_type, index_type>::magnetizing(T psi) const
    {
      return inv_Lm_ * psi + slope_ * (Math::ramp(psi - knee_) - Math::ramp(-psi - knee_));
    }

    /**
     * @brief Winding quantities of one phase of `winding`.
     *
     * @tparam T Scalar type of the kernel.
     * @param[in] winding Winding 1 or winding 2.
     * @param[in] u Per-unit voltage of the phase at the terminal.
     * @param[in] psi Magnetizing flux linkage of the node.
     * @param[in] series Series leakage current into the node.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] auto Transformer<scalar_type, index_type>::windingValues(const Winding& winding, T u, T psi, T series) const
        -> WindingValues<T>
    {
      const T m  = winding.beta * magnetizing(psi);
      const T iw = winding.gain * (series + m + winding.turns * winding.G * u);
      const T e  = winding.turns * (u - Rw_ * iw);
      return {e, m + winding.G * e, iw};
    }

    /**
     * @brief Write the series leakage and the magnetizing equations of each phase.
     *
     * The indices are constants, since Enzyme needs the layout at compile time.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    Transformer<scalar_type, index_type>::internalResidual(const T* y, const T* yp, const T* w, const T*, T* f) const
    {
#pragma clang loop unroll(full)
      for (size_t n = 0; n < PHASES; ++n)
      {
        const auto first  = windingValues(w1_, w1_.voltage(w + BUS1, n), y[PSI1 + n], y[I12 + n]);
        const auto second = windingValues(w2_, w2_.voltage(w + BUS2, n), y[PSI2 + n], -y[I12 + n]);

        f[I12 + n]  = X_ / omega_b_ * yp[I12 + n] + second.e - first.e;
        f[PSI1 + n] = yp[PSI1 + n] / omega_b_ - first.e;
        f[PSI2 + n] = yp[PSI2 + n] / omega_b_ - second.e;
      }
    }

    /**
     * @brief Write the winding currents to the residuals of the bus voltages.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    Transformer<scalar_type, index_type>::externalResidual(const T* y, const T*, const T* w, const T*, T* f) const
    {
      T iw1[PHASES];
      T iw2[PHASES];
#pragma clang loop unroll(full)
      for (size_t n = 0; n < PHASES; ++n)
      {
        iw1[n] = windingValues(w1_, w1_.voltage(w + BUS1, n), y[PSI1 + n], y[I12 + n]).iw;
        iw2[n] = windingValues(w2_, w2_.voltage(w + BUS2, n), y[PSI2 + n], -y[I12 + n]).iw;
      }
      w1_.inject(iw1, f + BUS1);
      w2_.inject(iw2, f + BUS2);
    }

    /**
     * @brief Write winding quantity `output` of each phase.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    Transformer<scalar_type, index_type>::output(OutputsT output, const T* y, const T*, const T* w, const T*, T* f) const
    {
#pragma clang loop unroll(full)
      for (size_t n = 0; n < PHASES; ++n)
      {
        const auto first  = windingValues(w1_, w1_.voltage(w + BUS1, n), y[PSI1 + n], y[I12 + n]);
        const auto second = windingValues(w2_, w2_.voltage(w + BUS2, n), y[PSI2 + n], -y[I12 + n]);

        // In TransformerOutputs order
        const std::array<T, Utilities::enum_size<TransformerOutputs>()> values{first.e, second.e, first.im, second.im, first.iw, second.iw};
        f[n] = values[static_cast<size_t>(output)];
      }
    }
  } // namespace EMT
} // namespace GridKit
