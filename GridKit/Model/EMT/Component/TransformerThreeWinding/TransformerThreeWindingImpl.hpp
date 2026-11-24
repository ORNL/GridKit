/**
 * @file TransformerThreeWindingImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT three-winding transformer.
 *
 * The model equations are in the `README.md` of this directory, and the
 * kernel arguments are defined in `ComponentImpl.hpp`.
 */

#pragma once

#include <array>
#include <cmath>
#include <numbers>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/Component/TransformerThreeWinding/TransformerThreeWinding.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct a three-winding transformer and its per-unit bases.
     *
     * Layout: internal variables `i13`, `i23`, `psi1`, `psi2` and `psi3` of
     * size 3, terminals `bus1`, `bus2` and `bus3` of size 3, and nine outputs
     * of size 3. `verify` reports invalid data.
     */
    template <typename scalar_type, typename index_type>
    TransformerThreeWinding<scalar_type, index_type>::TransformerThreeWinding(const ModelDataT& data)
      : ComponentT(std::type_identity<TransformerThreeWinding>{}),
        N_(data.N),
        S_(data.S),
        V1_(data.V1),
        V2_(data.V2),
        V3_(data.V3),
        f_b_(data.f),
        P1_(data.P1),
        P2_(data.P2),
        P3_(data.P3),
        tap_(data.tap),
        phase_(data.phase),
        R12_(data.R12),
        R23_(data.R23),
        R31_(data.R31),
        X12_(data.X12),
        X23_(data.X23),
        X31_(data.X31),
        I0_(data.I0),
        P0_(data.P0),
        knee_(data.knee),
        Lsat_(data.Lsat)
    {
      // An invalid map is replaced by the identity, which `verify` reports
      connectionMap(P1_, w1_.P);
      connectionMap(P2_, w2_.P);
      connectionMap(P3_, w3_.P);

      // A map without a phase shift keeps its exact entries
      if (phase_ != ZERO<RealT>)
      {
        shiftPhase(w1_.P, phase_);
      }

      omega_b_ = TWO<RealT> * std::numbers::pi_v<RealT> * f_b_;
      R3_      = HALF<RealT> * (R23_ + R31_ - R12_);
      X3_      = HALF<RealT> * (X23_ + X31_ - X12_);
      beta_    = ONE<RealT> / THREE<RealT>;
      if (S_ > ZERO<RealT>)
      {
        Gc_ = P0_ / S_;
      }
      G_ = beta_ * Gc_;

      const auto initializeWinding = [this](Winding& winding, RealT voltage, RealT turns)
      {
        const RealT Vw = windingVoltage(winding.P, voltage);
        if (Vw > ZERO<RealT>)
        {
          winding.Vpk = std::numbers::sqrt2_v<RealT> * Vw;
          winding.Ipk = std::numbers::sqrt2_v<RealT> * S_ / (THREE<RealT> * Vw);
        }
        winding.turns = turns;
      };
      initializeWinding(w1_, V1_, ONE<RealT>);
      initializeWinding(w2_, V2_, tap_);
      initializeWinding(w3_, V3_, tap_);

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
    Model::ConfigurationChecks TransformerThreeWinding<scalar_type, index_type>::verify() const
    {
      MapT map;

      Model::ConfigurationChecks checks("TransformerThreeWinding");
      checks.check(N_ == PHASES, "N must be 3");
      checks.check(S_ > ZERO<RealT> && V1_ > ZERO<RealT> && V2_ > ZERO<RealT> && V3_ > ZERO<RealT> && f_b_ > ZERO<RealT>,
                   "S, V1, V2, V3 and f must be positive");
      checks.check(tap_ > ZERO<RealT>, "tap must be positive");
      checks.check(R12_ >= ZERO<RealT> && R23_ >= ZERO<RealT> && R31_ >= ZERO<RealT>, "R12, R23 and R31 must be nonnegative");
      checks.check(X12_ > ZERO<RealT> && X23_ > ZERO<RealT> && X31_ > ZERO<RealT>, "X12, X23 and X31 must be positive");
      checks.check(R3_ * R3_ <= R23_ * R31_, "R12, R23 and R31 must give a positive semidefinite resistance matrix");
      checks.check(X3_ * X3_ < X23_ * X31_, "X12, X23 and X31 must give a positive definite leakage matrix");
      checks.check(P0_ >= ZERO<RealT>, "P0 must be nonnegative");
      checks.check(I0_ > Gc_, "I0 must exceed P0 / S");
      checks.check(knee_ > ZERO<RealT> && Lsat_ > ZERO<RealT>, "knee and Lsat must be positive");
      checks.check(connectionMap(P1_, map) && connectionMap(P2_, map) && connectionMap(P3_, map),
                   "each column of P1, P2 and P3 must have one 1, at most one -1 and zeros elsewhere");
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
    bool TransformerThreeWinding<scalar_type, index_type>::connectionMap(const Matrix<RealT>& matrix, MapT& map)
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
    void TransformerThreeWinding<scalar_type, index_type>::shiftPhase(MapT& map, RealT angle)
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
    typename TransformerThreeWinding<scalar_type, index_type>::RealT
    TransformerThreeWinding<scalar_type, index_type>::windingVoltage(const MapT& map, RealT V)
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
    [[gnu::always_inline]] T TransformerThreeWinding<scalar_type, index_type>::magnetizing(T psi) const
    {
      return inv_Lm_ * psi + slope_ * (Math::ramp(psi - knee_) - Math::ramp(-psi - knee_));
    }

    /**
     * @brief Winding quantities of one phase of `winding`.
     *
     * @tparam T Scalar type of the kernel.
     * @param[in] winding Winding 1, 2 or 3.
     * @param[in] u Per-unit voltage of the phase at the terminal.
     * @param[in] psi Magnetizing flux linkage of the node.
     * @param[in] series Leakage current out of the node.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] auto TransformerThreeWinding<scalar_type, index_type>::windingValues(const Winding& winding, T u, T psi, T series) const
        -> WindingValues<T>
    {
      const T e  = winding.turns * u;
      const T im = beta_ * magnetizing(psi) + G_ * e;
      return {e, im, winding.turns * (series + im)};
    }

    /**
     * @brief Write the two leakage loops and the magnetizing equations of each phase.
     *
     * The indices are constants, since Enzyme needs the layout at compile time.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    TransformerThreeWinding<scalar_type, index_type>::internalResidual(const T* y, const T* yp, const T* w, const T*, T* f) const
    {
#pragma clang loop unroll(full)
      for (size_t n = 0; n < PHASES; ++n)
      {
        const T e1 = w1_.turns * w1_.voltage(w + BUS1, n);
        const T e2 = w2_.turns * w2_.voltage(w + BUS2, n);
        const T e3 = w3_.turns * w3_.voltage(w + BUS3, n);

        f[I13 + n]  = (X31_ * yp[I13 + n] + X3_ * yp[I23 + n]) / omega_b_ + R31_ * y[I13 + n] + R3_ * y[I23 + n] + e3 - e1;
        f[I23 + n]  = (X3_ * yp[I13 + n] + X23_ * yp[I23 + n]) / omega_b_ + R3_ * y[I13 + n] + R23_ * y[I23 + n] + e3 - e2;
        f[PSI1 + n] = yp[PSI1 + n] / omega_b_ - e1;
        f[PSI2 + n] = yp[PSI2 + n] / omega_b_ - e2;
        f[PSI3 + n] = yp[PSI3 + n] / omega_b_ - e3;
      }
    }

    /**
     * @brief Write the winding currents to the residuals of the bus voltages.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    TransformerThreeWinding<scalar_type, index_type>::externalResidual(const T* y, const T*, const T* w, const T*, T* f) const
    {
      T iw1[PHASES];
      T iw2[PHASES];
      T iw3[PHASES];
#pragma clang loop unroll(full)
      for (size_t n = 0; n < PHASES; ++n)
      {
        const T i3 = -(y[I13 + n] + y[I23 + n]);

        iw1[n] = windingValues(w1_, w1_.voltage(w + BUS1, n), y[PSI1 + n], y[I13 + n]).iw;
        iw2[n] = windingValues(w2_, w2_.voltage(w + BUS2, n), y[PSI2 + n], y[I23 + n]).iw;
        iw3[n] = windingValues(w3_, w3_.voltage(w + BUS3, n), y[PSI3 + n], i3).iw;
      }
      w1_.inject(iw1, f + BUS1);
      w2_.inject(iw2, f + BUS2);
      w3_.inject(iw3, f + BUS3);
    }

    /**
     * @brief Write winding quantity `output` of each phase.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    TransformerThreeWinding<scalar_type, index_type>::output(OutputsT output, const T* y, const T*, const T* w, const T*, T* f) const
    {
#pragma clang loop unroll(full)
      for (size_t n = 0; n < PHASES; ++n)
      {
        const T i3 = -(y[I13 + n] + y[I23 + n]);

        const auto first  = windingValues(w1_, w1_.voltage(w + BUS1, n), y[PSI1 + n], y[I13 + n]);
        const auto second = windingValues(w2_, w2_.voltage(w + BUS2, n), y[PSI2 + n], y[I23 + n]);
        const auto third  = windingValues(w3_, w3_.voltage(w + BUS3, n), y[PSI3 + n], i3);

        // In TransformerThreeWindingOutputs order
        const std::array<T, Utilities::enum_size<TransformerThreeWindingOutputs>()> values{first.e, second.e, third.e, first.im, second.im, third.im, first.iw, second.iw, third.iw};
        f[n] = values[static_cast<size_t>(output)];
      }
    }
  } // namespace EMT
} // namespace GridKit
