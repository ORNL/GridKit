/**
 * @file LoadRLCImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT parallel RLC load.
 *
 * The model equations are in the `README.md` of this directory, and the
 * kernel arguments are defined in `ComponentImpl.hpp`.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <optional>

#include <GridKit/Model/EMT/Component/Load/LoadRLC/LoadRLC.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct a parallel RLC load with `N` phases.
     *
     * Layout: internal variable `i_l` with one entry per phase that has an
     * inductor, terminal `bus` of size `N`, and output `i` of size `N`.
     * `verify` reports invalid data.
     */
    template <typename scalar_type, typename index_type>
    LoadRLC<scalar_type, index_type>::LoadRLC(const ModelDataT& data)
      : N_(data.N),
        R_(data.R),
        L_(data.L),
        C_(data.C)
    {
      // An array that does not have N entries gives no branch past its end, which `verify` reports
      const auto branch = [](const std::optional<PhaseValues>& values, IdxT n) -> std::optional<RealT>
      {
        if (!values || n >= values->size())
        {
          return std::nullopt;
        }
        return (*values)[n];
      };

      phases_.resize(N_);
      IdxT inductors = 0;
      for (IdxT n = 0; n < N_; ++n)
      {
        Phase& phase = phases_[n];
        if (const auto r = branch(R_, n))
        {
          phase.G = ONE<RealT> / *r;
        }
        if (const auto l = branch(L_, n))
        {
          phase.L   = *l;
          phase.i_l = inductors++;
        }
        if (const auto c = branch(C_, n))
        {
          phase.C = *c;
        }
      }

      this->setLayout({inductors}, {N_}, {}, {N_});

      const IdxT I_L = this->local(LoadRLCInternalVariables::i_l);
      for (Phase& phase : phases_)
      {
        phase.i_l += I_L;
      }
    }

    /**
     * @brief Current injection of the phase into the bus.
     *
     * @tparam T Scalar type of the kernel.
     * @param[in] y Internal values of the load.
     * @param[in] v Phase voltage.
     * @param[in] vp Derivative of the phase voltage.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    T LoadRLC<scalar_type, index_type>::Phase::current(const T* y, const T& v, const T& vp) const
    {
      T current = -(G * v + C * vp);
      if (L > ZERO<RealT>)
      {
        current += y[i_l];
      }
      return current;
    }

    /**
     * @brief Check that `N` is positive, that `R`, `L` and `C` each have `N`
     *        entries when given, and that each entry given is positive and
     *        finite.
     *
     * @return The failed checks.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks LoadRLC<scalar_type, index_type>::verify() const
    {
      const auto sized = [this](const std::optional<PhaseValues>& values)
      {
        return !values || values->size() == N_;
      };
      const auto positive = [](const std::optional<PhaseValues>& values)
      {
        if (!values)
        {
          return true;
        }
        return std::ranges::all_of(*values, [](const std::optional<RealT>& value)
                                   { return !value || (*value > ZERO<RealT> && std::isfinite(*value) && std::isfinite(ONE<RealT> / *value)); });
      };

      Model::ConfigurationChecks checks("LoadRLC");
      checks.check(N_ > 0, "N must be positive");
      checks.check(sized(R_) && sized(L_) && sized(C_), "R, L and C must each have N entries");
      checks.check(positive(R_) && positive(L_) && positive(C_), "each entry of R, L and C must be positive and finite");
      return checks;
    }

    /**
     * @brief Set the tolerances of the inductor currents and of the load
     *        current to `rel_tol` times the system power base over the bus
     *        base voltage.
     *
     * @param[in] rel_tol Relative tolerance of the integrator.
     * @return 0.
     *
     * @pre The bus at the terminal has set its tolerances with the same `rel_tol`.
     */
    template <typename scalar_type, typename index_type>
    int LoadRLC<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      const IdxT BUS = this->terminalOffset(LoadRLCTerminals::bus);
      const auto i   = this->outputTolerance(LoadRLCOutputs::i);
      auto*      tol = this->abs_tol_.getData();
      for (IdxT n = 0; n < N_; ++n)
      {
        const Phase& phase  = phases_[n];
        const RealT  v_base = this->externalTolerance(BUS + n) / rel_tol;
        i[n]                = rel_tol * this->va_system_base_ / v_base;
        if (phase.L > ZERO<RealT>)
        {
          tol[phase.i_l] = i[n];
        }
      }
      return 0;
    }

    /**
     * @brief Write the inductor equation of each phase that has an inductor.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    void LoadRLC<scalar_type, index_type>::internalResidual(const T*, const T* yp, const T* w, const T*, T* f) const
    {
      const IdxT BUS = this->terminalOffset(LoadRLCTerminals::bus);
      for (IdxT n = 0; n < N_; ++n)
      {
        const Phase& phase = phases_[n];
        if (phase.L > ZERO<RealT>)
        {
          f[phase.i_l] = phase.L * yp[phase.i_l] + w[BUS + n];
        }
      }
    }

    /**
     * @brief Write the load current to the residual of the bus voltage.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    void LoadRLC<scalar_type, index_type>::externalResidual(const T* y, const T*, const T* w, const T* wp, T* f) const
    {
      const IdxT BUS = this->terminalOffset(LoadRLCTerminals::bus);
      for (IdxT n = 0; n < N_; ++n)
      {
        f[BUS + n] = phases_[n].current(y, w[BUS + n], wp[BUS + n]);
      }
    }

    /**
     * @brief Write the load current `i`.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    void LoadRLC<scalar_type, index_type>::output(OutputsT, const T* y, const T*, const T* w, const T* wp, T* f) const
    {
      const IdxT BUS = this->terminalOffset(LoadRLCTerminals::bus);
      for (IdxT n = 0; n < N_; ++n)
      {
        f[n] = phases_[n].current(y, w[BUS + n], wp[BUS + n]);
      }
    }
  } // namespace EMT
} // namespace GridKit
