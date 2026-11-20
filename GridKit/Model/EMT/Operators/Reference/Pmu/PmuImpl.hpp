/**
 * @file PmuImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT phasor measurement unit.
 *
 * The model equations are in the `README.md` of this directory.
 */

#pragma once

#include <cmath>
#include <numbers>
#include <numeric>
#include <type_traits>
#include <vector>

#include <GridKit/Model/EMT/ComponentImpl.hpp>
#include <GridKit/Model/EMT/Operators/Reference/Pll/PllImpl.hpp>
#include <GridKit/Model/EMT/Operators/Reference/Pmu/Pmu.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct a phasor measurement unit.
     *
     * Layout: internal variables `vms`, `p` and `q`, inputs `v` and `i` of
     * size 3, and outputs `vrms`, `theta` and `omega`. Submodel `pll` reads
     * the same voltage, and input `pll` reads its angle and integrator state.
     */
    template <typename scalar_type, typename index_type>
    Pmu<scalar_type, index_type>::Pmu(const ModelDataT& data)
      : ComponentT(std::type_identity<Pmu>{})
    {
      std::vector<IdxT> v(NUM_INPUT_VARIABLES[static_cast<size_t>(PmuInputs::v)]);
      std::iota(v.begin(), v.end(), this->size() + this->inputOffset(PmuInputs::v));

      auto& pll = this->template addSubmodel<PllT>(PmuSubmodels::pll, data.pll);
      pll.connectInput(PllInputs::v, *this, v);
      this->connectInput(PmuInputs::pll, pll, {static_cast<IdxT>(PllInternalVariables::theta), static_cast<IdxT>(PllInternalVariables::x)});
      pll_ = &pll;
    }

    /**
     * @brief Submodel `symbol`.
     */
    template <typename scalar_type, typename index_type>
    typename Pmu<scalar_type, index_type>::PllT& Pmu<scalar_type, index_type>::submodel(SubmodelsT symbol)
    {
      return ComponentT::template submodel<PllT>(symbol);
    }

    /**
     * @brief A unit has no parameter of its own to check.
     *
     * @return No failed check.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks Pmu<scalar_type, index_type>::verify() const
    {
      return Model::ConfigurationChecks("Pmu");
    }

    /**
     * @brief Set the tolerances from the base voltage of the measured bus and
     *        the system power base.
     *
     * `theta` and `omega` take the tolerances that `Pll` sets for its own.
     *
     * @param[in] rel_tol Relative tolerance of the integrator.
     * @return 0.
     *
     * @pre The system has bound this model and set the tolerance of its bus.
     */
    template <typename scalar_type, typename index_type>
    int Pmu<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      const auto  voltage    = this->externalTolerance(V);
      const RealT omega_base = TWO<RealT> * std::numbers::pi_v<RealT> * this->freq_system_base_;

      this->absoluteTolerance(PmuInternalVariables::vms)[0] = voltage * voltage;
      this->absoluteTolerance(PmuInternalVariables::p)[0]   = rel_tol * this->va_system_base_;
      this->absoluteTolerance(PmuInternalVariables::q)[0]   = rel_tol * this->va_system_base_;
      this->outputTolerance(PmuOutputs::vrms)[0]            = voltage;
      this->outputTolerance(PmuOutputs::theta)[0]           = rel_tol;
      this->outputTolerance(PmuOutputs::omega)[0]           = rel_tol * omega_base;
      return 0;
    }

    /**
     * @brief Write the residuals of the mean-square voltage and of the powers.
     *
     * Each follows its instantaneous value over half a cycle of the system
     * frequency base. The mean over the phases is removed, so on a balanced
     * bus `vms` is the square of the line-to-line RMS voltage.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    Pmu<scalar_type, index_type>::internalResidual(const T* y, const T* yp, const T* w, const T*, T* f) const
    {
      const auto VMS = static_cast<size_t>(PmuInternalVariables::vms);
      const auto P   = static_cast<size_t>(PmuInternalVariables::p);
      const auto Q   = static_cast<size_t>(PmuInternalVariables::q);

      const RealT rate = TWO<RealT> * this->freq_system_base_;

      const T v_mean = (w[V] + w[V + 1] + w[V + 2]) / THREE<RealT>;
      const T va     = w[V] - v_mean;
      const T vb     = w[V + 1] - v_mean;
      const T vc     = w[V + 2] - v_mean;

      const T squared  = va * va + vb * vb + vc * vc;
      const T power    = va * w[I] + vb * w[I + 1] + vc * w[I + 2];
      const T reactive = ((vb - vc) * w[I] + (vc - va) * w[I + 1] + (va - vb) * w[I + 2]) / std::numbers::sqrt3_v<RealT>;

      f[VMS] = -yp[VMS] + rate * (squared - y[VMS]);
      f[P]   = -yp[P] + rate * (power - y[P]);
      f[Q]   = -yp[Q] + rate * (reactive - y[Q]);
    }

    /**
     * @brief Write output `output`.
     *
     * `vrms` is the square root of `vms`. `theta` is the angle of submodel
     * `pll`, and `omega` its frequency, which `Pll` writes from its own
     * variables and the measured voltage.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    [[gnu::always_inline]] void
    Pmu<scalar_type, index_type>::output(OutputsT output, const T* y, const T*, const T* w, const T*, T* f) const
    {
      if (output == PmuOutputs::theta)
      {
        f[0] = w[PLL + static_cast<size_t>(PllInternalVariables::theta)];
        return;
      }
      if (output == PmuOutputs::omega)
      {
        pll_->output(PllOutputs::omega, w + PLL, static_cast<const T*>(nullptr), w + V, static_cast<const T*>(nullptr), f);
        return;
      }

      const auto VMS = static_cast<size_t>(PmuInternalVariables::vms);

      // The square root has no derivative on a dead bus, where the product keeps the dependency on `vms`
      T vrms = ZERO<RealT> * y[VMS];
      if (y[VMS] > ZERO<RealT>)
      {
        vrms = std::sqrt(y[VMS]);
      }
      f[0] = vrms;
    }
  } // namespace EMT
} // namespace GridKit
