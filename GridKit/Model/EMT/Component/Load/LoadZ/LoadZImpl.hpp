/**
 * @file LoadZImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT impedance load.
 *
 * The model equations are in the `README.md` of this directory, and the
 * kernel arguments are defined in `ComponentImpl.hpp`.
 */

#pragma once

#include <numeric>
#include <vector>

#include <GridKit/Model/EMT/Component/Load/LoadZ/LoadZ.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct an impedance load with `N` phases.
     *
     * Layout: internal variable `i` of size `N`, terminal `bus` of size `N`.
     * Submodel `z` reads `i` as its input, and its output is added to the
     * residual of `i`.
     */
    template <typename scalar_type, typename index_type>
    LoadZ<scalar_type, index_type>::LoadZ(const ModelDataT& data)
      : N_(data.N)
    {
      this->setLayout({N_}, {N_}, {}, {});

      std::vector<IdxT> current(N_);
      std::iota(current.begin(), current.end(), this->local(LoadZInternalVariables::i));

      // An impedance that is not N by N does not connect, which `verify` reports
      auto& z = this->template addSubmodel<VectorFitT>(LoadZSubmodels::z, data.Z);
      z.connectInput(VectorFitInputs::u, *this, current);
      z.connectOutput(VectorFitOutputs::y, *this, current);
    }

    /**
     * @brief Submodel `symbol`, the impedance operator.
     */
    template <typename scalar_type, typename index_type>
    typename LoadZ<scalar_type, index_type>::VectorFitT& LoadZ<scalar_type, index_type>::submodel(SubmodelsT symbol)
    {
      return ComponentT::template submodel<VectorFitT>(symbol);
    }

    /**
     * @brief Check that `N` is positive and that the impedance is `N` by `N`.
     *
     * @return The failed checks.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks LoadZ<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("LoadZ");
      checks.check(N_ > 0, "N must be positive");

      const auto& z = ComponentT::template submodel<VectorFitT>(LoadZSubmodels::z);
      checks.check(z.outputs() == N_ && z.inputs() == N_, "Z must be N by N");
      return checks;
    }

    /**
     * @brief Set the current tolerances to `rel_tol` times the system power
     *        base over the bus base voltage.
     *
     * @param[in] rel_tol Relative tolerance of the integrator.
     * @return 0.
     *
     * @pre The bus at the terminal has set its tolerances with the same `rel_tol`.
     */
    template <typename scalar_type, typename index_type>
    int LoadZ<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      const auto i   = this->absoluteTolerance(LoadZInternalVariables::i);
      const IdxT BUS = this->terminalOffset(LoadZTerminals::bus);
      for (IdxT n = 0; n < N_; ++n)
      {
        const RealT v_base = this->externalTolerance(BUS + n) / rel_tol;
        i[n]               = rel_tol * this->va_system_base_ / v_base;
      }
      return 0;
    }

    /**
     * @brief Write the bus voltage to the residual of `i`.
     *
     * The output of the impedance submodel is added to the same residual.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    void LoadZ<scalar_type, index_type>::internalResidual(const T*, const T*, const T* w, const T*, T* f) const
    {
      const IdxT I   = this->local(LoadZInternalVariables::i);
      const IdxT BUS = this->terminalOffset(LoadZTerminals::bus);
      for (IdxT n = 0; n < N_; ++n)
      {
        f[I + n] = w[BUS + n];
      }
    }

    /**
     * @brief Write `i` to the residual of the bus voltage.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    void LoadZ<scalar_type, index_type>::externalResidual(const T* y, const T*, const T*, const T*, T* f) const
    {
      const IdxT I   = this->local(LoadZInternalVariables::i);
      const IdxT BUS = this->terminalOffset(LoadZTerminals::bus);
      for (IdxT n = 0; n < N_; ++n)
      {
        f[BUS + n] = y[I + n];
      }
    }
  } // namespace EMT
} // namespace GridKit
