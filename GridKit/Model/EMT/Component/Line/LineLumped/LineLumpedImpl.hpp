/**
 * @file LineLumpedImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT lumped line.
 */

#pragma once

#include <algorithm>
#include <numeric>
#include <vector>

#include <GridKit/Model/EMT/Component/Line/LineLumped/LineLumped.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Construct a line of `K` conductors on `N` phases.
     *
     * Layout: internal variable `i12` of size `K`, terminals `bus1` and
     * `bus2` of size `N`. Submodel `z` reads `i12` as its input, and its
     * output is added to the residual of `i12`. When the case gives `Yp`,
     * submodels `y1` and `y2` read the conductor voltages at their terminals,
     * and their outputs are added to the residuals of those voltages.
     */
    template <typename scalar_type, typename index_type>
    LineLumped<scalar_type, index_type>::LineLumped(const ModelDataT& data)
      : N_(data.N),
        K_(data.K),
        dx_(data.dx),
        shunt_(!omitted(data.Yp))
    {
      for (const IdxT conductor : data.conductors)
      {
        c_.push_back(conductor - 1);
      }
      this->setLayout({K_}, {N_, N_}, {}, {});

      std::vector<IdxT> i12(K_);
      std::vector<IdxT> v1;
      std::vector<IdxT> v2;
      std::iota(i12.begin(), i12.end(), this->local(LineLumpedInternalVariables::i12));
      for (const IdxT phase : c_)
      {
        v1.push_back(this->local(LineLumpedTerminals::bus1, phase));
        v2.push_back(this->local(LineLumpedTerminals::bus2, phase));
      }

      // An operator that is not K by K does not connect, which `verify` reports
      auto& z = this->template addSubmodel<VectorFitT>(LineLumpedSubmodels::z, scaled(data.Zp, dx_));
      z.connectInput(VectorFitInputs::u, *this, i12);
      z.connectOutput(VectorFitOutputs::y, *this, i12);

      if (shunt_)
      {
        // Each shunt carries half of the line admittance, and its current leaves the bus
        const auto half = scaled(data.Yp, -dx_ / TWO<RealT>);
        auto&      y1   = this->template addSubmodel<VectorFitT>(LineLumpedSubmodels::y1, half);
        auto&      y2   = this->template addSubmodel<VectorFitT>(LineLumpedSubmodels::y2, half);
        y1.connectInput(VectorFitInputs::u, *this, v1);
        y1.connectOutput(VectorFitOutputs::y, *this, v1);
        y2.connectInput(VectorFitInputs::u, *this, v2);
        y2.connectOutput(VectorFitOutputs::y, *this, v2);
      }
    }

    /**
     * @brief Submodel `symbol`, an impedance or admittance operator.
     *
     * @throws std::invalid_argument for `y1` or `y2` of a line without `Yp`.
     */
    template <typename scalar_type, typename index_type>
    typename LineLumped<scalar_type, index_type>::VectorFitT& LineLumped<scalar_type, index_type>::submodel(SubmodelsT symbol)
    {
      return ComponentT::template submodel<VectorFitT>(symbol);
    }

    /**
     * @brief Check the phase and conductor counts, the phase of each
     *        conductor, the length, and that the operators are `K` by `K`.
     *
     * @return The failed checks.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks LineLumped<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("LineLumped");
      checks.check(N_ > 0 && K_ > 0, "N and K must be positive");
      checks.check(c_.size() == K_, "conductors must have K values");

      bool              named = true;
      std::vector<bool> covered(N_, false);
      for (const IdxT phase : c_)
      {
        if (phase < N_)
        {
          covered[phase] = true;
        }
        else
        {
          named = false;
        }
      }
      checks.check(named, "conductors must name phases 1 to N");
      checks.check(std::ranges::all_of(covered, [](bool phase)
                                       { return phase; }),
                   "every phase must have a conductor");
      checks.check(dx_ > ZERO<RealT>, "dx must be positive");

      const auto& z = ComponentT::template submodel<VectorFitT>(LineLumpedSubmodels::z);
      checks.check(z.outputs() == K_ && z.inputs() == K_, "Zp must be K by K");
      if (shunt_)
      {
        const auto& y1 = ComponentT::template submodel<VectorFitT>(LineLumpedSubmodels::y1);
        checks.check(y1.outputs() == K_ && y1.inputs() == K_, "Yp must be K by K");
      }
      return checks;
    }

    /**
     * @brief Set the current tolerances to `rel_tol` times the system power
     *        base over the base voltage of each conductor's phase at terminal 1.
     *
     * @param[in] rel_tol Relative tolerance of the integrator.
     * @return 0.
     *
     * @pre The bus at `bus1` has set its tolerances with the same `rel_tol`.
     */
    template <typename scalar_type, typename index_type>
    int LineLumped<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      const auto i12  = this->absoluteTolerance(LineLumpedInternalVariables::i12);
      const IdxT BUS1 = this->terminalOffset(LineLumpedTerminals::bus1);
      for (IdxT k = 0; k < K_; ++k)
      {
        const RealT v_base = this->externalTolerance(BUS1 + c_[k]) / rel_tol;
        i12[k]             = rel_tol * this->va_system_base_ / v_base;
      }
      return 0;
    }

    /**
     * @brief Write the voltage across each conductor to the residual of its
     *        series current, to which the output of the impedance is added.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    void LineLumped<scalar_type, index_type>::internalResidual(const T*, const T*, const T* w, const T*, T* f) const
    {
      const IdxT I12  = this->local(LineLumpedInternalVariables::i12);
      const IdxT BUS1 = this->terminalOffset(LineLumpedTerminals::bus1);
      const IdxT BUS2 = this->terminalOffset(LineLumpedTerminals::bus2);
      for (IdxT k = 0; k < K_; ++k)
      {
        f[I12 + k] = w[BUS2 + c_[k]] - w[BUS1 + c_[k]];
      }
    }

    /**
     * @brief Add the series current of each conductor to the residuals of the
     *        voltages of its phase.
     *
     * @pre The residuals at the terminals are zero on entry.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    void LineLumped<scalar_type, index_type>::externalResidual(const T* y, const T*, const T*, const T*, T* f) const
    {
      const IdxT I12  = this->local(LineLumpedInternalVariables::i12);
      const IdxT BUS1 = this->terminalOffset(LineLumpedTerminals::bus1);
      const IdxT BUS2 = this->terminalOffset(LineLumpedTerminals::bus2);
      for (IdxT k = 0; k < K_; ++k)
      {
        f[BUS1 + c_[k]] -= y[I12 + k];
        f[BUS2 + c_[k]] += y[I12 + k];
      }
    }
  } // namespace EMT
} // namespace GridKit
