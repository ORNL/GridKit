/**
 * @file LineLumpedImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT lumped line.
 */

#pragma once

#include <algorithm>
#include <stdexcept>
#include <vector>

#include <GridKit/Model/EMT/Component/Line/LineLumped/LineLumped.hpp>
#include <GridKit/Model/EMT/ComponentImpl.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    LineLumped<scalar_type, index_type>::LineLumped(const ModelDataT& data)
      : N_(data.N),
        K_(data.K),
        dx_(data.dx)
    {
      if (data.conductors.size() != K_)
      {
        throw std::invalid_argument("conductors must have K values");
      }
      for (const IdxT conductor : data.conductors)
      {
        c_.push_back(conductor - 1);
        if (c_.back() >= N_)
        {
          throw std::invalid_argument("conductors must name phases 1 to N");
        }
      }
      this->setLayout({K_, K_, K_}, {N_, N_});

      std::vector<IdxT> v1(K_);
      std::vector<IdxT> v2(K_);
      for (IdxT k = 0; k < K_; ++k)
      {
        v1[k] = this->local(LineLumpedPorts::bus1, c_[k]);
        v2[k] = this->local(LineLumpedPorts::bus2, c_[k]);
      }

      const IdxT I12  = this->local(LineLumpedInternalVariables::i12);
      const IdxT ISH1 = this->local(LineLumpedInternalVariables::i_sh1);
      const IdxT ISH2 = this->local(LineLumpedInternalVariables::i_sh2);

      auto& z  = this->template addSubmodel<VectorFitT>(LineLumpedSubmodels::z, data.Zp, K_, K_, dx_);
      auto& y1 = this->template addSubmodel<VectorFitT>(LineLumpedSubmodels::y1, data.Yp, K_, K_, dx_);
      auto& y2 = this->template addSubmodel<VectorFitT>(LineLumpedSubmodels::y2, data.Yp, K_, K_, dx_);
      if (z.connect(VectorFitPorts::input, *this, I12) != 0
          || z.connect(VectorFitPorts::out, *this, I12) != 0
          || y1.connect(VectorFitPorts::input, *this, v1) != 0
          || y1.connect(VectorFitPorts::out, *this, ISH1) != 0
          || y2.connect(VectorFitPorts::input, *this, v2) != 0
          || y2.connect(VectorFitPorts::out, *this, ISH2) != 0)
      {
        throw std::logic_error("LineLumped: cannot wire the impedance and admittances");
      }
    }

    template <typename scalar_type, typename index_type>
    typename LineLumped<scalar_type, index_type>::VectorFitT& LineLumped<scalar_type, index_type>::submodel(SubmodelsT symbol)
    {
      return ComponentT::template submodel<VectorFitT>(symbol);
    }

    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks LineLumped<scalar_type, index_type>::verify() const
    {
      std::vector<bool> covered(N_, false);
      for (const IdxT phase : c_)
      {
        covered[phase] = true;
      }

      Model::ConfigurationChecks checks("LineLumped");
      checks.check(N_ > 0 && K_ > 0, "N and K must be positive");
      checks.check(std::ranges::all_of(covered, [](bool phase)
                                       { return phase; }),
                   "every phase must have a conductor");
      checks.check(dx_ > ZERO<RealT>, "dx must be positive");
      return checks;
    }

    /// Currents at the system power base over the base voltage of each conductor's phase at terminal 1
    template <typename scalar_type, typename index_type>
    int LineLumped<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      const auto i12   = this->absoluteTolerance(LineLumpedInternalVariables::i12);
      const auto i_sh1 = this->absoluteTolerance(LineLumpedInternalVariables::i_sh1);
      const auto i_sh2 = this->absoluteTolerance(LineLumpedInternalVariables::i_sh2);
      for (IdxT k = 0; k < K_; ++k)
      {
        const RealT v_base  = this->absoluteTolerance(LineLumpedPorts::bus1, c_[k]) / rel_tol;
        const RealT current = rel_tol * this->va_system_base_ / v_base;
        i12[k]              = current;
        i_sh1[k]            = current;
        i_sh2[k]            = current;
      }
      return 0;
    }

    template <typename scalar_type, typename index_type>
    template <typename T>
    void LineLumped<scalar_type, index_type>::internalResidual(const T* y, const T*, const T* w, const T*, T* f) const
    {
      const IdxT I12  = this->local(LineLumpedInternalVariables::i12);
      const IdxT ISH1 = this->local(LineLumpedInternalVariables::i_sh1);
      const IdxT ISH2 = this->local(LineLumpedInternalVariables::i_sh2);
      const IdxT BUS1 = this->external(LineLumpedPorts::bus1);
      const IdxT BUS2 = this->external(LineLumpedPorts::bus2);
      for (IdxT k = 0; k < K_; ++k)
      {
        f[I12 + k]  = w[BUS2 + c_[k]] - w[BUS1 + c_[k]];
        f[ISH1 + k] = TWO<RealT> * y[ISH1 + k];
        f[ISH2 + k] = TWO<RealT> * y[ISH2 + k];
      }
    }

    template <typename scalar_type, typename index_type>
    template <typename T>
    void LineLumped<scalar_type, index_type>::externalResidual(const T* y, const T*, const T*, const T*, T* f) const
    {
      const IdxT I12  = this->local(LineLumpedInternalVariables::i12);
      const IdxT ISH1 = this->local(LineLumpedInternalVariables::i_sh1);
      const IdxT ISH2 = this->local(LineLumpedInternalVariables::i_sh2);
      const IdxT BUS1 = this->external(LineLumpedPorts::bus1);
      const IdxT BUS2 = this->external(LineLumpedPorts::bus2);
      for (IdxT k = 0; k < K_; ++k)
      {
        f[BUS1 + c_[k]] += y[ISH1 + k] - y[I12 + k];
        f[BUS2 + c_[k]] += y[ISH2 + k] + y[I12 + k];
      }
    }
  } // namespace EMT
} // namespace GridKit
