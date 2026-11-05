/**
 * @file VectorFitImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT VectorFit operator.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include <GridKit/Model/EMT/ComponentImpl.hpp>
#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFit.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    VectorFit<scalar_type, index_type>::VectorFit(const ModelDataT& data, IdxT N, IdxT K, RealT scale)
      : N_(N),
        K_(K),
        D_(coefficients(data.D, N, K, scale)),
        E_(coefficients(data.E, N, K, scale))
    {
      setDerivedParameters(data, scale);

      const auto Qr = static_cast<IdxT>(ar_.size());
      const auto Qc = static_cast<IdxT>(ac_.size());
      this->setLayout({K_ * (Qr + Qc), K_ * Qc}, {K_, N_});
    }

    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks VectorFit<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("VectorFit");
      checks.check(N_ > 0 && K_ > 0, "N and K must be positive");
      checks.check(D_.size() == N_ * K_ && E_.size() == N_ * K_, "D and E must be N by K");
      checks.check(residues_shaped_, "residues must be N by K by Q");
      checks.check(real_residues_, "real poles must have real residues");
      checks.check(conjugate_pairs_, "complex poles and their residues must come in conjugate pairs");
      return checks;
    }

    /// Memory states at the tolerance of their input over the pole magnitude
    template <typename scalar_type, typename index_type>
    int VectorFit<scalar_type, index_type>::setAbsoluteTolerance(RealT)
    {
      const auto scale = [](RealT input, RealT magnitude)
      {
        if (magnitude == ZERO<RealT>)
        {
          return input;
        }
        return input / magnitude;
      };

      const auto w  = this->absoluteTolerance(VectorFitInternalVariables::w);
      const auto v  = this->absoluteTolerance(VectorFitInternalVariables::v);
      const auto Qr = static_cast<IdxT>(ar_.size());
      const auto Qc = static_cast<IdxT>(ac_.size());
      for (IdxT k = 0; k < K_; ++k)
      {
        const RealT input = this->absoluteTolerance(VectorFitPorts::input, k);
        for (IdxT q = 0; q < Qr; ++q)
        {
          w[q * K_ + k] = scale(input, std::abs(ar_[q]));
        }
        for (IdxT q = 0; q < Qc; ++q)
        {
          w[(Qr + q) * K_ + k] = scale(input, std::hypot(ac_[q], wc_[q]));
          v[q * K_ + k]        = scale(input, std::hypot(ac_[q], wc_[q]));
        }
      }
      return 0;
    }

    template <typename scalar_type, typename index_type>
    template <typename T>
    void VectorFit<scalar_type, index_type>::internalResidual(const T* y, const T* yp, const T* u, const T*, T* f) const
    {
      const auto Qr    = static_cast<IdxT>(ar_.size());
      const auto Qc    = static_cast<IdxT>(ac_.size());
      const IdxT W     = this->local(VectorFitInternalVariables::w);
      const IdxT V     = this->local(VectorFitInternalVariables::v);
      const IdxT INPUT = this->external(VectorFitPorts::input);

      for (IdxT q = 0; q < Qr; ++q)
      {
        for (IdxT k = 0; k < K_; ++k)
        {
          const IdxT i = W + q * K_ + k;
          f[i]         = -yp[i] + ar_[q] * y[i] + u[INPUT + k];
        }
      }

      for (IdxT q = 0; q < Qc; ++q)
      {
        for (IdxT k = 0; k < K_; ++k)
        {
          const IdxT i = W + (Qr + q) * K_ + k;
          const IdxT j = V + q * K_ + k;
          f[i]         = -yp[i] + ac_[q] * y[i] - wc_[q] * y[j] + u[INPUT + k];
          f[j]         = -yp[j] + wc_[q] * y[i] + ac_[q] * y[j];
        }
      }
    }

    template <typename scalar_type, typename index_type>
    template <typename T>
    void VectorFit<scalar_type, index_type>::externalResidual(const T* y, const T*, const T* u, const T* up, T* f) const
    {
      const auto Qr    = static_cast<IdxT>(ar_.size());
      const auto Qc    = static_cast<IdxT>(ac_.size());
      const IdxT W     = this->local(VectorFitInternalVariables::w);
      const IdxT V     = this->local(VectorFitInternalVariables::v);
      const IdxT INPUT = this->external(VectorFitPorts::input);
      const IdxT OUT   = this->external(VectorFitPorts::out);
      const IdxT NK    = N_ * K_;

      for (IdxT n = 0; n < N_; ++n)
      {
        for (IdxT k = 0; k < K_; ++k)
        {
          const IdxT nk  = n * K_ + k;
          f[OUT + n]    += D_[nk] * u[INPUT + k] + E_[nk] * up[INPUT + k];
          for (IdxT q = 0; q < Qr; ++q)
          {
            f[OUT + n] += Ar_[q * NK + nk] * y[W + q * K_ + k];
          }
          for (IdxT q = 0; q < Qc; ++q)
          {
            f[OUT + n] += TWO<RealT> * (Ac_[q * NK + nk] * y[W + (Qr + q) * K_ + k] - Bc_[q * NK + nk] * y[V + q * K_ + k]);
          }
        }
      }
    }

    /// Row-major and scaled; an omitted coefficient is zero and a misshaped one is empty
    template <typename scalar_type, typename index_type>
    std::vector<typename VectorFit<scalar_type, index_type>::RealT>
    VectorFit<scalar_type, index_type>::coefficients(const Matrix<RealT>& matrix, IdxT N, IdxT K, RealT scale)
    {
      std::vector<RealT> values;
      if (matrix.empty())
      {
        values.assign(N * K, ZERO<RealT>);
        return values;
      }

      const bool shaped = matrix.size() == N && std::ranges::all_of(matrix, [K](const auto& row)
                                                                    { return row.size() == K; });
      if (!shaped)
      {
        return values;
      }

      for (const auto& row : matrix)
      {
        for (const RealT value : row)
        {
          values.push_back(scale * value);
        }
      }
      return values;
    }

    /// Real poles and the first pole of each conjugate pair, with their scaled residues
    template <typename scalar_type, typename index_type>
    void VectorFit<scalar_type, index_type>::setDerivedParameters(const ModelDataT& data, RealT scale)
    {
      const auto Q  = static_cast<IdxT>(data.poles.size());
      const IdxT NK = N_ * K_;

      // Residue (n, k) of pole q at nk Q + q
      std::vector<Pair<RealT>> residues;
      for (const auto& row : data.residues)
      {
        residues_shaped_ = residues_shaped_ && row.size() == K_;
        for (const auto& entry : row)
        {
          residues_shaped_ = residues_shaped_ && entry.size() == Q;
          residues.insert(residues.end(), entry.begin(), entry.end());
        }
      }
      residues_shaped_ = residues_shaped_ && residues.size() == NK * Q;
      if (!residues_shaped_)
      {
        return;
      }

      IdxT q = 0;
      while (q < Q)
      {
        const RealT a = data.poles[q][0];
        const RealT b = data.poles[q][1];
        if (b == ZERO<RealT>)
        {
          ar_.push_back(a);
          for (IdxT nk = 0; nk < NK; ++nk)
          {
            const Pair<RealT>& residue = residues[nk * Q + q];
            Ar_.push_back(scale * residue[0]);
            real_residues_ = real_residues_ && residue[1] == ZERO<RealT>;
          }
          q += 1;
          continue;
        }

        const IdxT partner = q + 1;
        bool       paired  = partner < Q && data.poles[partner][0] == a && data.poles[partner][1] == -b;
        ac_.push_back(a);
        wc_.push_back(b);
        for (IdxT nk = 0; nk < NK; ++nk)
        {
          const Pair<RealT>& residue = residues[nk * Q + q];
          Ac_.push_back(scale * residue[0]);
          Bc_.push_back(scale * residue[1]);
          if (paired)
          {
            const Pair<RealT>& conjugate = residues[nk * Q + partner];
            paired                       = conjugate[0] == residue[0] && conjugate[1] == -residue[1];
          }
        }
        conjugate_pairs_  = conjugate_pairs_ && paired;
        q                += 2;
      }
    }
  } // namespace EMT
} // namespace GridKit
