/**
 * @file VectorFitImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of an EMT VectorFit operator.
 *
 * The model equations are in the `README.md` of this directory, and the
 * kernel arguments are defined in `ComponentImpl.hpp`. Its kernels name the
 * external values and derivatives `u` and `up`.
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
    /**
     * @brief Construct the operator that the data describes.
     *
     * The output dimension `N` and the input dimension `K` are the rows and
     * columns of the first of `D`, `E` and the residues that the data gives.
     * Layout: internal variables `w` and `v`, input `u` of size `K`, output
     * `y` of size `N`. `verify` reports invalid data.
     *
     * @param[in] data Coefficients, poles and residues.
     */
    template <typename scalar_type, typename index_type>
    VectorFit<scalar_type, index_type>::VectorFit(const ModelDataT& data)
    {
      const auto shape = [this](const auto& matrix)
      {
        if (N_ == 0 && !matrix.empty())
        {
          N_ = static_cast<IdxT>(matrix.size());
          K_ = static_cast<IdxT>(matrix.front().size());
        }
      };
      shape(data.D);
      shape(data.E);
      shape(data.residues);

      D_ = coefficients(data.D);
      E_ = coefficients(data.E);
      setDerivedParameters(data);

      const auto Qr = static_cast<IdxT>(ar_.size());
      const auto Qc = static_cast<IdxT>(ac_.size());
      this->setLayout({K_ * (Qr + Qc), K_ * Qc}, {}, {K_}, {N_});
    }

    /**
     * @brief Check the sizes and shapes of the data and the pairing of
     *        complex poles.
     *
     * @return The failed checks.
     */
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

    /**
     * @brief Set the tolerance of each memory state to the tolerance of its
     *        input over the pole magnitude.
     *
     * The argument is unused. A pole at zero keeps the input tolerance.
     *
     * @return 0.
     *
     * @pre The component at input `u` has set its tolerances.
     */
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
        const RealT input = this->externalTolerance(this->inputOffset(VectorFitInputs::u) + k);
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

    /**
     * @brief Write the residuals of the memory states `w` and `v`.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    void VectorFit<scalar_type, index_type>::internalResidual(const T* y, const T* yp, const T* u, const T*, T* f) const
    {
      const auto Qr    = static_cast<IdxT>(ar_.size());
      const auto Qc    = static_cast<IdxT>(ac_.size());
      const IdxT W     = this->local(VectorFitInternalVariables::w);
      const IdxT V     = this->local(VectorFitInternalVariables::v);
      const IdxT INPUT = this->inputOffset(VectorFitInputs::u);

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

    /**
     * @brief Write the output `y`.
     *
     * @pre The entries of `y` are zero on entry.
     */
    template <typename scalar_type, typename index_type>
    template <typename T>
    void VectorFit<scalar_type, index_type>::output(OutputsT, const T* y, const T*, const T* u, const T* up, T* f) const
    {
      const auto Qr    = static_cast<IdxT>(ar_.size());
      const auto Qc    = static_cast<IdxT>(ac_.size());
      const IdxT W     = this->local(VectorFitInternalVariables::w);
      const IdxT V     = this->local(VectorFitInternalVariables::v);
      const IdxT INPUT = this->inputOffset(VectorFitInputs::u);
      const IdxT NK    = N_ * K_;

      for (IdxT n = 0; n < N_; ++n)
      {
        for (IdxT k = 0; k < K_; ++k)
        {
          const IdxT nk  = n * K_ + k;
          f[n]          += D_[nk] * u[INPUT + k] + E_[nk] * up[INPUT + k];
          for (IdxT q = 0; q < Qr; ++q)
          {
            f[n] += Ar_[q * NK + nk] * y[W + q * K_ + k];
          }
          for (IdxT q = 0; q < Qc; ++q)
          {
            f[n] += TWO<RealT> * (Ac_[q * NK + nk] * y[W + (Qr + q) * K_ + k] - Bc_[q * NK + nk] * y[V + q * K_ + k]);
          }
        }
      }
    }

    /**
     * @brief Row-major copy of `matrix`.
     *
     * @param[in] matrix Coefficient of `N` rows and `K` columns, empty if omitted.
     * @return `N * K` values, all zero if `matrix` is empty, or no values if
     *         `matrix` is misshaped.
     *
     * @pre `N_` and `K_` are set.
     */
    template <typename scalar_type, typename index_type>
    std::vector<typename VectorFit<scalar_type, index_type>::RealT>
    VectorFit<scalar_type, index_type>::coefficients(const Matrix<RealT>& matrix) const
    {
      std::vector<RealT> values;
      if (matrix.empty())
      {
        values.assign(N_ * K_, ZERO<RealT>);
        return values;
      }

      const bool shaped = matrix.size() == N_ && std::ranges::all_of(matrix, [this](const auto& row)
                                                                     { return row.size() == K_; });
      if (!shaped)
      {
        return values;
      }

      for (const auto& row : matrix)
      {
        for (const RealT value : row)
        {
          values.push_back(value);
        }
      }
      return values;
    }

    /**
     * @brief Store the real poles and the first pole of each conjugate pair,
     *        with their residues.
     *
     * @param[in] data Poles and residues.
     *
     * @pre `N_` and `K_` are set.
     * @post The flags that `verify` reports are set. Misshaped residues store
     *       no poles.
     */
    template <typename scalar_type, typename index_type>
    void VectorFit<scalar_type, index_type>::setDerivedParameters(const ModelDataT& data)
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
            Ar_.push_back(residue[0]);
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
          Ac_.push_back(residue[0]);
          Bc_.push_back(residue[1]);
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
