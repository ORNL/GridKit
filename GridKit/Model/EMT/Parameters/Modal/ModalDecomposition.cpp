/**
 * @file ModalDecomposition.cpp
 *
 * @brief Per-sample spectral factorization of Gamma with identity
 *        assignment, cluster alignment, and exact biorthogonal duals.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#include "ModalDecomposition.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace
{
  /// Guards relative comparisons against vanishing magnitudes.
  constexpr double MAGNITUDE_FLOOR = 1.0e-300;

  /**
   * @brief Greedy maximal assignment: labels claim raw columns in
   *        descending overlap order until every label owns one column.
   */
  std::vector<Eigen::Index> assignByOverlap(Eigen::MatrixXd overlap)
  {
    const auto                count = overlap.rows();
    std::vector<Eigen::Index> raw_for_label(static_cast<size_t>(count), 0);

    for (Eigen::Index pass = 0; pass < count; ++pass)
    {
      Eigen::Index best_label = 0;
      Eigen::Index best_raw   = 0;
      double       best       = -1.0;
      for (Eigen::Index label = 0; label < count; ++label)
      {
        for (Eigen::Index raw = 0; raw < count; ++raw)
        {
          if (overlap(label, raw) > best)
          {
            best       = overlap(label, raw);
            best_label = label;
            best_raw   = raw;
          }
        }
      }
      raw_for_label[static_cast<size_t>(best_label)] = best_raw;
      overlap.row(best_label).setConstant(-1.0);
      overlap.col(best_raw).setConstant(-1.0);
    }
    return raw_for_label;
  }

  /**
   * @brief Clusters of labels whose eigenvalues sit within the
   *        relative gap of one another; singletons included, so one
   *        alignment law covers every mode.
   */
  std::vector<std::vector<Eigen::Index>>
  clustersByGap(const Eigen::VectorXcd& lambda, double gap)
  {
    const auto                count = lambda.size();
    std::vector<Eigen::Index> parent(static_cast<size_t>(count));
    std::iota(parent.begin(), parent.end(), Eigen::Index{0});

    const auto find = [&parent](Eigen::Index node)
    {
      while (parent[static_cast<size_t>(node)] != node)
      {
        node = parent[static_cast<size_t>(node)];
      }
      return node;
    };

    for (Eigen::Index i = 0; i < count; ++i)
    {
      for (Eigen::Index j = i + 1; j < count; ++j)
      {
        const double distance = std::abs(lambda(i) - lambda(j));
        const double scale =
            std::abs(lambda(i)) + std::abs(lambda(j)) + MAGNITUDE_FLOOR;
        if (distance <= gap * scale)
        {
          parent[static_cast<size_t>(find(j))] = find(i);
        }
      }
    }

    std::vector<std::vector<Eigen::Index>> clusters;
    std::vector<Eigen::Index>              slot_of_root(static_cast<size_t>(count), -1);
    for (Eigen::Index i = 0; i < count; ++i)
    {
      auto& slot = slot_of_root[static_cast<size_t>(find(i))];
      if (slot < 0)
      {
        slot = static_cast<Eigen::Index>(clusters.size());
        clusters.emplace_back();
      }
      clusters[static_cast<size_t>(slot)].push_back(i);
    }
    return clusters;
  }

  /**
   * @brief Polar factor of a square matrix: the unitary closest to it.
   */
  Eigen::MatrixXcd polarFactor(const Eigen::MatrixXcd& matrix)
  {
    const Eigen::JacobiSVD<Eigen::MatrixXcd> svd(
        matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
    return svd.matrixU() * svd.matrixV().adjoint();
  }
} // namespace

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename scalar_type, typename index_type>
      ModalDecomposition<scalar_type, index_type>::ModalDecomposition(
          IdxT       conductor_count,
          RealT      line_length,
          Parameters parameters)
        : parameters_(parameters),
          layout_(conductor_count),
          length_(line_length),
          values_(static_cast<size_t>(layout_.n), RealT{0})
      {
        if (conductor_count < 1 || !(line_length > RealT{0})
            || !(parameters_.cluster_gap >= RealT{0}))
        {
          throw std::runtime_error(
              "ModalDecomposition requires a positive mode count, line "
              "length, and cluster gap");
        }
      }

      template <typename scalar_type, typename index_type>
      void ModalDecomposition<scalar_type, index_type>::reset()
      {
        tv_previous_.clear();
        ti_previous_.clear();
        lambda_previous_.clear();
      }

      template <typename scalar_type, typename index_type>
      int ModalDecomposition<scalar_type, index_type>::decompose(
          RealT omega, const ScalarT* alpha, const ScalarT* beta)
      {
        if (!(omega > RealT{0}))
        {
          return -1;
        }

        const auto k = static_cast<Eigen::Index>(layout_.K);

        Eigen::MatrixXcd gamma(k, k);
        for (Eigen::Index i = 0; i < k; ++i)
        {
          for (Eigen::Index j = 0; j < k; ++j)
          {
            const auto entry = static_cast<size_t>(i * k + j);
            gamma(i, j)      = ComplexT{alpha[entry], beta[entry]};
          }
        }

        const Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(gamma,
                                                                 true);
        if (solver.info() != Eigen::Success)
        {
          return -1;
        }
        Eigen::MatrixXcd v      = solver.eigenvectors();
        Eigen::VectorXcd lambda = solver.eigenvalues();

        if (tv_previous_.empty())
        {
          // Canonical first frame: modes ascend by (Im, Re) and each
          // column's largest entry is made real positive.
          std::vector<Eigen::Index> order(static_cast<size_t>(k));
          std::iota(order.begin(), order.end(), Eigen::Index{0});
          std::sort(order.begin(),
                    order.end(),
                    [&lambda](Eigen::Index lhs, Eigen::Index rhs)
                    {
                      const auto l = lambda(lhs);
                      const auto r = lambda(rhs);
                      return (l.imag() < r.imag())
                             || (l.imag() == r.imag()
                                 && l.real() < r.real());
                    });

          Eigen::MatrixXcd ordered_v(k, k);
          Eigen::VectorXcd ordered_lambda(k);
          for (Eigen::Index m = 0; m < k; ++m)
          {
            ordered_v.col(m)  = v.col(order[static_cast<size_t>(m)]);
            ordered_lambda(m) = lambda(order[static_cast<size_t>(m)]);
          }
          for (Eigen::Index m = 0; m < k; ++m)
          {
            Eigen::Index pivot = 0;
            ordered_v.col(m).cwiseAbs().maxCoeff(&pivot);
            const double norm = std::abs(ordered_v(pivot, m));
            if (norm > MAGNITUDE_FLOOR)
            {
              ordered_v.col(m) *= std::conj(ordered_v(pivot, m)) / norm;
            }
          }
          v      = ordered_v;
          lambda = ordered_lambda;
        }
        else
        {
          const Eigen::Map<const Eigen::MatrixXcd> tv_previous(
              tv_previous_.data(), k, k);
          const Eigen::Map<const Eigen::MatrixXcd> ti_previous(
              ti_previous_.data(), k, k);
          const Eigen::Map<const Eigen::VectorXcd> lambda_previous(
              lambda_previous_.data(), k);

          // Identity follows the eigenvectors: the previous duals
          // score every raw vector and each label claims its best
          // match, so eigenvalue crossings cannot swap modes.
          const auto raw_for_label =
              assignByOverlap((ti_previous.adjoint() * v).cwiseAbs());

          Eigen::MatrixXcd assigned_v(k, k);
          Eigen::VectorXcd assigned_lambda(k);
          for (Eigen::Index label = 0; label < k; ++label)
          {
            const auto raw         = raw_for_label[static_cast<size_t>(label)];
            assigned_v.col(label)  = v.col(raw);
            assigned_lambda(label) = lambda(raw);
          }

          // Align every cluster to the previous frame by the polar
          // factor of its overlap, the closest admissible frame. The
          // orthonormal span stays accurate where the member vectors
          // are not individually defined.
          for (const auto& cluster :
               clustersByGap(assigned_lambda, parameters_.cluster_gap))
          {
            const auto size = static_cast<Eigen::Index>(cluster.size());

            Eigen::MatrixXcd members(k, size);
            Eigen::MatrixXcd reference(k, size);
            for (Eigen::Index c = 0; c < size; ++c)
            {
              const auto label = cluster[static_cast<size_t>(c)];
              members.col(c)   = assigned_v.col(label);
              reference.col(c) = tv_previous.col(label);
            }

            const Eigen::HouseholderQR<Eigen::MatrixXcd> qr(members);
            const Eigen::MatrixXcd                       span =
                qr.householderQ() * Eigen::MatrixXcd::Identity(k, size);
            const Eigen::MatrixXcd aligned =
                span * polarFactor(span.adjoint() * reference);

            // Eigenvalues keep their own trajectories through the
            // cluster: match each label to its nearest predecessor.
            std::vector<bool> taken(cluster.size(), false);
            Eigen::VectorXcd  matched(size);
            for (Eigen::Index c = 0; c < size; ++c)
            {
              const auto   label   = cluster[static_cast<size_t>(c)];
              Eigen::Index winner  = 0;
              double       nearest = std::numeric_limits<double>::infinity();
              for (Eigen::Index d = 0; d < size; ++d)
              {
                if (taken[static_cast<size_t>(d)])
                {
                  continue;
                }
                const double distance = std::abs(
                    assigned_lambda(cluster[static_cast<size_t>(d)])
                    - lambda_previous(label));
                if (distance < nearest)
                {
                  nearest = distance;
                  winner  = d;
                }
              }
              taken[static_cast<size_t>(winner)] = true;
              matched(c) =
                  assigned_lambda(cluster[static_cast<size_t>(winner)]);
            }

            for (Eigen::Index c = 0; c < size; ++c)
            {
              const auto label       = cluster[static_cast<size_t>(c)];
              assigned_v.col(label)  = aligned.col(c);
              assigned_lambda(label) = matched(c);
            }
          }

          v      = assigned_v;
          lambda = assigned_lambda;
        }

        // Duals computed fresh, never transported: Ti^H Tv = I to
        // machine precision at every sample.
        const Eigen::MatrixXcd ti = v.inverse().adjoint();

        const auto tv_r = static_cast<size_t>(layout_.tv_r);
        const auto tv_i = static_cast<size_t>(layout_.tv_i);
        const auto ti_r = static_cast<size_t>(layout_.ti_r);
        const auto ti_i = static_cast<size_t>(layout_.ti_i);
        for (Eigen::Index i = 0; i < k; ++i)
        {
          for (Eigen::Index j = 0; j < k; ++j)
          {
            const auto entry      = static_cast<size_t>(i * k + j);
            values_[tv_r + entry] = v(i, j).real();
            values_[tv_i + entry] = v(i, j).imag();
            values_[ti_r + entry] = ti(i, j).real();
            values_[ti_i + entry] = ti(i, j).imag();
          }
        }

        const auto a   = static_cast<size_t>(layout_.a);
        const auto b   = static_cast<size_t>(layout_.b);
        const auto h_r = static_cast<size_t>(layout_.h_r);
        const auto h_i = static_cast<size_t>(layout_.h_i);
        const auto tau = static_cast<size_t>(layout_.tau);
        for (Eigen::Index m = 0; m < k; ++m)
        {
          const auto     entry  = static_cast<size_t>(m);
          const ComplexT value  = lambda(m);
          const ComplexT factor = std::exp(-length_ * value);

          values_[a + entry]   = value.real();
          values_[b + entry]   = value.imag();
          values_[h_r + entry] = factor.real();
          values_[h_i + entry] = factor.imag();
          values_[tau + entry] = length_ * value.imag() / omega;
        }

        tv_previous_.assign(v.data(), v.data() + k * k);
        ti_previous_.assign(ti.data(), ti.data() + k * k);
        lambda_previous_.assign(lambda.data(), lambda.data() + k);
        return 0;
      }

      template class ModalDecomposition<double, long int>;
      template class ModalDecomposition<double, size_t>;
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
