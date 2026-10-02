// /**
//  * @file Interpolation.cpp
//  * @brief Implementation of barycentric Lagrange polynomial interpolation.
//  *
//  * This file implements the LagrangeInterpolant class, which evaluates
//  * polynomial interpolants over a time interval using either Chebyshev or
//  * equidistant interpolation nodes. More interpolation nodes can be added as
//  * need. The interpolation is performed on the reference interval [-1, 1], while
//  * physical times are mapped to the reference interval using an affine
//  * transformation.
//  * @author Abdourahman Barry
//  */
// #pragma once
// #include <cassert>
// #include <cmath>
// #include <cstddef>
// #include <memory>
// #include <numeric>
// #include <vector>

// #include <GridKit/LinearAlgebra/Vector/Vector.hpp>
// #include <GridKit/LinearAlgebra/Vector/VectorHandler.hpp>
// #include <GridKit/MemoryUtilities/MemoryUtils.hpp>
// #include <GridKit/ScalarTraits.hpp>

// namespace GridKit
// {

//   enum class NodeType
//   {
//     Chebyshev,
//     Equidistant
//   };

//   template <class ScalarT, typename IdxT>
//   class LagrangeInterpolant
//   {

//   public:
//     using Vector      = GridKit::LinearAlgebra::Vector<ScalarT, IdxT>;
//     using MemorySpace = GridKit::memory::MemorySpace;
//     using Handler     = GridKit::LinearAlgebra::VectorHandler<ScalarT, IdxT>;
//     using Matrix      = std::vector<std::shared_ptr<Vector>>;
//     using RealT       = typename GridKit::ScalarTraits<ScalarT>::RealT;

//     /**
//      * @brief Constructs a Lagrange interpolant from a collection of node values.
//      *
//      * The constructor stores a reference to the provided node values, sets the
//      * physical time interval, creates interpolation nodes on the reference interval
//      * [-1, 1], and computes the corresponding barycentric weights.
//      *
//      * @param node_values Values at the interpolation nodes.
//      * @param t0 Start time of the physical interval.
//      * @param dt Length of the physical time interval.
//      * @param node_type Type of interpolation nodes to use.
//      * @param memspace Memory space used for vector allocation and operations.
//      */
//     LagrangeInterpolant(Matrix&     node_values,
//                         RealT       t0,
//                         RealT       dt,
//                         NodeType    node_type = NodeType::Chebyshev,
//                         MemorySpace memspace  = MemorySpace::HOST)
//       : node_values_(node_values),
//         node_type_(node_type),
//         memspace_(memspace),
//         t0_(t0),
//         dt_(dt)
//     {
//       y_ = std::make_unique<Vector>(node_values_[0]->getSize());
//       y_->allocate(memspace_);
//       y_->setToZero(memspace_);

//       setTimeInterval(t0, dt);
//       setNodes(node_type);
//     }

//     /**
//      * @brief Sets the interpolation nodes and computes barycentric weights.
//      *
//      * The number of interpolation nodes is determined from the number of entries
//      * in node_values_. The nodes are created on the reference interval [-1, 1].
//      *
//      * @param node_type Type of nodes to generate.
//      *
//      * @throws std::invalid_argument if node_values_ is empty.
//      */
//     void setNodes(NodeType node_type)
//     {
//       if (node_values_.size() <= 0)
//       {
//         throw std::invalid_argument("node_values must have at least one column.");
//       }

//       node_type_           = node_type;
//       node_times_          = makeNodes(node_values_.size(), node_type_);
//       barycentric_weights_ = computeBarycentricWeights(node_times_);
//     }

//     /**
//      * @brief Sets the physical time interval for interpolation.
//      *
//      * The physical interval is [t0, t0 + dt]. Evaluation times are mapped from
//      * this interval to the reference interval [-1, 1].
//      *
//      * @param t0 Start time of the interval.
//      * @param dt Length of the interval.
//      *
//      * @throws std::invalid_argument if dt is negative.
//      */
//     void setTimeInterval(RealT t0, RealT dt)
//     {
//       if (dt < 0.0)
//       {
//         throw std::invalid_argument("dt must be non-negative");
//       }

//       t0_ = t0;
//       dt_ = dt;
//     }

//     /**
//      * @brief Evaluates the interpolant at a physical time t.
//      *
//      * The time t is first mapped to the reference interval [-1, 1]. Then the
//      * barycentric Lagrange basis values are computed and used to form a weighted
//      * linear combination of the node values.
//      *
//      * @param t Physical time at which to evaluate the interpolant.
//      *
//      * @return Pointer to the interpolated vector value.
//      *
//      * @note If there is only one node value, this function returns the existing
//      * node value pointer directly. Otherwise, it allocates a new Vector using
//      * `new`. The caller is responsible for deleting the returned pointer in that
//      * case, unless ownership is managed elsewhere.
//      *
//      * @warning Returning a raw pointer allocated with `new` can easily cause memory
//      * ownership issues. Consider returning `std::unique_ptr<Vector>` instead.
//      */
//     void evaluate(RealT t, Vector& y)
//     {
//       checkReady();

//       assert(y.getSize() == node_values_[0]->getSize());

//       const size_t n = node_values_.size();

//       if (n == 1)
//       {
//         y.copyFromExternal(*node_values_[0], memspace_, memspace_);
//       }

//       const double              x     = toReferenceTime(t);
//       const std::vector<double> basis = interpolationBasis(x);

//       y.setToConst(0.0, memspace_);

//       for (size_t i = 0; i < n; ++i)
//       {
//         vector_handler_.axpy(
//             basis[i],
//             node_values_[i].get(),
//             &y,
//             memspace_);
//       }

//       y.setDataUpdated(memspace_);
//     }

//     static std::vector<RealT> getPhysicalOutputTimes(size_t num_nodes, RealT t0, RealT dt, NodeType node_type = NodeType::Chebyshev)
//     {

//       if (num_nodes == 1)
//       {
//         return {t0 + dt};
//       }

//       const auto output_nodes = makeNodes(num_nodes, node_type);

//       std::vector<RealT> output_times(num_nodes);

//       for (size_t i = 0; i < num_nodes; ++i)
//       {
//         output_times[i] = t0 + static_cast<RealT>(0.5) * dt * (output_nodes[i] + static_cast<RealT>(1.0));
//       }

//       return output_times;
//     }

//     /**
//      * @brief Returns the stored interpolation node values.
//      *
//      * @return Constant reference to the node value matrix.
//      */
//     const LagrangeInterpolant::Matrix& getNodeValues() const
//     {
//       return node_values_;
//     }

//     /**
//      * @brief Returns the interpolation node locations on [-1, 1].
//      *
//      * @return Constant reference to the vector of node locations.
//      */
//     const std::vector<RealT>& getNodeTimes() const
//     {
//       return node_times_;
//     }

//     /**
//      * @brief Returns the barycentric interpolation weights.
//      *
//      * @return Constant reference to the barycentric weights.
//      */
//     const std::vector<RealT>& getBarycentricWeights() const
//     {
//       return barycentric_weights_;
//     }

//     /**
//      * @brief Returns the start time of the physical interval.
//      *
//      * @return Start time t0.
//      */
//     RealT getStartTime() const
//     {
//       return t0_;
//     }

//     /**
//      * @brief Returns the physical time step length.
//      *
//      * @return Time interval length dt.
//      */
//     RealT getTimeStep() const
//     {
//       return dt_;
//     }

//     /**
//      * @brief Returns the interpolation node type.
//      *
//      * @return Current node type.
//      */
//     NodeType getNodeType() const
//     {
//       return node_type_;
//     }

//     /**
//      * @brief Generates interpolation nodes on the reference interval [-1, 1].
//      *
//      * Supports Chebyshev and equidistant nodes. If only one node is requested, the
//      * single node is placed at 0.
//      *
//      * @param num_nodes Number of interpolation nodes.
//      * @param type Type of nodes to generate.
//      *
//      * @return Vector containing the generated node locations.
//      *
//      * @throws std::invalid_argument if num_nodes is not positive.
//      */
//     static std::vector<RealT> makeNodes(size_t num_nodes, NodeType type)
//     {
//       if (num_nodes <= 0)
//       {
//         throw std::invalid_argument("num_nodes must be positive.");
//       }

//       std::vector<RealT> nodes(num_nodes);

//       if (num_nodes == 1)
//       {
//         nodes[0] = 0.0;
//         return nodes;
//       }

//       switch (type)
//       {
//       case NodeType::Chebyshev:
//         for (size_t k = 0; k < num_nodes; ++k)
//         {
//           const RealT k_real = static_cast<RealT>(num_nodes - k - 1);
//           const RealT n_real = static_cast<RealT>(num_nodes - 1);

//           nodes[k] = std::cos(M_PI * k_real / n_real);
//         }
//         break;

//       case NodeType::Equidistant:
//         for (size_t k = 0; k < num_nodes; ++k)
//         {
//           nodes[k] = static_cast<RealT>(-1.0) + static_cast<RealT>(2.0) * static_cast<RealT>(k) / static_cast<RealT>(num_nodes - 1);
//         }
//         break;
//       }

//       return nodes;
//     }

//   private:
//     /**
//      * @brief Checks whether the interpolator has valid node values.
//      *
//      * @throws std::runtime_error if no node values are available.
//      */
//     void checkReady() const
//     {
//       if (node_values_.size() == 0)
//       {
//         throw std::runtime_error("Interpolator has no node values.");
//       }
//     }

//     /**
//      * @brief Maps a physical time to the reference interval [-1, 1].
//      *
//      * The affine transformation used is:
//      *
//      * x = 2 * (t - t0) / dt - 1
//      *
//      * where [t0, t0 + dt] is the physical interval.
//      *
//      * @param t Physical time.
//      * @return Corresponding reference coordinate.
//      */
//     RealT toReferenceTime(RealT t) const
//     {
//       return 2.0 * (t - t0_) / dt_ - 1.0;
//     }

//     /**
//      * @brief Computes barycentric weights for Lagrange interpolation.
//      *
//      * For each node x_k, the barycentric weight is computed as:
//      *
//      * w_k = 1 / product_{j != k} (x_k - x_j)
//      *
//      * These weights are later used to efficiently evaluate the Lagrange basis.
//      *
//      * @param nodes Interpolation node locations.
//      * @return Vector of barycentric weights.
//      */
//     static std::vector<RealT> computeBarycentricWeights(const std::vector<RealT>& nodes)
//     {
//       size_t             n = nodes.size();
//       std::vector<RealT> weights(n);

//       for (size_t k = 0; k < n; ++k)
//       {
//         RealT w = 1.0;
//         for (size_t j = 0; j < n; ++j)
//         {
//           if (j == k)
//           {
//             continue;
//           }
//           w /= (nodes[k] - nodes[j]);
//         }
//         weights[k] = w;
//       }

//       return weights;
//     }

//     /**
//      * @brief Evaluates all Lagrange basis functions at a reference coordinate.
//      *
//      * This function uses the barycentric interpolation formula. If the evaluation
//      * point coincides with one of the interpolation nodes, the corresponding basis
//      * value is set to 1 and all others are set to 0.
//      *
//      * @param x Reference coordinate where the basis should be evaluated.
//      * @return Vector containing the value of each basis function at x.
//      */
//     std::vector<RealT> interpolationBasis(RealT x) const
//     {
//       size_t             n = node_times_.size();
//       std::vector<RealT> result(n);

//       for (size_t k = 0; k < n; ++k)
//       {
//         RealT diff = x - node_times_[k];

//         if (std::abs(diff) < 1e-14)
//         {
//           std::fill(result.begin(), result.end(), 0.0);
//           result[k] = 1.0;
//           return result;
//         }

//         result[k] = barycentric_weights_[k] / diff;
//       }

//       RealT s = std::accumulate(result.begin(), result.end(), 0.0);

//       for (auto& val : result)
//       {
//         val *= (1.0 / s);
//       }

//       return result;
//     }

//     Matrix&            node_values_;
//     std::vector<RealT> node_times_;
//     std::vector<RealT> barycentric_weights_;

//     NodeType    node_type_;
//     MemorySpace memspace_;
//     Handler     vector_handler_;

//     std::unique_ptr<Vector> y_;

//     RealT t0_;
//     RealT dt_;
//   };

// } // namespace GridKit

/**
 * @file Interpolation.hpp
 * @brief Implementation of barycentric Lagrange polynomial interpolation.
 *
 * This file implements the LagrangeInterpolant class, which evaluates
 * polynomial interpolants over a time interval using either Chebyshev or
 * equidistant interpolation nodes. More interpolation nodes can be added as
 * need. The interpolation is performed on the reference interval [-1, 1], while
 * physical times are mapped to the reference interval using an affine
 * transformation.
 * @author Abdourahman Barry
 */
#pragma once
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <GridKit/LinearAlgebra/Vector/Vector.hpp>
#include <GridKit/LinearAlgebra/Vector/VectorHandler.hpp>
#include <GridKit/MemoryUtilities/MemoryUtils.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{

  enum class NodeType
  {
    Chebyshev,
    Equidistant
  };

  template <class ScalarT, typename IdxT>
  class LagrangeInterpolant
  {

  public:
    using Vector      = GridKit::LinearAlgebra::Vector<ScalarT, IdxT>;
    using MemorySpace = GridKit::memory::MemorySpace;
    using Handler     = GridKit::LinearAlgebra::VectorHandler<ScalarT, IdxT>;
    using Matrix      = std::vector<std::shared_ptr<Vector>>;
    using RealT       = typename GridKit::ScalarTraits<ScalarT>::RealT;

    /**
     * @brief Constructs a Lagrange interpolant from a collection of node values.
     *
     * The constructor copies the list of node-value pointers (the underlying
     * vectors are shared, not copied), sets the physical time interval, creates
     * interpolation nodes on the reference interval [-1, 1], and computes the
     * corresponding barycentric weights.
     *
     * Because the pointer list is owned by the interpolant, the caller's
     * `node_values` container may be a temporary. Updates written into the
     * shared vectors are seen by the interpolant without rebuilding it.
     *
     * @param node_values Values at the interpolation nodes.
     * @param t0 Start time of the physical interval.
     * @param dt Length of the physical time interval.
     * @param node_type Type of interpolation nodes to use.
     * @param memspace Memory space used for vector allocation and operations.
     */
    LagrangeInterpolant(const Matrix& node_values,
                        RealT         t0,
                        RealT         dt,
                        NodeType      node_type = NodeType::Chebyshev,
                        MemorySpace   memspace  = MemorySpace::HOST)
      : node_values_(node_values),
        node_type_(node_type),
        memspace_(memspace),
        t0_(t0),
        dt_(dt)
    {
      setTimeInterval(t0, dt);
      setNodes(node_type);
    }

    /**
     * @brief Sets the interpolation nodes and computes barycentric weights.
     *
     * The number of interpolation nodes is determined from the number of entries
     * in node_values_. The nodes are created on the reference interval [-1, 1].
     *
     * @param node_type Type of nodes to generate.
     *
     * @throws std::invalid_argument if node_values_ is empty.
     */
    void setNodes(NodeType node_type)
    {
      if (node_values_.empty())
      {
        throw std::invalid_argument("node_values must have at least one column.");
      }

      node_type_           = node_type;
      node_times_          = makeNodes(node_values_.size(), node_type_);
      barycentric_weights_ = computeBarycentricWeights(node_times_);
      basis_.assign(node_values_.size(), static_cast<RealT>(0));
    }

    /**
     * @brief Sets the physical time interval for interpolation.
     *
     * The physical interval is [t0, t0 + dt]. Evaluation times are mapped from
     * this interval to the reference interval [-1, 1].
     *
     * @param t0 Start time of the interval.
     * @param dt Length of the interval.
     *
     * @throws std::invalid_argument if dt is negative.
     */
    void setTimeInterval(RealT t0, RealT dt)
    {
      if (dt < 0.0)
      {
        throw std::invalid_argument("dt must be non-negative");
      }

      t0_ = t0;
      dt_ = dt;
    }

    /**
     * @brief Evaluates the interpolant at a physical time t into y.
     *
     * The time t is first mapped to the reference interval [-1, 1]. Then the
     * barycentric Lagrange basis values are computed and used to form a weighted
     * linear combination of the node values.
     *
     * @param t Physical time at which to evaluate the interpolant.
     * @param y Output vector; must have the same size as the node values.
     */
    void evaluate(RealT t, Vector& y)
    {
      checkReady();

      assert(y.getSize() == node_values_[0]->getSize());

      const size_t n = node_values_.size();

      if (n == 1)
      {
        y.copyFromExternal(*node_values_[0], memspace_, memspace_);
        return;
      }

      const RealT x = toReferenceTime(t);
      interpolationBasis(x, basis_);

      y.setToConst(0.0, memspace_);

      for (size_t i = 0; i < n; ++i)
      {
        vector_handler_.axpy(
            basis_[i],
            node_values_[i].get(),
            &y,
            memspace_);
      }

      y.setDataUpdated(memspace_);
    }

    static std::vector<RealT> getPhysicalOutputTimes(size_t num_nodes, RealT t0, RealT dt, NodeType node_type = NodeType::Chebyshev)
    {
      if (num_nodes == 1)
      {
        return {t0 + dt};
      }

      const auto output_nodes = makeNodes(num_nodes, node_type);

      std::vector<RealT> output_times(num_nodes);

      for (size_t i = 0; i < num_nodes; ++i)
      {
        output_times[i] = t0 + static_cast<RealT>(0.5) * dt * (output_nodes[i] + static_cast<RealT>(1.0));
      }

      return output_times;
    }

    /**
     * @brief Returns the stored interpolation node values.
     *
     * @return Constant reference to the node value matrix.
     */
    const Matrix& getNodeValues() const
    {
      return node_values_;
    }

    /**
     * @brief Returns the interpolation node locations on [-1, 1].
     *
     * @return Constant reference to the vector of node locations.
     */
    const std::vector<RealT>& getNodeTimes() const
    {
      return node_times_;
    }

    /**
     * @brief Returns the barycentric interpolation weights.
     *
     * @return Constant reference to the barycentric weights.
     */
    const std::vector<RealT>& getBarycentricWeights() const
    {
      return barycentric_weights_;
    }

    /**
     * @brief Returns the start time of the physical interval.
     *
     * @return Start time t0.
     */
    RealT getStartTime() const
    {
      return t0_;
    }

    /**
     * @brief Returns the physical time step length.
     *
     * @return Time interval length dt.
     */
    RealT getTimeStep() const
    {
      return dt_;
    }

    /**
     * @brief Returns the interpolation node type.
     *
     * @return Current node type.
     */
    NodeType getNodeType() const
    {
      return node_type_;
    }

    /**
     * @brief Generates interpolation nodes on the reference interval [-1, 1].
     *
     * Supports Chebyshev and equidistant nodes. If only one node is requested, the
     * single node is placed at 0.
     *
     * @param num_nodes Number of interpolation nodes.
     * @param type Type of nodes to generate.
     *
     * @return Vector containing the generated node locations.
     *
     * @throws std::invalid_argument if num_nodes is not positive.
     */
    static std::vector<RealT> makeNodes(size_t num_nodes, NodeType type)
    {
      if (num_nodes == 0)
      {
        throw std::invalid_argument("num_nodes must be positive.");
      }

      std::vector<RealT> nodes(num_nodes);

      if (num_nodes == 1)
      {
        nodes[0] = 0.0;
        return nodes;
      }

      switch (type)
      {
      case NodeType::Chebyshev:
        for (size_t k = 0; k < num_nodes; ++k)
        {
          const RealT k_real = static_cast<RealT>(num_nodes - k - 1);
          const RealT n_real = static_cast<RealT>(num_nodes - 1);

          nodes[k] = std::cos(M_PI * k_real / n_real);
        }
        break;

      case NodeType::Equidistant:
        for (size_t k = 0; k < num_nodes; ++k)
        {
          nodes[k] = static_cast<RealT>(-1.0) + static_cast<RealT>(2.0) * static_cast<RealT>(k) / static_cast<RealT>(num_nodes - 1);
        }
        break;
      }

      return nodes;
    }

  private:
    /**
     * @brief Checks whether the interpolator has valid node values.
     *
     * @throws std::runtime_error if no node values are available.
     */
    void checkReady() const
    {
      if (node_values_.empty())
      {
        throw std::runtime_error("Interpolator has no node values.");
      }
    }

    /**
     * @brief Maps a physical time to the reference interval [-1, 1].
     *
     * x = 2 * (t - t0) / dt - 1, where [t0, t0 + dt] is the physical interval.
     */
    RealT toReferenceTime(RealT t) const
    {
      return static_cast<RealT>(2.0) * (t - t0_) / dt_ - static_cast<RealT>(1.0);
    }

    /**
     * @brief Computes barycentric weights for Lagrange interpolation.
     *
     * w_k = 1 / product_{j != k} (x_k - x_j)
     */
    static std::vector<RealT> computeBarycentricWeights(const std::vector<RealT>& nodes)
    {
      const size_t       n = nodes.size();
      std::vector<RealT> weights(n);

      for (size_t k = 0; k < n; ++k)
      {
        RealT w = 1.0;
        for (size_t j = 0; j < n; ++j)
        {
          if (j == k)
          {
            continue;
          }
          w /= (nodes[k] - nodes[j]);
        }
        weights[k] = w;
      }

      return weights;
    }

    /**
     * @brief Evaluates all Lagrange basis functions at a reference coordinate.
     *
     * Uses the barycentric interpolation formula. If the evaluation point
     * coincides with a node, that basis value is 1 and all others are 0.
     * Writes into the caller-provided buffer, which must have size n.
     */
    void interpolationBasis(RealT x, std::vector<RealT>& result) const
    {
      const size_t n = node_times_.size();
      assert(result.size() == n);

      for (size_t k = 0; k < n; ++k)
      {
        const RealT diff = x - node_times_[k];

        if (std::abs(diff) < static_cast<RealT>(1e-14))
        {
          std::fill(result.begin(), result.end(), static_cast<RealT>(0));
          result[k] = static_cast<RealT>(1);
          return;
        }

        result[k] = barycentric_weights_[k] / diff;
      }

      const RealT s     = std::accumulate(result.begin(), result.end(), static_cast<RealT>(0));
      const RealT inv_s = static_cast<RealT>(1) / s;

      for (auto& val : result)
      {
        val *= inv_s;
      }
    }

    Matrix             node_values_; ///< Owned list of shared node-value vectors.
    std::vector<RealT> node_times_;
    std::vector<RealT> barycentric_weights_;
    std::vector<RealT> basis_; ///< Scratch buffer for basis values; sized in setNodes().

    NodeType    node_type_;
    MemorySpace memspace_;
    Handler     vector_handler_;

    RealT t0_;
    RealT dt_;
  };

} // namespace GridKit