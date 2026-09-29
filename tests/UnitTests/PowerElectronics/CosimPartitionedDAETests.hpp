#pragma once

/**
 * @file PartitionedDaeComponents.hpp
 *
 * The 4-variable partitioned DAE test problem split into two GridKit
 * components, following the same pattern as HiresComponents.hpp.
 *
 *   x' = y - x^3 + z                  (Component: internal, differential)
 *   y' = -x + y(1 - y^2) - w          (Bus + Component, differential)
 *   0  = sin(x - y) - z               (Component: internal, algebraic)
 *   0  = w + y - z - pi/2             (Bus + Component, algebraic)
 *
 * Residuals are written as F = g(y) - M y', matching the HIRES convention,
 * and Jacobians as dF/dy + alpha * dF/dy'.
 *
 * Variable ownership:
 *   - y and w live on the bus node (external variables). The bus component
 *     supplies every term of their equations that depends only on y and w.
 *   - x and z are internal to PartitionedDaeComponent, which also adds the
 *     x- and z-dependent terms (-x and -z) to the bus equations through its
 *     external residual.
 *
 * Local index layout (same as HIRES): externals first, then internals.
 *   PartitionedDaeComponent: 0 = y, 1 = w, 2 = x, 3 = z
 *   PartitionedDaeBus:       0 = y, 1 = w
 *
 * Initial condition (consistent, so IDACalcIC only has to find y'):
 *   x = pi, y = pi/2, z = 1, w = 1
 *   check: sin(pi - pi/2) - 1 = 0 and 1 + pi/2 - 1 - pi/2 = 0
 *
 * Time span: [0, 5]. Reference solution at t = 5:
 *   x = -1.017881864730135e+00
 *   y = -1.641238691143580e+00
 *   z =  5.837639304756466e-01
 *   w =  3.795798948414124e+00
 */

#include <cmath>
#include <vector>

#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>
#include <GridKit/Model/PowerElectronics/NodeBase.hpp>

namespace GridKit
{
  namespace PartitionedDae
  {
    constexpr double PI = 3.14159265358979323846;
  } // namespace PartitionedDae

  /// Owns x (differential) and z (algebraic); adds -x to the y equation and
  /// -z to the w equation on the bus.
  template <class ScalarT, typename IdxT>
  class PartitionedDaeComponent : public CircuitComponent<ScalarT, IdxT>
  {
    using RealT = typename CircuitComponent<ScalarT, IdxT>::RealT;
    using NodeT = typename PowerElectronics::NodeBase<ScalarT, IdxT>;

    using CircuitComponent<ScalarT, IdxT>::size_;
    using CircuitComponent<ScalarT, IdxT>::nnz_;
    using CircuitComponent<ScalarT, IdxT>::alpha_;
    using CircuitComponent<ScalarT, IdxT>::y_ext_;
    using CircuitComponent<ScalarT, IdxT>::y_int_;
    using CircuitComponent<ScalarT, IdxT>::yp_ext_;
    using CircuitComponent<ScalarT, IdxT>::yp_int_;
    using CircuitComponent<ScalarT, IdxT>::tag_;
    using CircuitComponent<ScalarT, IdxT>::abs_tol_;
    using CircuitComponent<ScalarT, IdxT>::f_ext_;
    using CircuitComponent<ScalarT, IdxT>::f_int_;
    using CircuitComponent<ScalarT, IdxT>::idc_;
    using CircuitComponent<ScalarT, IdxT>::extern_indices_;
    using CircuitComponent<ScalarT, IdxT>::n_extern_;
    using CircuitComponent<ScalarT, IdxT>::n_intern_;

  public:
    PartitionedDaeComponent(NodeT* bus, IdxT id)
      : node_ref_(bus)
    {
      size_           = 4;
      n_intern_       = 2;
      n_extern_       = 2;
      extern_indices_ = {0, 1};
      idc_            = id;
      nnz_            = 8;
    }

    int allocate() final
    {
      CircuitComponent<ScalarT, IdxT>::allocate();
      this->setExternalConnectionNodes(0, node_ref_->getNodeConnection(0));
      this->setExternalConnectionNodes(1, node_ref_->getNodeConnection(1));
      return 0;
    }

    int initialize() final
    {
      return 0;
    }

    /// Local order: 0 = y (diff), 1 = w (alg), 2 = x (diff), 3 = z (alg).
    int tagDifferentiable() final
    {
      tag_[0] = true;
      tag_[1] = false;
      tag_[2] = true;
      tag_[3] = false;
      return 0;
    }

    int evaluateInternalResidual() final
    {
      const ScalarT x = y_int_[0];
      const ScalarT z = y_int_[1];
      const ScalarT y = *y_ext_[0];

      f_int_[0] = -yp_int_[0] + y - x * x * x + z; // x equation
      f_int_[1] = std::sin(x - y) - z;             // z equation (algebraic)
      return 0;
    }

    int evaluateExternalResidual() final
    {
      *f_ext_[0] += -y_int_[0]; // -x into the y equation
      *f_ext_[1] += -y_int_[1]; // -z into the w equation
      return 0;
    }

    int evaluateJacobian() final
    {
      this->zeroJacMatrix();

      const RealT x = y_int_[0];
      const RealT y = *y_ext_[0];
      const RealT c = std::cos(x - y);

      // x equation: d/dx, d/dz, d/dy
      std::vector<IdxT>  row = {2, 2, 2};
      std::vector<IdxT>  col = {2, 3, 0};
      std::vector<RealT> val = {-3.0 * x * x - alpha_, 1.0, 1.0};
      this->setJacValues(row, col, val);

      // z equation: d/dx, d/dz, d/dy
      row = {3, 3, 3};
      col = {2, 3, 0};
      val = {c, -1.0, -c};
      this->setJacValues(row, col, val);

      // Contributions to the bus equations
      row = {0, 1};
      col = {2, 3};
      val = {-1.0, -1.0};
      this->setJacValues(row, col, val);

      return 0;
    }

    int evaluateIntegrand() final
    {
      return 0;
    }

    int initializeAdjoint() final
    {
      return 0;
    }

    int evaluateAdjointResidual() final
    {
      return 0;
    }

    int evaluateAdjointIntegrand() final
    {
      return 0;
    }

    int setAbsoluteTolerance(RealT rel_tol) final
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return 0;
    }

    CircuitComponent<ScalarT, IdxT>* clone() const
    {
      return new PartitionedDaeComponent<ScalarT, IdxT>(*this);
    }

  private:
    NodeT* node_ref_;
  };

  /// Bus self terms:
  ///   y equation: -y' + y(1 - y^2) - w
  ///   w equation:  w + y - pi/2
  template <class ScalarT, typename IdxT>
  class PartitionedDaeBus : public CircuitComponent<ScalarT, IdxT>
  {
    using RealT = typename CircuitComponent<ScalarT, IdxT>::RealT;
    using NodeT = typename PowerElectronics::NodeBase<ScalarT, IdxT>;

    using CircuitComponent<ScalarT, IdxT>::size_;
    using CircuitComponent<ScalarT, IdxT>::nnz_;
    using CircuitComponent<ScalarT, IdxT>::alpha_;
    using CircuitComponent<ScalarT, IdxT>::y_ext_;
    using CircuitComponent<ScalarT, IdxT>::yp_ext_;
    using CircuitComponent<ScalarT, IdxT>::tag_;
    using CircuitComponent<ScalarT, IdxT>::abs_tol_;
    using CircuitComponent<ScalarT, IdxT>::f_ext_;
    using CircuitComponent<ScalarT, IdxT>::idc_;
    using CircuitComponent<ScalarT, IdxT>::extern_indices_;
    using CircuitComponent<ScalarT, IdxT>::n_extern_;
    using CircuitComponent<ScalarT, IdxT>::n_intern_;

  public:
    PartitionedDaeBus(NodeT* bus, IdxT id)
      : node_ref_(bus)
    {
      size_           = 2;
      n_intern_       = 0;
      n_extern_       = 2;
      extern_indices_ = {0, 1};
      idc_            = id;
      nnz_            = 4;
    }

    int allocate() final
    {
      CircuitComponent<ScalarT, IdxT>::allocate();
      this->setExternalConnectionNodes(0, node_ref_->getNodeConnection(0));
      this->setExternalConnectionNodes(1, node_ref_->getNodeConnection(1));
      return 0;
    }

    int initialize() final
    {
      return 0;
    }

    /// Local order: 0 = y (diff), 1 = w (alg).
    int tagDifferentiable() final
    {
      tag_[0] = true;
      tag_[1] = false;
      return 0;
    }

    int evaluateInternalResidual() final
    {
      return 0;
    }

    int evaluateExternalResidual() final
    {
      const ScalarT y = *y_ext_[0];
      const ScalarT w = *y_ext_[1];

      *f_ext_[0] += -*yp_ext_[0] + y * (1.0 - y * y) - w;
      *f_ext_[1] += w + y - PartitionedDae::PI / 2.0;
      return 0;
    }

    int evaluateJacobian() final
    {
      this->zeroJacMatrix();

      const RealT y = *y_ext_[0];

      std::vector<IdxT>  row = {0, 0, 1, 1};
      std::vector<IdxT>  col = {0, 1, 0, 1};
      std::vector<RealT> val = {1.0 - 3.0 * y * y - alpha_, -1.0, 1.0, 1.0};
      this->setJacValues(row, col, val);

      return 0;
    }

    int evaluateIntegrand() final
    {
      return 0;
    }

    int initializeAdjoint() final
    {
      return 0;
    }

    int evaluateAdjointResidual() final
    {
      return 0;
    }

    int evaluateAdjointIntegrand() final
    {
      return 0;
    }

    int setAbsoluteTolerance(RealT rel_tol) final
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return 0;
    }

    CircuitComponent<ScalarT, IdxT>* clone() const
    {
      return new PartitionedDaeBus<ScalarT, IdxT>(*this);
    }

  private:
    NodeT* node_ref_;
  };

} // namespace GridKit