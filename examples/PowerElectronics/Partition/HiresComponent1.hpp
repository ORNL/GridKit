

#pragma once

#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>
#include <GridKit/Model/PowerElectronics/NodeBase.hpp>

namespace GridKit
{
  /*!
   * @brief Declaration of a Hires Component 1 class.
   *
   */
  template <class ScalarT, typename IdxT>
  class HiresComponent1 : public CircuitComponent<ScalarT, IdxT>
  {
    using RealT = typename CircuitComponent<ScalarT, IdxT>::RealT;
    using NodeT = typename PowerElectronics::NodeBase<ScalarT, IdxT>;

    using CircuitComponent<ScalarT, IdxT>::size_;
    using CircuitComponent<ScalarT, IdxT>::nnz_;
    using CircuitComponent<ScalarT, IdxT>::time_;
    using CircuitComponent<ScalarT, IdxT>::alpha_;
    using CircuitComponent<ScalarT, IdxT>::y_;
    using CircuitComponent<ScalarT, IdxT>::y_int_;
    using CircuitComponent<ScalarT, IdxT>::yp_;
    using CircuitComponent<ScalarT, IdxT>::yp_int_;
    using CircuitComponent<ScalarT, IdxT>::tag_;
    using CircuitComponent<ScalarT, IdxT>::abs_tol_;
    using CircuitComponent<ScalarT, IdxT>::f_;
    using CircuitComponent<ScalarT, IdxT>::f_int_;
    using CircuitComponent<ScalarT, IdxT>::g_;
    using CircuitComponent<ScalarT, IdxT>::yB_;
    using CircuitComponent<ScalarT, IdxT>::ypB_;
    using CircuitComponent<ScalarT, IdxT>::fB_;
    using CircuitComponent<ScalarT, IdxT>::gB_;
    using CircuitComponent<ScalarT, IdxT>::param_;
    using CircuitComponent<ScalarT, IdxT>::idc_;

    using CircuitComponent<ScalarT, IdxT>::extern_indices_;
    using CircuitComponent<ScalarT, IdxT>::n_extern_;
    using CircuitComponent<ScalarT, IdxT>::n_intern_;

  public:
    HiresComponent1(NodeT* bus, IdxT id)
      : node_ref_(bus)
    {
      size_           = 5;
      n_intern_       = 3;
      n_extern_       = 2;
      extern_indices_ = {0, 1};
      idc_            = id;
      nnz_            = 12;
    }

    ~HiresComponent1()
    {
    }

    int allocate()
    {
      CircuitComponent<ScalarT, IdxT>::allocate();

      this->setExternalConnectionNodes(0, node_ref_->getNodeConnection(0));
      this->setExternalConnectionNodes(1, node_ref_->getNodeConnection(1));

      return 0;
    }

    int initialize()
    {
      return 0;
    }

    int tagDifferentiable()
    {
      return 0;
    }

    int evaluateInternalResidual()
    {
      auto* y = y_.getData();

      // Internals
      f_int_[0] = -yp_int_[0] - 1.71 * y_int_[0] + 0.43 * y_int_[1] + 8.32 * y_int_[2] + 0.0007;
      f_int_[1] = -yp_int_[1] + 1.71 * y_int_[0] - 8.75 * y_int_[1];
      f_int_[2] = -yp_int_[2] - 10.03 * y_int_[2] + 0.43 * y[0] + 0.035 * y[1];

      return 0;
    }

    int evaluateExternalResidual()
    {
      auto* y = y_.getData();
      auto* f = f_.getData();

      // outputs
      f[0] = 8.32 * y_int_[1] + 1.71 * y_int_[2] - 0.1 * y[0];
      f[1] = -0.7 * y[1];

      f_.setDataUpdated();

      return 0;
    }

    int evaluateJacobian()
    {

      this->zeroJacMatrix();

      // Internal Jacobian Entries
      std::vector<IdxT>  row = {2, 2, 2, 3, 3, 4, 4, 4};
      std::vector<IdxT>  col = {2, 3, 4, 2, 3, 4, 0, 1};
      std::vector<RealT> val = {-1.71 - alpha_, 0.43, 8.32, 1.71, -8.75 - alpha_, -10.03 - alpha_, 0.43, 0.035};

      this->setJacValues(row, col, val);

      // External Jacobian Entries
      row = {0, 0, 0, 1};
      col = {3, 4, 0, 1};
      val = {8.32, 1.71, -0.1, -0.7};

      this->setJacValues(row, col, val);

      return 0;
    }

    int evaluateIntegrand()
    {
      return 0;
    }

    int initializeAdjoint()
    {
      return 0;
    }

    int evaluateAdjointResidual()
    {
      return 0;
    }

    int evaluateAdjointIntegrand()
    {
      return 0;
    }

    /**
     * @brief Compute the absolute tolerance for each variable in the model
     *
     * @param rel_tol The relative tolerance which can be used to pick the
     *        absolute tolerance.
     * @tparam ScalarT Scalar data type
     * @tparam IdxT Index data type
     * @return int 0 if successful, non-zero otherwise.
     *
     * This represents a "noise" level close to zero for which pure relative
     * error cannot be used.
     */
    int setAbsoluteTolerance(RealT rel_tol)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return 0;
    }

  private:
    NodeT* node_ref_;
  };
} // namespace GridKit
