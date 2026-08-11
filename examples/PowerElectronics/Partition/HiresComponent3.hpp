#pragma once

#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>
#include <GridKit/Model/PowerElectronics/NodeBase.hpp>

namespace GridKit
{
  /*!
   * @brief Hires Component 3 class.
   *
   */
  template <class ScalarT, typename IdxT>
  class HiresComponent3 : public CircuitComponent<ScalarT, IdxT>
  {
    using RealT = typename CircuitComponent<ScalarT, IdxT>::RealT;
    using NodeT = typename PowerElectronics::NodeBase<ScalarT, IdxT>;

    using CircuitComponent<ScalarT, IdxT>::size_;
    using CircuitComponent<ScalarT, IdxT>::nnz_;
    using CircuitComponent<ScalarT, IdxT>::time_;
    using CircuitComponent<ScalarT, IdxT>::alpha_;
    using CircuitComponent<ScalarT, IdxT>::y_ext_;
    using CircuitComponent<ScalarT, IdxT>::y_int_;
    using CircuitComponent<ScalarT, IdxT>::yp_ext_;
    using CircuitComponent<ScalarT, IdxT>::yp_int_;
    using CircuitComponent<ScalarT, IdxT>::tag_;
    using CircuitComponent<ScalarT, IdxT>::abs_tol_;
    using CircuitComponent<ScalarT, IdxT>::f_ext_;
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
    HiresComponent3(NodeT* bus, IdxT id)
      : node_ref_(bus)
    {
      size_           = 5;
      n_intern_       = 3;
      n_extern_       = 2;
      extern_indices_ = {0, 1};
      idc_            = id;
      nnz_            = 15;
    }

    ~HiresComponent3()
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

      // Internals
      f_int_[0] = -yp_int_[0] - 280 * y_int_[0] * y_int_[2] + 0.69 * *y_ext_[0] + 1.71 * *y_ext_[1] - 0.43 * y_int_[0] + 0.69 * y_int_[1];
      f_int_[1] = -yp_int_[1] + 280 * y_int_[0] * y_int_[2] - 1.81 * y_int_[1];
      f_int_[2] = -yp_int_[2] - 280 * y_int_[0] * y_int_[2] + 1.81 * y_int_[1];

      return 0;
    }

    int evaluateExternalResidual()
    {
      // Externals
      *f_ext_[0] += -0.02 * *y_ext_[0];
      *f_ext_[1] += -0.045 * *y_ext_[1] + 0.43 * y_int_[0] + 0.43 * y_int_[1];

      return 0;
    }

    int evaluateJacobian()
    {
      this->zeroJacMatrix();

      // Internal Jacobian Entries [row 1]
      std::vector<IdxT>  row = {2, 2, 2, 2, 2};
      std::vector<IdxT>  col = {2, 3, 4, 0, 1};
      std::vector<RealT> val = {-280 * y_int_[2] - 0.43 - alpha_, 0.69, -280 * y_int_[0], 0.69, 1.71};

      this->setJacValues(row, col, val);

      // Internal Jacobian Entries [row 2]
      row = {3, 3, 3};
      col = {2, 3, 4};
      val = {280 * y_int_[2], -1.81 - alpha_, 280 * y_int_[0]};

      this->setJacValues(row, col, val);

      // Internal Jacobian Entries [row 3]
      row = {4, 4, 4};
      col = {2, 3, 4};
      val = {-280 * y_int_[2], 1.81, -280 * y_int_[0] - alpha_};

      this->setJacValues(row, col, val);

      // External Jacobian Entries
      row = {0, 1, 1, 1};
      col = {0, 2, 3, 1};
      val = {-0.02, 0.43, 0.43, -0.045};

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
