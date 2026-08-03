#pragma once

#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>
#include <GridKit/Model/PowerElectronics/NodeBase.hpp>

namespace GridKit
{
  /*!
   * @brief Declaration of a Hires Bus Component.
   *
   */
  template <class ScalarT, typename IdxT>
  class HiresBus : public CircuitComponent<ScalarT, IdxT>
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
    HiresBus(NodeT* bus, IdxT id)
      : node_ref_(bus)
    {
      size_           = 2;
      n_intern_       = 0;
      n_extern_       = 2;
      extern_indices_ = {0, 1};
      idc_            = id;
      nnz_            = 2;
    }

    ~HiresBus()
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
      return 0;
    }

    int evaluateExternalResidual()
    {
      *f_ext_[0] += -*yp_ext_[0] - *y_ext_[0];
      *f_ext_[1] += -*yp_ext_[1] - *y_ext_[1];

      this->getResidual().setDataUpdated();

      return 0;
    }

    int evaluateJacobian()
    {
      this->zeroJacMatrix();

      std::vector<IdxT>  row = {0, 1};
      std::vector<IdxT>  col = {0, 1};
      std::vector<RealT> val = {-alpha_ - 1.0, -alpha_ - 1.0};

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