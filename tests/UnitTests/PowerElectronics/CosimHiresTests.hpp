#pragma once

/**
 * @file HiresComponents.hpp
 *
 * HIRES problem split into three GridKit components (copied from the
 * SubsystemModel HIRES unit test so it can be used outside the test suite).
 *
 * Assembled equations (standard HIRES, written as y' = g(y)):
 *
 *   y1' = -1.71 y1 + 0.43 y2 + 8.32 y3 + 0.0007                  (Component1)
 *   y2' =  1.71 y1 - 8.75 y2                                      (Component1)
 *   y3' = -10.03 y3 + 0.43 y4 + 0.035 y5                          (Component1)
 *   y4' =  8.32 y2 + 1.71 y3 - 1.12 y4                            (Bus + C1 + C3)
 *   y5' = -1.745 y5 + 0.43 y6 + 0.43 y7                           (Bus + C1 + C3)
 *   y6' = -280 y6 y8 + 0.69 y4 + 1.71 y5 - 0.43 y6 + 0.69 y7      (Component3)
 *   y7' =  280 y6 y8 - 1.81 y7                                    (Component3)
 *   y8' = -280 y6 y8 + 1.81 y7                                    (Component3)
 *
 * Every residual has the form F_i = g_i(y) - y_i', so all variables are
 * differential. y4 and y5 live on the bus node; each component adds its
 * contribution to the bus equations through its external residual.
 */

#include <vector>

#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>
#include <GridKit/Model/PowerElectronics/NodeBase.hpp>

namespace GridKit
{
  /// Equations y1, y2, y3 plus its contribution to the bus equations y4, y5.
  template <class ScalarT, typename IdxT>
  class HiresComponent1 : public CircuitComponent<ScalarT, IdxT>
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
    using CircuitComponent<ScalarT, IdxT>::abs_tol_;
    using CircuitComponent<ScalarT, IdxT>::f_ext_;
    using CircuitComponent<ScalarT, IdxT>::f_int_;
    using CircuitComponent<ScalarT, IdxT>::idc_;
    using CircuitComponent<ScalarT, IdxT>::extern_indices_;
    using CircuitComponent<ScalarT, IdxT>::n_extern_;
    using CircuitComponent<ScalarT, IdxT>::n_intern_;
    using CircuitComponent<ScalarT, IdxT>::tag_;

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

    int tagDifferentiable() final
    {
      std::fill(tag_.begin(), tag_.end(), true);
      return 0;
    }

    int evaluateInternalResidual() final
    {
      f_int_[0] = -yp_int_[0] - 1.71 * y_int_[0] + 0.43 * y_int_[1] + 8.32 * y_int_[2] + 0.0007;
      f_int_[1] = -yp_int_[1] + 1.71 * y_int_[0] - 8.75 * y_int_[1];
      f_int_[2] = -yp_int_[2] - 10.03 * y_int_[2] + 0.43 * *y_ext_[0] + 0.035 * *y_ext_[1];
      return 0;
    }

    int evaluateExternalResidual() final
    {
      *f_ext_[0] += 8.32 * y_int_[1] + 1.71 * y_int_[2] - 0.1 * *y_ext_[0];
      *f_ext_[1] += -0.7 * *y_ext_[1];
      return 0;
    }

    int evaluateJacobian() final
    {
      this->zeroJacMatrix();

      std::vector<IdxT>  row = {2, 2, 2, 3, 3, 4, 4, 4};
      std::vector<IdxT>  col = {2, 3, 4, 2, 3, 4, 0, 1};
      std::vector<RealT> val = {-1.71 - alpha_, 0.43, 8.32, 1.71, -8.75 - alpha_, -10.03 - alpha_, 0.43, 0.035};
      this->setJacValues(row, col, val);

      row = {0, 0, 0, 1};
      col = {3, 4, 0, 1};
      val = {8.32, 1.71, -0.1, -0.7};
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
      return new HiresComponent1<ScalarT, IdxT>(*this);
    }

  private:
    NodeT* node_ref_;
  };

  /// Self terms of the bus equations: y4' + y4 and y5' + y5.
  template <class ScalarT, typename IdxT>
  class HiresBus : public CircuitComponent<ScalarT, IdxT>
  {
    using RealT = typename CircuitComponent<ScalarT, IdxT>::RealT;
    using NodeT = typename PowerElectronics::NodeBase<ScalarT, IdxT>;

    using CircuitComponent<ScalarT, IdxT>::size_;
    using CircuitComponent<ScalarT, IdxT>::nnz_;
    using CircuitComponent<ScalarT, IdxT>::alpha_;
    using CircuitComponent<ScalarT, IdxT>::y_ext_;
    using CircuitComponent<ScalarT, IdxT>::yp_ext_;
    using CircuitComponent<ScalarT, IdxT>::abs_tol_;
    using CircuitComponent<ScalarT, IdxT>::f_ext_;
    using CircuitComponent<ScalarT, IdxT>::idc_;
    using CircuitComponent<ScalarT, IdxT>::extern_indices_;
    using CircuitComponent<ScalarT, IdxT>::n_extern_;
    using CircuitComponent<ScalarT, IdxT>::n_intern_;
    using CircuitComponent<ScalarT, IdxT>::tag_;

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

    int tagDifferentiable() final
    {
      std::fill(tag_.begin(), tag_.end(), true);
      return 0;
    }

    int evaluateInternalResidual() final
    {
      return 0;
    }

    int evaluateExternalResidual() final
    {
      *f_ext_[0] += -*yp_ext_[0] - *y_ext_[0];
      *f_ext_[1] += -*yp_ext_[1] - *y_ext_[1];
      return 0;
    }

    int evaluateJacobian() final
    {
      this->zeroJacMatrix();
      std::vector<IdxT>  row = {0, 1};
      std::vector<IdxT>  col = {0, 1};
      std::vector<RealT> val = {-alpha_ - 1.0, -alpha_ - 1.0};
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
      return new HiresBus<ScalarT, IdxT>(*this);
    }

  private:
    NodeT* node_ref_;
  };

  /// Equations y6, y7, y8 plus its contribution to the bus equations y4, y5.
  template <class ScalarT, typename IdxT>
  class HiresComponent3 : public CircuitComponent<ScalarT, IdxT>
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
    using CircuitComponent<ScalarT, IdxT>::abs_tol_;
    using CircuitComponent<ScalarT, IdxT>::f_ext_;
    using CircuitComponent<ScalarT, IdxT>::f_int_;
    using CircuitComponent<ScalarT, IdxT>::idc_;
    using CircuitComponent<ScalarT, IdxT>::extern_indices_;
    using CircuitComponent<ScalarT, IdxT>::n_extern_;
    using CircuitComponent<ScalarT, IdxT>::n_intern_;
    using CircuitComponent<ScalarT, IdxT>::tag_;

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

    int tagDifferentiable() final
    {
      std::fill(tag_.begin(), tag_.end(), true);
      return 0;
    }

    int evaluateInternalResidual() final
    {
      f_int_[0] = -yp_int_[0] - 280 * y_int_[0] * y_int_[2] + 0.69 * *y_ext_[0] + 1.71 * *y_ext_[1]
                  - 0.43 * y_int_[0] + 0.69 * y_int_[1];
      f_int_[1] = -yp_int_[1] + 280 * y_int_[0] * y_int_[2] - 1.81 * y_int_[1];
      f_int_[2] = -yp_int_[2] - 280 * y_int_[0] * y_int_[2] + 1.81 * y_int_[1];
      return 0;
    }

    int evaluateExternalResidual() final
    {
      *f_ext_[0] += -0.02 * *y_ext_[0];
      *f_ext_[1] += -0.045 * *y_ext_[1] + 0.43 * y_int_[0] + 0.43 * y_int_[1];
      return 0;
    }

    int evaluateJacobian() final
    {
      this->zeroJacMatrix();

      std::vector<IdxT>  row = {2, 2, 2, 2, 2};
      std::vector<IdxT>  col = {2, 3, 4, 0, 1};
      std::vector<RealT> val = {-280 * y_int_[2] - 0.43 - alpha_, 0.69, -280 * y_int_[0], 0.69, 1.71};
      this->setJacValues(row, col, val);

      row = {3, 3, 3};
      col = {2, 3, 4};
      val = {280 * y_int_[2], -1.81 - alpha_, 280 * y_int_[0]};
      this->setJacValues(row, col, val);

      row = {4, 4, 4};
      col = {2, 3, 4};
      val = {-280 * y_int_[2], 1.81, -280 * y_int_[0] - alpha_};
      this->setJacValues(row, col, val);

      row = {0, 1, 1, 1};
      col = {0, 2, 3, 1};
      val = {-0.02, 0.43, 0.43, -0.045};
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
      return new HiresComponent3<ScalarT, IdxT>(*this);
    }

  private:
    NodeT* node_ref_;
  };

} // namespace GridKit