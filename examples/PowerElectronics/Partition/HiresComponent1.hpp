

#pragma once

#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>

namespace GridKit
{
  /*!
   * @brief Declaration of a Hires Component 1 class.
   *
   */
  template <class ScalarT, typename IdxT>
  class HiresComponent1 : public CircuitComponent<ScalarT, IdxT>
  {
    using RealT   = typename CircuitComponent<ScalarT, IdxT>::RealT;
    using MatrixT = typename CircuitComponent<RealT, IdxT>::MatrixT;

    using CircuitComponent<ScalarT, IdxT>::size_;
    using CircuitComponent<ScalarT, IdxT>::nnz_;
    using CircuitComponent<ScalarT, IdxT>::time_;
    using CircuitComponent<ScalarT, IdxT>::alpha_;
    using CircuitComponent<ScalarT, IdxT>::y_;
    using CircuitComponent<ScalarT, IdxT>::yp_;
    using CircuitComponent<ScalarT, IdxT>::tag_;
    using CircuitComponent<ScalarT, IdxT>::f_;
    using CircuitComponent<ScalarT, IdxT>::g_;
    using CircuitComponent<ScalarT, IdxT>::yB_;
    using CircuitComponent<ScalarT, IdxT>::ypB_;
    using CircuitComponent<ScalarT, IdxT>::fB_;
    using CircuitComponent<ScalarT, IdxT>::gB_;
    using CircuitComponent<ScalarT, IdxT>::jac_;
    using CircuitComponent<ScalarT, IdxT>::param_;
    using CircuitComponent<ScalarT, IdxT>::idc_;

    using CircuitComponent<ScalarT, IdxT>::extern_indices_;
    using CircuitComponent<ScalarT, IdxT>::n_extern_;
    using CircuitComponent<ScalarT, IdxT>::n_intern_;

  public:
    HiresComponent1(IdxT id)
    {
      size_           = 5;
      n_intern_       = 3;
      n_extern_       = 2;
      extern_indices_ = {3, 4};
      idc_            = id;
    }

    ~HiresComponent1()
    {
    }

    int allocate()
    {
      y_.resize(static_cast<size_t>(size_));
      yp_.resize(static_cast<size_t>(size_));
      f_.resize(static_cast<size_t>(size_));

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

    int evaluateResidual()
    {
      f_[0] = yp_[0] + 1.71 * y_[0] - 0.43 * y_[1] - 8.32 * y_[2] - 0.0007;
      f_[1] = yp_[1] - 1.71 * y_[0] + 8.75 * y_[1];
      f_[2] = yp_[2] + 10.03 * y_[2] - 0.43 * y_[3] - 0.035 * y_[4];
      f_[3] = -8.32 * y_[1] - 1.71 * y_[2] + 0.1 * y_[3];
      f_[4] = +0.7 * y_[4];

      return 0;
    }

    int evaluateJacobian()
    {
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
  };
} // namespace GridKit
