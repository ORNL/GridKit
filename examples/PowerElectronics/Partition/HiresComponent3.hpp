

#pragma once

#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>

namespace GridKit
{
  /*!
   * @brief Declaration of a Hires Component 3 class.
   *
   */
  template <class ScalarT, typename IdxT>
  class HiresComponent3 : public CircuitComponent<ScalarT, IdxT>
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
    HiresComponent3(IdxT id)
    {
      size_           = 5;
      n_intern_       = 3;
      n_extern_       = 2;
      extern_indices_ = {0, 1};
      idc_            = id;
    }

    ~HiresComponent3()
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
      // Outputs
      f_[0] = 0.02 * y_[0];
      f_[1] = 0.045 * y_[1] - 0.43 * y_[2] - 0.43 * y_[3];

      // Internals
      f_[2] = yp_[2] + 280 * y_[2] * y_[4] - 0.69 * y_[0] - 1.71 * y_[1] + 0.43 * y_[2] - 0.69 * y_[3];
      f_[3] = yp_[3] - 280 * y_[2] * y_[4] + 1.81 * y_[3];
      f_[4] = yp_[4] + 280 * y_[2] * y_[4] - 1.81 * y_[3];

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
