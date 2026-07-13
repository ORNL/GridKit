#pragma once

#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>

namespace GridKit
{
  /*!
   * @brief Declaration of a Hires Bus Component.
   *
   */
  template <class ScalarT, typename IdxT>
  class HiresBus : public CircuitComponent<ScalarT, IdxT>
  {
    using RealT   = typename CircuitComponent<ScalarT, IdxT>::RealT;
    using MatrixT = typename CircuitComponent<RealT, IdxT>::MatrixT;

    using CircuitComponent<ScalarT, IdxT>::size_;
    using CircuitComponent<ScalarT, IdxT>::nnz_;
    using CircuitComponent<ScalarT, IdxT>::time_;
    using CircuitComponent<ScalarT, IdxT>::alpha_;
    using CircuitComponent<ScalarT, IdxT>::y_;
    using CircuitComponent<ScalarT, IdxT>::y_int_;
    using CircuitComponent<ScalarT, IdxT>::yp_;
    using CircuitComponent<ScalarT, IdxT>::yp_int_;
    using CircuitComponent<ScalarT, IdxT>::tag_;
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
    HiresBus(IdxT id)
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
      f_[0] = -yp_[0] - y_[0];
      f_[1] = -yp_[1] - y_[1];
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
  };
} // namespace GridKit
