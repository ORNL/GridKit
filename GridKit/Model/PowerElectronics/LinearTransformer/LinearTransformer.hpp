

#pragma once

#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>

namespace GridKit
{
  template <class ScalarT, typename IdxT>
  class BaseBus;
}

namespace GridKit
{
  /*!
   * @brief Declaration of a LinearTransformer class.
   *
   */
  template <class ScalarT, typename IdxT>
  class LinearTransformer : public CircuitComponent<ScalarT, IdxT>
  {
    using RealT = typename CircuitComponent<ScalarT, IdxT>::RealT;

    using CircuitComponent<ScalarT, IdxT>::size_;
    using CircuitComponent<ScalarT, IdxT>::nnz_;
    using CircuitComponent<ScalarT, IdxT>::time_;
    using CircuitComponent<ScalarT, IdxT>::alpha_;
    using CircuitComponent<ScalarT, IdxT>::y_ext_;
    using CircuitComponent<ScalarT, IdxT>::y_int_;
    using CircuitComponent<ScalarT, IdxT>::yp_ext_;
    using CircuitComponent<ScalarT, IdxT>::yp_int_;
    using CircuitComponent<ScalarT, IdxT>::tag_;
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
    LinearTransformer(IdxT id, RealT L0, RealT L1, RealT R0, RealT R1, RealT M);
    virtual ~LinearTransformer();

    int initialize();
    int tagDifferentiable();
    int evaluateInternalResidual() final;
    int evaluateExternalResidual() final;
    int evaluateJacobian();
    int evaluateIntegrand();

    int initializeAdjoint();
    int evaluateAdjointResidual();
    // int evaluateAdjointJacobian();
    int evaluateAdjointIntegrand();

  private:
    RealT L0_;
    RealT L1_;
    RealT R0_;
    RealT R1_;
    RealT M_;
  };
} // namespace GridKit
