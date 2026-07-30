

#pragma once

#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>

namespace GridKit
{
  template <typename scalar_type, typename index_type>
  class BaseBus;
}

namespace GridKit
{
  /*!
   * @brief Declaration of a LinearTransformer class.
   *
   */
  template <typename scalar_type, typename index_type>
  class LinearTransformer : public CircuitComponent<scalar_type, index_type>
  {
    using CircuitComponent<scalar_type, index_type>::size_;
    using CircuitComponent<scalar_type, index_type>::nnz_;
    using CircuitComponent<scalar_type, index_type>::time_;
    using CircuitComponent<scalar_type, index_type>::alpha_;
    using CircuitComponent<scalar_type, index_type>::y_;
    using CircuitComponent<scalar_type, index_type>::y_int_;
    using CircuitComponent<scalar_type, index_type>::yp_;
    using CircuitComponent<scalar_type, index_type>::yp_int_;
    using CircuitComponent<scalar_type, index_type>::tag_;
    using CircuitComponent<scalar_type, index_type>::f_;
    using CircuitComponent<scalar_type, index_type>::f_int_;
    using CircuitComponent<scalar_type, index_type>::g_;
    using CircuitComponent<scalar_type, index_type>::yB_;
    using CircuitComponent<scalar_type, index_type>::ypB_;
    using CircuitComponent<scalar_type, index_type>::fB_;
    using CircuitComponent<scalar_type, index_type>::gB_;
    using CircuitComponent<scalar_type, index_type>::param_;
    using CircuitComponent<scalar_type, index_type>::idc_;

    using CircuitComponent<scalar_type, index_type>::extern_indices_;
    using CircuitComponent<scalar_type, index_type>::n_extern_;
    using CircuitComponent<scalar_type, index_type>::n_intern_;

  public:
    using ScalarT = scalar_type;
    using IdxT    = index_type;
    using RealT   = typename CircuitComponent<ScalarT, IdxT>::RealT;

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
