

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
   * @brief Declaration of a TransmissionLine class.
   *
   * Model from Adam Birchfield paper (medium distances < 2km).
   * See also textbooks "Power System Analysis" by Grainger and "Power System Dynamics and Stability" by Sauer & Pai
   *
   * @note Not used in the Microgrid model.
   */
  template <typename scalar_type, typename index_type>
  class TransmissionLine : public CircuitComponent<scalar_type, index_type>
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

    TransmissionLine(IdxT id, RealT R, RealT X, RealT B);
    virtual ~TransmissionLine();

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
    RealT R_;
    RealT X_;
    RealT B_;
    RealT YReMat_;
    RealT YImMatDi_;
    RealT YImMatOff_;
  };
} // namespace GridKit
