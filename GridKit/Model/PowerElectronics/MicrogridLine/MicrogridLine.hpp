

#pragma once

#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>
#include <GridKit/Model/PowerElectronics/NodeBase.hpp>

namespace GridKit
{
  template <typename scalar_type, typename index_type>
  class BaseBus;
}

namespace GridKit
{
  /*!
   * @brief Declaration of a MicrogridLine class.
   *
   */
  template <typename scalar_type, typename index_type>
  class MicrogridLine : public CircuitComponent<scalar_type, index_type>
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
    using CircuitComponent<scalar_type, index_type>::abs_tol_;
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
    using NodeT   = typename PowerElectronics::NodeBase<ScalarT, IdxT>;

    MicrogridLine(IdxT id, RealT R, RealT L, NodeT* node_ref, NodeT* bus1, NodeT* bus2);
    virtual ~MicrogridLine();

    int initialize();
    int allocate() final;
    int tagDifferentiable();
    int setAbsoluteTolerance(RealT);
    int evaluateInternalResidual() final;
    int evaluateExternalResidual() final;
    int evaluateJacobian();
    int evaluateIntegrand();

    int initializeAdjoint();
    int evaluateAdjointResidual();
    // int evaluateAdjointJacobian();
    int evaluateAdjointIntegrand();

  private:
    RealT  R_;
    RealT  L_;
    NodeT* node_ref_;
    NodeT* bus1_;
    NodeT* bus2_;
  };
} // namespace GridKit
