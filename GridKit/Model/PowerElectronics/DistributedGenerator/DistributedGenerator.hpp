

#pragma once

#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>

namespace GridKit
{
  template <class ScalarT, typename IdxT>
  class BaseBus;

  template <class RealT, typename IdxT>
  struct DistributedGeneratorParameters
  {
    RealT wb_;
    RealT wc_;
    RealT mp_;
    RealT Vn_;
    RealT nq_;
    RealT F_;
    RealT Kiv_;
    RealT Kpv_;
    RealT Kic_;
    RealT Kpc_;
    RealT Cf_;
    RealT rLf_;
    RealT Lf_;
    RealT rLc_;
    RealT Lc_;
  };
} // namespace GridKit

namespace GridKit
{
  /*!
   * @brief Declaration of a DistributedGenerator class.
   *
   */
  template <class ScalarT, typename IdxT>
  class DistributedGenerator : public CircuitComponent<ScalarT, IdxT>
  {
    using real_type = typename CircuitComponent<ScalarT, IdxT>::real_type;

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
    DistributedGenerator(IdxT                                            id,
                         DistributedGeneratorParameters<real_type, IdxT> parm,
                         bool                                            reference_frame);
    virtual ~DistributedGenerator();

    int allocate();
    int initialize();
    int tagDifferentiable();
    int evaluateResidual();
    int evaluateJacobian();
    int evaluateIntegrand();
    int initializeAdjoint();
    int evaluateAdjointResidual();
    // int evaluateAdjointJacobian();
    int evaluateAdjointIntegrand();

  private:
    ScalarT wb_;
    ScalarT wc_;
    ScalarT mp_;
    ScalarT Vn_;
    ScalarT nq_;
    ScalarT F_;
    ScalarT Kiv_;
    ScalarT Kpv_;
    ScalarT Kic_;
    ScalarT Kpc_;
    ScalarT Cf_;
    ScalarT rLf_;
    ScalarT Lf_;
    ScalarT rLc_;
    ScalarT Lc_;
    bool    refframe_;
  };
} // namespace GridKit
