

#pragma once

#include <tuple>

#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>

namespace GridKit
{
  template <typename scalar_type, typename index_type>
  class BaseBus;
}

namespace GridKit
{
  /*!
   * @brief Declaration of a SynchronousMachine class.
   *
   */
  template <typename scalar_type, typename index_type>
  class SynchronousMachine : public CircuitComponent<scalar_type, index_type>
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

    SynchronousMachine(IdxT id, RealT Lls, std::tuple<RealT, RealT> Llkq, RealT Llfd, RealT Llkd, RealT Lmq, RealT Lmd, RealT Rs, std::tuple<RealT, RealT> Rkq, RealT Rfd, RealT Rkd, RealT RJ, RealT P, RealT mub);
    virtual ~SynchronousMachine();

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
    RealT                    Lls_;
    std::tuple<RealT, RealT> Llkq_;
    RealT                    Llfd_;
    RealT                    Llkd_;
    RealT                    Lmq_;
    RealT                    Lmd_;
    RealT                    Rs_;
    std::tuple<RealT, RealT> Rkq_;
    RealT                    Rfd_;
    RealT                    Rkd_;
    RealT                    RJ_;
    RealT                    P_;
    RealT                    mub_;
  };
} // namespace GridKit
