

#pragma once

#include <tuple>

#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>

namespace GridKit
{
  template <class ScalarT, typename IdxT>
  class BaseBus;
}

namespace GridKit
{
  /*!
   * @brief Declaration of a SynchronousMachine class.
   *
   */
  template <class ScalarT, typename IdxT>
  class SynchronousMachine : public CircuitComponent<ScalarT, IdxT>
  {
    using RealT = typename CircuitComponent<ScalarT, IdxT>::RealT;

    using CircuitComponent<ScalarT, IdxT>::size_;
    using CircuitComponent<ScalarT, IdxT>::nnz_;
    using CircuitComponent<ScalarT, IdxT>::time_;
    using CircuitComponent<ScalarT, IdxT>::alpha_;
    using CircuitComponent<ScalarT, IdxT>::y_;
    using CircuitComponent<ScalarT, IdxT>::int_;
    using CircuitComponent<ScalarT, IdxT>::yp_;
    using CircuitComponent<ScalarT, IdxT>::intp_;
    using CircuitComponent<ScalarT, IdxT>::tag_;
    using CircuitComponent<ScalarT, IdxT>::f_;
    using CircuitComponent<ScalarT, IdxT>::int_f_;
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
