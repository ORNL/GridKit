

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
    SynchronousMachine(IdxT id, real_type Lls, std::tuple<real_type, real_type> Llkq, real_type Llfd, real_type Llkd, real_type Lmq, real_type Lmd, real_type Rs, std::tuple<real_type, real_type> Rkq, real_type Rfd, real_type Rkd, real_type RJ, real_type P, real_type mub);
    virtual ~SynchronousMachine();

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
    real_type                        Lls_;
    std::tuple<real_type, real_type> Llkq_;
    real_type                        Llfd_;
    real_type                        Llkd_;
    real_type                        Lmq_;
    real_type                        Lmd_;
    real_type                        Rs_;
    std::tuple<real_type, real_type> Rkq_;
    real_type                        Rfd_;
    real_type                        Rkd_;
    real_type                        RJ_;
    real_type                        P_;
    real_type                        mub_;
  };
} // namespace GridKit
