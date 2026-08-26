

#pragma once

#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>
#include <GridKit/Model/PowerElectronics/NodeBase.hpp>

namespace GridKit
{
  template <class ScalarT, typename IdxT>
  class BaseBus;
}

namespace GridKit
{
  /*!
   * @brief Declaration of a Resistor class.
   *
   */
  template <class ScalarT, typename IdxT>
  class Resistor : public CircuitComponent<ScalarT, IdxT>
  {
    using RealT = typename CircuitComponent<ScalarT, IdxT>::RealT;
    using NodeT = typename PowerElectronics::NodeBase<ScalarT, IdxT>;

    using CircuitComponent<ScalarT, IdxT>::size_;
    using CircuitComponent<ScalarT, IdxT>::nnz_;
    using CircuitComponent<ScalarT, IdxT>::time_;
    using CircuitComponent<ScalarT, IdxT>::alpha_;
    using CircuitComponent<ScalarT, IdxT>::y_ext_;
    using CircuitComponent<ScalarT, IdxT>::y_int_;
    using CircuitComponent<ScalarT, IdxT>::yp_ext_;
    using CircuitComponent<ScalarT, IdxT>::yp_int_;
    using CircuitComponent<ScalarT, IdxT>::tag_;
    using CircuitComponent<ScalarT, IdxT>::abs_tol_;
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
    Resistor(IdxT id, RealT R, NodeT* node1, NodeT* node2);
    virtual ~Resistor();

    int initialize();
    int allocate() final;
    int tagDifferentiable();
    int setAbsoluteTolerance(RealT);
    int evaluateInternalResidual() final;
    int evaluateExternalResidual() final;
    int evaluateJacobian();
    int evaluateIntegrand();

    int  initializeAdjoint();
    int  evaluateAdjointResidual();
    // int evaluateAdjointJacobian();
    int  evaluateAdjointIntegrand();
    bool isCloneable() const;

    CircuitComponent<ScalarT, IdxT>* clone() const;

  private:
    RealT  R_;
    NodeT* node1_;
    NodeT* node2_;
  };
} // namespace GridKit
