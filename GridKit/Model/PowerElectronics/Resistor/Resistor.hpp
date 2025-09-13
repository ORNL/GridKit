

#pragma once

#include <GridKit/Model/EvaluatorMixins.hpp>
#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>

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
  class Resistor : public CircuitComponent<ScalarT, IdxT>, public Mixin::Evaluator::CsrJacobian<ScalarT, IdxT, Resistor>
  {
    using RealT = typename CircuitComponent<ScalarT, IdxT>::RealT;

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

    using typename Model::Evaluator<ScalarT, IdxT>::CsrJacobian;

  public:
    const static size_t SIZE = 2;

    Resistor(IdxT id, RealT R);
    virtual ~Resistor();

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

    template <bool INCLUDE_DIAGONALS, bool KEEP_SORTED, bool USE_TEMPLATE>
    CsrJacobian buildCsrJacobian(LinearAlgebra::CsrBuilder<RealT, IdxT, INCLUDE_DIAGONALS, KEEP_SORTED, USE_TEMPLATE> builder);

  private:
    RealT R_;
  };
} // namespace GridKit
