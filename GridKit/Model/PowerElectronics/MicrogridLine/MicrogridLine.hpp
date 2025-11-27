

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
   * @brief Declaration of a MicrogridLine class.
   *
   */
  template <class ScalarT, typename IdxT>
  class MicrogridLine : public CircuitComponent<ScalarT, IdxT>
  {
    using RealT   = typename CircuitComponent<ScalarT, IdxT>::RealT;
    using MatrixT = typename CircuitComponent<RealT, IdxT>::MatrixT;

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
    MicrogridLine(IdxT id, RealT R, RealT L);
    MicrogridLine(const MicrogridLine<ScalarT, IdxT>& other);
    CircuitComponent<ScalarT, IdxT>* clone() const override
    {
      return new MicrogridLine<ScalarT, IdxT>(*this);
    }
    
    CircuitComponent<ScalarT, IdxT>& operator=(const CircuitComponent<ScalarT, IdxT>& other);
    virtual ~MicrogridLine();

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

    RealT getResistance() const
    {
      return R_;
    }

    template<class S>
    RealT setResistance(S r)
    {
      R_ = r;
    }

    RealT getInductance() const
    {
      return L_;
    }

    template<class S>
    RealT setInductance(S l)
    {
      L_ = l;
    }

  private:
    RealT R_;
    RealT L_;
  };
} // namespace GridKit
