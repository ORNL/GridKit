

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
   * @brief Declaration of a MicrogridBusDQ class.
   *
   */
  template <class ScalarT, typename IdxT>
  class MicrogridBusDQ : public CircuitComponent<ScalarT, IdxT>
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

  public:
    MicrogridBusDQ(IdxT id, RealT RN);
    MicrogridBusDQ(const MicrogridBusDQ& other);

    CircuitComponent<ScalarT, IdxT>* clone() const override
    {
      return new MicrogridBusDQ<ScalarT, IdxT>(*this);
    }
    CircuitComponent<ScalarT, IdxT>& operator=(const CircuitComponent<ScalarT, IdxT>& other);

    virtual ~MicrogridBusDQ();

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

    RealT getVirtualResistance() const
    {
      return RN_;
    }

    template<class S>
    void setVirtualResistance(S rn)
    {
      RN_ = rn;
    }



  private:
    RealT RN_;
  };
} // namespace GridKit
