
#pragma once

#include <GridKit/Model/PowerElectronics/Component.hpp>
#include <GridKit/Model/PowerElectronics/NodeBase.hpp>

namespace GridKit
{
  namespace PowerElectronics
  {
    /*!
     * @brief Declaration of a MicrogridBusDQ class.
     *
     */
    template <class ScalarT, typename IdxT>
    class MicrogridBusDQ : public Component<ScalarT, IdxT>
    {
      using RealT = typename Component<ScalarT, IdxT>::RealT;
      using NodeT = typename PowerElectronics::NodeBase<ScalarT, IdxT>;

      using Component<ScalarT, IdxT>::size_;
      using Component<ScalarT, IdxT>::nnz_;
      using Component<ScalarT, IdxT>::time_;
      using Component<ScalarT, IdxT>::alpha_;
      using Component<ScalarT, IdxT>::y_ext_;
      using Component<ScalarT, IdxT>::y_int_;
      using Component<ScalarT, IdxT>::yp_ext_;
      using Component<ScalarT, IdxT>::yp_int_;
      using Component<ScalarT, IdxT>::tag_;
      using Component<ScalarT, IdxT>::abs_tol_;
      using Component<ScalarT, IdxT>::f_ext_;
      using Component<ScalarT, IdxT>::f_int_;
      using Component<ScalarT, IdxT>::g_;
      using Component<ScalarT, IdxT>::yB_;
      using Component<ScalarT, IdxT>::ypB_;
      using Component<ScalarT, IdxT>::fB_;
      using Component<ScalarT, IdxT>::gB_;
      using Component<ScalarT, IdxT>::param_;
      using Component<ScalarT, IdxT>::idc_;

      using Component<ScalarT, IdxT>::extern_indices_;
      using Component<ScalarT, IdxT>::n_extern_;
      using Component<ScalarT, IdxT>::n_intern_;

    public:
      MicrogridBusDQ(IdxT id, RealT RN, NodeT* node1);
      virtual ~MicrogridBusDQ();

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

      Component<ScalarT, IdxT>* clone() const;

    private:
      RealT  RN_;
      NodeT* node1_;
    };
  } // namespace PowerElectronics
} // namespace GridKit
