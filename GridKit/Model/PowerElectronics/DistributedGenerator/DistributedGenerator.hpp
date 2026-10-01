
#pragma once

#include <GridKit/Model/PowerElectronics/Component.hpp>
#include <GridKit/Model/PowerElectronics/NodeBase.hpp>

namespace GridKit
{
  namespace PowerElectronics
  {
    /*!
     * @brief Declaration of a DistributedGenerator parameter struct.
     *
     */
    template <typename RealT, typename IdxT>
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

    /*!
     * @brief Declaration of a DistributedGenerator class.
     *
     */
    template <class ScalarT, typename IdxT>
    class DistributedGenerator : public Component<ScalarT, IdxT>
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
      using Component<ScalarT, IdxT>::abs_tol_;
      using Component<ScalarT, IdxT>::tag_;
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
      DistributedGenerator(IdxT                                        id,
                           DistributedGeneratorParameters<RealT, IdxT> parm,
                           bool                                        reference_frame,
                           NodeT*                                      node_ref,
                           NodeT*                                      node_bus);
      virtual ~DistributedGenerator();

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
      bool  refframe_;

      NodeT* node_ref_;
      NodeT* node_bus_;
    };
  } // namespace PowerElectronics
} // namespace GridKit
