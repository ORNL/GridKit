/**
 * @file SexsPti.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the SEXS-PTI exciter model.
 */

#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ConnectedElement.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      template <typename RealP, typename IdxP>
      struct SexsPtiData;
    } // namespace Exciter

    template <typename ScalarP, typename IdxP>
    class BusBase;

    template <typename ScalarP, typename IdxP>
    class SignalNode;

  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      template <typename, typename>
      class SexsPti;

      /// Internal variables of a `SexsPti`.
      enum class SexsPtiInternalVariables : size_t
      {
        VR,  ///< Lead-lag block state
        EFD, ///< Exciter field voltage output
        VTR, ///< Terminal voltage error signal
        MAXIMUM,
      };

      /// External variables of a `SexsPti`.
      enum class SexsPtiExternalVariables : size_t
      {
        VS, ///< Stabilizer output signal
        MAXIMUM,
      };
    } // namespace Exciter

    template <typename ScalarP, typename IdxP>
    struct ConnectedElementTraits<Exciter::SexsPti<ScalarP, IdxP>>
    {
      using SexsPtiT = Exciter::SexsPti<ScalarP, IdxP>;

      using ElementT           = SexsPtiT;
      using ScalarT            = ScalarP;
      using IdxT               = IdxP;
      using RealT              = typename ScalarTraits<ScalarT>::RealT;
      using ModelDataT         = Exciter::SexsPtiData<RealT, IdxT>;
      using InternalVariablesT = Exciter::SexsPtiInternalVariables;
      using ExternalVariablesT = Exciter::SexsPtiExternalVariables;
      using InterfaceT         = Component<ScalarT, IdxT>;
    };

    namespace Exciter
    {
      template <typename ScalarP, typename IdxP>
      class SexsPti : public ConnectedElement<SexsPti<ScalarP, IdxP>>
      {
        using ConnectedElement<SexsPti>::gridkit_component_id_;
        using ConnectedElement<SexsPti>::alpha_;
        using ConnectedElement<SexsPti>::f_;
        using ConnectedElement<SexsPti>::nnz_;
        using ConnectedElement<SexsPti>::size_;
        using ConnectedElement<SexsPti>::tag_;
        using ConnectedElement<SexsPti>::time_;
        using ConnectedElement<SexsPti>::y_;
        using ConnectedElement<SexsPti>::yp_;
        using ConnectedElement<SexsPti>::wb_;
        using ConnectedElement<SexsPti>::J_;
        using ConnectedElement<SexsPti>::J_rows_buffer_;
        using ConnectedElement<SexsPti>::J_cols_buffer_;
        using ConnectedElement<SexsPti>::J_vals_buffer_;
        using ConnectedElement<SexsPti>::variable_indices_;
        using ConnectedElement<SexsPti>::residual_indices_;
        using ConnectedElement<SexsPti>::signals_;
        using ConnectedElement<SexsPti>::monitor_;

      public:
        using ScalarT    = typename ConnectedElement<SexsPti>::ScalarT;
        using IdxT       = typename ConnectedElement<SexsPti>::IdxT;
        using RealT      = typename ConnectedElement<SexsPti>::RealT;
        using ModelDataT = SexsPtiData<RealT, IdxT>;
        using SignalT    = SignalNode<ScalarT, IdxT>;
        using BusT       = BusBase<ScalarT, IdxT>;
        using MonitorT   = typename ConnectedElement<SexsPti>::MonitorT;

        SexsPti(BusT* bus);
        SexsPti(BusT* bus, const ModelDataT& data);
        ~SexsPti();

        int setGridKitComponentID(IdxT) override final;
        int allocate() override final;
        int verify() const override final;
        int initialize() override final;
        int tagDifferentiable() override final;
        int evaluateResidual() override final;
        int evaluateJacobian() override final;

        __attribute__((always_inline)) inline int evaluateInternalResidual(
            ScalarT*, ScalarT*, ScalarT*, ScalarT*, ScalarT*);

      private:
        BusT* bus_{nullptr};

        RealT Ta_{0};
        RealT Tb_{0};
        RealT Te_{0};
        RealT K_{0};
        RealT Efdmax_{0};
        RealT Efdmin_{0};

        int missing_param_count_{0};

        ScalarT vref_{0};
        ScalarT vOEL_{0};
        ScalarT vUEL_{0};

        void initModelParams(const ModelDataT& data);
        void initializeMonitor();

        std::vector<ScalarT> ws_;
        std::vector<IdxT>    ws_indices_;
      };

    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
