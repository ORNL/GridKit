#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ConnectedElement.hpp>
#include <GridKit/Model/PhasorDynamics/LoadZIP/LoadZIPData.hpp>

// Forward declarations.
namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename ScalarP, typename IdxP>
    class BusBase;

    template <typename RealP, typename IdxP>
    struct LoadZIPData;
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename, typename>
    class LoadZIP;

    enum class LoadZIPInternalVariables : size_t
    {
      MAXIMUM
    };

    enum class LoadZIPExternalVariables : size_t
    {
      MAXIMUM
    };

    template <typename ScalarP, typename IdxP>
    struct ConnectedElementTraits<LoadZIP<ScalarP, IdxP>>
    {
      using LoadZIPT = LoadZIP<ScalarP, IdxP>;

      using ElementT           = LoadZIPT;
      using ScalarT            = ScalarP;
      using IdxT               = IdxP;
      using RealT              = typename ScalarTraits<ScalarT>::RealT;
      using ModelDataT         = LoadZIPData<RealT, IdxT>;
      using InternalVariablesT = LoadZIPInternalVariables;
      using ExternalVariablesT = LoadZIPExternalVariables;
      using InterfaceT         = Component<ScalarT, IdxT>;
    };

    /*!
     * @brief Implementation of a ZIP load.
     *
     */
    template <typename ScalarP, typename IdxP>
    class LoadZIP : public ConnectedElement<LoadZIP<ScalarP, IdxP>>
    {
      using ConnectedElement<LoadZIP>::gridkit_component_id_;
      using ConnectedElement<LoadZIP>::size_;
      using ConnectedElement<LoadZIP>::nnz_;
      using ConnectedElement<LoadZIP>::time_;
      using ConnectedElement<LoadZIP>::alpha_;
      using ConnectedElement<LoadZIP>::y_;
      using ConnectedElement<LoadZIP>::yp_;
      using ConnectedElement<LoadZIP>::tag_;
      using ConnectedElement<LoadZIP>::wb_;
      using ConnectedElement<LoadZIP>::h_;
      using ConnectedElement<LoadZIP>::f_;
      using ConnectedElement<LoadZIP>::J_;
      using ConnectedElement<LoadZIP>::J_rows_buffer_;
      using ConnectedElement<LoadZIP>::J_cols_buffer_;
      using ConnectedElement<LoadZIP>::J_vals_buffer_;
      using ConnectedElement<LoadZIP>::variable_indices_;
      using ConnectedElement<LoadZIP>::residual_indices_;

    public:
      using ScalarT    = typename ConnectedElement<LoadZIP>::ScalarT;
      using IdxT       = typename ConnectedElement<LoadZIP>::IdxT;
      using RealT      = typename ConnectedElement<LoadZIP>::RealT;
      using BusT       = BusBase<ScalarT, IdxT>;
      using ModelDataT = LoadZIPData<RealT, IdxT>;

      LoadZIP(BusT* bus);
      LoadZIP(BusT* bus, RealT P0, RealT Q0, RealT V0, RealT alphaI, RealT alphaP);
      LoadZIP(BusT* bus, const ModelDataT& data);
      ~LoadZIP();

      int setGridKitComponentID(IdxT) override final;
      int allocate() override final;
      int initialize() override final;
      int tagDifferentiable() override final;
      int evaluateResidual() override final;
      int evaluateJacobian() override final;

      int verify() const override final
      {
        return 0;
      }

    public:
      void setP0(RealT P0)
      {
        P0_ = P0;
      }

      void setQ0(RealT Q0)
      {
        Q0_ = Q0;
      }

      void setV0(RealT V0)
      {
        V0_ = V0;
      }

      void setalphaI(RealT alphaI)
      {
        alphaI_ = alphaI;
      }

      void setalphaP(RealT alphaP)
      {
        alphaP_ = alphaP;
      }

    private:
      void setDerivedParams();

      ScalarT& Vr()
      {
        return bus_->Vr();
      }

      ScalarT& Vi()
      {
        return bus_->Vi();
      }

      ScalarT& Ir()
      {
        return bus_->Ir();
      }

      ScalarT& Ii()
      {
        return bus_->Ii();
      }

    public:
      __attribute__((always_inline)) inline int evaluateBusResidual(ScalarT*, ScalarT*, ScalarT*, ScalarT*);
      __attribute__((always_inline)) inline int evaluateInternalResidual(ScalarT*, ScalarT*, ScalarT*, ScalarT*);

    private:
      BusT* bus_{nullptr};
      RealT P0_{0};
      RealT Q0_{0};
      RealT V0_{1.0};
      RealT alphaI_{0};
      RealT alphaP_{0};
    };

  } // namespace PhasorDynamics
} // namespace GridKit
