#pragma once

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ConnectedElement.hpp>
#include <GridKit/Model/PhasorDynamics/Load/LoadData.hpp>

// Forward declarations.
namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename ScalarP, typename IdxP>
    class BusBase;

    template <typename RealP, typename IdxP>
    struct LoadData;
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename, typename>
    class Load;

    enum class LoadInternalVariables : size_t
    {
      MAXIMUM
    };

    enum class LoadExternalVariables : size_t
    {
      MAXIMUM
    };

    template <typename ScalarP, typename IdxP>
    struct ConnectedElementTraits<Load<ScalarP, IdxP>>
    {
      using LoadT = Load<ScalarP, IdxP>;

      using ElementT           = LoadT;
      using ScalarT            = ScalarP;
      using IdxT               = IdxP;
      using RealT              = typename ScalarTraits<ScalarT>::RealT;
      using ModelDataT         = LoadData<RealT, IdxT>;
      using InternalVariablesT = LoadInternalVariables;
      using ExternalVariablesT = LoadExternalVariables;
      using InterfaceT         = Component<ScalarT, IdxT>;
    };

    /*!
     * @brief Implementation of a constant load.
     *
     */
    template <typename ScalarP, typename IdxP>
    class Load : public ConnectedElement<Load<ScalarP, IdxP>>
    {
      using ConnectedElement<Load>::gridkit_component_id_;
      using ConnectedElement<Load>::size_;
      using ConnectedElement<Load>::nnz_;
      using ConnectedElement<Load>::time_;
      using ConnectedElement<Load>::alpha_;
      using ConnectedElement<Load>::y_;
      using ConnectedElement<Load>::yp_;
      using ConnectedElement<Load>::tag_;
      using ConnectedElement<Load>::wb_;
      using ConnectedElement<Load>::h_;
      using ConnectedElement<Load>::J_rows_buffer_;
      using ConnectedElement<Load>::J_cols_buffer_;
      using ConnectedElement<Load>::J_vals_buffer_;
      using ConnectedElement<Load>::variable_indices_;
      using ConnectedElement<Load>::residual_indices_;

    public:
      using ScalarT    = typename ConnectedElement<Load>::ScalarT;
      using IdxT       = typename ConnectedElement<Load>::IdxT;
      using RealT      = typename ConnectedElement<Load>::RealT;
      using BusT       = BusBase<ScalarT, IdxT>;
      using ModelDataT = LoadData<RealT, IdxT>;
      using MonitorT   = typename ConnectedElement<Load>::MonitorT;

      Load(BusT* bus);
      Load(BusT* bus, RealT R, RealT X);
      Load(BusT* bus, const ModelDataT& data);
      virtual ~Load();

      virtual int setGridKitComponentID(IdxT) override final;
      virtual int allocate() override final;
      virtual int initialize() override final;
      virtual int tagDifferentiable() override final;
      virtual int evaluateResidual() override final;
      virtual int evaluateJacobian() override final;

      virtual int verify() const override final
      {
        return 0;
      }

    public:
      void setR(RealT R)
      {
        R_ = R;
      }

      void setX(RealT X)
      {
        // std::cout << "Setting X ...\n";
        X_ = X;
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

    private:
      BusT* bus_{nullptr};
      RealT R_{0.1};
      RealT X_{0.01};

      /* Derivied parameters */
      RealT b_;
      RealT g_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
