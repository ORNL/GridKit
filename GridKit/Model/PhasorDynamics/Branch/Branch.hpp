/**
 * @file Branch.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Declaration of a phasor dynamics branch model.
 *
 * The model uses Cartesian coordinates.
 *
 */
#pragma once

#include <GridKit/Model/PhasorDynamics/Branch/BranchData.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ConnectedElement.hpp>

// Forward declarations.
namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename ScalarP, typename IdxP>
    class BusBase;

    template <typename RealP, typename IdxP>
    struct BranchData;
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename, typename>
    class Branch;

    enum class BranchInternalVariables : size_t
    {
      MAXIMUM
    };

    enum class BranchExternalVariables : size_t
    {
      MAXIMUM
    };

    template <typename ScalarP, typename IdxP>
    struct ConnectedElementTraits<Branch<ScalarP, IdxP>>
    {
      using BranchT = Branch<ScalarP, IdxP>;

      using ElementT           = BranchT;
      using ScalarT            = ScalarP;
      using IdxT               = IdxP;
      using RealT              = typename ScalarTraits<ScalarT>::RealT;
      using ModelDataT         = BranchData<RealT, IdxT>;
      using InternalVariablesT = BranchInternalVariables;
      using ExternalVariablesT = BranchExternalVariables;
      using InterfaceT         = Component<ScalarT, IdxT>;
    };

    /**
     * @brief Implementation of a pi-model branch between two buses.
     *
     * The model is implemented in Cartesian coordinates. Positive current
     * direction is into the busses.
     *
     */
    template <typename ScalarP, typename IdxP>
    class Branch : public ConnectedElement<Branch<ScalarP, IdxP>>
    {
      using ConnectedElement<Branch>::gridkit_component_id_;
      using ConnectedElement<Branch>::size_;
      using ConnectedElement<Branch>::nnz_;
      using ConnectedElement<Branch>::time_;
      using ConnectedElement<Branch>::alpha_;
      using ConnectedElement<Branch>::y_;
      using ConnectedElement<Branch>::yp_;
      using ConnectedElement<Branch>::tag_;
      using ConnectedElement<Branch>::f_;
      using ConnectedElement<Branch>::wb_;
      using ConnectedElement<Branch>::h_;
      using ConnectedElement<Branch>::J_;
      using ConnectedElement<Branch>::J_rows_buffer_;
      using ConnectedElement<Branch>::J_cols_buffer_;
      using ConnectedElement<Branch>::J_vals_buffer_;
      using ConnectedElement<Branch>::variable_indices_;
      using ConnectedElement<Branch>::residual_indices_;
      using ConnectedElement<Branch>::monitor_;

    public:
      using ScalarT    = typename ConnectedElement<Branch>::ScalarT;
      using IdxT       = typename ConnectedElement<Branch>::IdxT;
      using RealT      = typename ConnectedElement<Branch>::RealT;
      using BusT       = BusBase<ScalarT, IdxT>;
      using ModelDataT = BranchData<RealT, IdxT>;
      using MonitorT   = typename ConnectedElement<Branch>::MonitorT;

      Branch(BusT* bus1, BusT* bus2);
      Branch(BusT* bus1, BusT* bus2, RealT R, RealT X, RealT G, RealT B);
      Branch(BusT* bus1, BusT* bus2, const ModelDataT& data);
      virtual ~Branch();

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

      void setR(RealT R)
      {
        R_ = R;
        setDerivedParams();
      }

      void setX(RealT X)
      {
        // std::cout << "Setting X ...\n";
        X_ = X;
        setDerivedParams();
      }

      void setG(RealT G)
      {
        G_ = G;
        setDerivedParams();
      }

      void setB(RealT B)
      {
        B_ = B;
        setDerivedParams();
      }

      const Model::VariableMonitorBase* getMonitor() const override;

    private:
      void initializeParameters(const ModelDataT& data);
      void initializeMonitor();
      void setDerivedParams();

      ScalarT& Vr1()
      {
        return bus1_->Vr();
      }

      ScalarT& Vi1()
      {
        return bus1_->Vi();
      }

      ScalarT& Ir1()
      {
        return bus1_->Ir();
      }

      ScalarT& Ii1()
      {
        return bus1_->Ii();
      }

      ScalarT& Vr2()
      {
        return bus2_->Vr();
      }

      ScalarT& Vi2()
      {
        return bus2_->Vi();
      }

      ScalarT& Ir2()
      {
        return bus2_->Ir();
      }

      ScalarT& Ii2()
      {
        return bus2_->Ii();
      }

    public:
      __attribute__((always_inline)) inline int evaluateBusResidual11(ScalarT*, ScalarT*, ScalarT*, ScalarT*);
      __attribute__((always_inline)) inline int evaluateBusResidual12(ScalarT*, ScalarT*, ScalarT*, ScalarT*);
      __attribute__((always_inline)) inline int evaluateBusResidual21(ScalarT*, ScalarT*, ScalarT*, ScalarT*);
      __attribute__((always_inline)) inline int evaluateBusResidual22(ScalarT*, ScalarT*, ScalarT*, ScalarT*);

    private:
      BusT* bus1_;
      BusT* bus2_;
      RealT R_{0.0};
      RealT X_{0.0};
      RealT G_{0.0};
      RealT B_{0.0};
      IdxT  bus1_id_{0};
      IdxT  bus2_id_{0};

      /* Derived parameters */
      RealT b_;
      RealT g_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
