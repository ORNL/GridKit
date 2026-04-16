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
#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

// Forward declarations.
namespace GridKit
{
  namespace PhasorDynamics
  {
    template <class ScalarT, typename IdxT>
    class BusBase;

    template <typename RealT, typename IdxT>
    struct BranchData;
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Implementation of a pi-model branch between two buses.
     *
     * The model is implemented in Cartesian coordinates. Positive current
     * direction is into the busses.
     *
     */
    template <class ScalarT, typename IdxT>
    class Branch : public Component<ScalarT, IdxT>
    {
      using Component<ScalarT, IdxT>::gridkit_component_id_;
      using Component<ScalarT, IdxT>::size_;
      using Component<ScalarT, IdxT>::nnz_;
      using Component<ScalarT, IdxT>::time_;
      using Component<ScalarT, IdxT>::alpha_;
      using Component<ScalarT, IdxT>::y_;
      using Component<ScalarT, IdxT>::yp_;
      using Component<ScalarT, IdxT>::tag_;
      using Component<ScalarT, IdxT>::f_;
      using Component<ScalarT, IdxT>::wb_;
      using Component<ScalarT, IdxT>::h_;
      using Component<ScalarT, IdxT>::J_;
      using Component<ScalarT, IdxT>::J_rows_buffer_;
      using Component<ScalarT, IdxT>::J_cols_buffer_;
      using Component<ScalarT, IdxT>::J_vals_buffer_;
      using Component<ScalarT, IdxT>::variable_indices_;
      using Component<ScalarT, IdxT>::residual_indices_;

    public:
      using RealT           = typename Component<ScalarT, IdxT>::RealT;
      using bus_type        = BusBase<ScalarT, IdxT>;
      using model_data_type = BranchData<RealT, IdxT>;
      using MonitorT        = Model::VariableMonitor<Branch, BranchData>;

      Branch(bus_type* bus1, bus_type* bus2);
      Branch(bus_type* bus1, bus_type* bus2, RealT R, RealT X, RealT G, RealT B);
      Branch(bus_type* bus1, bus_type* bus2, const model_data_type& data);
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
      void initializeParameters(const model_data_type& data);
      void initializeMonitor();
      void setDerivedParams();

      ScalarT& Vr1()
      {
        return bus1_->v(0);
      }

      ScalarT& Vi1()
      {
        return bus1_->v(1);
      }

      ScalarT& Ir1()
      {
        return bus1_->i(0);
      }

      ScalarT& Ii1()
      {
        return bus1_->i(1);
      }

      ScalarT& Vr2()
      {
        return bus2_->v(0);
      }

      ScalarT& Vi2()
      {
        return bus2_->v(1);
      }

      ScalarT& Ir2()
      {
        return bus2_->i(0);
      }

      ScalarT& Ii2()
      {
        return bus2_->i(1);
      }

    public:
      __attribute__((always_inline)) inline int evaluateBusResidual11(ScalarT*, ScalarT*, ScalarT*, ScalarT*);
      __attribute__((always_inline)) inline int evaluateBusResidual12(ScalarT*, ScalarT*, ScalarT*, ScalarT*);
      __attribute__((always_inline)) inline int evaluateBusResidual21(ScalarT*, ScalarT*, ScalarT*, ScalarT*);
      __attribute__((always_inline)) inline int evaluateBusResidual22(ScalarT*, ScalarT*, ScalarT*, ScalarT*);

    private:
      bus_type* bus1_;
      bus_type* bus2_;
      RealT     R_{0.0};
      RealT     X_{0.0};
      RealT     G_{0.0};
      RealT     B_{0.0};
      IdxT      bus1_id_{0};
      IdxT      bus2_id_{0};

      /* Derived parameters */
      RealT b_;
      RealT g_;

      /// Variable monitor
      std::unique_ptr<MonitorT> monitor_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
