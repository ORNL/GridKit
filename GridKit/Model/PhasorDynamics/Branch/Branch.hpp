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
#include <GridKit/Definitions.hpp>

// Forward declarations.
namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename scalar_type, typename index_type>
    class BusBase;

    template <typename real_type, typename index_type>
    struct BranchData;
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Implementation of a line or off-nominal transformer branch between two buses.
     *
     * The model is implemented in Cartesian coordinates. Positive current
     * direction is into the busses.
     *
     */
    template <typename scalar_type, typename index_type>
    class Branch : public Component<scalar_type, index_type>
    {
      using Component<scalar_type, index_type>::gridkit_component_id_;
      using Component<scalar_type, index_type>::size_;
      using Component<scalar_type, index_type>::nnz_;
      using Component<scalar_type, index_type>::time_;
      using Component<scalar_type, index_type>::alpha_;
      using Component<scalar_type, index_type>::y_;
      using Component<scalar_type, index_type>::yp_;
      using Component<scalar_type, index_type>::tag_;
      using Component<scalar_type, index_type>::f_;
      using Component<scalar_type, index_type>::wb_;
      using Component<scalar_type, index_type>::h_;
      using Component<scalar_type, index_type>::J_rows_buffer_;
      using Component<scalar_type, index_type>::J_cols_buffer_;
      using Component<scalar_type, index_type>::J_vals_buffer_;
      using Component<scalar_type, index_type>::variable_indices_;
      using Component<scalar_type, index_type>::residual_indices_;
      using Component<scalar_type, index_type>::allocated_;

    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename Component<ScalarT, IdxT>::RealT;
      using BusT       = BusBase<ScalarT, IdxT>;
      using ModelDataT = BranchData<RealT, IdxT>;
      using MonitorT   = Model::VariableMonitor<Branch, BranchData>;

      Branch(BusT* bus1, BusT* bus2);
      Branch(BusT* bus1,
             BusT* bus2,
             RealT R,
             RealT X,
             RealT G,
             RealT B,
             RealT tap   = 1.0,
             RealT phase = 0.0);
      Branch(BusT* bus1, BusT* bus2, const ModelDataT& data);
      virtual ~Branch();

      virtual int setGridKitComponentID(IdxT) override final;
      virtual int allocate() override final;
      virtual int initialize() override final;
      virtual int tagDifferentiable() override final;
      virtual int setAbsoluteTolerance(RealT rel_tol) override final;
      virtual int evaluateResidual() override final;
      virtual int evaluateJacobian() override final;
      virtual int verify() const override final;

      void setR(RealT R)
      {
        R_ = R;
        setDerivedParams();
      }

      void setX(RealT X)
      {
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

      void setTap(RealT tap)
      {
        tap_ = tap;
        setDerivedParams();
      }

      void setPhase(RealT phase)
      {
        phase_ = phase;
        setDerivedParams();
      }

      const Model::VariableMonitorBase* getMonitor() const override;

    private:
      void initializeParameters(const ModelDataT& data);
      void initializeMonitor();
      void setDerivedParams();
      void terminalCurrent1(ScalarT& Ir, ScalarT& Ii);
      void terminalCurrent2(ScalarT& Ir, ScalarT& Ii);
      bool readRealParameter(const ModelDataT&               data,
                             typename ModelDataT::Parameters parameter,
                             RealT&                          target);

      FORCE_INLINE static void addAdmittanceContribution(const RealT   G,
                                                         const RealT   B,
                                                         const ScalarT Vr,
                                                         const ScalarT Vi,
                                                         ScalarT&      Ir,
                                                         ScalarT&      Ii);

      FORCE_INLINE static void evaluateAdmittanceBlock(const RealT    G,
                                                       const RealT    B,
                                                       const ScalarT* wb,
                                                       ScalarT*       h);

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
      FORCE_INLINE int evaluateBusResidual11(
          const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);
      FORCE_INLINE int evaluateBusResidual12(
          const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);
      FORCE_INLINE int evaluateBusResidual21(
          const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);
      FORCE_INLINE int evaluateBusResidual22(
          const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

    private:
      BusT* bus1_;
      BusT* bus2_;
      RealT R_{0.0};
      RealT X_{0.0};
      RealT G_{0.0};
      RealT B_{0.0};
      RealT Gmag_{0.0};
      RealT Bmag_{0.0};
      RealT tap_{1.0};
      RealT phase_{0.0};
      IdxT  bus1_id_{0};
      IdxT  bus2_id_{0};

      RealT g11_{0.0};
      RealT b11_{0.0};
      RealT g12_{0.0};
      RealT b12_{0.0};
      RealT g21_{0.0};
      RealT b21_{0.0};
      RealT g22_{0.0};
      RealT b22_{0.0};

      int parameter_error_count_{0};

      /// Variable monitor
      std::unique_ptr<MonitorT> monitor_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
