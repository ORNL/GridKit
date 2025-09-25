/**
 * @file Branch.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Declaration of a phasor dynamics branch model.
 *
 * The model uses Cartesian coordinates.
 *
 */
#pragma once

#include <Model/PhasorDynamics/Component.hpp>
#include <Model/PhasorDynamics/ComponentSignals.hpp>

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

      using real_type       = typename Component<ScalarT, IdxT>::real_type;
      using bus_type        = BusBase<ScalarT, IdxT>;
      using model_data_type = BranchData<real_type, IdxT>;

    public:
      Branch(bus_type* bus1, bus_type* bus2);
      Branch(bus_type* bus1, bus_type* bus2, real_type R, real_type X, real_type G, real_type B);
      Branch(bus_type* bus1, bus_type* bus2, const model_data_type& data);
      virtual ~Branch();

      virtual int setGridKitComponentID(IdxT) override;
      virtual int allocate() override;
      virtual int initialize() override;
      virtual int tagDifferentiable() override;
      virtual int evaluateResidual() override;
      virtual int evaluateJacobian() override;

      virtual void updateTime(real_type /* t */, real_type /* a */) override
      {
      }

    public:
      void setR(real_type R)
      {
        R_ = R;
      }

      void setX(real_type X)
      {
        // std::cout << "Setting X ...\n";
        X_ = X;
      }

      void setG(real_type G)
      {
        G_ = G;
      }

      void setB(real_type B)
      {
        B_ = B;
      }

    private:
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

    private:
      bus_type* bus1_;
      bus_type* bus2_;
      real_type R_{0.0};
      real_type X_{0.0};
      real_type G_{0.0};
      real_type B_{0.0};
      IdxT      bus1_id_{0};
      IdxT      bus2_id_{0};
    };

  } // namespace PhasorDynamics
} // namespace GridKit
