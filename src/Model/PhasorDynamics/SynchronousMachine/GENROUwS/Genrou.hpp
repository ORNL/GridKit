/**
 * @file Genrou.cpp
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Declaration of a GENROU generator model.
 *
 */

#pragma once

#include <Model/PhasorDynamics/Component.hpp>

// Forward declarations.
namespace GridKit
{
  namespace PhasorDynamics
  {
    template <class ScalarT, typename IdxT>
    class BusBase;

    template <typename RealT, typename IdxT>
    struct GenrouData;
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {

    template <class ScalarT, typename IdxT>
    class Genrou : public Component<ScalarT, IdxT>
    {
      using Component<ScalarT, IdxT>::alpha_;
      using Component<ScalarT, IdxT>::f_;
      using Component<ScalarT, IdxT>::nnz_;
      using Component<ScalarT, IdxT>::size_;
      using Component<ScalarT, IdxT>::tag_;
      using Component<ScalarT, IdxT>::time_;
      using Component<ScalarT, IdxT>::y_;
      using Component<ScalarT, IdxT>::yp_;

      using real_type       = typename Component<ScalarT, IdxT>::real_type;
      using bus_type        = BusBase<ScalarT, IdxT>;
      using model_data_type = GenrouData<real_type, IdxT>;

    public:
      Genrou(bus_type* bus, IdxT unit_id);
      Genrou(bus_type*              bus,
             bus_type*              pmech_signal,
             bus_type*              omega_signal,
             const model_data_type& data);

      // TODO remove, becoming hard to manage
      Genrou(bus_type* bus,
             IdxT      unit_id,
             ScalarT   p0,
             ScalarT   q0,
             real_type H,
             real_type D,
             real_type Ra,
             real_type Tdop,
             real_type Tdopp,
             real_type Tqopp,
             real_type Tqop,
             real_type Xd,
             real_type Xdp,
             real_type Xdpp,
             real_type Xq,
             real_type Xqp,
             real_type Xqpp,
             real_type Xl,
             real_type S10,
             real_type S12);
      ~Genrou() = default;

      int allocate() override;
      int initialize() override;
      int tagDifferentiable() override;
      int evaluateResidual() override;

      // Still to be implemented
      int evaluateJacobian() override;

      void updateTime(real_type /* t */, real_type /* a */) override
      {
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

    private:
      /* Identification */
      bus_type* bus_;
      bus_type* pmech_signal_;
      bus_type* omega_signal_;

      IdxT busID_{0};
      IdxT unit_id_; //< @todo this should be removed

      /* Initial terminal conditions */
      ScalarT p0_{0.0};
      ScalarT q0_{0.0};

      /* Input parameters */
      real_type H_{0.0};
      real_type D_{0.0};
      real_type Ra_{0.0};
      real_type Tdop_{0.0};
      real_type Tdopp_{0.0};
      real_type Tqopp_{0.0};
      real_type Tqop_{0.0};
      real_type Xd_{0.0};
      real_type Xdp_{0.0};
      real_type Xdpp_{0.0};
      real_type Xq_{0.0};
      real_type Xqp_{0.0};
      real_type Xqpp_{0.0};
      real_type Xl_{0.0};
      real_type S10_{0.0};
      real_type S12_{0.0};

      /* Derivied parameters */
      real_type SA_;
      real_type SB_;
      real_type Xd1_;
      real_type Xd2_;
      real_type Xd3_;
      real_type Xd4_;
      real_type Xd5_;
      real_type Xq1_;
      real_type Xq2_;
      real_type Xq3_;
      real_type Xq4_;
      real_type Xq5_;
      real_type Xqd_;
      real_type G_;
      real_type B_;

      /* Setpoints for control variables (determined at initialization) */
      real_type pmech_set_;
      real_type efd_set_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
