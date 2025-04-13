/* GENROU Component - Adam Birchfield */
#pragma once

#define _USE_MATH_DEFINES
#include <Model/PhasorDynamics/BusBase.hpp>
#include <Model/PhasorDynamics/Component.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using ComponentT = Component<double, size_t>;
    using BaseBusT   = BusBase<double, size_t>;

    class ClassicalGen : public ComponentT
    {
      using ComponentT::alpha_;
      using ComponentT::f_;
      using ComponentT::fB_;
      using ComponentT::g_;
      using ComponentT::gB_;
      using ComponentT::nnz_;
      using ComponentT::param_;
      using ComponentT::size_;
      using ComponentT::tag_;
      using ComponentT::time_;
      using ComponentT::y_;
      using ComponentT::yB_;
      using ComponentT::yp_;
      using ComponentT::ypB_;

    public:
      ClassicalGen(BaseBusT* bus, int unit_id, double pmech, double ep)
        : bus_(bus),
          unit_id_(unit_id),
          busID_(0),
          H_(3),
          D_(0),
          Ra_(0),
          Xdp_(0.2),
          pmech_(pmech),
          ep_(ep)
      {
        size_ = 5;
        set_derived_params();
      }

      ClassicalGen(BaseBusT* bus,
             int       unit_id,
             double    H,
             double    D,
             double    Ra,
             double    Xdp,
             double    pmech,
             double    ep)
        : bus_(bus),
          unit_id_(unit_id),
          busID_(0),
          H_(H),
          D_(D),
          Ra_(Ra),
          Xdp_(Xdp),
          pmech_(pmech),
          ep_(ep)
      {
        size_ = 5;
        set_derived_params();
      }

      void set_derived_params()
      {
        g  = Ra_ / (Ra_ * Ra_ + Xdp_ * Xdp_);
        b   = Xdp_ / (Ra_ * Ra_ + Xdp_ * Xdp_);
      }

      ~ClassicalGen()
      {
      }

      int allocate() override
      {
        f_.resize(size_);
        y_.resize(size_);
        yp_.resize(size_);
        tag_.resize(size_);
        fB_.resize(size_);
        yB_.resize(size_);
        ypB_.resize(size_);
        return 0;
      }

      int initialize() override
      {
        return  0;
      }

      int tagDifferentiable() override
      {
        return 0;
      }

      int evaluateResidual() override
      {
        /* Read variables */
        double delta  = y_[0];
        double omega  = y_[1];
        double telec  = y_[2];
        double ir     = y_[3];
        double ii     = y_[4];

        /* Read derivatives */
        double delta_dot = yp_[0];
        double omega_dot = yp_[1];

        /* 6 GENROU differential equations */
        f_[0] = delta_dot - omega * (2 * M_PI * 60);
        f_[1] = omega_dot - (1.0 / (2 * H_)) * ((pmech_ - D_ * omega) / (1 + omega) - telec);
        
        /* 11 GENROU algebraic equations */
        f_[2] = telec - (1.0/(1.0 + omega))*(g*ep_*ep_ - ep_*(cos(delta)*(g*Vr() - b*Vi()) + sin(delta)*(b*Vr() + g*Vi())));

        f_[3] = ir + g*Vr() - b * Vi()  - ep_*(g*cos(delta) -b*sin(delta));
        f_[4] = ii + b*Vr() +  g * Vi() - ep_*(b*cos(delta) + g*sin(delta));

    
        /* Current balance */
        Ir() += - (g*Vr() - b * Vi()  - ep_*(g*cos(delta) -b*sin(delta)));
        Ii() += - (b*Vr() +  g * Vi() - ep_*(b*cos(delta) + g*sin(delta)));

        // printf("GENROU residual\n");
        // for (int i = 0 ; i < 21; ++i) printf("%d: %g\n", i, f_[i]);

        // printf("GENROU inr %g Vr %g B %g Vi %g G %g\n", inr, Vr(), B_, Vi(), G_);
        // printf("GENROU Ii = %g\n", inr - Vr()*B_ - Vi()*G_);

        return 0;
      }

      int evaluateJacobian() override
      {
        /* TODO */
        return 0;
      }

      /* Don't know what to do with any of these */
      int evaluateIntegrand() override
      {
        return 0;
      }

      int initializeAdjoint() override
      {
        return 0;
      }

      int evaluateAdjointResidual() override
      {
        return 0;
      }

      int evaluateAdjointIntegrand() override
      {
        return 0;
      }

      void updateTime(double t, double a) override
      {
      }

    private:
      double& Vr()
      {
        return bus_->Vr();
      }

      double& Vi()
      {
        return bus_->Vi();
      }

      double& Ir()
      {
        return bus_->Ir();
      }

      double& Ii()
      {
        return bus_->Ii();
      }

    private:
      /* Identification */
      BaseBusT* bus_;
      const int busID_;
      int       unit_id_;

      /* Input parameters */ 
      double H_;
      double D_;
      double Ra_;
      double Xdp_;

      /* Derivied parameters */
      double g;
      double b;

      /* Setpoints for control variables (determined at initialization) */
      double pmech_;
      double ep_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit