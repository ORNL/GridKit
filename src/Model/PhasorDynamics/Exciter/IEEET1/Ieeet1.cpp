/**
 * @file Ieeet1.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * 
 * @brief Definition of a IEEET1 Exciter.
 *
 */

#include <cmath>
#include <iostream>


#include <Model/PhasorDynamics/Bus/Bus.hpp>

#include <Model/PhasorDynamics/Exciter/IEEET1/Ieeet1Data.hpp>
#include <Model/PhasorDynamics/SignalNode/SignalNode.hpp>

#include "Ieeet1.hpp"

#define _USE_MATH_DEFINES

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {

      /**
       * @brief  Constructor for IEEET1 Exciter
       *
       * @param data  Data object to store parameters
       * @param bus   Signal used for terminal reference vmag
       * @param speed Signal used for machine relative speed
       * @param efd   Signal used for E field
       * @tparam ScalarT Scalar data type
       * @tparam IdxT Index data type
       */
      template <class ScalarT, typename IdxT>
      Ieeet1<ScalarT, IdxT>::Ieeet1(          
          signal_type* efd_signal,
          signal_type* speed_signal,
          bus_type* bus,
          const model_data_type& data)
      :   efd_signal_(efd_signal),
          speed_signal_(speed_signal),
          bus_(bus),
          Tr_(data.Tr),  
          Ka_(data.Ka),  
          Ta_(data.Ta),
          Ke_(data.Ke),  
          Te_(data.Te), 
          Kf_(data.Kf),
          Tf_(data.Tf),  
          Vrmin_(data.Vrmin),   
          Vrmax_(data.Vrmax), 
          E1_(data.E1), 
          E2_(data.E2),  
          Se1_(data.Se1),  
          Se2_(data.Se2),  
          Ispdlim_(data.Ispdlim)
      {

        // 9 Internal Variables
        size_ = 9;

        // Derived Parameters
        ScalarT SR = std::sqrt( Se2_ / Se1_);

        // Technically two solutions
        
        // Solution 1 (Aligned with PW)
        SA_ = (SR * E1_ - E2_) / (SR - 1);
        SB_ = Se1_ / (E1_ - SA_) / (E1_ - SA_);

        // Solution 2
        // SA_ = (SR * E1_ + E2_) / (SR + 1);
        // SB_ = Se1_ / (E1_ - SA_) / (E1_ - SA_);
      }


      /**
       * @brief Allocate memory for model
       *
       */
      template <class ScalarT, typename IdxT>
      int Ieeet1<ScalarT, IdxT>::allocate()
      {
        auto size = static_cast<size_t>(size_); // avoid compiler warnings
        f_.resize(size);
        y_.resize(size);
        yp_.resize(size);
        tag_.resize(size);

        // Set output signal after allocation
        // The signal is accessible to the generator
        // y_[7] is output efd state
        if (efd_signal_)
        {
          efd_signal_->set(&y_[7]);
        }
        return 0;
      }

      /**
       * @brief Initialization of the Exciter
       * 
       * Sets/configures all of the initial values of the exciter
       * by assuming no saturation and steady-state.
       *
       */
      template <class ScalarT, typename IdxT>
      int Ieeet1<ScalarT, IdxT>::initialize()
      {

        // Read External Variables

        ScalarT efd0{0};
        if (efd_signal_)
        {
          efd0 = y_[7]; //<- TODO generator needs to set efd initial value
        }

        // Terminal Voltage
        ScalarT vreal = bus_->Vr();
        ScalarT vimag = bus_->Vi();
        Ec_  = std::sqrt(vreal*vreal + vimag*vimag);


        // Derived from External initial values
        ScalarT vr    = Ke_       * efd0;        
        ScalarT vfx   = Kf_ / Tf_ * efd0;  
        ScalarT vtr   = Ke_ / Ka_ * efd0;          

        // Vref (setpoint = terminal + error)
        vref_ = Ec_ + vtr;

        // IVP for Internal Variables
        y_[0] = Ec_;  // y0 - vts  - Sensed term volt 
        y_[1] = vr;   // y1 - vr   - Votlage reg
        y_[2] = efd0; // y2 - efdp - Efd pre mult
        y_[3] = vfx;  // y3 - vfx  - Exciter feedback
        y_[4] = vtr;  // y4 - vtr  - Term Volt Err
        y_[5] = 0;    // y5 - vf   - Feedback volt
        y_[6] = 0;    // y6 - ve   - Excit. Cntrl Volt
        y_[7] = efd0; // y7 - efd  - Efd
        y_[8] = 0;    // y8 - ksat - Saturation

        // Steady State Conditions
        yp_[0] = 0.0; 
        yp_[1] = 0.0; 
        yp_[2] = 0.0; 
        yp_[3] = 0.0; 
        yp_[4] = 0.0; 
        yp_[5] = 0.0;
        yp_[6] = 0.0; 
        yp_[7] = 0.0;
        yp_[8] = 0.0;

        return 0;
      }

      /**
       * @brief  Identify differential variables.
       *
       * @tparam ScalarT Scalar data type
       * @tparam IdxT Index data type
       * @return int 0
       */
      template <class ScalarT, typename IdxT>
      int Ieeet1<ScalarT, IdxT>::tagDifferentiable()
      {

        tag_[0] = true;  // y0 - vts  - Sensed term volt 
        tag_[1] = true;  // y1 - vr   - Votlage reg
        tag_[2] = true;  // y2 - efdp - Efd pre mult
        tag_[3] = true;  // y3 - vfx  - Exciter feedback
        tag_[4] = false; // y4 - vtr  - Term Volt Err
        tag_[5] = false; // y5 - vf   - Feedback volt
        tag_[6] = false; // y6 - ve   - Excit. Cntrl Volt
        tag_[7] = false; // y7 - efd  - Efd
        tag_[8] = false; // y8 - ksat - Saturation

        return 0;
      }

      /**
       * @brief  Scaled sigmoid activation function
       *
       * @param[in] x State variable
       * @tparam ScalarT Scalar data type
       * @tparam IdxT Index data type
       * @return Sigmoid approximation evaluated at x
       *
       * @warning This needs to be abstracted to be used
       *          across phasor dynamics. Identical pattern is
       *          being used in TGOV1 model.
       * 
       *   Algebraic approximation of transcendental sigmoid.
       */
      template <class ScalarT, typename IdxT>
      ScalarT Ieeet1<ScalarT, IdxT>::sigmoid(ScalarT x)
      {
        return ((0.5 * mu_ * x) / (1.0 + std::abs(mu_ * x))) + 0.5;
      }

      /**
       * @brief  Indicator function for lower valve limit violation.
       *
       * @param[in] x State variable
       * @param[in] f Conditional derivative of state variable
       * @tparam ScalarT Scalar data type
       * @tparam IdxT Index data type
       * @return Lower Smooth Indicator value at coordinates (x, f(x))
       *
       * @warning This needs to be abstracted to be used
       *          across phasor dynamics. Identical pattern is
       *          being used in TGOV1 model.
       */
      template <class ScalarT, typename IdxT>
      ScalarT Ieeet1<ScalarT, IdxT>::indicator_low(ScalarT x, ScalarT f)
      {
        return (this->sigmoid(Vrmin_ - x)) * (this->sigmoid(-f));
      }

      /**
       * @brief Indicator function for high valve limit violation.
       *
       * @param[in] x State variable
       * @param[in] f Conditional derivative of state variable
       * @tparam ScalarT Scalar data type
       * @tparam IdxT Index data type
       * @return Higher Smooth Indicator value at coordinates (x, f(x))
       *
       * @warning This needs to be abstracted to be used
       *          across phasor dynamics. Identical pattern is
       *          being used in TGOV1 model.
       */
      template <class ScalarT, typename IdxT>
      ScalarT Ieeet1<ScalarT, IdxT>::indicator_high(ScalarT x, ScalarT f)
      {
        return (this->sigmoid(x - Vrmax_)) * (this->sigmoid(f));
      }

      /**
       * @brief Net Indicator function for regulator limits
       *
       * @param[in] x State variable
       * @param[in] f Conditional derivative of state variable
       * @tparam ScalarT Scalar data type
       * @tparam IdxT Index data type
       * @return <DESCRIPTION OF RETURN VALUE>
       *
       * @warning This needs to be abstracted to be used
       *          across phasor dynamics. Identical pattern is
       *          being used in TGOV1 model.
       */
      template <class ScalarT, typename IdxT>
      ScalarT Ieeet1<ScalarT, IdxT>::indicator(ScalarT x, ScalarT f)
      {
        return (1 - this->indicator_low(x, f)) * (1 - this->indicator_high(x, f));
      }

      /**
       * @brief Residuals of system equations
       *
       */
      template <class ScalarT, typename IdxT>
      int Ieeet1<ScalarT, IdxT>::evaluateResidual()
      {

        // Input Variables
        ScalarT omega{0};
        if (speed_signal_)
        {
          omega = speed_signal_->read();
        }

        // NOTE
        // No external variable 'write' needed,
        // as it is available to others dierctly throguh 
        // pointer.

        // Read E comp (terminal voltage, unless compensation impedence)
        ScalarT vreal = bus_->Vr();
        ScalarT vimag = bus_->Vi();
        Ec_  = std::sqrt(vreal*vreal + vimag*vimag);

        // Read Internal Variables
        ScalarT vts   = y_[0]; // y0 - Sensed term volt
        ScalarT vr    = y_[1]; // y1 - Votlage reg
        ScalarT efdp  = y_[2]; // y2 - Efd pre mult
        ScalarT vfx   = y_[3]; // y3 - Exciter feedback
        ScalarT vtr   = y_[4]; // y4 - Term Volt Err
        ScalarT vf    = y_[5]; // y5 - Feedback volt
        ScalarT ve    = y_[6]; // y6 - Excit. Cntrl Volt
        ScalarT efd   = y_[7]; // y7 - Efd (EXTERNALLY ACCESSABLE)
        ScalarT ksat  = y_[8]; // y8 - Saturation

        // Read Internal Derivatives
        ScalarT vts_dot  = yp_[0];
        ScalarT vr_dot   = yp_[1];
        ScalarT efdp_dot = yp_[2];
        ScalarT vfx_dot  = yp_[3];

        // The 'pre-limit' derivative of Pv
        ScalarT f        = -vr + Ka_ * vtr;
        ScalarT vr_ind   = this->indicator(vr, f);

        // Internal Differential Equations
        f_[0] = -vts_dot  + (Ec_ - vts) / Tr_; 
        f_[1] = -vr_dot   + vr_ind * f  / Ta_;
        f_[2] = -efdp_dot + (vr - ve - Ke_ * efdp) / Te_;
        f_[3] = -vfx_dot  + vf / Tf_;

        // Internal Algebraic Equations
        f_[4] = -vtr      + vref_ + vUEL_ + vOEL_ + vS_ - vf - vts;
        f_[5] = -vf       + (efdp * Kf_) / Tf_ - vfx;
        f_[6] = -ve       + ksat * efdp;
        f_[7] = -efd      + efdp + omega * efdp * Ispdlim_;
        
        // TODO smooth approaximation for Enzyme
        // NOTE seems about double PW saturation.
        f_[8] = -ksat;
        if (efdp > SA_)
          f_[8] += SB_ * (efdp - SA_) * (efdp - SA_); //  / 2 ?


        return 0;
      }

      /**
       * @brief Jacobian evaluation not implemented yet
       *
       * @tparam ScalarT - Scalar data type
       * @tparam IdxT    - Index data type
       * @return int - error code, 0 = success
       */
      template <class ScalarT, typename IdxT>
      int Ieeet1<ScalarT, IdxT>::evaluateJacobian()
      {
        std::cout << "Jacobian evaluation not implemented!" << std::endl;
        return 0;
      }

      // Available template instantiations
      template class Ieeet1<double, long int>;
      template class Ieeet1<double, size_t>;

    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit