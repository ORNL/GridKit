/**
 * @file Ieeet1Impl.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 *
 * @brief Definition of a IEEET1 Exciter.
 *
 */

#include <cmath>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/IEEET1/Ieeet1.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/IEEET1/Ieeet1Data.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

#define _USE_MATH_DEFINES

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      using Log = ::GridKit::Utilities::Logger;

      /**
       * @brief  Constructor for IEEET1 Exciter
       *
       * @tparam ScalarT Scalar data type
       * @tparam IdxT Index data type
       */
      template <class ScalarT, typename IdxT>
      Ieeet1<ScalarT, IdxT>::Ieeet1(bus_type* bus)
        : bus_(bus)
      {
        size_ = 9;
      }

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
      Ieeet1<ScalarT, IdxT>::Ieeet1(signal_type*           efd_signal,
                                    signal_type*           speed_signal,
                                    bus_type*              bus,
                                    const model_data_type& data)
        : efd_signal_(efd_signal),
          speed_signal_(speed_signal),
          bus_(bus),
          monitor_(std::make_unique<MonitorT>(data))
      {

        // Parse data struct into model
        this->initModelParams(data);

        initializeMonitor();

        // 9 Internal Variables
        size_ = 9;
      }

      /**
       * @brief  Constructor for IEEET1 Exciter
       *
       * @param bus   Signal used for terminal reference vmag
       * @param data  Data object to store parameters
       * @tparam ScalarT Scalar data type
       * @tparam IdxT Index data type
       */
      template <class ScalarT, typename IdxT>
      Ieeet1<ScalarT, IdxT>::Ieeet1(bus_type*              bus,
                                    const model_data_type& data)
        : bus_(bus),
          monitor_(std::make_unique<MonitorT>(data))
      {

        // Parse data struct into model
        this->initModelParams(data);

        initializeMonitor();

        // 9 Internal Variables
        size_ = 9;
      }

      template <class ScalarT, typename IdxT>
      Ieeet1<ScalarT, IdxT>::~Ieeet1()
      {
      }

      /**
       * @brief Set the component ID
       */
      template <class ScalarT, typename IdxT>
      int Ieeet1<ScalarT, IdxT>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      /**
       * @brief Allocate memory for model
       *
       */
      template <class ScalarT, typename IdxT>
      int Ieeet1<ScalarT, IdxT>::allocate()
      {
        // Resize component model data
        auto size = static_cast<size_t>(size_); // avoid compiler warnings
        f_.resize(size);
        y_.resize(size);
        yp_.resize(size);
        tag_.resize(size);
        abs_tol_.resize(size);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        // Resize bus data
        wb_.resize(2);

        // Resize signal variable data
        ws_.resize(1);
        ws_indices_.resize(1);
        ws_[0]         = 0.0;
        ws_indices_[0] = INVALID_INDEX<IdxT>;

        // Default variable and residual index mapping to local index
        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        // Set output signals
        if (signals_.template isAssigned<Ieeet1InternalVariables::EFD>())
        {
          signals_.template getSignalNode<Ieeet1InternalVariables::EFD>()->set(&y_[7], &(this->getVariableIndex(7)));
        }

        return 0;
      }

      /**
       * @brief verify method checks that attached signals are also linked
       */
      template <class ScalarT, typename IdxT>
      int Ieeet1<ScalarT, IdxT>::verify() const
      {
        static constexpr auto OMEGA = Ieeet1ExternalVariables::OMEGA;

        int ret = 0;

        if (signals_.template isAttached<OMEGA>())
        {
          if (!signals_.template isLinked<OMEGA>())
          {
            Log::error() << "Ieeet1: omega signal attached with no linked generator\n";
            ret += 1;
          }
        }

        return ret;
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

        // External Variables
        ScalarT efd0{0};

        // Initial Efd set by generator
        // The exciter object has no way of knowing if the generator
        // has set the initial value for Efd.
        // TODO: Build protections in system initialization call to
        // ensure Efd is initialized externally before the exciter initializes
        // other variables.
        if (signals_.template isAssigned<Ieeet1InternalVariables::EFD>())
        {
          efd0 = y_[7]; ///<- generator needs to be initialized first
        }

        // Terminal Voltage
        ScalarT vreal = bus_->Vr();
        ScalarT vimag = bus_->Vi();
        ScalarT Ec    = std::sqrt(vreal * vreal + vimag * vimag);

        // Derived from External initial values
        ScalarT vr  = Ke_ * efd0;
        ScalarT vfx = Kf_ / Tf_ * efd0;
        ScalarT vtr = Ke_ / Ka_ * efd0;

        // Vref (setpoint = terminal + error)
        vref_ = Ec + vtr;

        // IVP for Internal Variables
        y_[0] = Ec;   // y0 - vts  - Sensed term volt
        y_[1] = vr;   // y1 - vr   - Voltage reg
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
        tag_[1] = true;  // y1 - vr   - Voltage reg
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
       * @brief Set absolute tolerance
       */
      template <class ScalarT, typename IdxT>
      int Ieeet1<ScalarT, IdxT>::setAbsoluteTolerance(RealT rel_tol)
      {
        std::fill(abs_tol_.begin(), abs_tol_.end(), rel_tol);
        return 0;
      }

      /**
       * @brief Internal Residual
       *
       */
      template <class ScalarT, typename IdxT>
      __attribute__((always_inline)) inline int Ieeet1<ScalarT, IdxT>::evaluateInternalResidual(
          ScalarT* y,
          ScalarT* yp,
          ScalarT* wb,
          ScalarT* ws,
          ScalarT* f)
      {
        // Read E comp (terminal voltage, unless compensation impedance)
        ScalarT vreal = wb[0];
        ScalarT vimag = wb[1];
        ScalarT Ec    = std::sqrt(vreal * vreal + vimag * vimag);

        // Read Internal Variables
        ScalarT vts  = y[0]; // y0 - Sensed term volt
        ScalarT vr   = y[1]; // y1 - Voltage reg
        ScalarT efdp = y[2]; // y2 - Efd pre mult
        ScalarT vfx  = y[3]; // y3 - Exciter feedback
        ScalarT vtr  = y[4]; // y4 - Term Volt Err
        ScalarT vf   = y[5]; // y5 - Feedback volt
        ScalarT ve   = y[6]; // y6 - Excit. Cntrl Volt
        ScalarT efd  = y[7]; // y7 - Efd
        ScalarT ksat = y[8]; // y8 - Saturation

        // Read Internal Derivatives
        ScalarT vts_dot  = yp[0];
        ScalarT vr_dot   = yp[1];
        ScalarT efdp_dot = yp[2];
        ScalarT vfx_dot  = yp[3];

        // Set signal variable aliases
        ScalarT omega = ws[0];

        // The 'pre-limit' derivative of Pv
        ScalarT func            = -vr + Ka_ * vtr;
        ScalarT func_normalized = func * 12.0 / 240.0; // Arbitrary normalization to use standardized sigmoid
        ScalarT vr_ind          = Math::indicator(Vrmin_, Vrmax_, vr, func_normalized);

        // Internal Differential Equations
        f[0] = -vts_dot + (Ec - vts) / Tr_;
        f[1] = -vr_dot + vr_ind * func / Ta_;
        f[2] = -efdp_dot + (vr - ve - Ke_ * efdp) / Te_;
        f[3] = -vfx_dot + vf / Tf_;

        // Internal Algebraic Equations
        f[4] = -vts + vref_ + vUEL_ + vOEL_ + vS_ - vtr - vf;
        f[5] = -vf + (efdp * Kf_) / Tf_ - vfx;
        f[6] = -ve + ksat * efdp;
        f[7] = -efd + efdp + omega * efdp * Ispdlim_;

        ScalarT efd_sat = (efdp - SA_) * (Math::sigmoid(efdp - SA_));
        f[8]            = -ksat + SB_ * efd_sat * efd_sat;

        return 0;
      }

      /**
       * @brief Residual evaluation
       *
       */
      template <class ScalarT, typename IdxT>
      int Ieeet1<ScalarT, IdxT>::evaluateResidual()
      {
        // Set Input Variables
        // Meta PR Note: This seems to be very slow,
        // but I see how read/write ownership may require this
        //
        // I believe implementing the equivalent to signal->read()
        // at the system level would address this, by routing
        // external signals into a generic inputs_ vector
        // at the same time as the internal state values y_
        // are recieved from IDA.
        if (signals_.template isAttached<Ieeet1ExternalVariables::OMEGA>())
        {
          ws_[0]         = signals_.template readExternalVariable<Ieeet1ExternalVariables::OMEGA>();
          ws_indices_[0] = signals_.template readExternalVariableIndex<Ieeet1ExternalVariables::OMEGA>();
        }

        // Bus voltages
        wb_[0] = bus_->Vr();
        wb_[1] = bus_->Vi();

        // Residual evaluation
        evaluateInternalResidual(y_.data(), yp_.data(), wb_.data(), ws_.data(), f_.data());

        return 0;
      }

      /**
       * @brief Initialization Exciter Parameters from data structure
       */
      template <class ScalarT, typename IdxT>
      void Ieeet1<ScalarT, IdxT>::initModelParams(const model_data_type& data)
      {

        if (data.parameters.contains(model_data_type::Parameters::Tr))
        {
          Tr_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Tr));
        }
        if (data.parameters.contains(model_data_type::Parameters::Ka))
        {
          Ka_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Ka));
        }
        if (data.parameters.contains(model_data_type::Parameters::Ta))
        {
          Ta_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Ta));
        }
        if (data.parameters.contains(model_data_type::Parameters::Ke))
        {
          Ke_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Ke));
        }
        if (data.parameters.contains(model_data_type::Parameters::Te))
        {
          Te_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Te));
        }
        if (data.parameters.contains(model_data_type::Parameters::Kf))
        {
          Kf_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Kf));
        }
        if (data.parameters.contains(model_data_type::Parameters::Tf))
        {
          Tf_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Tf));
        }
        if (data.parameters.contains(model_data_type::Parameters::Vrmin))
        {
          Vrmin_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Vrmin));
        }
        if (data.parameters.contains(model_data_type::Parameters::Vrmax))
        {
          Vrmax_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Vrmax));
        }
        if (data.parameters.contains(model_data_type::Parameters::E1))
        {
          E1_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::E1));
        }
        if (data.parameters.contains(model_data_type::Parameters::E2))
        {
          E2_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::E2));
        }
        if (data.parameters.contains(model_data_type::Parameters::Se1))
        {
          Se1_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Se1));
        }
        if (data.parameters.contains(model_data_type::Parameters::Se2))
        {
          Se2_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Se2));
        }
        if (data.parameters.contains(model_data_type::Parameters::Ispdlim))
        {
          Ispdlim_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::Ispdlim));
        }

        // Derived Parameters
        RealT SR = std::sqrt(Se2_ / Se1_);

        // Solution 1 (Aligned with PW)
        SA_ = (SR * E1_ - E2_) / (SR - 1);
        SB_ = Se1_ / (E1_ - SA_) / (E1_ - SA_);
      }

      template <class ScalarT, typename IdxT>
      const Model::VariableMonitorBase* Ieeet1<ScalarT, IdxT>::getMonitor() const
      {
        return monitor_.get();
      }

      template <class ScalarT, typename IdxT>
      void Ieeet1<ScalarT, IdxT>::initializeMonitor()
      {
        using Variable = model_data_type::MonitorableVariables;
        monitor_->set(Variable::efd, [this]
                      { return efd_signal_->read(); });
        // monitor_->set(Variable::ksat, [this] { return ?; });
      }
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
