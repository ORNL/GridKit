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
       */
      template <typename scalar_type, typename index_type>
      Ieeet1<scalar_type, index_type>::Ieeet1(BusT* bus)
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
       */
      template <typename scalar_type, typename index_type>
      Ieeet1<scalar_type, index_type>::Ieeet1(SignalT*          efd_signal,
                                              SignalT*          speed_signal,
                                              BusT*             bus,
                                              const ModelDataT& data)
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
       */
      template <typename scalar_type, typename index_type>
      Ieeet1<scalar_type, index_type>::Ieeet1(BusT*             bus,
                                              const ModelDataT& data)
        : bus_(bus),
          monitor_(std::make_unique<MonitorT>(data))
      {

        // Parse data struct into model
        this->initModelParams(data);

        initializeMonitor();

        // 9 Internal Variables
        size_ = 9;
      }

      template <typename scalar_type, typename index_type>
      Ieeet1<scalar_type, index_type>::~Ieeet1()
      {
      }

      /**
       * @brief Set the component ID
       */
      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      /**
       * @brief Allocate memory for model
       *
       */
      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::allocate()
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
        ws_.resize(2);
        ws_indices_.resize(2);
        ws_[0]         = 0.0;
        ws_indices_[0] = INVALID_INDEX<IdxT>;
        ws_[1]         = 0.0;
        ws_indices_[1] = INVALID_INDEX<IdxT>;

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
      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::verify() const
      {
        static constexpr auto OMEGA = Ieeet1ExternalVariables::OMEGA;
        static constexpr auto VS    = Ieeet1ExternalVariables::VS;

        int ret = 0;

        if (signals_.template isAttached<OMEGA>())
        {
          if (!signals_.template isLinked<OMEGA>())
          {
            Log::error() << "Ieeet1: omega signal attached with no linked generator\n";
            ret += 1;
          }
        }

        if (signals_.template isAttached<VS>())
        {
          if (!signals_.template isLinked<VS>())
          {
            Log::error() << "Ieeet1: VS signal attached with no linked source\n";
            ret += 1;
          }
        }

        return ret;
      }

      /**
       * @brief Initialization of the Exciter
       *
       * Solves for a steady-state initial condition that satisfies
       * F(y, yp=0, t=0) = 0 exactly for every residual equation.
       *
       * Inputs:
       *   - EFD assigned by the generator (read from y_[7]).
       *   - Bus voltage, used to form the sensed terminal voltage Ec.
       *   - Attached external signals (omega, V_S)
       *
       * Saturation is included via ksat computed from efdp and SA, SB.
       */
      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::initialize()
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

        ScalarT omega{0};
        ScalarT vs{0};
        if (signals_.template isAttached<Ieeet1ExternalVariables::OMEGA>())
        {
          omega = signals_.template readExternalVariable<Ieeet1ExternalVariables::OMEGA>();
        }
        if (signals_.template isAttached<Ieeet1ExternalVariables::VS>())
        {
          vs = signals_.template readExternalVariable<Ieeet1ExternalVariables::VS>();
        }

        // Terminal Voltage
        ScalarT vreal = bus_->Vr();
        ScalarT vimag = bus_->Vi();
        ScalarT Ec    = std::sqrt(vreal * vreal + vimag * vimag);

        ScalarT efdp = efd0 / (ONE<RealT> + omega * Ispdlim_);
        ScalarT ksat = SB_ * Math::qramp(efdp - SA_);
        ScalarT ve   = ksat * efdp;
        ScalarT vr   = Ke_ * efdp + ve;
        ScalarT vtr  = vr / Ka_;
        ScalarT vf{0};
        ScalarT vfx = (Kf_ / Tf_) * efdp;

        vref_ = Ec + vtr + vf - vUEL_ - vOEL_ - vs;

        y_[0] = Ec;   // y0 - vts  - Sensed term volt
        y_[1] = vr;   // y1 - vr   - Voltage reg
        y_[2] = efdp; // y2 - efdp - Efd pre mult
        y_[3] = vfx;  // y3 - vfx  - Exciter feedback
        y_[4] = vtr;  // y4 - vtr  - Term Volt Err
        y_[5] = vf;   // y5 - vf   - Feedback volt
        y_[6] = ve;   // y6 - ve   - Excit. Cntrl Volt
        y_[7] = efd0; // y7 - efd  - Efd
        y_[8] = ksat; // y8 - ksat - Saturation

        for (size_t i = 0; i < yp_.size(); ++i)
        {
          yp_[i] = 0.0;
        }

        return 0;
      }

      /**
       * @brief  Identify differential variables.
       *
       * @return int 0
       */
      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::tagDifferentiable()
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
       * @brief Compute the absolute tolerance for each variable in the model
       *
       * @param rel_tol The relative tolerance which can be used to pick the
       *        absolute tolerance.
       * @tparam ScalarT Scalar data type
       * @tparam IdxT Index data type
       * @return int 0 if successful, non-zero otherwise.
       *
       * This represents a "noise" level close to zero for which pure relative
       * error cannot be used.
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
      template <typename scalar_type, typename index_type>
      FORCE_INLINE int Ieeet1<scalar_type, index_type>::evaluateInternalResidual(
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
        ScalarT omega     = ws[0];
        ScalarT vs_signal = ws[1];

        // The 'pre-limit' derivative of Vr.
        ScalarT func            = (-vr + Ka_ * vtr) / Ta_;
        ScalarT func_normalized = func / static_cast<RealT>(500.0); // TODO This is arbitrary, need more general conditioning method that is fast
        ScalarT vr_ind          = Math::indicator(vr, func_normalized, Vrmin_, Vrmax_);

        // Internal Differential Equations
        f[0] = -vts_dot + (Ec - vts) / Tr_;
        f[1] = -vr_dot + vr_ind * func;
        f[2] = -efdp_dot + (vr - ve - Ke_ * efdp) / Te_;
        f[3] = -vfx_dot + vf / Tf_;

        // Internal Algebraic Equations
        f[4] = -vts + vref_ + vUEL_ + vOEL_ + vs_signal - vtr - vf;
        f[5] = -vf + (efdp * Kf_) / Tf_ - vfx;
        f[6] = -ve + ksat * efdp;
        f[7] = -efd + efdp + omega * efdp * Ispdlim_;

        f[8] = -ksat + SB_ * Math::qramp(efdp - SA_);

        return 0;
      }

      /**
       * @brief Residual evaluation
       *
       */
      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::evaluateResidual()
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

        // VS signal (stabilizer output, optional)
        ws_[1]         = 0.0;
        ws_indices_[1] = INVALID_INDEX<IdxT>;
        if (signals_.template isAttached<Ieeet1ExternalVariables::VS>())
        {
          ws_[1]         = signals_.template readExternalVariable<Ieeet1ExternalVariables::VS>();
          ws_indices_[1] = signals_.template readExternalVariableIndex<Ieeet1ExternalVariables::VS>();
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
      template <typename scalar_type, typename index_type>
      void Ieeet1<scalar_type, index_type>::initModelParams(const ModelDataT& data)
      {

        if (data.parameters.contains(ModelDataT::Parameters::Tr))
        {
          Tr_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::Tr));
        }
        if (data.parameters.contains(ModelDataT::Parameters::Ka))
        {
          Ka_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::Ka));
        }
        if (data.parameters.contains(ModelDataT::Parameters::Ta))
        {
          Ta_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::Ta));
        }
        if (data.parameters.contains(ModelDataT::Parameters::Ke))
        {
          Ke_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::Ke));
        }
        if (data.parameters.contains(ModelDataT::Parameters::Te))
        {
          Te_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::Te));
        }
        if (data.parameters.contains(ModelDataT::Parameters::Kf))
        {
          Kf_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::Kf));
        }
        if (data.parameters.contains(ModelDataT::Parameters::Tf))
        {
          Tf_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::Tf));
        }
        if (data.parameters.contains(ModelDataT::Parameters::Vrmin))
        {
          Vrmin_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::Vrmin));
        }
        if (data.parameters.contains(ModelDataT::Parameters::Vrmax))
        {
          Vrmax_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::Vrmax));
        }
        if (data.parameters.contains(ModelDataT::Parameters::E1))
        {
          E1_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::E1));
        }
        if (data.parameters.contains(ModelDataT::Parameters::E2))
        {
          E2_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::E2));
        }
        if (data.parameters.contains(ModelDataT::Parameters::Se1))
        {
          Se1_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::Se1));
        }
        if (data.parameters.contains(ModelDataT::Parameters::Se2))
        {
          Se2_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::Se2));
        }
        if (data.parameters.contains(ModelDataT::Parameters::Ispdlim))
        {
          Ispdlim_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::Ispdlim));
        }

        // Derived Parameters
        RealT SR = std::sqrt(Se2_ / Se1_);

        // Solution 1 (Aligned with PW)
        SA_ = (SR * E1_ - E2_) / (SR - 1);
        SB_ = Se1_ / (E1_ - SA_) / (E1_ - SA_);
      }

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* Ieeet1<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      template <typename scalar_type, typename index_type>
      void Ieeet1<scalar_type, index_type>::initializeMonitor()
      {
        using Variable = ModelDataT::MonitorableVariables;
        monitor_->set(Variable::efd, [this]
                      { return efd_signal_->read(); });
        monitor_->set(Variable::ksat, [this]
                      { return SB_ * Math::qramp(y_[2] - SA_); });
      }
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
