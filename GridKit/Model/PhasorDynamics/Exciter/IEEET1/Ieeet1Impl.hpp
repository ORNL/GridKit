/**
 * @file Ieeet1Impl.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 *
 * @brief Definition of a IEEET1 Exciter.
 *
 */

#include <algorithm>
#include <cmath>
#include <iostream>
#include <mutex>

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/IEEET1/Ieeet1.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/IEEET1/Ieeet1Data.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/ConfigurationChecks.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>
#include <GridKit/Utilities/ParameterReader.hpp>

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
        : Ieeet1(bus, ModelDataT{})
      {
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
        if (!allocated_)
        {
          this->allocateVectors(size_);
        }
        auto size = static_cast<size_t>(size_); // avoid compiler warnings

        tag_.resize(size);

        variable_indices_.resize(size);
        residual_indices_.resize(size);

        // Default variable and residual index mapping to local index
        for (IdxT j = 0; j < size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }

        // Resize bus data
        wb_.resize(2);

        // Resize signal variable data
        const auto signal_size = static_cast<size_t>(Ieeet1ExternalVariables::MAXIMUM);
        ws_.resize(static_cast<IdxT>(signal_size));
        ws_.setToZero();
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        // Set output signals
        if (signals_.template isAssigned<Ieeet1InternalVariables::EFD>())
        {
          auto* y = y_.getData();
          signals_.template getSignalNode<Ieeet1InternalVariables::EFD>()->set(&y[7], &(this->getVariableIndex(7)));
        }

        allocated_ = true;
        return 0;
      }

      /**
       * @brief Verify parameter values and attached signal links
       */
      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::verify() const
      {
        Utilities::ConfigurationChecks checks("Ieeet1");

        checks.check(Ka_ > ZERO<RealT>, "Ka must be positive");
        checks.check(Vrmin_ <= Vrmax_, "Vrmin must be less than or equal to Vrmax");
        checks.check(Ispdlim_ == ZERO<RealT> || Ispdlim_ == ONE<RealT>,
                     "Ispdlim must be 0 or 1");

        const bool saturation_disabled =
            Se1_ == ZERO<RealT> && Se2_ == ZERO<RealT>;

        if (!saturation_disabled)
        {
          checks.check(E1_ > ZERO<RealT>, "E1 must be positive when saturation is enabled");
          checks.check(E2_ > ZERO<RealT>, "E2 must be positive when saturation is enabled");
          checks.check(Se1_ >= ZERO<RealT>, "Se1 must be non-negative when saturation is enabled");
          checks.check(Se2_ >= ZERO<RealT>, "Se2 must be non-negative when saturation is enabled");

          const bool sat_ordered = (E2_ > E1_ && Se2_ > Se1_) || (E2_ < E1_ && Se2_ < Se1_);
          checks.check(sat_ordered, "E1/E2 and Se1/Se2 must be ordered consistently");
        }

        signals_.template checkOptional<Ieeet1ExternalVariables::OMEGA>(checks, "speed");
        signals_.template checkOptional<Ieeet1ExternalVariables::VREF>(checks, "vref");
        signals_.template checkOptional<Ieeet1ExternalVariables::VS>(checks, "vs");
        signals_.template checkOptional<Ieeet1ExternalVariables::VUEL>(checks, "vuel");
        signals_.template checkOptional<Ieeet1ExternalVariables::VOEL>(checks, "voel");

        return static_cast<int>(parameter_error_count_) + checks.errorCount();
      }

      /**
       * @brief Initialization of the Exciter
       *
       * Solves for a steady-state initial condition that satisfies
       * F(y, yp=0, t=0) = 0 exactly for every residual equation.
       *
       * Inputs:
       *   - EFD assigned by the generator.
       *   - Bus voltage, used to form the sensed terminal voltage magnitude.
       *   - Attached external signals (omega, V_S, V_UEL, V_OEL)
       *
       * Enabled saturation is included via ksat computed from efdp and SA, SB.
       * The resolved V_ref is written to an attached vref signal.
       *
       * @warning IEEE Std 421.5-2016 states: “In some programs, if
       *          \f$K_{E}\f$ is entered as zero, \f$K_{E}\f$ is automatically
       *          calculated by the program to represent a self-excited shunt
       *          field and a trimmed rheostat as its initial condition.” GridKit
       *          preserves the configured \f$K_{E}\f$ and resolves
       *          \f$K_{E}^{\mathrm{eff}}\f$ using the PSS/E-compatible
       *          \f$V_R = V_R^{\max}/10 = 0.1 V_R^{\max}\f$ rule. The divisor
       *          10 is unitless and represents 10% of the maximum regulator
       *          output.
       */
      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::initialize()
      {
        if (verify() != 0)
        {
          Log::error() << "Ieeet1: cannot initialize with invalid configuration\n";
          return 1;
        }

        // External Variables
        ScalarT efd0{0};
        auto*   y  = y_.getData();
        auto*   yp = yp_.getData();

        // Initial Efd set by generator
        // The exciter object has no way of knowing if the generator
        // has set the initial value for Efd.
        // TODO: Build protections in system initialization call to
        // ensure Efd is initialized externally before the exciter initializes
        // other variables.
        if (signals_.template isAssigned<Ieeet1InternalVariables::EFD>())
        {
          efd0 = y[7]; ///<- generator needs to be initialized first
        }

        // Setpoint members provide the defaults for unattached signals.
        const ScalarT omega = signals_.template readOrDefault<Ieeet1ExternalVariables::OMEGA>(omega_set_);
        const ScalarT vs    = signals_.template readOrDefault<Ieeet1ExternalVariables::VS>(vs_set_);
        const ScalarT vuel  = signals_.template readOrDefault<Ieeet1ExternalVariables::VUEL>(vuel_set_);
        const ScalarT voel  = signals_.template readOrDefault<Ieeet1ExternalVariables::VOEL>(voel_set_);

        uel_on_ = ZERO<RealT>;
        if (signals_.template isAttached<Ieeet1ExternalVariables::VUEL>())
        {
          uel_on_ = ONE<RealT>;
        }

        oel_on_ = ZERO<RealT>;
        if (signals_.template isAttached<Ieeet1ExternalVariables::VOEL>())
        {
          oel_on_ = ONE<RealT>;
        }

        // Terminal Voltage
        ScalarT vreal = bus_->Vr();
        ScalarT vimag = bus_->Vi();
        ScalarT Ec    = std::sqrt(vreal * vreal + vimag * vimag);

        ScalarT efdp = efd0 / (ONE<RealT> + omega * Ispdlim_);
        ScalarT ksat = SB_ * Math::qramp(efdp - SA_);
        if (Ke_ == ZERO<RealT>)
        {
          Ke_eff_ = (Vrmax_ / 10.0 - static_cast<RealT>(ksat))
                    / static_cast<RealT>(efdp);
          if (!std::isfinite(Ke_eff_))
          {
            Log::error() << "Ieeet1: derived effective Ke must be finite\n";
            return 1;
          }
          Log::misc() << "Ieeet1: Ke is zero so effective Ke is derived during initialization\n";
        }
        else
        {
          Ke_eff_ = Ke_;
        }

        ScalarT ve  = ksat;
        ScalarT vr  = Ke_eff_ * efdp + ve;
        ScalarT vtr = vr / Ka_;
        ScalarT vf{0};
        ScalarT vfx = (Kf_ / Tf_) * efdp;

        const ScalarT vref = Ec + vtr + vf - vs - uel_on_ * vuel - oel_on_ * voel;

        y[0] = Ec;   // y0 - vts  - Sensed term volt
        y[1] = vr;   // y1 - vr   - Voltage reg
        y[2] = efdp; // y2 - efdp - Efd pre mult
        y[3] = vfx;  // y3 - vfx  - Exciter feedback
        y[4] = vtr;  // y4 - vtr  - Term Volt Err
        y[5] = vf;   // y5 - vf   - Feedback volt
        y[6] = ve;   // y6 - ve   - Excit. Cntrl Volt
        y[7] = efd0; // y7 - efd  - Efd
        y[8] = ksat; // y8 - ksat - Saturation

        for (IdxT i = 0; i < yp_.getSize(); ++i)
        {
          yp[i] = 0.0;
        }

        omega_set_ = omega;
        vref_set_  = vref;
        vs_set_    = vs;
        vuel_set_  = vuel;
        voel_set_  = voel;

        if (signals_.template isAttached<Ieeet1ExternalVariables::VREF>())
        {
          signals_.template writeExternalVariable<Ieeet1ExternalVariables::VREF>(vref_set_);
        }

        y_.setDataUpdated();
        yp_.setDataUpdated();

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
       * @tparam scalar_type Scalar data type
       * @tparam index_type Index data type
       * @return int 0 if successful, non-zero otherwise.
       *
       * This represents a "noise" level close to zero for which pure relative
       * error cannot be used.
       */
      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
      {
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
        return 0;
      }

      /**
       * @brief Internal Residual
       *
       */
      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline int Ieeet1<scalar_type, index_type>::evaluateInternalResidual(
          const ScalarT* y,
          const ScalarT* yp,
          const ScalarT* wb,
          const ScalarT* ws,
          ScalarT*       f)
      {
        const auto OMEGA = static_cast<size_t>(Ieeet1ExternalVariables::OMEGA);
        const auto VREF  = static_cast<size_t>(Ieeet1ExternalVariables::VREF);
        const auto VS    = static_cast<size_t>(Ieeet1ExternalVariables::VS);
        const auto VUEL  = static_cast<size_t>(Ieeet1ExternalVariables::VUEL);
        const auto VOEL  = static_cast<size_t>(Ieeet1ExternalVariables::VOEL);

        // Read bus voltage components
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
        ScalarT omega = ws[OMEGA];
        ScalarT vref  = ws[VREF];
        ScalarT vs    = ws[VS];
        ScalarT vuel  = ws[VUEL];
        ScalarT voel  = ws[VOEL];

        // The 'pre-limit' derivative of Vr.
        ScalarT func = (-vr + Ka_ * vtr) / Ta_;

        // Internal Differential Equations
        f[0] = -vts_dot + (Ec - vts) / Tr_;
        f[1] = -vr_dot + Math::antiwindup(vr, func, Vrmin_, Vrmax_);
        f[2] = -efdp_dot + (vr - ve - Ke_eff_ * efdp) / Te_;
        f[3] = -vfx_dot + vf / Tf_;

        // Internal Algebraic Equations
        f[4] = -vts + vref + vs + uel_on_ * vuel + oel_on_ * voel - vtr - vf;
        f[5] = -Tf_ * (vf + vfx) + Kf_ * efdp;
        f[6] = -ve + ksat;
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
        auto* ws = ws_.getData();

        // Attached signals are read live; unattached ones keep the latched value.
        auto* ws_indices = ws_indices_.data();
        signals_.template refreshWorkspace<Ieeet1ExternalVariables::OMEGA>(omega_set_, ws, ws_indices);
        signals_.template refreshWorkspace<Ieeet1ExternalVariables::VREF>(vref_set_, ws, ws_indices);
        signals_.template refreshWorkspace<Ieeet1ExternalVariables::VS>(vs_set_, ws, ws_indices);
        signals_.template refreshWorkspace<Ieeet1ExternalVariables::VUEL>(vuel_set_, ws, ws_indices);
        signals_.template refreshWorkspace<Ieeet1ExternalVariables::VOEL>(voel_set_, ws, ws_indices);

        // Bus voltages
        auto* wb = wb_.getData();
        wb[0]    = bus_->Vr();
        wb[1]    = bus_->Vi();

        // Residual evaluation
        const auto* y  = y_.getData();
        const auto* yp = yp_.getData();
        auto*       f  = f_.getData();
        evaluateInternalResidual(y, yp, wb, ws, f);

        f_.setDataUpdated();

        return 0;
      }

      /**
       * @brief Initialization Exciter Parameters from data structure
       */
      template <typename scalar_type, typename index_type>
      void Ieeet1<scalar_type, index_type>::initModelParams(const ModelDataT& data)
      {
        using Parameter = typename ModelDataT::Parameters;

        parameter_error_count_ = 0;

        Utilities::ConfigurationChecks checks("Ieeet1");
        Utilities::ParameterReader     reader(data, checks);

        reader.loadReal(Parameter::Tr, Tr_);
        reader.loadReal(Parameter::Ka, Ka_);
        reader.loadReal(Parameter::Ta, Ta_);
        reader.loadReal(Parameter::Ke, Ke_);
        reader.loadReal(Parameter::Te, Te_);
        reader.loadReal(Parameter::Kf, Kf_);
        reader.loadReal(Parameter::Tf, Tf_);
        reader.loadReal(Parameter::Vrmin, Vrmin_);
        reader.loadReal(Parameter::Vrmax, Vrmax_);
        reader.loadReal(Parameter::E1, E1_);
        reader.loadReal(Parameter::E2, E2_);
        reader.loadReal(Parameter::Se1, Se1_);
        reader.loadReal(Parameter::Se2, Se2_);
        reader.loadReal(Parameter::Ispdlim, Ispdlim_);

        parameter_error_count_ = static_cast<IdxT>(checks.errorCount());

        setDerivedParameters();
      }

      /**
       * @brief Static method to log time constant warnings
       *
       * @note Used in combination with static std:once_flag and std:call_once,
       *       to reduce the number of times the warning is printed.
       */
      template <typename scalar_type, typename index_type>
      void Ieeet1<scalar_type, index_type>::logTimeConstantWarning()
      {
        Log::warning() << "Ieeet1: Tr, Ta, Te, and Tf below "
                       << TIME_CONSTANT_MINIMUM
                       << " s are raised to that floor\n";
      }

      /**
       * @brief Resolve the parameter-derived constants
       */
      template <typename scalar_type, typename index_type>
      void Ieeet1<scalar_type, index_type>::setDerivedParameters()
      {
        if (Tr_ < TIME_CONSTANT_MINIMUM || Ta_ < TIME_CONSTANT_MINIMUM
            || Te_ < TIME_CONSTANT_MINIMUM || Tf_ < TIME_CONSTANT_MINIMUM)
        {
          static std::once_flag time_constant_warning_flag_;
          std::call_once(time_constant_warning_flag_,
                         &logTimeConstantWarning);
        }

        Tr_ = std::max(Tr_, TIME_CONSTANT_MINIMUM);
        Ta_ = std::max(Ta_, TIME_CONSTANT_MINIMUM);
        Te_ = std::max(Te_, TIME_CONSTANT_MINIMUM);
        Tf_ = std::max(Tf_, TIME_CONSTANT_MINIMUM);

        SA_ = ZERO<RealT>;
        SB_ = ZERO<RealT>;

        const bool saturation_disabled =
            Se1_ == ZERO<RealT> && Se2_ == ZERO<RealT>;

        if (saturation_disabled)
        {
          return;
        }

        const bool sat_ordered = (E2_ > E1_ && Se2_ > Se1_) || (E2_ < E1_ && Se2_ < Se1_);
        if (E1_ <= ZERO<RealT> || E2_ <= ZERO<RealT>
            || Se1_ < ZERO<RealT> || Se2_ < ZERO<RealT>
            || !sat_ordered)
        {
          return;
        }

        if (Se1_ == ZERO<RealT>)
        {
          const RealT dE = E2_ - E1_;
          SA_            = E1_;
          SB_            = Se2_ * E2_ / (dE * dE);
          return;
        }

        if (Se2_ == ZERO<RealT>)
        {
          const RealT dE = E1_ - E2_;
          SA_            = E2_;
          SB_            = Se1_ * E1_ / (dE * dE);
          return;
        }

        const RealT C = std::sqrt(Se2_ * E2_ / (Se1_ * E1_));

        // Solution 1 (Aligned with PW)
        SA_ = (C * E1_ - E2_) / (C - ONE<RealT>);
        SB_ = Se1_ * E1_ / ((E1_ - SA_) * (E1_ - SA_));
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
                      { return y_.getData()[7]; });
        monitor_->set(Variable::ksat, [this]
                      { return SB_ * Math::qramp(y_.getData()[2] - SA_); });
      }
    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
