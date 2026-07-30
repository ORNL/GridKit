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

#include <GridKit/Model/PhasorDynamics/Exciter/IEEET1/Ieeet1.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/IEEET1/Ieeet1Data.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

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
      Ieeet1<scalar_type, index_type>::Ieeet1()
        : Ieeet1(ModelDataT{})
      {
      }

      /**
       * @brief  Constructor for IEEET1 Exciter
       *
       * @param data  Data object to store parameters
       */
      template <typename scalar_type, typename index_type>
      Ieeet1<scalar_type, index_type>::Ieeet1(const ModelDataT& data)
        : monitor_(std::make_unique<MonitorT>(data))
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

        // Resize coupling data
        this->allocateExternalVectors(static_cast<IdxT>(Ieeet1ExternalVariables::MAXIMUM));

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
        static constexpr auto OMEGA = Ieeet1ExternalVariables::OMEGA;
        static constexpr auto VREAL = Ieeet1ExternalVariables::VREAL;
        static constexpr auto VIMAG = Ieeet1ExternalVariables::VIMAG;
        static constexpr auto VS    = Ieeet1ExternalVariables::VS;

        int ret = 0;

        auto check = [&](bool condition, const char* message)
        {
          if (!condition)
          {
            Log::error() << "Ieeet1: " << message << '\n';
            ret += 1;
          }
        };

        check(signals_.template isAttached<VREAL>(), "VREAL signal is not attached");
        check(signals_.template isAttached<VIMAG>(), "VIMAG signal is not attached");

        check(Ka_ > ZERO<RealT>, "Ka must be positive");
        check(Vrmin_ <= Vrmax_, "Vrmin must be less than or equal to Vrmax");
        check(Ispdlim_ == ZERO<RealT> || Ispdlim_ == ONE<RealT>,
              "Ispdlim must be 0 or 1");

        const bool saturation_disabled =
            Se1_ == ZERO<RealT> && Se2_ == ZERO<RealT>;

        if (!saturation_disabled)
        {
          check(E1_ > ZERO<RealT>, "E1 must be positive when saturation is enabled");
          check(E2_ > ZERO<RealT>, "E2 must be positive when saturation is enabled");
          check(Se1_ > ZERO<RealT>, "Se1 must be positive when saturation is enabled");
          check(Se2_ > ZERO<RealT>, "Se2 must be positive when saturation is enabled");
          check(E1_ != E2_, "E1 and E2 must differ when saturation is enabled");
          check(Se1_ != Se2_, "Se1 and Se2 must differ when saturation is enabled");
        }

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
       *   - EFD assigned by the generator.
       *   - Bus voltage, used to form the sensed terminal voltage magnitude.
       *   - Attached external signals (omega, V_S)
       *
       * Enabled saturation is included via ksat computed from efdp and SA, SB.
       */
      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::initialize()
      {

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
        ScalarT vreal = signals_.template readExternalVariable<Ieeet1ExternalVariables::VREAL>();
        ScalarT vimag = signals_.template readExternalVariable<Ieeet1ExternalVariables::VIMAG>();
        ScalarT Ec    = std::sqrt(vreal * vreal + vimag * vimag);

        ScalarT efdp = efd0 / (ONE<RealT> + omega * Ispdlim_);
        ScalarT ksat = SB_ * Math::qramp(efdp - SA_);
        ScalarT ve   = ksat * efdp;
        ScalarT vr   = Ke_ * efdp + ve;
        ScalarT vtr  = vr / Ka_;
        ScalarT vf{0};
        ScalarT vfx = (Kf_ / Tf_) * efdp;

        vref_ = Ec + vtr + vf - vUEL_ - vOEL_ - vs;

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
          const ScalarT* y_ext,
          ScalarT*       f)
      {
        // Read bus voltage components
        ScalarT vreal = y_ext[1];
        ScalarT vimag = y_ext[2];
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
        ScalarT omega     = y_ext[0];
        ScalarT vs_signal = y_ext[3];

        // The 'pre-limit' derivative of Vr.
        ScalarT func = (-vr + Ka_ * vtr) / Ta_;

        // Internal Differential Equations
        f[0] = -vts_dot + (Ec - vts) / Tr_;
        f[1] = -vr_dot + Math::antiwindup(vr, func, Vrmin_, Vrmax_);
        f[2] = -efdp_dot + (vr - ve - Ke_ * efdp) / Te_;
        f[3] = -vfx_dot + vf / Tf_;

        // Internal Algebraic Equations
        f[4] = -vts + vref_ + vUEL_ + vOEL_ + vs_signal - vtr - vf;
        f[5] = -Tf_ * (vf + vfx) + Kf_ * efdp;
        f[6] = -ve + ksat * efdp;
        f[7] = -efd + efdp + omega * efdp * Ispdlim_;
        f[8] = -ksat + SB_ * Math::qramp(efdp - SA_);

        return 0;
      }

      /**
       * @brief Gather external variables and index maps.
       *
       */
      template <typename scalar_type, typename index_type>
      void Ieeet1<scalar_type, index_type>::gatherExternalVariables()
      {
        // Set input variables.
        if (signals_.template isAttached<Ieeet1ExternalVariables::OMEGA>())
        {
          y_ext_[0]                = signals_.template readExternalVariable<Ieeet1ExternalVariables::OMEGA>();
          variable_indices_ext_[0] = signals_.template readExternalVariableIndex<Ieeet1ExternalVariables::OMEGA>();
        }

        // Bus voltages
        static constexpr auto VREAL = Ieeet1ExternalVariables::VREAL;
        static constexpr auto VIMAG = Ieeet1ExternalVariables::VIMAG;

        y_ext_[1]                = signals_.template readExternalVariable<VREAL>();
        y_ext_[2]                = signals_.template readExternalVariable<VIMAG>();
        variable_indices_ext_[1] = signals_.template readExternalVariableIndex<VREAL>();
        variable_indices_ext_[2] = signals_.template readExternalVariableIndex<VIMAG>();

        // VS signal (stabilizer output, optional)
        y_ext_[3]                = 0.0;
        variable_indices_ext_[3] = INVALID_INDEX<IdxT>;
        if (signals_.template isAttached<Ieeet1ExternalVariables::VS>())
        {
          y_ext_[3]                = signals_.template readExternalVariable<Ieeet1ExternalVariables::VS>();
          variable_indices_ext_[3] = signals_.template readExternalVariableIndex<Ieeet1ExternalVariables::VS>();
        }
      }

      /**
       * @brief Residual evaluation
       *
       */
      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::evaluateInternalResidual()
      {
        gatherExternalVariables();

        // Residual evaluation
        const auto* y  = y_.getData();
        const auto* yp = yp_.getData();
        auto*       f  = f_.getData();
        evaluateInternalResidual(y, yp, y_ext_.data(), f);

        f_.setDataUpdated();

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::evaluateResidual()
      {
        evaluateInternalResidual();
        return this->evaluateExternalResidual();
      }

      /**
       * @brief Initialization Exciter Parameters from data structure
       */
      template <typename scalar_type, typename index_type>
      void Ieeet1<scalar_type, index_type>::initModelParams(const ModelDataT& data)
      {
        using Parameter = typename ModelDataT::Parameters;

        if (data.parameters.contains(Parameter::Tr))
        {
          Tr_ = std::get<RealT>(data.parameters.at(Parameter::Tr));
        }
        if (data.parameters.contains(Parameter::Ka))
        {
          Ka_ = std::get<RealT>(data.parameters.at(Parameter::Ka));
        }
        if (data.parameters.contains(Parameter::Ta))
        {
          Ta_ = std::get<RealT>(data.parameters.at(Parameter::Ta));
        }
        if (data.parameters.contains(Parameter::Ke))
        {
          Ke_ = std::get<RealT>(data.parameters.at(Parameter::Ke));
        }
        if (data.parameters.contains(Parameter::Te))
        {
          Te_ = std::get<RealT>(data.parameters.at(Parameter::Te));
        }
        if (data.parameters.contains(Parameter::Kf))
        {
          Kf_ = std::get<RealT>(data.parameters.at(Parameter::Kf));
        }
        if (data.parameters.contains(Parameter::Tf))
        {
          Tf_ = std::get<RealT>(data.parameters.at(Parameter::Tf));
        }
        if (data.parameters.contains(Parameter::Vrmin))
        {
          Vrmin_ = std::get<RealT>(data.parameters.at(Parameter::Vrmin));
        }
        if (data.parameters.contains(Parameter::Vrmax))
        {
          Vrmax_ = std::get<RealT>(data.parameters.at(Parameter::Vrmax));
        }
        if (data.parameters.contains(Parameter::E1))
        {
          E1_ = std::get<RealT>(data.parameters.at(Parameter::E1));
        }
        if (data.parameters.contains(Parameter::E2))
        {
          E2_ = std::get<RealT>(data.parameters.at(Parameter::E2));
        }
        if (data.parameters.contains(Parameter::Se1))
        {
          Se1_ = std::get<RealT>(data.parameters.at(Parameter::Se1));
        }
        if (data.parameters.contains(Parameter::Se2))
        {
          Se2_ = std::get<RealT>(data.parameters.at(Parameter::Se2));
        }
        if (data.parameters.contains(Parameter::Ispdlim))
        {
          Ispdlim_ = std::get<RealT>(data.parameters.at(Parameter::Ispdlim));
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

        if (E1_ <= ZERO<RealT> || E2_ <= ZERO<RealT> || E1_ == E2_
            || Se1_ <= ZERO<RealT> || Se2_ <= ZERO<RealT> || Se1_ == Se2_)
        {
          return;
        }

        const RealT C = std::sqrt(Se2_ / Se1_);

        // Solution 1 (Aligned with PW)
        SA_ = (C * E1_ - E2_) / (C - ONE<RealT>);
        SB_ = Se1_ / ((E1_ - SA_) * (E1_ - SA_));
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
