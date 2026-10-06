/**
 * @file BranchBreakersImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the breaker-terminated branch model.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <numbers>
#include <variant>

#include <GridKit/Model/PhasorDynamics/Branch/BranchBreakers/BranchBreakers.hpp>
#include <GridKit/Model/PhasorDynamics/Branch/BranchBreakers/BranchBreakersData.hpp>
#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Enum.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Construct a breaker-terminated branch without parameters
     *
     * The model is sized and every parameter keeps its documented default.
     * No monitor is created.
     *
     * @param[in] bus1 Bus-1 terminal, the tapped side.
     * @param[in] bus2 Bus-2 terminal.
     */
    template <typename scalar_type, typename index_type>
    BranchBreakers<scalar_type, index_type>::BranchBreakers(BusT* bus1, BusT* bus2)
      : bus1_(bus1),
        bus2_(bus2)
    {
      size_ = static_cast<IdxT>(Utilities::enum_size<BranchBreakersInternalVariables>());
      setDerivedParameters();
    }

    /**
     * @brief Construct a breaker-terminated branch from model data
     *
     * @param[in] bus1 Bus-1 terminal, the tapped side.
     * @param[in] bus2 Bus-2 terminal.
     * @param[in] data Parameters, ports, and monitored-variable selections.
     */
    template <typename scalar_type, typename index_type>
    BranchBreakers<scalar_type, index_type>::BranchBreakers(BusT* bus1, BusT* bus2, const ModelDataT& data)
      : bus1_(bus1),
        bus2_(bus2),
        monitor_(std::make_unique<MonitorT>(data))
    {
      initializeParameters(data);
      initializeMonitor();
      size_ = static_cast<IdxT>(Utilities::enum_size<BranchBreakersInternalVariables>());
    }

    /**
     * @brief Destroy the breaker-terminated branch.
     */
    template <typename scalar_type, typename index_type>
    BranchBreakers<scalar_type, index_type>::~BranchBreakers()
    {
    }

    /**
     * @brief Set the component ID
     *
     * @param[in] component_id Identifier assigned by the system model.
     * @return Zero on success.
     */
    template <typename scalar_type, typename index_type>
    int BranchBreakers<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /**
     * @brief Allocate the model vectors and wire the terminal-current outputs
     *
     * Sizes the state, residual, bus-interface, and signal-interface buffers,
     * seeds the identity index maps, and points assigned current outputs at
     * the published terminal currents.
     *
     * @return Zero on success.
     */
    template <typename scalar_type, typename index_type>
    int BranchBreakers<scalar_type, index_type>::allocate()
    {
      if (!allocated_)
      {
        this->allocateVectors(size_);
      }
      auto size = static_cast<size_t>(size_);

      tag_.assign(size, false);
      variable_indices_.resize(size);
      residual_indices_.resize(size);

      wb_.resize(4);
      wb_.setToZero();
      h_.resize(4);
      h_.setToZero();

      const auto signal_size = Utilities::enum_size<BranchBreakersExternalVariables>();
      ws_.resize(static_cast<IdxT>(signal_size));
      ws_.setToZero();
      ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, j);
        this->setResidualIndex(j, j);
      }

      if (auto port = ports_.out.template port<BranchBreakersSignalOutputs::ir1>())
      {
        port.link(&ir1_, &current_index_);
      }
      if (auto port = ports_.out.template port<BranchBreakersSignalOutputs::ii1>())
      {
        port.link(&ii1_, &current_index_);
      }
      if (auto port = ports_.out.template port<BranchBreakersSignalOutputs::ir2>())
      {
        port.link(&ir2_, &current_index_);
      }
      if (auto port = ports_.out.template port<BranchBreakersSignalOutputs::ii2>())
      {
        port.link(&ii2_, &current_index_);
      }

      allocated_ = true;
      return 0;
    }

    /**
     * @brief Validate the breaker-terminated branch configuration
     *
     * @return Number of configuration errors; zero when valid.
     */
    template <typename scalar_type, typename index_type>
    int BranchBreakers<scalar_type, index_type>::verify() const
    {
      int ret = static_cast<int>(parameter_error_count_);

      auto check = [&](bool condition, const char* message)
      {
        if (!condition)
        {
          Log::error() << "BranchBreakers: " << message << '\n';
          ret += 1;
        }
      };

      check(bus1_ != nullptr, "bus1 pointer is null");
      check(bus2_ != nullptr, "bus2 pointer is null");
      check(R_ * R_ + X_ * X_ > ZERO<RealT>, "R and X cannot both be zero");
      check(tap_ > ZERO<RealT>, "tap must be positive");

      // An open side is Kron-reduced through the opposite diagonal entry.
      check(g11_ * g11_ + b11_ * b11_ > ZERO<RealT>, "bus-1 diagonal admittance must be nonzero");
      check(g22_ * g22_ + b22_ * b22_ > ZERO<RealT>, "bus-2 diagonal admittance must be nonzero");

      auto check_attached_signal =
          [&]<BranchBreakersSignalInputs variable>(const char* name)
      {
        if (ports_.in.template port<variable>().connected()
            && !ports_.in.template port<variable>().linked())
        {
          Log::error() << "BranchBreakers: " << name << " signal attached with no linked source\n";
          ret += 1;
        }
      };

      check_attached_signal.template operator()<BranchBreakersSignalInputs::trip1>("trip1");
      check_attached_signal.template operator()<BranchBreakersSignalInputs::reset1>("reset1");
      check_attached_signal.template operator()<BranchBreakersSignalInputs::trip2>("trip2");
      check_attached_signal.template operator()<BranchBreakersSignalInputs::reset2>("reset2");

      return ret;
    }

    /**
     * @brief Initialize both breakers closed
     *
     * Publishes the closed-branch terminal currents at the initial bus voltages.
     *
     * @return Zero on success; nonzero when the configuration is rejected.
     */
    template <typename scalar_type, typename index_type>
    int BranchBreakers<scalar_type, index_type>::initialize()
    {
      const auto Z1 = static_cast<size_t>(BranchBreakersInternalVariables::Z1);
      const auto Z2 = static_cast<size_t>(BranchBreakersInternalVariables::Z2);

      if (!allocated_)
      {
        Log::error() << "BranchBreakers: allocate must complete before initialize\n";
        return 1;
      }

      if (verify() > 0)
      {
        Log::error() << "BranchBreakers: cannot initialize with invalid configuration\n";
        return 1;
      }

      auto* y = y_.getData();
      y[Z1]   = ONE<RealT>;
      y[Z2]   = ONE<RealT>;

      y_.setDataUpdated();
      yp_.setToConst(static_cast<ScalarT>(ZERO<RealT>));

      ScalarT current[4];
      terminalCurrents(current);
      ir1_ = current[0];
      ii1_ = current[1];
      ir2_ = current[2];
      ii2_ = current[3];

      if constexpr (std::is_same_v<scalar_type, DependencyTracking::Variable>)
      {
        this->initializeDependencyTrackingVariableNumbers();
      }

      return 0;
    }

    /**
     * @brief Identify the differential variables
     *
     * Both breaker latches carry derivatives.
     *
     * @return Zero on success.
     */
    template <typename scalar_type, typename index_type>
    int BranchBreakers<scalar_type, index_type>::tagDifferentiable()
    {
      std::fill(tag_.begin(), tag_.end(), true);
      return 0;
    }

    /**
     * @brief Compute the absolute tolerance for each variable in the model
     *
     * Every internal variable receives @p rel_tol as its absolute tolerance.
     *
     * @param[in] rel_tol Solver relative tolerance.
     * @return Zero on success.
     */
    template <typename scalar_type, typename index_type>
    int BranchBreakers<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return 0;
    }

    /**
     * @brief Residuals of system equations
     *
     * Refreshes the bus and signal interface buffers, evaluates the breaker
     * latches, adds the terminal currents to the bus residuals, and publishes
     * them to the current outputs. An unattached trip or reset port reads zero.
     *
     * @return Zero on success.
     */
    template <typename scalar_type, typename index_type>
    int BranchBreakers<scalar_type, index_type>::evaluateResidual()
    {
      const auto TRIP1  = static_cast<size_t>(BranchBreakersExternalVariables::TRIP1);
      const auto RESET1 = static_cast<size_t>(BranchBreakersExternalVariables::RESET1);
      const auto TRIP2  = static_cast<size_t>(BranchBreakersExternalVariables::TRIP2);
      const auto RESET2 = static_cast<size_t>(BranchBreakersExternalVariables::RESET2);

      auto* ws = ws_.getData();

      ws[TRIP1]  = ZERO<RealT>;
      ws[RESET1] = ZERO<RealT>;
      ws[TRIP2]  = ZERO<RealT>;
      ws[RESET2] = ZERO<RealT>;
      std::fill(ws_indices_.begin(), ws_indices_.end(), INVALID_INDEX<IdxT>);

      if (auto trip1_port = ports_.in.template port<BranchBreakersSignalInputs::trip1>())
      {
        ws[TRIP1]          = trip1_port.readSignal();
        ws_indices_[TRIP1] = trip1_port.signalVariableIndex();
      }
      if (auto reset1_port = ports_.in.template port<BranchBreakersSignalInputs::reset1>())
      {
        ws[RESET1]          = reset1_port.readSignal();
        ws_indices_[RESET1] = reset1_port.signalVariableIndex();
      }
      if (auto trip2_port = ports_.in.template port<BranchBreakersSignalInputs::trip2>())
      {
        ws[TRIP2]          = trip2_port.readSignal();
        ws_indices_[TRIP2] = trip2_port.signalVariableIndex();
      }
      if (auto reset2_port = ports_.in.template port<BranchBreakersSignalInputs::reset2>())
      {
        ws[RESET2]          = reset2_port.readSignal();
        ws_indices_[RESET2] = reset2_port.signalVariableIndex();
      }

      auto* wb = wb_.getData();
      wb[0]    = Vr1();
      wb[1]    = Vi1();
      wb[2]    = Vr2();
      wb[3]    = Vi2();

      const auto* y  = y_.getData();
      const auto* yp = yp_.getData();
      auto*       f  = f_.getData();
      auto*       h  = h_.getData();

      evaluateInternalResidual(y, yp, wb, ws, f);
      evaluateBusResidual(y, yp, wb, h);

      Ir1() += h[0];
      Ii1() += h[1];
      Ir2() += h[2];
      Ii2() += h[3];

      ir1_ = h[0];
      ii1_ = h[1];
      ir2_ = h[2];
      ii2_ = h[3];

      if (bus1_->size() > 0)
      {
        bus1_->getResidual().setDataUpdated();
      }
      if (bus2_->size() > 0)
      {
        bus2_->getResidual().setDataUpdated();
      }
      f_.setDataUpdated();
      return 0;
    }

    /**
     * @brief Access the monitor
     *
     * @return Monitor for this model, or nullptr when the model was
     *         constructed without data.
     */
    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* BranchBreakers<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    /**
     * @brief Evaluate the breaker latches
     *
     * The body is kept free of branches and loops so sparse automatic
     * differentiation resolves a fixed structure.
     *
     * @param[in] y Internal variables in BranchBreakersInternalVariables order.
     * @param[in] yp Internal derivatives in the same enum order.
     * @param[in] wb Bus-1 then bus-2 \f$(V_{\mathrm{r}},V_{\mathrm{i}})\f$ voltage components.
     * @param[in] ws Signal values in BranchBreakersExternalVariables order.
     * @param[out] f Residuals in BranchBreakersInternalVariables order.
     * @return Zero on success.
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline int
    BranchBreakers<scalar_type, index_type>::evaluateInternalResidual(
        const ScalarT*                  y,
        const ScalarT*                  yp,
        [[maybe_unused]] const ScalarT* wb,
        const ScalarT*                  ws,
        ScalarT*                        f)
    {
      const auto Z1 = static_cast<size_t>(BranchBreakersInternalVariables::Z1);
      const auto Z2 = static_cast<size_t>(BranchBreakersInternalVariables::Z2);

      const auto TRIP1  = static_cast<size_t>(BranchBreakersExternalVariables::TRIP1);
      const auto RESET1 = static_cast<size_t>(BranchBreakersExternalVariables::RESET1);
      const auto TRIP2  = static_cast<size_t>(BranchBreakersExternalVariables::TRIP2);
      const auto RESET2 = static_cast<size_t>(BranchBreakersExternalVariables::RESET2);

      const ScalarT z1 = y[Z1];
      const ScalarT z2 = y[Z2];

      const ScalarT z1_dot = yp[Z1];
      const ScalarT z2_dot = yp[Z2];

      const ScalarT trip1  = ws[TRIP1];
      const ScalarT reset1 = ws[RESET1];
      const ScalarT trip2  = ws[TRIP2];
      const ScalarT reset2 = ws[RESET2];

      // Closed-state latches: trip is the latch reset, so trip has priority.
      f[Z1] = -z1_dot + Math::latch(z1, reset1, trip1) / Tlatch_;
      f[Z2] = -z2_dot + Math::latch(z2, reset2, trip2) / Tlatch_;

      return 0;
    }

    /**
     * @brief Evaluate the terminal currents into the buses
     *
     * @param[in] y Internal variables in BranchBreakersInternalVariables order.
     * @param[in] yp Internal derivatives in the same enum order.
     * @param[in] wb Bus-1 then bus-2 \f$(V_{\mathrm{r}},V_{\mathrm{i}})\f$ voltage components.
     * @param[out] h Bus-1 then bus-2 current injections.
     * @return Zero on success.
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline int
    BranchBreakers<scalar_type, index_type>::evaluateBusResidual(
        const ScalarT*                  y,
        [[maybe_unused]] const ScalarT* yp,
        const ScalarT*                  wb,
        ScalarT*                        h)
    {
      const auto Z1 = static_cast<size_t>(BranchBreakersInternalVariables::Z1);
      const auto Z2 = static_cast<size_t>(BranchBreakersInternalVariables::Z2);

      const ScalarT z1 = y[Z1];
      const ScalarT z2 = y[Z2];

      const ScalarT vr1 = wb[0];
      const ScalarT vi1 = wb[1];
      const ScalarT vr2 = wb[2];
      const ScalarT vi2 = wb[3];

      // Bilinear in the closed fractions; each corner is the Kron reduction of the open sides.
      const ScalarT u1  = Math::above(z1, HALF<RealT>);
      const ScalarT u2  = Math::above(z2, HALF<RealT>);
      const ScalarT g11 = u1 * (g11_ - (ONE<RealT> - u2) * gk1_);
      const ScalarT b11 = u1 * (b11_ - (ONE<RealT> - u2) * bk1_);
      const ScalarT g12 = u1 * u2 * g12_;
      const ScalarT b12 = u1 * u2 * b12_;
      const ScalarT g21 = u1 * u2 * g21_;
      const ScalarT b21 = u1 * u2 * b21_;
      const ScalarT g22 = u2 * (g22_ - (ONE<RealT> - u1) * gk2_);
      const ScalarT b22 = u2 * (b22_ - (ONE<RealT> - u1) * bk2_);

      h[0] = g11 * vr1 - b11 * vi1 + g12 * vr2 - b12 * vi2;
      h[1] = b11 * vr1 + g11 * vi1 + b12 * vr2 + g12 * vi2;
      h[2] = g21 * vr1 - b21 * vi1 + g22 * vr2 - b22 * vi2;
      h[3] = b21 * vr1 + g21 * vi1 + b22 * vr2 + g22 * vi2;

      return 0;
    }

    //
    //  Private methods
    //

    /**
     * @brief Read the parameters out of the model data
     *
     * No parameter is required; every parameter keeps the default documented
     * in the model README when omitted. A non-numeric or non-finite value is
     * counted and reported by verify() rather than throwing. Integer JSON
     * values are accepted for real parameters.
     *
     * @param[in] data Parameters and monitored-variable selections.
     */
    template <typename scalar_type, typename index_type>
    void BranchBreakers<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Params = typename ModelDataT::Parameters;

      parameter_error_count_ = 0;

      auto load_real = [&](auto key, RealT& target, const char* name)
      {
        if (!data.parameters.contains(key))
        {
          return;
        }

        const auto& value = data.parameters.at(key);
        RealT       parsed_value{};
        if (const auto* real_value = std::get_if<RealT>(&value))
        {
          parsed_value = *real_value;
        }
        else if (const auto* index_value = std::get_if<IdxT>(&value))
        {
          parsed_value = static_cast<RealT>(*index_value);
        }
        else
        {
          Log::error() << "BranchBreakers: parameter '" << name << "' must be numeric\n";
          ++parameter_error_count_;
          return;
        }

        if (!std::isfinite(parsed_value))
        {
          Log::error() << "BranchBreakers: parameter '" << name << "' must be finite\n";
          ++parameter_error_count_;
          return;
        }

        target = parsed_value;
      };

      load_real(Params::R, R_, "R");
      load_real(Params::X, X_, "X");
      load_real(Params::G, G_, "G");
      load_real(Params::B, B_, "B");
      load_real(Params::Gmag, Gmag_, "Gmag");
      load_real(Params::Bmag, Bmag_, "Bmag");
      load_real(Params::tap, tap_, "tap");
      load_real(Params::phase, phase_, "phase");
      load_real(Params::Tbrk, Tbrk_, "Tbrk");
      setDerivedParameters();
    }

    /**
     * @brief Bind the monitorable variables to the terminal currents
     */
    template <typename scalar_type, typename index_type>
    void BranchBreakers<scalar_type, index_type>::initializeMonitor()
    {
      using Variable = typename ModelDataT::MonitorableVariables;

      monitor_->set(Variable::ir1, [this]
                    {
                      ScalarT current[4];
                      terminalCurrents(current);
                      return current[0]; });
      monitor_->set(Variable::ii1, [this]
                    {
                      ScalarT current[4];
                      terminalCurrents(current);
                      return current[1]; });
      monitor_->set(Variable::im1, [this]
                    {
                      ScalarT current[4];
                      terminalCurrents(current);
                      return std::sqrt(current[0] * current[0] + current[1] * current[1]); });
      monitor_->set(Variable::p1, [this]
                    {
                      ScalarT current[4];
                      terminalCurrents(current);
                      return Vr1() * current[0] + Vi1() * current[1]; });
      monitor_->set(Variable::q1, [this]
                    {
                      ScalarT current[4];
                      terminalCurrents(current);
                      return Vi1() * current[0] - Vr1() * current[1]; });
      monitor_->set(Variable::ir2, [this]
                    {
                      ScalarT current[4];
                      terminalCurrents(current);
                      return current[2]; });
      monitor_->set(Variable::ii2, [this]
                    {
                      ScalarT current[4];
                      terminalCurrents(current);
                      return current[3]; });
      monitor_->set(Variable::im2, [this]
                    {
                      ScalarT current[4];
                      terminalCurrents(current);
                      return std::sqrt(current[2] * current[2] + current[3] * current[3]); });
      monitor_->set(Variable::p2, [this]
                    {
                      ScalarT current[4];
                      terminalCurrents(current);
                      return Vr2() * current[2] + Vi2() * current[3]; });
      monitor_->set(Variable::q2, [this]
                    {
                      ScalarT current[4];
                      terminalCurrents(current);
                      return Vi2() * current[2] - Vr2() * current[3]; });
    }

    /**
     * @brief Static method to log time constant warnings
     */
    template <typename scalar_type, typename index_type>
    void BranchBreakers<scalar_type, index_type>::logTimeConstantWarning()
    {
      Log::warning() << "BranchBreakers: Tbrk below "
                     << TIME_CONSTANT_MINIMUM
                     << " s is raised to that floor to keep the breaker latches well posed\n";
    }

    /**
     * @brief Resolve the closed admittance, the open-side Kron terms, and the latch time constant
     */
    template <typename scalar_type, typename index_type>
    void BranchBreakers<scalar_type, index_type>::setDerivedParameters()
    {
      if (Tbrk_ < ZERO<RealT>)
      {
        Log::error() << "BranchBreakers: Tbrk must be non-negative\n";
        ++parameter_error_count_;
      }
      if (Tbrk_ < TIME_CONSTANT_MINIMUM)
      {
        logTimeConstantWarning();
        Tbrk_ = TIME_CONSTANT_MINIMUM;
      }

      // A full command drives a latch through half travel in Tbrk.
      Tlatch_ = Tbrk_ / std::numbers::ln2_v<RealT>;

      g11_ = ZERO<RealT>;
      b11_ = ZERO<RealT>;
      g12_ = ZERO<RealT>;
      b12_ = ZERO<RealT>;
      g21_ = ZERO<RealT>;
      b21_ = ZERO<RealT>;
      g22_ = ZERO<RealT>;
      b22_ = ZERO<RealT>;
      gk1_ = ZERO<RealT>;
      bk1_ = ZERO<RealT>;
      gk2_ = ZERO<RealT>;
      bk2_ = ZERO<RealT>;

      const RealT denom = R_ * R_ + X_ * X_;
      if (denom == ZERO<RealT> || tap_ == ZERO<RealT>)
      {
        return;
      }

      const RealT g_br    = R_ / denom;
      const RealT b_br    = -X_ / denom;
      const RealT inv_tap = ONE<RealT> / tap_;
      const RealT cos_ph  = std::cos(phase_);
      const RealT sin_ph  = std::sin(phase_);

      const RealT g_diag = -g_br;
      const RealT b_diag = -b_br;

      g11_ = g_diag * inv_tap * inv_tap - HALF<RealT> * G_ - Gmag_;
      b11_ = b_diag * inv_tap * inv_tap - HALF<RealT> * B_ - Bmag_;

      g12_ = (g_br * cos_ph - b_br * sin_ph) * inv_tap;
      b12_ = (b_br * cos_ph + g_br * sin_ph) * inv_tap;

      g21_ = (g_br * cos_ph + b_br * sin_ph) * inv_tap;
      b21_ = (b_br * cos_ph - g_br * sin_ph) * inv_tap;

      g22_ = g_diag - HALF<RealT> * G_;
      b22_ = b_diag - HALF<RealT> * B_;

      // Open-side Kron terms; verify() rejects a zero diagonal entry.
      const RealT g_product = g12_ * g21_ - b12_ * b21_;
      const RealT b_product = g12_ * b21_ + b12_ * g21_;
      const RealT denom1    = g11_ * g11_ + b11_ * b11_;
      const RealT denom2    = g22_ * g22_ + b22_ * b22_;

      if (denom2 > ZERO<RealT>)
      {
        gk1_ = (g_product * g22_ + b_product * b22_) / denom2;
        bk1_ = (b_product * g22_ - g_product * b22_) / denom2;
      }
      if (denom1 > ZERO<RealT>)
      {
        gk2_ = (g_product * g11_ + b_product * b11_) / denom1;
        bk2_ = (b_product * g11_ - g_product * b11_) / denom1;
      }
    }

    /**
     * @brief Terminal currents into the buses at the present state and voltages
     *
     * @param[out] current Bus-1 then bus-2 current components.
     */
    template <typename scalar_type, typename index_type>
    void BranchBreakers<scalar_type, index_type>::terminalCurrents(ScalarT* current)
    {
      const ScalarT wb[4] = {Vr1(), Vi1(), Vr2(), Vi2()};
      evaluateBusResidual(y_.getData(), yp_.getData(), wb, current);
    }

    template <typename scalar_type, typename index_type>
    scalar_type& BranchBreakers<scalar_type, index_type>::Vr1()
    {
      return bus1_->Vr();
    }

    template <typename scalar_type, typename index_type>
    scalar_type& BranchBreakers<scalar_type, index_type>::Vi1()
    {
      return bus1_->Vi();
    }

    template <typename scalar_type, typename index_type>
    scalar_type& BranchBreakers<scalar_type, index_type>::Vr2()
    {
      return bus2_->Vr();
    }

    template <typename scalar_type, typename index_type>
    scalar_type& BranchBreakers<scalar_type, index_type>::Vi2()
    {
      return bus2_->Vi();
    }

    template <typename scalar_type, typename index_type>
    scalar_type& BranchBreakers<scalar_type, index_type>::Ir1()
    {
      return bus1_->Ir();
    }

    template <typename scalar_type, typename index_type>
    scalar_type& BranchBreakers<scalar_type, index_type>::Ii1()
    {
      return bus1_->Ii();
    }

    template <typename scalar_type, typename index_type>
    scalar_type& BranchBreakers<scalar_type, index_type>::Ir2()
    {
      return bus2_->Ir();
    }

    template <typename scalar_type, typename index_type>
    scalar_type& BranchBreakers<scalar_type, index_type>::Ii2()
    {
      return bus2_->Ii();
    }
  } // namespace PhasorDynamics
} // namespace GridKit
