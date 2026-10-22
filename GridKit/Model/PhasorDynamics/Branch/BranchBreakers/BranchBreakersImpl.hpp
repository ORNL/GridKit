/**
 * @file BranchBreakersImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the BranchBreakers phasor-dynamics branch model.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <format>
#include <numbers>

#include <GridKit/Model/ConfigurationChecks.hpp>
#include <GridKit/Model/ParameterReader.hpp>
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
     * @brief Construct a BranchBreakers from model data
     *
     * @param[in] bus1 Bus-1 terminal, the tapped side.
     * @param[in] bus2 Bus-2 terminal.
     * @param[in] data Parameters, signal assignments, and monitored-variable
     *                 selections.
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
     * @brief Destroy the BranchBreakers.
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
     * @brief Allocate the latches and wire the terminal-current outputs
     *
     * The terminal currents are the bus-residual contributions in h_. Each
     * attached output port publishes its entry of h_ as a computed output;
     * SystemModel assigns output_indices_ after the system variables.
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

      wb_.resize(4); // Vr1, Vi1, Vr2, Vi2
      wb_.setToZero();
      wb_indices_.assign(static_cast<size_t>(wb_.getSize()), INVALID_INDEX<IdxT>);

      const auto signal_size = Utilities::enum_size<BranchBreakersExternalVariables>();
      ws_.resize(static_cast<IdxT>(signal_size));
      ws_.setToZero();
      ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, j);
        this->setResidualIndex(j, j);
      }

      const auto output_size = Utilities::enum_size<BranchBreakersSignalOutputs>();
      h_.resize(static_cast<IdxT>(output_size));
      h_.setToZero();
      output_indices_.assign(output_size, INVALID_INDEX<IdxT>);

      auto* h = h_.getData();
      for (const auto output : Utilities::enum_values<BranchBreakersSignalOutputs>())
      {
        const auto k = static_cast<size_t>(output);
        if (auto port = ports_.out[output])
        {
          port.link(&h[k], &output_indices_[k]);
        }
      }

      allocated_ = true;
      return 0;
    }

    /**
     * @brief Validate the BranchBreakers configuration
     *
     * Checks the terminal buses, the branch impedance, the tap, the closed
     * self-admittances divided by in the open-end reductions, and attached
     * breaker commands.
     *
     * @return The configuration checks; passed() when valid.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks BranchBreakers<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks("BranchBreakers");

      checks.check(bus1_ != nullptr, "bus1 pointer is null");
      checks.check(bus2_ != nullptr, "bus2 pointer is null");
      checks.check(R_ * R_ + X_ * X_ > ZERO<RealT>, "R and X cannot both be zero");
      checks.check(tap_ > ZERO<RealT>, "tap must be positive");
      checks.check(g11_ * g11_ + b11_ * b11_ > ZERO<RealT>, "Y11 must be nonzero");
      checks.check(g22_ * g22_ + b22_ * b22_ > ZERO<RealT>, "Y22 must be nonzero");

      // An attached port must resolve to writable signal storage. The
      // enumerator is a template argument, so each port names itself once.
      auto check_attached_signal =
          [&]<BranchBreakersSignalInputs variable>(const char* name)
      {
        if (ports_.in.template port<variable>().connected()
            && !ports_.in.template port<variable>().linked())
        {
          checks.fail(std::format("{} signal attached with no linked source", name));
        }
      };

      check_attached_signal.template operator()<BranchBreakersSignalInputs::trip1>("trip1");
      check_attached_signal.template operator()<BranchBreakersSignalInputs::close1>("close1");
      check_attached_signal.template operator()<BranchBreakersSignalInputs::trip2>("trip2");
      check_attached_signal.template operator()<BranchBreakersSignalInputs::close2>("close2");

      return checks;
    }

    /**
     * @brief Initialize both breakers closed
     *
     * @return Zero on success; nonzero when the configuration is rejected.
     */
    template <typename scalar_type, typename index_type>
    int BranchBreakers<scalar_type, index_type>::initialize()
    {
      const auto Z1 = static_cast<size_t>(BranchBreakersInternalVariables::Z1);
      const auto Z2 = static_cast<size_t>(BranchBreakersInternalVariables::Z2);

      if (!this->verifyAndLog())
      {
        Log::error() << "BranchBreakers: cannot initialize with invalid configuration\n";
        return 1;
      }

      auto* y = y_.getData();
      y[Z1]   = ONE<RealT>;
      y[Z2]   = ONE<RealT>;

      y_.setDataUpdated();
      yp_.setToConst(static_cast<ScalarT>(ZERO<RealT>));

      // For DependencyTracking::Variable, set variable numbers
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
     * @brief Residuals of the breaker latches and the bus currents
     *
     * Refreshes the command buffer, evaluates the latch residuals, and adds
     * the terminal currents to the bus residuals. Unattached commands are
     * zero.
     *
     * @return Zero on success.
     */
    template <typename scalar_type, typename index_type>
    int BranchBreakers<scalar_type, index_type>::evaluateResidual()
    {
      const auto TRIP1  = static_cast<size_t>(BranchBreakersExternalVariables::TRIP1);
      const auto CLOSE1 = static_cast<size_t>(BranchBreakersExternalVariables::CLOSE1);
      const auto TRIP2  = static_cast<size_t>(BranchBreakersExternalVariables::TRIP2);
      const auto CLOSE2 = static_cast<size_t>(BranchBreakersExternalVariables::CLOSE2);

      const auto IR1 = static_cast<size_t>(BranchBreakersSignalOutputs::ir1);
      const auto II1 = static_cast<size_t>(BranchBreakersSignalOutputs::ii1);
      const auto IR2 = static_cast<size_t>(BranchBreakersSignalOutputs::ir2);
      const auto II2 = static_cast<size_t>(BranchBreakersSignalOutputs::ii2);

      ws_.setToZero();
      std::fill(ws_indices_.begin(), ws_indices_.end(), INVALID_INDEX<IdxT>);
      auto* ws = ws_.getData();

      if (auto trip1_port = ports_.in.template port<BranchBreakersSignalInputs::trip1>())
      {
        ws[TRIP1]          = trip1_port.readSignal();
        ws_indices_[TRIP1] = trip1_port.signalVariableIndex();
      }
      if (auto close1_port = ports_.in.template port<BranchBreakersSignalInputs::close1>())
      {
        ws[CLOSE1]          = close1_port.readSignal();
        ws_indices_[CLOSE1] = close1_port.signalVariableIndex();
      }
      if (auto trip2_port = ports_.in.template port<BranchBreakersSignalInputs::trip2>())
      {
        ws[TRIP2]          = trip2_port.readSignal();
        ws_indices_[TRIP2] = trip2_port.signalVariableIndex();
      }
      if (auto close2_port = ports_.in.template port<BranchBreakersSignalInputs::close2>())
      {
        ws[CLOSE2]          = close2_port.readSignal();
        ws_indices_[CLOSE2] = close2_port.signalVariableIndex();
      }

      evaluateInternalResidual(y_.getData(), yp_.getData(), nullptr, ws, f_.getData());
      f_.setDataUpdated();

      evaluateOutputs();

      const auto* h  = h_.getData();
      bus1_->Ir()   += h[IR1];
      bus1_->Ii()   += h[II1];
      bus2_->Ir()   += h[IR2];
      bus2_->Ii()   += h[II2];

      if (bus1_->size() > 0)
      {
        bus1_->getResidual().setDataUpdated();
      }
      if (bus2_->size() > 0)
      {
        bus2_->getResidual().setDataUpdated();
      }

      return 0;
    }

    /**
     * @brief Refresh the terminal voltages and the terminal currents in h_.
     *
     * @return Zero on success.
     */
    template <typename scalar_type, typename index_type>
    int BranchBreakers<scalar_type, index_type>::evaluateOutputs()
    {
      auto* wb = wb_.getData();
      wb[0]    = bus1_->Vr();
      wb[1]    = bus1_->Vi();
      wb[2]    = bus2_->Vr();
      wb[3]    = bus2_->Vi();

      return evaluateOutputs(y_.getData(), yp_.getData(), wb, ws_.getData(), h_.getData());
    }

    /**
     * @brief Access the monitor
     *
     * @return Monitor for this model.
     */
    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* BranchBreakers<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    /**
     * @brief Evaluate the breaker latch residuals.
     *
     * @param[in] y Breaker latches in BranchBreakersInternalVariables order.
     * @param[in] yp Latch derivatives in the same order.
     * @param[in] wb Unused; the latches do not read the buses.
     * @param[in] ws Commands in BranchBreakersExternalVariables order.
     * @param[out] f Residuals in BranchBreakersInternalVariables order.
     * @return Zero on success.
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline int BranchBreakers<scalar_type, index_type>::evaluateInternalResidual(
        const ScalarT*                  y,
        const ScalarT*                  yp,
        [[maybe_unused]] const ScalarT* wb,
        const ScalarT*                  ws,
        ScalarT*                        f)
    {
      const auto Z1 = static_cast<size_t>(BranchBreakersInternalVariables::Z1);
      const auto Z2 = static_cast<size_t>(BranchBreakersInternalVariables::Z2);

      const auto TRIP1  = static_cast<size_t>(BranchBreakersExternalVariables::TRIP1);
      const auto CLOSE1 = static_cast<size_t>(BranchBreakersExternalVariables::CLOSE1);
      const auto TRIP2  = static_cast<size_t>(BranchBreakersExternalVariables::TRIP2);
      const auto CLOSE2 = static_cast<size_t>(BranchBreakersExternalVariables::CLOSE2);

      f[Z1] = -yp[Z1] + Math::latch(y[Z1], ws[CLOSE1], ws[TRIP1]) / Tlatch_;
      f[Z2] = -yp[Z2] + Math::latch(y[Z2], ws[CLOSE2], ws[TRIP2]) / Tlatch_;

      return 0;
    }

    /**
     * @brief Evaluate the terminal currents entering both buses.
     *
     * @param[in] y Breaker latches.
     * @param[in] yp Latch derivatives, unused.
     * @param[in] wb Terminal voltages (Vr1, Vi1, Vr2, Vi2).
     * @param[in] ws Commands, unused.
     * @param[out] o Terminal currents in BranchBreakersSignalOutputs order.
     * @return Zero on success.
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline int BranchBreakers<scalar_type, index_type>::evaluateOutputs(
        const ScalarT*                  y,
        const ScalarT*                  yp,
        const ScalarT*                  wb,
        [[maybe_unused]] const ScalarT* ws,
        ScalarT*                        o)
    {
      const auto IR1 = static_cast<size_t>(BranchBreakersSignalOutputs::ir1);
      const auto II1 = static_cast<size_t>(BranchBreakersSignalOutputs::ii1);
      const auto IR2 = static_cast<size_t>(BranchBreakersSignalOutputs::ir2);
      const auto II2 = static_cast<size_t>(BranchBreakersSignalOutputs::ii2);

      ScalarT i11[2];
      ScalarT i12[2];
      ScalarT i21[2];
      ScalarT i22[2];

      evaluateBusResidual11(y, yp, wb, i11);
      evaluateBusResidual12(y, yp, wb + 2, i12);
      evaluateBusResidual21(y, yp, wb, i21);
      evaluateBusResidual22(y, yp, wb + 2, i22);

      o[IR1] = i11[0] + i12[0];
      o[II1] = i11[1] + i12[1];
      o[IR2] = i21[0] + i22[0];
      o[II2] = i21[1] + i22[1];

      return 0;
    }

    /**
     * @brief Bus-1 current from the bus-1 voltage
     *
     * With the bus-2 breaker open, bus 1 sees the open-end admittance.
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline int BranchBreakers<scalar_type, index_type>::evaluateBusResidual11(
        const ScalarT*                  y,
        [[maybe_unused]] const ScalarT* yp,
        const ScalarT*                  wb,
        ScalarT*                        h)
    {
      const auto Z1 = static_cast<size_t>(BranchBreakersInternalVariables::Z1);
      const auto Z2 = static_cast<size_t>(BranchBreakersInternalVariables::Z2);

      // 1 - above(z; 1/2) is below(z; 1/2). Reading z2 through both above()
      // and below() makes Enzyme auto-sparsity drop the z2 column.
      const ScalarT u1    = Math::above(y[Z1], HALF<RealT>);
      const ScalarT u2    = Math::above(y[Z2], HALF<RealT>);
      const ScalarT ubar2 = ONE<RealT> - u2;

      evaluateAdmittanceBlock(u1 * (u2 * g11_ + ubar2 * g1o_),
                              u1 * (u2 * b11_ + ubar2 * b1o_),
                              wb,
                              h);
      return 0;
    }

    /**
     * @brief Bus-1 current from the bus-2 voltage
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline int BranchBreakers<scalar_type, index_type>::evaluateBusResidual12(
        const ScalarT*                  y,
        [[maybe_unused]] const ScalarT* yp,
        const ScalarT*                  wb,
        ScalarT*                        h)
    {
      const auto Z1 = static_cast<size_t>(BranchBreakersInternalVariables::Z1);
      const auto Z2 = static_cast<size_t>(BranchBreakersInternalVariables::Z2);

      const ScalarT u12 = Math::above(y[Z1], HALF<RealT>) * Math::above(y[Z2], HALF<RealT>);

      evaluateAdmittanceBlock(u12 * g12_, u12 * b12_, wb, h);
      return 0;
    }

    /**
     * @brief Bus-2 current from the bus-1 voltage
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline int BranchBreakers<scalar_type, index_type>::evaluateBusResidual21(
        const ScalarT*                  y,
        [[maybe_unused]] const ScalarT* yp,
        const ScalarT*                  wb,
        ScalarT*                        h)
    {
      const auto Z1 = static_cast<size_t>(BranchBreakersInternalVariables::Z1);
      const auto Z2 = static_cast<size_t>(BranchBreakersInternalVariables::Z2);

      const ScalarT u12 = Math::above(y[Z1], HALF<RealT>) * Math::above(y[Z2], HALF<RealT>);

      evaluateAdmittanceBlock(u12 * g21_, u12 * b21_, wb, h);
      return 0;
    }

    /**
     * @brief Bus-2 current from the bus-2 voltage
     *
     * With the bus-1 breaker open, bus 2 sees the open-end admittance.
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline int BranchBreakers<scalar_type, index_type>::evaluateBusResidual22(
        const ScalarT*                  y,
        [[maybe_unused]] const ScalarT* yp,
        const ScalarT*                  wb,
        ScalarT*                        h)
    {
      const auto Z1 = static_cast<size_t>(BranchBreakersInternalVariables::Z1);
      const auto Z2 = static_cast<size_t>(BranchBreakersInternalVariables::Z2);

      // 1 - above(z; 1/2) is below(z; 1/2), as in evaluateBusResidual11().
      const ScalarT u1    = Math::above(y[Z1], HALF<RealT>);
      const ScalarT ubar1 = ONE<RealT> - u1;
      const ScalarT u2    = Math::above(y[Z2], HALF<RealT>);

      evaluateAdmittanceBlock(u2 * (u1 * g22_ + ubar1 * g2o_),
                              u2 * (u1 * b22_ + ubar1 * b2o_),
                              wb,
                              h);
      return 0;
    }

    /**
     * @brief Current through an admittance g + jb at one terminal voltage.
     *
     * @param[in] g Conductance.
     * @param[in] b Susceptance.
     * @param[in] wb Terminal voltage (Vr, Vi).
     * @param[out] h Current (Ir, Ii).
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline void BranchBreakers<scalar_type, index_type>::evaluateAdmittanceBlock(
        const ScalarT  g,
        const ScalarT  b,
        const ScalarT* wb,
        ScalarT*       h)
    {
      h[0] = g * wb[0] - b * wb[1];
      h[1] = b * wb[0] + g * wb[1];
    }

    /**
     * @brief Read the parameters out of the model data
     *
     * No parameter is required; every parameter keeps its default when
     * omitted. A non-numeric or non-finite real or a negative breaker time
     * throws std::invalid_argument. Integer JSON values are accepted for real
     * parameters.
     *
     * @param[in] data Parameters and monitored-variable selections.
     */
    template <typename scalar_type, typename index_type>
    void BranchBreakers<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Params = typename ModelDataT::Parameters;

      const Model::ParameterReader reader(data, "BranchBreakers");

      R_     = reader.readReal(Params::R).value_or(R_);
      X_     = reader.readReal(Params::X).value_or(X_);
      G_     = reader.readReal(Params::G).value_or(G_);
      B_     = reader.readReal(Params::B).value_or(B_);
      Gmag_  = reader.readReal(Params::Gmag).value_or(Gmag_);
      Bmag_  = reader.readReal(Params::Bmag).value_or(Bmag_);
      tap_   = reader.readReal(Params::tap).value_or(tap_);
      phase_ = reader.readReal(Params::phase).value_or(phase_);
      Tbrk_  = reader.readNonNegative(Params::Tbrk).value_or(Tbrk_);
      setDerivedParameters();
    }

    /**
     * @brief Bind the monitorable variables to the terminal currents
     */
    template <typename scalar_type, typename index_type>
    void BranchBreakers<scalar_type, index_type>::initializeMonitor()
    {
      using Variable = typename ModelDataT::MonitorableVariables;

      const auto IR1 = static_cast<size_t>(BranchBreakersSignalOutputs::ir1);
      const auto II1 = static_cast<size_t>(BranchBreakersSignalOutputs::ii1);
      const auto IR2 = static_cast<size_t>(BranchBreakersSignalOutputs::ir2);
      const auto II2 = static_cast<size_t>(BranchBreakersSignalOutputs::ii2);

      monitor_->set(Variable::ir1, [this]
                    { return terminalCurrents()[IR1]; });
      monitor_->set(Variable::ii1, [this]
                    { return terminalCurrents()[II1]; });
      monitor_->set(Variable::im1, [this]
                    {
                      const auto i = terminalCurrents();
                      return std::sqrt(i[IR1] * i[IR1] + i[II1] * i[II1]); });
      monitor_->set(Variable::p1, [this]
                    {
                      const auto i = terminalCurrents();
                      return bus1_->Vr() * i[IR1] + bus1_->Vi() * i[II1]; });
      monitor_->set(Variable::q1, [this]
                    {
                      const auto i = terminalCurrents();
                      return bus1_->Vi() * i[IR1] - bus1_->Vr() * i[II1]; });
      monitor_->set(Variable::ir2, [this]
                    { return terminalCurrents()[IR2]; });
      monitor_->set(Variable::ii2, [this]
                    { return terminalCurrents()[II2]; });
      monitor_->set(Variable::im2, [this]
                    {
                      const auto i = terminalCurrents();
                      return std::sqrt(i[IR2] * i[IR2] + i[II2] * i[II2]); });
      monitor_->set(Variable::p2, [this]
                    {
                      const auto i = terminalCurrents();
                      return bus2_->Vr() * i[IR2] + bus2_->Vi() * i[II2]; });
      monitor_->set(Variable::q2, [this]
                    {
                      const auto i = terminalCurrents();
                      return bus2_->Vi() * i[IR2] - bus2_->Vr() * i[II2]; });
    }

    /**
     * @brief Terminal currents at the current state
     *
     * Leaves the interface buffers untouched, so monitors do not disturb a
     * residual evaluation.
     *
     * @return Terminal currents in BranchBreakersSignalOutputs order.
     */
    template <typename scalar_type, typename index_type>
    auto BranchBreakers<scalar_type, index_type>::terminalCurrents() -> CurrentsT
    {
      const ScalarT wb[] = {bus1_->Vr(), bus1_->Vi(), bus2_->Vr(), bus2_->Vi()};
      CurrentsT     currents{};
      evaluateOutputs(y_.getData(), yp_.getData(), wb, nullptr, currents.data());
      return currents;
    }

    /**
     * @brief Resolve the breaker latch time and the closed and open-end admittances
     *
     * Degenerate data keeps zero admittances so the coefficients stay finite;
     * verify() reports it.
     */
    template <typename scalar_type, typename index_type>
    void BranchBreakers<scalar_type, index_type>::setDerivedParameters()
    {
      Tbrk_   = std::max(Tbrk_, TIME_CONSTANT_MINIMUM);
      Tlatch_ = Tbrk_ / std::numbers::ln2_v<RealT>;

      g11_ = ZERO<RealT>;
      b11_ = ZERO<RealT>;
      g12_ = ZERO<RealT>;
      b12_ = ZERO<RealT>;
      g21_ = ZERO<RealT>;
      b21_ = ZERO<RealT>;
      g22_ = ZERO<RealT>;
      b22_ = ZERO<RealT>;
      g1o_ = ZERO<RealT>;
      b1o_ = ZERO<RealT>;
      g2o_ = ZERO<RealT>;
      b2o_ = ZERO<RealT>;

      const RealT denom = R_ * R_ + X_ * X_;
      if (denom == ZERO<RealT> || tap_ == ZERO<RealT>)
      {
        return;
      }

      // Closed admittance, as in Branch
      const RealT g_br    = R_ / denom;
      const RealT b_br    = -X_ / denom;
      const RealT inv_tap = ONE<RealT> / tap_;
      const RealT cos_ph  = std::cos(phase_);
      const RealT sin_ph  = std::sin(phase_);

      g11_ = -g_br * inv_tap * inv_tap - HALF<RealT> * G_ - Gmag_;
      b11_ = -b_br * inv_tap * inv_tap - HALF<RealT> * B_ - Bmag_;
      g12_ = (g_br * cos_ph - b_br * sin_ph) * inv_tap;
      b12_ = (b_br * cos_ph + g_br * sin_ph) * inv_tap;
      g21_ = (g_br * cos_ph + b_br * sin_ph) * inv_tap;
      b21_ = (b_br * cos_ph - g_br * sin_ph) * inv_tap;
      g22_ = -g_br - HALF<RealT> * G_;
      b22_ = -b_br - HALF<RealT> * B_;

      const RealT y11_squared = g11_ * g11_ + b11_ * b11_;
      const RealT y22_squared = g22_ * g22_ + b22_ * b22_;
      if (y11_squared == ZERO<RealT> || y22_squared == ZERO<RealT>)
      {
        return;
      }

      // Y12 Y21, then Y11 - Y12 Y21 / Y22 and Y22 - Y21 Y12 / Y11
      const RealT g_product = g12_ * g21_ - b12_ * b21_;
      const RealT b_product = g12_ * b21_ + b12_ * g21_;

      g1o_ = g11_ - (g_product * g22_ + b_product * b22_) / y22_squared;
      b1o_ = b11_ - (b_product * g22_ - g_product * b22_) / y22_squared;
      g2o_ = g22_ - (g_product * g11_ + b_product * b11_) / y11_squared;
      b2o_ = b22_ - (b_product * g11_ - g_product * b11_) / y11_squared;
    }
  } // namespace PhasorDynamics
} // namespace GridKit
