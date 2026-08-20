#pragma once

#include <array>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Controller/REECB/Reecb.hpp>
#include <GridKit/Model/PhasorDynamics/Controller/REECB/ReecbData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/VariableMonitorController.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <typename scalar_type, typename index_type>
    class ControllerReecbTests
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      ControllerReecbTests()  = default;
      ~ControllerReecbTests() = default;

      static constexpr RealT kTol =
          static_cast<RealT>(100.0) * std::numeric_limits<RealT>::epsilon();

      /// Validate construction, row layout, defaults, parameters, buses,
      /// signal links, and the time-constant floor.
      TestOutcome validation()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 0.0);

        PhasorDynamics::Controller::Reecb<ScalarT, IdxT> empty(&bus);
        success *= (empty.size() == static_cast<IdxT>(Utilities::enum_size<Vars>()));
        success *= (empty.getMonitor() == nullptr);

        const std::array<Vars, Utilities::enum_size<Vars>()> row_order{{
            Vars::VMEAS,
            Vars::PMEAS,
            Vars::XPIQ,
            Vars::XPIV,
            Vars::QV,
            Vars::PORD,
            Vars::VT,
            Vars::VSAFE,
            Vars::SDIP,
            Vars::IQV,
            Vars::QREF,
            Vars::EQ,
            Vars::VPIQ,
            Vars::EPIV,
            Vars::RPORD,
            Vars::ILMAX,
            Vars::ILCAP,
            Vars::IQMAX,
            Vars::IPMAX,
            Vars::IQBASE,
            Vars::IQRAW,
            Vars::IQCMD,
            Vars::IPCMD,
        }};
        for (size_t row = 0; row < row_order.size(); ++row)
        {
          success *= (index(row_order[row]) == row);
        }

        Fixture<ScalarT> configured(makeData());
        success *= (configured.reecb.size() == static_cast<IdxT>(Utilities::enum_size<Vars>()));
        success *= (configured.reecb.getMonitor() != nullptr);
        success *= (configured.reecb.verify() == 0);
        success *= (configured.reecb.initialize() != 0);
        success *= (configured.reecb.allocate() == 0);
        success *= (configured.reecb.tagDifferentiable() == 0);
        success *= (static_cast<size_t>(configured.reecb.getResidual().getSize())
                    == Utilities::enum_size<Vars>());

        for (size_t row = 0; row < Utilities::enum_size<Vars>(); ++row)
        {
          const bool expected = row <= index(Vars::PORD);
          if (configured.reecb.tag()[row] != expected)
          {
            std::cout << "REECB differentiability tag " << row << " mismatch\n";
            success = false;
          }
        }

        Fixture<ScalarT> documented_defaults(makeMinimalData());
        success *= (documented_defaults.reecb.verify() == 0);
        success *= defaultsMatchDocumentedValues();

        // Integer JSON values are accepted for real parameters; booleans are
        // not numeric.
        auto integer_numeric                    = makeData();
        integer_numeric.parameters[Params::mva] = static_cast<IdxT>(100);
        Fixture<ScalarT> integer_parameter(integer_numeric);
        success *= (integer_parameter.reecb.verify() == 0);
        success *= invalidParameterCase(Params::mva, true);

        const RealT                  nan      = std::numeric_limits<RealT>::quiet_NaN();
        const RealT                  infinity = std::numeric_limits<RealT>::infinity();
        const std::array<Params, 26> real_parameters{{
            Params::mva,
            Params::Trv,
            Params::Tp,
            Params::Vref0,
            Params::Vdip,
            Params::Vup,
            Params::dbd1,
            Params::dbd2,
            Params::kqv,
            Params::Iql1,
            Params::Iqh1,
            Params::Qmax,
            Params::Qmin,
            Params::Kqp,
            Params::Kqi,
            Params::Vmax,
            Params::Vmin,
            Params::Kvp,
            Params::Kvi,
            Params::Tiq,
            Params::Tpord,
            Params::dPmax,
            Params::dPmin,
            Params::Pmax,
            Params::Pmin,
            Params::Imax,
        }};
        for (const Params parameter : real_parameters)
        {
          success *= invalidParameterCase(parameter, nan);
          success *= invalidParameterCase(parameter, infinity);
          success *= invalidParameterCase(parameter, -infinity);
        }

        success *= invalidParameterCase(Params::mva, 0.0);
        success *= invalidParameterCase(Params::Trv, -0.1);
        success *= invalidParameterCase(Params::Vdip, 2.0);
        success *= invalidParameterCase(Params::dbd1, 0.1);
        success *= invalidParameterCase(Params::dbd2, -0.1);
        success *= invalidParameterCase(Params::Iql1, 2.0);
        success *= invalidParameterCase(Params::Qmin, 3.0);
        success *= invalidParameterCase(Params::Vmin, 2.0);
        success *= invalidParameterCase(Params::dPmin, 0.0);
        success *= invalidParameterCase(Params::dPmax, 0.0);
        success *= invalidParameterCase(Params::Pmin, 3.0);
        success *= invalidParameterCase(Params::Imax, 0.0);

        const std::array<Params, 5> nonnegative_gains{{
            Params::kqv,
            Params::Kqp,
            Params::Kqi,
            Params::Kvp,
            Params::Kvi,
        }};
        for (const Params gain : nonnegative_gains)
        {
          success *= invalidParameterCase(gain, -0.1);
        }

        const std::array<Params, 4> flag_parameters{{
            Params::PfFlag,
            Params::VFlag,
            Params::QFlag,
            Params::Pqflag,
        }};
        const std::array<bool, 2>   valid_flag_values{{false, true}};
        const std::array<IdxT, 3>   invalid_integral_flag_values{{
            static_cast<IdxT>(0),
            static_cast<IdxT>(1),
            static_cast<IdxT>(2),
        }};
        const std::array<RealT, 5>  invalid_real_flag_values{{
            static_cast<RealT>(0.0),
            static_cast<RealT>(0.5),
            static_cast<RealT>(1.0),
            nan,
            infinity,
        }};
        for (const Params flag : flag_parameters)
        {
          for (const bool value : valid_flag_values)
          {
            auto data             = makeData();
            data.parameters[flag] = value;
            Fixture<ScalarT> model(data);
            success *= (model.reecb.verify() == 0);
          }

          for (const IdxT value : invalid_integral_flag_values)
          {
            success *= invalidParameterCase(flag, value);
          }

          for (const RealT value : invalid_real_flag_values)
          {
            success *= invalidParameterCase(flag, value);
          }
        }

        PhasorDynamics::Controller::Reecb<ScalarT, IdxT> busless(nullptr, makeData());
        busless.setSystemBase(kNominalFrequency, kSystemBaseVa);
        success *= (busless.verify() > 0);

        success *= unlinkedSignalRejected<Ext::pe>();
        success *= unlinkedSignalRejected<Ext::qgen>();
        success *= unlinkedSignalRejected<Ext::qext>();
        success *= unlinkedSignalRejected<Ext::pfaref>();
        success *= unlinkedSignalRejected<Ext::pref>();

        auto floor_data                      = makeData();
        floor_data.parameters[Params::Trv]   = 0.0;
        floor_data.parameters[Params::Tp]    = 0.0;
        floor_data.parameters[Params::Tiq]   = 0.0;
        floor_data.parameters[Params::Tpord] = 0.0;

        Fixture<ScalarT> floored(floor_data);
        success *= floored.initialize(kInitialIqcmd, kInitialIpcmd);
        success *= (floored.evaluate() == 0);
        success *= allResidualsWithinInitTolerance(floored.reecb);

        // Each floored lag turns a half-unit state offset into a rate of 500,
        // and saturates the active-power ramp limiter.
        setState(floored.reecb, {{Vars::VMEAS, 0.5}});
        success *= (floored.evaluate() == 0);
        success *= residualsMatch(floored.reecb,
                                  {{Vars::VMEAS, 500.0}},
                                  "floored voltage filter");

        setState(floored.reecb,
                 {{Vars::VMEAS, 1.0},
                  {Vars::PMEAS, 1.0},
                  {Vars::QV, 1.0},
                  {Vars::PORD, 1.0}});
        success *= (floored.evaluate() == 0);
        success *= residualsMatch(floored.reecb,
                                  {{Vars::PMEAS, 500.0},
                                   {Vars::QV, 500.0},
                                   {Vars::RPORD, 1.0}},
                                  "floored time constants");

        // The algebraic slew row supplies the differential order row.
        setState(floored.reecb, {{Vars::RPORD, 1.0}});
        success *= (floored.evaluate() == 0);
        success *= residualsMatch(floored.reecb,
                                  {{Vars::PORD, 1.0}, {Vars::RPORD, 0.0}},
                                  "floored active-power order rate");

        return success.report(__func__);
      }

      /// Check initialization state, known-input preservation, unknown-reference
      /// publication, command aliases, latches, monitors, and the power base.
      TestOutcome initializationAndSignals()
      {
        TestStatus success = true;

        Fixture<ScalarT> fixture(makeData(), 0.8, 0.6);
        fixture.attachAllInputs(99.0);
        fixture.input(Ext::pe)    = kInitialIpcmd;
        fixture.input(Ext::qgen)  = kInitialIqcmd;
        success                  *= fixture.initialize(kInitialIqcmd, kInitialIpcmd);
        success                  *= (fixture.evaluate() == 0);

        const std::array<VariableValue, Utilities::enum_size<Vars>()> initial_state{{
            {Vars::VMEAS, 1.0},
            {Vars::PMEAS, 1.5},
            {Vars::XPIQ, 0.0},
            {Vars::XPIV, 0.0},
            {Vars::QV, 1.5},
            {Vars::PORD, 1.5},
            {Vars::VT, 1.0},
            {Vars::VSAFE, 1.0},
            {Vars::SDIP, 1.0},
            {Vars::IQV, 0.0},
            {Vars::QREF, 1.5},
            {Vars::EQ, 0.0},
            {Vars::VPIQ, 0.5},
            {Vars::EPIV, 0.0},
            {Vars::RPORD, 0.0},
            {Vars::ILMAX, 2.0},
            {Vars::ILCAP, 2.0},
            {Vars::IQMAX, 2.0},
            {Vars::IPMAX, 2.5},
            {Vars::IQBASE, 0.0},
            {Vars::IQRAW, 1.5},
            {Vars::IQCMD, 0.75},
            {Vars::IPCMD, 0.75},
        }};
        success *= stateMatches(fixture.reecb, initial_state, "initialization");

        success *= scalarPreserved(fixture.iqcmd(), kInitialIqcmd, "preserved iqcmd");
        success *= scalarPreserved(fixture.ipcmd(), kInitialIpcmd, "preserved ipcmd");
        success *= scalarPreserved(fixture.input(Ext::pe), kInitialIpcmd, "preserved pe");
        success *= scalarPreserved(fixture.input(Ext::qgen), kInitialIqcmd, "preserved qgen");
        success *= scalarMatches(fixture.input(Ext::qext), 0.75, "published qext");
        success *= scalarMatches(fixture.input(Ext::pfaref), 0.0, "published pfaref");
        success *= scalarMatches(fixture.input(Ext::pref), 0.75, "published pref");
        success *= allResidualsWithinInitTolerance(fixture.reecb);

        success *= monitorMatches(fixture.reecb,
                                  {{kInitialIqcmd, kInitialIpcmd, 1.0, 1.5}},
                                  "initialization");

        constexpr RealT absolute_tolerance  = 2.5e-7;
        success                            *= (fixture.reecb.setAbsoluteTolerance(absolute_tolerance) == 0);
        const auto* tolerances              = fixture.reecb.absoluteTolerance().getData();
        for (size_t row = 0; row < Utilities::enum_size<Vars>(); ++row)
        {
          success *= valueUnchanged(tolerances[row], absolute_tolerance, "absolute tolerance", row);
        }

        // Unassigned command outputs keep the commands in the model vector.
        Fixture<ScalarT> latched(makeData(), 1.0, 0.0, kSystemBaseVa, false);
        success *= latched.initialize(kInitialIqcmd, kInitialIpcmd);
        success *= (latched.evaluate() == 0);
        success *= scalarPreserved(latched.iqcmd(), kInitialIqcmd, "unassigned iqcmd");
        success *= scalarPreserved(latched.ipcmd(), kInitialIpcmd, "unassigned ipcmd");
        success *= allResidualsWithinInitTolerance(latched.reecb);

        // Initialization preserves AD metadata on the owned command inputs.
        Fixture<DependencyTracking::Variable> tracked_commands(
            makeData(), 1.0, 0.0, kSystemBaseVa, false);
        success                  *= tracked_commands.prepare(kInitialIqcmd, kInitialIpcmd);
        auto*      tracked_state  = tracked_commands.reecb.y().getData();
        const auto iqcmd_index    = index(Vars::IQCMD);
        const auto ipcmd_index    = index(Vars::IPCMD);
        tracked_state[iqcmd_index].setVariableNumber(iqcmd_index);
        tracked_state[ipcmd_index].setVariableNumber(ipcmd_index);
        const auto iqcmd_dependencies  = tracked_state[iqcmd_index].getDependencies();
        const auto ipcmd_dependencies  = tracked_state[ipcmd_index].getDependencies();
        success                       *= (tracked_commands.reecb.initialize() == 0);
        success                       *= isEqual(tracked_state[iqcmd_index].getDependencies(), iqcmd_dependencies);
        success                       *= isEqual(tracked_state[ipcmd_index].getDependencies(), ipcmd_dependencies);

        // An omitted component rating falls back to the system power base, so
        // the same commands land on a different measured power.
        auto system_base_data = makeData();
        system_base_data.parameters.erase(Params::mva);
        Fixture<ScalarT> system_base(system_base_data, 1.0, 0.0, static_cast<RealT>(50.0e6));
        system_base.attachAllInputs();
        system_base.input(Ext::pe)  = 0.75;
        success                    *= system_base.initialize(kInitialIqcmd, 1.5);
        success                    *= (system_base.evaluate() == 0);
        success                    *= stateMatches(system_base.reecb,
                                                   {{Vars::PMEAS, 0.75},
                                                    {Vars::PORD, 1.5},
                                                    {Vars::ILMAX, 2.0}},
                                                   "omitted component rating");
        success                    *= allResidualsWithinInitTolerance(system_base.reecb);

        return success.report(__func__);
      }

      /// Check adjusted limits, initialization rejection, and atomicity.
      TestOutcome initializationDomain()
      {
        TestStatus success = true;

        const auto data = makeData();

        success *= initializationRejectedAtomically(data, 0.75, -0.1, "negative active-current command");
        success *= initializationRejectedAtomically(data, 0.75, std::numeric_limits<RealT>::infinity(), "nonfinite active-current command");
        success *= initializationRejectedAtomically(
            data, std::numeric_limits<RealT>::infinity(), 0.75, "nonfinite reactive-current command");

        auto pord_above                     = data;
        pord_above.parameters[Params::Pmax] = 1.0;
        Fixture<ScalarT> adjusted_pmax(pord_above);
        success *= adjusted_pmax.initialize(0.75, 0.75);
        success *= (adjusted_pmax.evaluate() == 0);
        success *= stateMatches(adjusted_pmax.reecb, {{Vars::PORD, 1.5}}, "adjusted Pmax");
        success *= allResidualsWithinInitTolerance(adjusted_pmax.reecb);
        setState(adjusted_pmax.reecb, {{Vars::PORD, 1.25}, {Vars::RPORD, 1.0}});
        success *= (adjusted_pmax.evaluate() == 0);
        success *= residualsMatch(adjusted_pmax.reecb, {{Vars::PORD, 1.0}}, "adjusted Pmax");

        auto pord_below                     = data;
        pord_below.parameters[Params::Pmin] = 2.0;
        Fixture<ScalarT> adjusted_pmin(pord_below);
        success *= adjusted_pmin.initialize(0.75, 0.75);
        success *= (adjusted_pmin.evaluate() == 0);
        success *= stateMatches(adjusted_pmin.reecb, {{Vars::PORD, 1.5}}, "adjusted Pmin");
        success *= allResidualsWithinInitTolerance(adjusted_pmin.reecb);
        setState(adjusted_pmin.reecb, {{Vars::PORD, 1.75}, {Vars::RPORD, -1.0}});
        success *= (adjusted_pmin.evaluate() == 0);
        success *= residualsMatch(adjusted_pmin.reecb, {{Vars::PORD, -1.0}}, "adjusted Pmin");

        auto expanded_current                     = data;
        expanded_current.parameters[Params::Imax] = 1.0;
        Fixture<ScalarT> adjusted_imax(expanded_current);
        success *= adjusted_imax.initialize(0.75, 0.75);
        success *= (adjusted_imax.evaluate() == 0);
        success *= stateMatches(adjusted_imax.reecb, {{Vars::ILMAX, 1.5}}, "adjusted Imax");
        success *= allResidualsWithinInitTolerance(adjusted_imax.reecb);

        auto reactive_pi                      = data;
        reactive_pi.parameters[Params::QFlag] = true;
        reactive_pi.parameters[Params::VFlag] = true;
        reactive_pi.parameters[Params::Kqi]   = 5.0;
        const std::array<RealT, 2> q_limits{{-1.25, 1.25}};
        for (const RealT qgen : q_limits)
        {
          Fixture<ScalarT> adjusted_q(reactive_pi);
          adjusted_q.attachAllInputs();
          adjusted_q.input(Ext::pe)    = 0.75;
          adjusted_q.input(Ext::qgen)  = qgen;
          success                     *= adjusted_q.initialize(0.75, 0.75);
          success                     *= (adjusted_q.evaluate() == 0);
          success                     *= allResidualsWithinInitTolerance(adjusted_q.reecb);
        }

        auto voltage_pi                    = reactive_pi;
        voltage_pi.parameters[Params::Kqi] = 0.0;
        voltage_pi.parameters[Params::Kvi] = 5.0;

        struct VoltageLimitCase
        {
          RealT voltage;
          RealT pe;
          RealT qgen;
        };

        const std::array<VoltageLimitCase, 2> voltage_limits{{
            {0.3, 0.3, 0.3},
            {1.6, 0.96, 0.8},
        }};
        for (const auto& test_case : voltage_limits)
        {
          Fixture<ScalarT> adjusted_v(voltage_pi, test_case.voltage);
          adjusted_v.attachAllInputs();
          adjusted_v.input(Ext::pe)    = test_case.pe;
          adjusted_v.input(Ext::qgen)  = test_case.qgen;
          success                     *= adjusted_v.initialize(0.75, 0.75);
          success                     *= (adjusted_v.evaluate() == 0);
          success                     *= allResidualsWithinInitTolerance(adjusted_v.reecb);
        }

        // Power-factor control needs a representable angle.
        auto power_factor                       = data;
        power_factor.parameters[Params::PfFlag] = true;

        auto late_data                     = power_factor;
        late_data.parameters[Params::Pmax] = 1.0;
        Fixture<ScalarT> late(late_data);
        late.attachAllInputs();
        late.input(Ext::pe)  = 0.0;
        success             *= late.prepare(0.75, 0.75);
        if (late.reecb.initialize() == 0)
        {
          std::cout << "Expected REECB initialization rejection after Pmax adjustment\n";
          success = false;
        }
        late.input(Ext::pref) = 0.75;
        setState(late.reecb, {{Vars::PORD, 1.25}, {Vars::VT, 1.0}});
        setDerivative(late.reecb, {{Vars::PORD, 0.0}});
        success *= (late.evaluate() == 0);
        success *= residualsMatch(late.reecb, {{Vars::PORD, 0.0}}, "rejected Pmax adjustment");
        success *= initializationRejectedAtomically(
            power_factor, 0.75, 0.75, "unrepresentable power-factor reference", 1.0e-8, 0.75);

        success *= initializationRejectedAtomically(data, 0.75, 0.75, "zero terminal voltage", 0.75, 0.75, 0.0);
        success *= initializationRejectedAtomically(
            data, 0.75, 0.75, "nonfinite active-power feedback", std::numeric_limits<RealT>::infinity(), 0.75);

        // Collapsed limits remain pinned at their initial output or expand to
        // include it.
        auto collapsed                     = reactive_pi;
        collapsed.parameters[Params::Kvi]  = 0.5;
        collapsed.parameters[Params::Qmin] = 1.2;
        collapsed.parameters[Params::Qmax] = 1.2;
        collapsed.parameters[Params::Vmin] = 1.0;
        collapsed.parameters[Params::Vmax] = 1.0;
        Fixture<ScalarT> collapsed_limits(collapsed);
        collapsed_limits.attachAllInputs();
        collapsed_limits.input(Ext::pe)    = 0.6;
        collapsed_limits.input(Ext::qgen)  = 0.6;
        success                           *= collapsed_limits.initialize(0.75, 0.75);
        success                           *= (collapsed_limits.evaluate() == 0);
        success                           *= allResidualsWithinInitTolerance(collapsed_limits.reecb);

        auto collapsed_reactive                     = collapsed;
        collapsed_reactive.parameters[Params::Vmin] = 0.5;
        collapsed_reactive.parameters[Params::Vmax] = 1.5;
        collapsed_reactive.parameters[Params::Kvi]  = 0.0;
        Fixture<ScalarT> expanded_reactive(collapsed_reactive);
        expanded_reactive.attachAllInputs();
        expanded_reactive.input(Ext::pe)    = 0.75;
        expanded_reactive.input(Ext::qgen)  = 0.3;
        success                            *= expanded_reactive.initialize(0.75, 0.75);
        success                            *= (expanded_reactive.evaluate() == 0);
        success                            *= allResidualsWithinInitTolerance(expanded_reactive.reecb);

        auto collapsed_voltage                     = collapsed;
        collapsed_voltage.parameters[Params::Qmin] = -2.0;
        collapsed_voltage.parameters[Params::Qmax] = 2.0;
        collapsed_voltage.parameters[Params::Kqi]  = 0.0;
        collapsed_voltage.parameters[Params::Vmin] = 1.4;
        collapsed_voltage.parameters[Params::Vmax] = 1.4;
        Fixture<ScalarT> expanded_voltage(collapsed_voltage);
        expanded_voltage.attachAllInputs();
        expanded_voltage.input(Ext::pe)    = 0.75;
        expanded_voltage.input(Ext::qgen)  = 0.75;
        success                           *= expanded_voltage.initialize(0.75, 0.75);
        success                           *= (expanded_voltage.evaluate() == 0);
        success                           *= allResidualsWithinInitTolerance(expanded_voltage.reecb);

        // Zero integral gains leave both controllers unconstrained, so any
        // reactive feedback initializes.
        auto zero_gains                    = reactive_pi;
        zero_gains.parameters[Params::Kqi] = 0.0;
        zero_gains.parameters[Params::Kvi] = 0.0;
        Fixture<ScalarT> unconstrained(zero_gains);
        unconstrained.attachAllInputs();
        unconstrained.input(Ext::pe)    = 0.75;
        unconstrained.input(Ext::qgen)  = 4.0;
        success                        *= unconstrained.initialize(0.75, 0.75);
        success                        *= (unconstrained.evaluate() == 0);
        success                        *= allResidualsWithinInitTolerance(unconstrained.reecb);

        // An invalid configuration is rejected before any state is written.
        auto invalid_data                     = data;
        invalid_data.parameters[Params::Imax] = 0.0;
        Fixture<ScalarT> invalid_fixture(invalid_data);
        invalid_fixture.attachAllInputs();
        success *= (invalid_fixture.reecb.allocate() == 0);
        poisonState(invalid_fixture, 0.75, 0.75);
        const auto invalid_y  = copyVector(invalid_fixture.reecb.y());
        const auto invalid_yp = copyVector(invalid_fixture.reecb.yp());
        if (invalid_fixture.reecb.initialize() == 0)
        {
          std::cout << "Expected REECB initialization rejection: invalid configuration\n";
          success = false;
        }
        success *= vectorUnchanged(invalid_fixture.reecb.y(), invalid_y, "state");
        success *= vectorUnchanged(invalid_fixture.reecb.yp(), invalid_yp, "derivative");

        return success.report(__func__);
      }

      /// The smooth-limiter inverse reproduces interior and boundary commands,
      /// including points that expand the current circle.
      TestOutcome initializationExactness()
      {
        TestStatus success = true;

        struct ExactnessCase
        {
          RealT       ipcmd;
          const char* label;
        };

        const std::array<ExactnessCase, 3> active_cases{{
            {1.0e-6, "near the lower active-current limit"},
            {0.75, "interior active-current command"},
            {1.249999, "near the upper active-current limit"},
        }};

        // The recovered order limits are widened so the reconstruction, not
        // the order limit, decides admissibility at the command endpoints.
        auto exactness_data                     = makeData();
        exactness_data.parameters[Params::Pmin] = -1.0;
        exactness_data.parameters[Params::Pmax] = 3.0;

        for (const auto& test_case : active_cases)
        {
          Fixture<ScalarT> fixture(exactness_data);
          success *= fixture.initialize(0.0, test_case.ipcmd);
          success *= (fixture.evaluate() == 0);
          success *= scalarPreserved(fixture.ipcmd(), test_case.ipcmd, test_case.label);
          success *= allResidualsWithinInitTolerance(fixture.reecb);
        }

        // An interior command recovers the ideal active-power order exactly.
        Fixture<ScalarT> interior(makeData());
        success *= interior.initialize(0.0, 0.75);
        success *= (interior.evaluate() == 0);
        success *= stateMatches(interior.reecb,
                                {{Vars::PORD, 1.5}},
                                "interior active-power order");

        // A command whose recovered order lands exactly on the order limit is
        // admitted rather than rejected.
        auto limit_data                     = makeData();
        limit_data.parameters[Params::Pmax] = 1.5;
        Fixture<ScalarT> at_limit(limit_data);
        success *= at_limit.initialize(0.0, 0.75);
        success *= (at_limit.evaluate() == 0);
        success *= stateMatches(at_limit.reecb, {{Vars::PORD, 1.5}}, "order at Pmax");
        success *= allResidualsWithinInitTolerance(at_limit.reecb);

        struct AsymmetricSlewCase
        {
          RealT       minimum;
          RealT       maximum;
          const char* label;
        };

        const std::array<AsymmetricSlewCase, 2> asymmetric_slews{{
            {-0.001, 0.1, "narrow negative ramp limit"},
            {-0.1, 0.001, "narrow positive ramp limit"},
        }};

        for (const auto& test_case : asymmetric_slews)
        {
          auto asymmetric_data                      = makeData();
          asymmetric_data.parameters[Params::dPmin] = test_case.minimum;
          asymmetric_data.parameters[Params::dPmax] = test_case.maximum;

          Fixture<ScalarT> asymmetric(asymmetric_data);
          asymmetric.attachAllInputs();
          success *= asymmetric.initialize(0.0, 0.75);
          success *= (asymmetric.evaluate() == 0);
          success *= stateMatches(asymmetric.reecb,
                                  {{Vars::PORD, 1.5}},
                                  test_case.label);
          success *= allResidualsWithinInitTolerance(asymmetric.reecb);
          success *= scalarMatches(asymmetric.input(Ext::pref),
                                   0.75,
                                   test_case.label);
        }

        // The reactive command shares the inverse, at both signs.
        const std::array<RealT, 2> reactive_commands{{
            static_cast<RealT>(0.999999),
            static_cast<RealT>(-0.999999),
        }};
        for (const RealT iqcmd : reactive_commands)
        {
          Fixture<ScalarT> reactive(exactness_data);
          success *= reactive.initialize(iqcmd, 0.75);
          success *= (reactive.evaluate() == 0);
          success *= scalarPreserved(reactive.iqcmd(), iqcmd, "near-limit reactive command");
          success *= allResidualsWithinInitTolerance(reactive.reecb);
        }

        struct BoundaryCase
        {
          bool        p_priority;
          RealT       iqcmd;
          RealT       ipcmd;
          RealT       ilmax;
          const char* label;
        };

        const std::array<BoundaryCase, 7> boundary_cases{{
            {true, 0.75, 0.0, 2.5, "zero active-current command"},
            {true, 0.0, 1.25, 0.0, "zero reactive-current capacity"},
            {true, 1.0, 0.75, 2.0, "upper reactive-current command"},
            {true, -1.0, 0.75, 2.0, "lower reactive-current command"},
            {false, 0.75, 1.0, 2.0, "upper active-current command"},
            {false, 1.25, 0.0, 0.0, "zero active-current capacity"},
            {true, 0.75, 1.5, 1.5, "expanded current circle"},
        }};
        for (const auto& test_case : boundary_cases)
        {
          auto boundary_data                       = exactness_data;
          boundary_data.parameters[Params::Pqflag] = test_case.p_priority;
          Fixture<ScalarT> boundary(boundary_data);
          success *= boundary.initialize(test_case.iqcmd, test_case.ipcmd);
          success *= (boundary.evaluate() == 0);
          success *= scalarPreserved(boundary.iqcmd(), test_case.iqcmd, test_case.label);
          success *= scalarPreserved(boundary.ipcmd(), test_case.ipcmd, test_case.label);
          success *= stateMatches(boundary.reecb, {{Vars::ILMAX, test_case.ilmax}}, test_case.label);
          success *= allResidualsWithinInitTolerance(boundary.reecb);
        }

        auto separated_data                     = exactness_data;
        separated_data.parameters[Params::Imax] = 1.0;
        Fixture<ScalarT> separated(separated_data);
        success *= separated.initialize(5.0e-13, 0.5);
        success *= (separated.evaluate() == 0);
        success *= scalarPreserved(separated.iqcmd(), 5.0e-13, "scale-separated current command");
        success *= allResidualsWithinInitTolerance(separated.reecb);

        // A low configured Imax requires representable bisection to preserve
        // a strict low-priority command.
        auto capacity_data                       = exactness_data;
        capacity_data.parameters[Params::mva]    = 100.0;
        capacity_data.parameters[Params::Pqflag] = true;
        capacity_data.parameters[Params::Imax]   = 0.1;

        const std::array<std::pair<RealT, RealT>, 2> capacity_cases{{
            {1.7, 0.3},
            {1.5, 0.5},
        }};
        for (const auto& [iqcmd, ipcmd] : capacity_cases)
        {
          Fixture<ScalarT> capacity_fixture(capacity_data);
          success           *= capacity_fixture.initialize(iqcmd, ipcmd);
          success           *= (capacity_fixture.evaluate() == 0);
          success           *= scalarPreserved(capacity_fixture.iqcmd(), iqcmd, "low-priority command");
          success           *= scalarPreserved(capacity_fixture.ipcmd(), ipcmd, "high-priority command");
          const RealT ilmax  = static_cast<RealT>(capacity_fixture.reecb.y().getData()[index(Vars::ILMAX)]);
          const RealT ilcap  = ilmax * ilmax
                               / std::sqrt(ilmax * ilmax + ReecbT::INITIALIZATION_TOLERANCE);
          if (ilcap < iqcmd)
          {
            std::cout << "REECB low-priority capacity does not include its initial command\n";
            success = false;
          }
          success *= allResidualsWithinInitTolerance(capacity_fixture.reecb);
        }

        auto nested_data                       = exactness_data;
        nested_data.parameters[Params::mva]    = 100.0;
        nested_data.parameters[Params::Pqflag] = true;
        nested_data.parameters[Params::QFlag]  = true;
        nested_data.parameters[Params::VFlag]  = false;
        nested_data.parameters[Params::Imax]   = 1.0;
        Fixture<ScalarT> nested(nested_data);
        nested.attachAllInputs();
        nested.input(Ext::pe)    = 0.6;
        nested.input(Ext::qgen)  = 0.8;
        success                 *= nested.initialize(0.8, 0.6);
        success                 *= (nested.evaluate() == 0);
        success                 *= scalarPreserved(nested.iqcmd(), 0.8, "nested-clamp reactive command");
        success                 *= scalarPreserved(nested.ipcmd(), 0.6, "nested-clamp active command");
        success                 *= allResidualsWithinInitTolerance(nested.reecb);

        auto exhausted_data                      = exactness_data;
        exhausted_data.parameters[Params::QFlag] = true;
        exhausted_data.parameters[Params::kqv]   = 1.0;
        exhausted_data.parameters[Params::Iql1]  = -0.4;
        exhausted_data.parameters[Params::Iqh1]  = 1.2;
        exhausted_data.parameters[Params::Vref0] = 2.2;
        Fixture<ScalarT> exhausted(exhausted_data);
        exhausted.attachAllInputs();
        exhausted.input(Ext::pe)  = 1.25;
        success                  *= exhausted.initialize(0.0, 1.25);
        success                  *= (exhausted.evaluate() == 0);
        success                  *= scalarPreserved(exhausted.iqcmd(), 0.0, "exhausted reactive-current capacity");
        success                  *= stateMatches(exhausted.reecb, {{Vars::ILMAX, 0.0}}, "injection does not expand current circle");
        success                  *= allResidualsWithinInitTolerance(exhausted.reecb);

        return success.report(__func__);
      }

      /// Check every residual row against an independent numerical answer key.
      /// The expected values are literals, not a second implementation of REECB.
      TestOutcome residualEquations()
      {
        TestStatus success = true;

        Fixture<ScalarT> fixture(makeResidualData(), kStateVr, kStateVi);
        fixture.attachAllInputs();
        setAnswerKeyInputs(fixture);
        success *= fixture.prepare(0.25, 0.35);
        setAnswerKeyState(fixture.reecb);
        success *= (fixture.evaluate() == 0);

        const std::array<VariableValue, Utilities::enum_size<Vars>()> expected_residuals{{
            {Vars::VMEAS, 0.99},
            {Vars::PMEAS, 0.145},
            {Vars::XPIQ, 0.21},
            {Vars::XPIV, 0.13},
            {Vars::QV, -0.05},
            {Vars::PORD, 0.26},
            {Vars::VT, -0.03},
            {Vars::VSAFE, 0.0},
            {Vars::SDIP, 0.0},
            {Vars::IQV, 0.0},
            {Vars::QREF, 0.0},
            {Vars::EQ, 0.0},
            {Vars::VPIQ, 0.0},
            {Vars::EPIV, 0.0},
            {Vars::RPORD, 0.0},
            {Vars::ILMAX, 0.32},
            {Vars::ILCAP, 0.0},
            {Vars::IQMAX, 0.0},
            {Vars::IPMAX, 0.0},
            {Vars::IQBASE, 0.0},
            {Vars::IQRAW, 0.0},
            {Vars::IQCMD, 0.19},
            {Vars::IPCMD, 0.05},
        }};

        success *= (static_cast<size_t>(fixture.reecb.getResidual().getSize())
                    == expected_residuals.size());
        for (size_t row = 0; row < expected_residuals.size(); ++row)
        {
          if (index(expected_residuals[row].variable) != row)
          {
            std::cout << "REECB residual key position " << row << " names row "
                      << variableName(expected_residuals[row].variable) << '\n';
            success = false;
          }
        }
        success *= residualsMatch(fixture.reecb,
                                  expected_residuals,
                                  "independent numerical answer key");

        return success.report(__func__);
      }

      /// Every selector combination initializes attached and unattached signals
      /// to a zero-residual state.
      TestOutcome selectorConfigurations()
      {
        TestStatus success = true;

        const std::array<bool, 2> selector_values{{false, true}};
        for (const bool pf : selector_values)
        {
          for (const bool voltage : selector_values)
          {
            for (const bool reactive : selector_values)
            {
              for (const bool p_priority : selector_values)
              {
                auto data                       = makeData();
                data.parameters[Params::PfFlag] = pf;
                data.parameters[Params::VFlag]  = voltage;
                data.parameters[Params::QFlag]  = reactive;
                data.parameters[Params::Pqflag] = p_priority;
                data.parameters[Params::Kqi]    = 0.0;
                data.parameters[Params::Kvi]    = 0.0;
                if (reactive)
                {
                  data.parameters[Params::Kvi] = 0.5;
                }
                if (reactive && voltage)
                {
                  data.parameters[Params::Kqi] = 0.4;
                }

                for (const bool attached : selector_values)
                {
                  Fixture<ScalarT> fixture(data);
                  if (attached)
                  {
                    fixture.attachAllInputs(7.0);
                    fixture.input(Ext::pe)   = 0.75;
                    fixture.input(Ext::qgen) = 0.75;
                  }

                  success *= fixture.initialize(0.75, 0.75);
                  success *= (fixture.evaluate() == 0);
                  success *= allResidualsWithinInitTolerance(fixture.reecb);
                  success *= scalarPreserved(fixture.iqcmd(), 0.75, "selector iqcmd");
                  success *= scalarPreserved(fixture.ipcmd(), 0.75, "selector ipcmd");
                  success *= stateMatches(fixture.reecb, {{Vars::ILMAX, 2.0}}, "selector ILMAX");

                  // Exactly one reactive path carries the operating point.
                  const auto* y = fixture.reecb.y().getData();
                  if (reactive)
                  {
                    success *= (y[index(Vars::XPIV)] != ZERO<RealT>);
                    if (voltage)
                    {
                      success *= (y[index(Vars::XPIQ)] != ZERO<RealT>);
                    }
                  }
                  else
                  {
                    success *= (y[index(Vars::QV)] != ZERO<RealT>);
                  }

                  if (attached)
                  {
                    success *= scalarPreserved(fixture.input(Ext::pe), 0.75, "selector pe");
                    success *= scalarPreserved(fixture.input(Ext::qgen), 0.75, "selector qgen");

                    RealT expected_qext = 0.75;
                    if (reactive && !voltage)
                    {
                      expected_qext = 1.0;
                    }
                    else if (pf)
                    {
                      expected_qext = 0.0;
                    }

                    RealT expected_pfaref = 0.0;
                    if (pf && (!reactive || voltage))
                    {
                      expected_pfaref = kUnitSlopeAngle;
                    }

                    success *= scalarMatches(fixture.input(Ext::qext), expected_qext, "published qext");
                    success *= scalarMatches(fixture.input(Ext::pfaref), expected_pfaref, "published pfaref");
                    success *= scalarMatches(fixture.input(Ext::pref), 0.75, "published pref");
                  }
                }
              }
            }
          }
        }

        return success.report(__func__);
      }

      /// Direct-voltage mode consumes and publishes the Volt/VAr reference
      /// without power-base conversion; the reactive modes convert it.
      TestOutcome voltVarReferenceBase()
      {
        TestStatus success = true;

        {
          auto data                      = makeData();
          data.parameters[Params::QFlag] = true;
          data.parameters[Params::VFlag] = false;
          data.parameters[Params::Kvi]   = 0.5;

          Fixture<ScalarT> fixture(data);
          fixture.attachAllInputs();
          success *= fixture.initialize(0.75, 0.75);
          success *= scalarMatches(fixture.input(Ext::qext), 1.0, "published voltage reference");
          success *= (fixture.evaluate() == 0);
          success *= allResidualsWithinInitTolerance(fixture.reecb);

          // A raised external voltage reference first enters the algebraic
          // voltage error without power-base conversion.
          fixture.input(Ext::qext)  = 1.2;
          success                  *= (fixture.evaluate() == 0);
          success                  *= residualsMatch(fixture.reecb,
                                                     {{Vars::EPIV, 0.2}},
                                                     "unconverted voltage reference");

          setState(fixture.reecb, {{Vars::EPIV, 0.2}});
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.reecb,
                                    {{Vars::XPIV, 0.1}},
                                    "unconverted voltage-reference rate");
        }

        {
          Fixture<ScalarT> fixture(makeData());
          fixture.attachAllInputs();
          success *= fixture.initialize(0.75, 0.75);
          success *= scalarMatches(fixture.input(Ext::qext), 0.75, "published system-base reactive power");
          success *= (fixture.evaluate() == 0);
          success *= allResidualsWithinInitTolerance(fixture.reecb);

          // The reactive reference keeps the power-base conversion before the
          // converted value drives the current-command lag.
          fixture.input(Ext::qext)  = 0.85;
          success                  *= (fixture.evaluate() == 0);
          success                  *= residualsMatch(fixture.reecb,
                                                     {{Vars::QREF, 0.2}},
                                                     "converted reactive reference");

          setState(fixture.reecb, {{Vars::QREF, 1.7}});
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.reecb,
                                    {{Vars::QV, 10.0}},
                                    "converted reactive-reference rate");
        }

        return success.report(__func__);
      }

      /// Check the reactive selector paths, the voltage-band gate, the
      /// reactive limits, both anti-windup gates, and the injection curve.
      TestOutcome reactiveControl()
      {
        TestStatus success = true;

        {
          // The constant-reactive path drives the current-command lag.
          auto data                      = makeResidualData();
          data.parameters[Params::QFlag] = false;
          data.parameters[Params::kqv]   = 0.0;
          Fixture<ScalarT> fixture(data);
          fixture.attachAllInputs();
          fixture.input(Ext::qext)  = 0.4;
          success                  *= fixture.prepare(0.0, 0.2);
          setControlState(fixture.reecb);
          setState(fixture.reecb, {{Vars::QV, 0.1}});
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.reecb,
                                    {{Vars::QREF, 0.8}},
                                    "constant-reactive reference");

          setState(fixture.reecb, {{Vars::QREF, 0.8}});
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.reecb, {{Vars::QV, 1.4}}, "constant-reactive lag");

          setState(fixture.reecb, {{Vars::VT, 0.0}});
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.reecb,
                                    {{Vars::SDIP, -1.0}},
                                    "constant-reactive voltage gate");

          setState(fixture.reecb, {{Vars::SDIP, 0.0}});
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.reecb, {{Vars::QV, 0.0}}, "gated constant-reactive lag");
        }

        {
          // Cascaded Volt/VAr control runs the reactive PI and bypasses the lag.
          auto data                      = makeResidualData();
          data.parameters[Params::QFlag] = true;
          data.parameters[Params::VFlag] = true;
          data.parameters[Params::kqv]   = 0.0;
          Fixture<ScalarT> fixture(data);
          fixture.attachAllInputs();
          fixture.input(Ext::qext)  = 0.1;
          fixture.input(Ext::qgen)  = -0.05;
          success                  *= fixture.prepare(0.0, 0.2);
          setControlState(fixture.reecb);
          setState(fixture.reecb, {{Vars::XPIQ, 0.82}});
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.reecb,
                                    {{Vars::QREF, 0.2}},
                                    "cascaded reactive-power reference");

          setState(fixture.reecb, {{Vars::QREF, 0.2}});
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.reecb,
                                    {{Vars::EQ, 0.3}},
                                    "cascaded reactive-power error");

          setState(fixture.reecb, {{Vars::EQ, 0.3}});
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.reecb,
                                    {{Vars::XPIQ, 0.12}, {Vars::QV, 0.0}},
                                    "reactive-power integral rate");

          setState(fixture.reecb, {{Vars::VT, 0.0}});
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.reecb,
                                    {{Vars::SDIP, -1.0}},
                                    "reactive-power voltage gate");

          setState(fixture.reecb, {{Vars::SDIP, 0.0}});
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.reecb, {{Vars::XPIQ, 0.0}}, "gated reactive-power integrator");
        }

        {
          // The direct-voltage reference bypasses the reactive PI.
          auto data                      = makeResidualData();
          data.parameters[Params::QFlag] = true;
          data.parameters[Params::VFlag] = false;
          data.parameters[Params::kqv]   = 0.0;
          Fixture<ScalarT> fixture(data);
          fixture.attachAllInputs();
          fixture.input(Ext::qext)  = 1.05;
          success                  *= fixture.prepare(0.0, 0.2);
          setControlState(fixture.reecb);
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.reecb,
                                    {{Vars::EPIV, 0.05}},
                                    "direct-voltage error");

          setState(fixture.reecb, {{Vars::EPIV, 0.05}});
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.reecb,
                                    {{Vars::XPIQ, 0.0}, {Vars::XPIV, 0.025}},
                                    "voltage-control integral rate");

          setState(fixture.reecb, {{Vars::VT, 2.0}});
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.reecb,
                                    {{Vars::SDIP, -1.0}},
                                    "voltage-control voltage gate");

          setState(fixture.reecb, {{Vars::SDIP, 0.0}});
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.reecb, {{Vars::XPIV, 0.0}}, "gated voltage-control integrator");
        }

        {
          // The reactive-power reference is limited before the error forms.
          auto data                      = makeResidualData();
          data.parameters[Params::QFlag] = true;
          data.parameters[Params::VFlag] = true;
          data.parameters[Params::kqv]   = 0.0;
          data.parameters[Params::Kqp]   = 0.0;

          struct ReactiveReferenceCase
          {
            RealT input;
            RealT limited_error;
            RealT expected_rate;
          };

          const std::array<ReactiveReferenceCase, 3> reference_cases{{
              {-0.6, -0.7, -0.28},
              {0.05, 0.1, 0.04},
              {0.6, 0.8, 0.32},
          }};
          for (const auto& test_case : reference_cases)
          {
            Fixture<ScalarT> fixture(data);
            fixture.attachAllInputs();
            fixture.input(Ext::qext)  = test_case.input;
            success                  *= fixture.prepare(0.0, 0.2);
            setControlState(fixture.reecb);
            setState(fixture.reecb,
                     {{Vars::XPIQ, 1.0}, {Vars::QREF, 2.0 * test_case.input}});
            success *= (fixture.evaluate() == 0);
            success *= residualsMatch(fixture.reecb,
                                      {{Vars::EQ, test_case.limited_error}},
                                      "reactive-power reference limit");

            setState(fixture.reecb, {{Vars::EQ, test_case.limited_error}});
            success *= (fixture.evaluate() == 0);
            success *= residualsMatch(fixture.reecb,
                                      {{Vars::XPIQ, test_case.expected_rate}},
                                      "reactive-power limited-error rate");
          }
        }

        {
          // Saturated probes sit beyond their limit by a margin, so a blocked
          // gate contributes nothing and an admitted gate passes the full rate.
          auto data                      = makeResidualData();
          data.parameters[Params::QFlag] = true;
          data.parameters[Params::VFlag] = true;
          data.parameters[Params::kqv]   = 0.0;
          data.parameters[Params::Kqp]   = 0.0;
          data.parameters[Params::Qmin]  = -3.0;
          data.parameters[Params::Qmax]  = 3.0;

          const std::array<AntiWindupCase, 5> antiwindup_cases{{
              {2.2, 1.0, 0.0},
              {2.2, -1.0, -0.8},
              {-0.2, -1.0, 0.0},
              {-0.2, 1.0, 0.8},
              {1.0, 1.0, 0.8},
          }};
          for (const auto& test_case : antiwindup_cases)
          {
            Fixture<ScalarT> fixture(data);
            fixture.attachAllInputs();
            fixture.input(Ext::qext)  = test_case.reference;
            success                  *= fixture.prepare(0.0, 0.2);
            setControlState(fixture.reecb);
            setState(fixture.reecb,
                     {{Vars::XPIQ, test_case.state},
                      {Vars::QREF, 2.0 * test_case.reference},
                      {Vars::EQ, 2.0 * test_case.reference}});
            success *= (fixture.evaluate() == 0);
            success *= residualsMatch(fixture.reecb,
                                      {{Vars::XPIQ, test_case.expected}},
                                      "reactive-power antiwindup");
          }
        }

        {
          // The voltage-control integrator saturates on the reactive-current
          // limit carried by the current circle.
          auto data                      = makeResidualData();
          data.parameters[Params::QFlag] = true;
          data.parameters[Params::VFlag] = false;
          data.parameters[Params::kqv]   = 0.0;
          data.parameters[Params::Kvp]   = 0.0;
          data.parameters[Params::Kvi]   = 1.0;

          const std::array<AntiWindupCase, 5> antiwindup_cases{{
              {1.0, 1.6, 0.0},
              {1.0, 0.4, -0.6},
              {-1.0, 0.4, 0.0},
              {-1.0, 1.6, 0.6},
              {0.0, 1.6, 0.6},
          }};
          for (const auto& test_case : antiwindup_cases)
          {
            Fixture<ScalarT> fixture(data);
            fixture.attachAllInputs();
            fixture.input(Ext::qext)  = test_case.reference;
            success                  *= fixture.prepare(0.0, 0.2);
            setControlState(fixture.reecb);
            setState(fixture.reecb,
                     {{Vars::XPIV, test_case.state},
                      {Vars::EPIV, test_case.reference - 1.0},
                      {Vars::ILMAX, 0.5},
                      {Vars::ILCAP, 0.5},
                      {Vars::IQMAX, 0.5}});
            success *= (fixture.evaluate() == 0);
            success *= residualsMatch(fixture.reecb,
                                      {{Vars::XPIV, test_case.expected}},
                                      "voltage-control antiwindup");
          }
        }

        {
          // With the reactive lag at zero the command row reads the injection
          // curve directly: a deadbanded voltage error scaled and limited.
          auto data                      = makeResidualData();
          data.parameters[Params::QFlag] = false;
          data.parameters[Params::kqv]   = 1.0;
          data.parameters[Params::dbd1]  = -0.6;
          data.parameters[Params::dbd2]  = 0.6;
          data.parameters[Params::Iql1]  = -1.2;
          data.parameters[Params::Iqh1]  = 1.5;
          data.parameters[Params::Vref0] = 3.0;

          const std::array<DrivenCase, 5> injection_cases{{
              {5.5, -1.2},
              {4.2, -0.6},
              {3.0, 0.0},
              {1.8, 0.6},
              {0.4, 1.5},
          }};
          for (const auto& test_case : injection_cases)
          {
            Fixture<ScalarT> fixture(data);
            success *= fixture.prepare(0.0, 0.2);
            setControlState(fixture.reecb);
            setState(fixture.reecb,
                     {{Vars::VMEAS, test_case.input},
                      {Vars::IQV, 0.0},
                      {Vars::IQRAW, 0.0},
                      {Vars::IQCMD, 0.0},
                      {Vars::ILMAX, 3.0},
                      {Vars::ILCAP, 3.0},
                      {Vars::IQMAX, 3.0}});
            success *= (fixture.evaluate() == 0);
            success *= residualsMatch(fixture.reecb,
                                      {{Vars::IQV, test_case.expected}},
                                      "reactive-current injection curve");

            setState(fixture.reecb, {{Vars::IQV, test_case.expected}});
            success *= (fixture.evaluate() == 0);
            success *= residualsMatch(fixture.reecb,
                                      {{Vars::IQRAW, test_case.expected}},
                                      "reactive-current injection sum");

            setState(fixture.reecb, {{Vars::IQRAW, test_case.expected}});
            success *= (fixture.evaluate() == 0);
            success *= residualsMatch(fixture.reecb,
                                      {{Vars::IQCMD, test_case.expected}},
                                      "reactive-current injection command");
          }
        }

        {
          // Power-factor control resolves the reactive reference from the
          // measured active power and the commanded angle.
          auto data                       = makeResidualData();
          data.parameters[Params::PfFlag] = true;
          data.parameters[Params::QFlag]  = false;
          data.parameters[Params::kqv]    = 0.0;
          Fixture<ScalarT> fixture(data);
          fixture.attachAllInputs();
          fixture.input(Ext::pfaref)  = std::atan(HALF<RealT>);
          success                    *= fixture.prepare(0.0, 0.2);
          setControlState(fixture.reecb);
          setState(fixture.reecb, {{Vars::PMEAS, 0.6}, {Vars::QV, 0.1}});
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.reecb,
                                    {{Vars::QREF, 0.3}},
                                    "power-factor reactive reference");

          setState(fixture.reecb, {{Vars::QREF, 0.3}});
          success *= (fixture.evaluate() == 0);
          success *= residualsMatch(fixture.reecb,
                                    {{Vars::QV, 0.4}},
                                    "power-factor reference");
        }

        return success.report(__func__);
      }

      /// Check the active-power ramp, its voltage gate and anti-windup, both
      /// command limits, the priority circle, and the signed continuation.
      TestOutcome activeCurrentControl()
      {
        TestStatus success = true;

        {
          // The ramp-rate limiter bounds the active-power order rate.
          const std::array<DrivenCase, 3> rate_cases{{
              {-1.0, -0.5},
              {0.2, 0.2},
              {1.0, 0.6},
          }};
          for (const auto& test_case : rate_cases)
          {
            Fixture<ScalarT> fixture(makeResidualData());
            fixture.attachAllInputs();
            fixture.input(Ext::pref)  = rampReference(0.5, test_case.input);
            success                  *= fixture.prepare(0.0, 0.2);
            setControlState(fixture.reecb);
            success *= (fixture.evaluate() == 0);
            success *= residualsMatch(fixture.reecb,
                                      {{Vars::RPORD, test_case.expected}},
                                      "active-power ramp limit");

            setState(fixture.reecb, {{Vars::RPORD, test_case.expected}});
            success *= (fixture.evaluate() == 0);
            success *= residualsMatch(fixture.reecb,
                                      {{Vars::PORD, test_case.expected}},
                                      "active-power order rate");
          }
        }

        {
          // Strongly asymmetric limits retain interior rates and bound each
          // direction independently.
          struct AsymmetricRateCase
          {
            RealT minimum;
            RealT maximum;
            RealT input;
            RealT expected;
          };

          const std::array<AsymmetricRateCase, 4> rate_cases{{
              {-0.001, 0.1, -0.002, -0.001},
              {-0.001, 0.1, 0.05, 0.05},
              {-0.1, 0.001, -0.05, -0.05},
              {-0.1, 0.001, 0.002, 0.001},
          }};

          for (const auto& test_case : rate_cases)
          {
            auto data                      = makeResidualData();
            data.parameters[Params::dPmin] = test_case.minimum;
            data.parameters[Params::dPmax] = test_case.maximum;

            Fixture<ScalarT> fixture(data);
            fixture.attachAllInputs();
            fixture.input(Ext::pref)  = rampReference(0.5, test_case.input);
            success                  *= fixture.prepare(0.0, 0.2);
            setControlState(fixture.reecb);
            success *= (fixture.evaluate() == 0);
            success *= residualsMatch(fixture.reecb,
                                      {{Vars::RPORD, test_case.expected}},
                                      "asymmetric active-power ramp limit");

            setState(fixture.reecb, {{Vars::RPORD, test_case.expected}});
            success *= (fixture.evaluate() == 0);
            success *= residualsMatch(fixture.reecb,
                                      {{Vars::PORD, test_case.expected}},
                                      "asymmetric active-power order rate");
          }
        }

        {
          // The voltage band gates the active-power order outside it.
          const std::array<DrivenCase, 3> gate_cases{{
              {0.0, 0.0},
              {1.0, 0.2},
              {2.0, 0.0},
          }};
          for (const auto& test_case : gate_cases)
          {
            Fixture<ScalarT> fixture(makeResidualData());
            fixture.attachAllInputs();
            fixture.input(Ext::pref)  = rampReference(0.5, 0.2);
            success                  *= fixture.prepare(0.0, 0.2);
            setControlState(fixture.reecb);
            setState(fixture.reecb, {{Vars::VT, test_case.input}});
            success          *= (fixture.evaluate() == 0);
            const RealT gate  = test_case.expected / static_cast<RealT>(0.2);
            success          *= residualsMatch(fixture.reecb,
                                               {{Vars::SDIP, gate - 1.0}},
                                               "active-power voltage gate");

            setState(fixture.reecb, {{Vars::SDIP, gate}, {Vars::RPORD, 0.2}});
            success *= (fixture.evaluate() == 0);
            success *= residualsMatch(fixture.reecb,
                                      {{Vars::PORD, test_case.expected}},
                                      "gated active-power order rate");
          }
        }

        {
          // Saturated probes sit beyond their limit by a margin, so a blocked
          // gate contributes nothing and an admitted gate passes the full rate.
          const std::array<AntiWindupCase, 5> antiwindup_cases{{
              {2.0, 1.0, 0.0},
              {2.0, -1.0, -0.5},
              {-0.4, -1.0, 0.0},
              {-0.4, 1.0, 0.6},
              {0.5, 1.0, 0.6},
          }};
          for (const auto& test_case : antiwindup_cases)
          {
            Fixture<ScalarT> fixture(makeResidualData());
            fixture.attachAllInputs();
            fixture.input(Ext::pref)  = rampReference(test_case.state, test_case.reference);
            success                  *= fixture.prepare(0.0, 0.2);
            setControlState(fixture.reecb);
            RealT limited_rate = 0.6;
            if (test_case.reference < ZERO<RealT>)
            {
              limited_rate = -0.5;
            }
            setState(fixture.reecb,
                     {{Vars::PORD, test_case.state}, {Vars::RPORD, limited_rate}});
            success *= (fixture.evaluate() == 0);
            success *= residualsMatch(fixture.reecb,
                                      {{Vars::PORD, test_case.expected}},
                                      "active-power antiwindup");
          }
        }

        {
          // P priority leaves the reactive command on the residual capacity.
          const std::array<DrivenCase, 5> reactive_limit_cases{{
              {-3.0, -2.0},
              {-1.0, -1.0},
              {0.0, 0.0},
              {1.0, 1.0},
              {3.0, 2.0},
          }};
          for (const auto& test_case : reactive_limit_cases)
          {
            Fixture<ScalarT> fixture(makeData());
            success *= fixture.prepare(0.0, 0.2);
            setControlState(fixture.reecb);
            setState(fixture.reecb,
                     {{Vars::ILMAX, 2.0},
                      {Vars::ILCAP, 2.0},
                      {Vars::IQMAX, 2.0},
                      {Vars::IQRAW, test_case.input},
                      {Vars::IQCMD, 0.0},
                      {Vars::QV, test_case.input}});
            success *= (fixture.evaluate() == 0);
            success *= residualsMatch(fixture.reecb,
                                      {{Vars::IQCMD, test_case.expected}},
                                      "reactive-command limit");
          }
        }

        {
          // Q priority leaves the active command on the residual capacity, and
          // the active command is one-sided.
          auto data                       = makeData();
          data.parameters[Params::Pqflag] = false;

          const std::array<DrivenCase, 4> active_limit_cases{{
              {-1.0, 0.0},
              {0.6, 0.6},
              {1.0, 1.0},
              {3.0, 2.0},
          }};
          for (const auto& test_case : active_limit_cases)
          {
            Fixture<ScalarT> fixture(data);
            success *= fixture.prepare(0.2, 0.0);
            setControlState(fixture.reecb);
            setState(fixture.reecb,
                     {{Vars::ILMAX, 2.0},
                      {Vars::ILCAP, 2.0},
                      {Vars::IPMAX, 2.0},
                      {Vars::IPCMD, 0.0},
                      {Vars::PORD, test_case.input}});
            success *= (fixture.evaluate() == 0);
            success *= residualsMatch(fixture.reecb,
                                      {{Vars::IPCMD, test_case.expected}},
                                      "active-command limit");
          }
        }

        {
          // The priority selector chooses which command consumes the circle.
          const std::array<std::pair<bool, RealT>, 2> priority_cases{{
              {true, 0.32},
              {false, 0.56},
          }};
          for (const auto& [p_priority, expected] : priority_cases)
          {
            auto data                       = makeResidualData();
            data.parameters[Params::Pqflag] = p_priority;
            Fixture<ScalarT> fixture(data, kStateVr, kStateVi);
            fixture.attachAllInputs();
            setAnswerKeyInputs(fixture);
            success *= fixture.prepare(0.25, 0.35);
            setAnswerKeyState(fixture.reecb);
            success           *= (fixture.evaluate() == 0);
            const char* label  = "Q-priority current circle";
            if (p_priority)
            {
              label = "P-priority current circle";
            }
            success *= residualsMatch(fixture.reecb,
                                      {{Vars::ILMAX, expected}},
                                      label);
          }
        }

        {
          // Selecting the priority command before forming the circle keeps an
          // overflowing inactive-command factor from producing NaN.
          const RealT               maximum = std::numeric_limits<RealT>::max();
          const RealT               limit   = maximum / 1024.0;
          const std::array<bool, 2> priorities{{false, true}};
          for (const bool p_priority : priorities)
          {
            auto data                       = makeData();
            data.parameters[Params::mva]    = 100.0;
            data.parameters[Params::Imax]   = limit;
            data.parameters[Params::Pqflag] = p_priority;

            Fixture<ScalarT> fixture(data);
            success *= fixture.prepare(0.0, 0.0);
            setControlState(fixture.reecb);
            RealT iqcmd = limit;
            RealT ipcmd = maximum;
            if (p_priority)
            {
              iqcmd = maximum;
              ipcmd = limit;
            }
            setState(fixture.reecb,
                     {{Vars::ILMAX, 0.0},
                      {Vars::IQCMD, iqcmd},
                      {Vars::IPCMD, ipcmd}});
            success *= (fixture.evaluate() == 0);
            success *= residualsMatch(fixture.reecb,
                                      {{Vars::ILMAX, 0.0}},
                                      "finite selected current circle");
            success *= allResidualsFinite(fixture.reecb);
          }
        }

        {
          // The signed-square continuation keeps a negative capacity iterate
          // finite, and its magnitude still bounds the low-priority command.
          auto data                       = makeData();
          data.parameters[Params::Imax]   = 1.0;
          data.parameters[Params::Pqflag] = true;

          const std::array<DrivenCase, 3> continuation_cases{{
              {-0.5, 1.0},
              {0.5, 0.5},
              {0.0, 0.75},
          }};
          for (const auto& test_case : continuation_cases)
          {
            Fixture<ScalarT> fixture(data);
            success *= fixture.prepare(0.0, 0.25);
            setControlState(fixture.reecb);
            setState(fixture.reecb,
                     {{Vars::ILMAX, test_case.input},
                      {Vars::ILCAP, 0.0},
                      {Vars::IQMAX, 0.0},
                      {Vars::IQRAW, 0.0},
                      {Vars::IPCMD, 0.25},
                      {Vars::IQCMD, 0.0},
                      {Vars::QV, 1.0}});
            success *= (fixture.evaluate() == 0);
            success *= residualsMatch(fixture.reecb,
                                      {{Vars::ILMAX, test_case.expected}},
                                      "signed capacity continuation");

            RealT expected_capacity = 0.0;
            if (test_case.input != ZERO<RealT>)
            {
              expected_capacity = 0.5;
            }
            success *= residualsMatch(fixture.reecb,
                                      {{Vars::ILCAP, expected_capacity}},
                                      "signed off-axis capacity");

            setState(fixture.reecb, {{Vars::ILCAP, expected_capacity}});
            success *= (fixture.evaluate() == 0);
            success *= residualsMatch(fixture.reecb,
                                      {{Vars::IQMAX, expected_capacity}},
                                      "reactive-current capacity");

            setState(fixture.reecb,
                     {{Vars::IQMAX, expected_capacity}, {Vars::IQRAW, 1.0}});
            success *= (fixture.evaluate() == 0);
            success *= residualsMatch(fixture.reecb,
                                      {{Vars::IQCMD, expected_capacity}},
                                      "capacity magnitude bound");
            success *= allResidualsFinite(fixture.reecb);
          }
        }

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// Enzyme and dependency tracking agree in every effective control mode
      /// and across the distinct nonlinear regions used by REECB.
      TestOutcome jacobian()
      {
        TestStatus success = true;

        struct ControlMode
        {
          const char* label;
          bool        pf_flag;
          bool        v_flag;
          bool        q_flag;
          bool        p_priority;
        };

        const std::array<ControlMode, 5> control_modes{{
            {"reactive-power reference with Q priority", false, false, false, false},
            {"power-factor reference with P priority", true, false, false, true},
            {"direct-voltage control with Q priority", false, false, true, false},
            {"cascaded reactive-power control with P priority", false, true, true, true},
            {"cascaded power-factor control with Q priority", true, true, true, false},
        }};
        for (const auto& mode : control_modes)
        {
          auto data                       = makeJacobianData();
          data.parameters[Params::PfFlag] = mode.pf_flag;
          data.parameters[Params::VFlag]  = mode.v_flag;
          data.parameters[Params::QFlag]  = mode.q_flag;
          data.parameters[Params::Pqflag] = mode.p_priority;

          success *= jacobiansMatch(
              dependencyTrackingJacobian(data, kNonunitAlpha, success),
              enzymeJacobian(data, kNonunitAlpha, success),
              mode.label);
        }

        struct CurrentCircleProbe
        {
          const char* label;
          bool        p_priority;
          RealT       ilmax;
        };

        const std::array<CurrentCircleProbe, 2> current_circle_probes{{
            {"negative signed current-circle capacity", true, -2.0},
            {"zero current-circle capacity", false, 0.0},
        }};
        for (const auto& probe : current_circle_probes)
        {
          auto data                        = makeJacobianData();
          data.parameters[Params::Pqflag]  = probe.p_priority;
          success                         *= jacobiansMatch(
              dependencyTrackingJacobian(
                  data, kNonunitAlpha, success, probe.ilmax),
              enzymeJacobian(data, kNonunitAlpha, success, probe.ilmax),
              probe.label);
        }

        struct NonlinearProbe
        {
          const char* label;
          RealT       mva;
          RealT       ilmax;
          RealT       epiv;
        };

        // The component base selects the active-order slew branch. Explicit
        // voltage-error states select either side of the moving current band.
        const std::array<NonlinearProbe, 4> nonlinear_probes{{
            {"upper active-order slew", 25.0, 2.0, 0.2},
            {"lower active-order slew", 100.0, 2.0, 0.2},
            {"below voltage-PI moving band", 50.0, 0.5, -1.0},
            {"above voltage-PI moving band", 50.0, 0.5, 1.0},
        }};
        for (const auto& probe : nonlinear_probes)
        {
          auto data                     = makeJacobianData();
          data.parameters[Params::mva]  = probe.mva;
          success                      *= jacobiansMatch(
              dependencyTrackingJacobian(
                  data, kNonunitAlpha, success, probe.ilmax, probe.epiv),
              enzymeJacobian(data, kNonunitAlpha, success, probe.ilmax, probe.epiv),
              probe.label);
        }

        auto injection_data                       = makeJacobianData();
        injection_data.parameters[Params::QFlag]  = false;
        injection_data.parameters[Params::kqv]    = 1.0;
        injection_data.parameters[Params::dbd1]   = -0.6;
        injection_data.parameters[Params::dbd2]   = 0.6;
        injection_data.parameters[Params::Iql1]   = -1.2;
        injection_data.parameters[Params::Iqh1]   = 1.5;
        injection_data.parameters[Params::Vref0]  = 2.2;
        success                                  *= jacobiansMatch(
            dependencyTrackingJacobian(injection_data, kNonunitAlpha, success),
            enzymeJacobian(injection_data, kNonunitAlpha, success),
            "voltage-error current injection");

        return success.report(__func__);
      }
#endif

    private:
      using Params = PhasorDynamics::Controller::ReecbParameters;
      using Vars   = PhasorDynamics::Controller::ReecbInternalVariables;
      using Mon    = PhasorDynamics::Controller::ReecbMonitorableVariables;
      using Data   = PhasorDynamics::Controller::ReecbData<RealT, IdxT>;
      using Ext    = typename Data::SignalInputs;
      using ReecbT = PhasorDynamics::Controller::Reecb<ScalarT, IdxT>;

      static constexpr size_t index(Vars variable)
      {
        return static_cast<size_t>(variable);
      }

      static constexpr size_t index(Ext variable)
      {
        return static_cast<size_t>(variable);
      }

      struct VariableValue
      {
        Vars  variable;
        RealT value;
      };

      struct DrivenCase
      {
        RealT input;
        RealT expected;
      };

      struct AntiWindupCase
      {
        RealT state;
        RealT reference;
        RealT expected;
      };

      /// Owns the terminal bus, REECB, the assigned command nodes, and the
      /// attached input nodes. Signal storage precedes the model so every
      /// referenced node outlives REECB; copying would invalidate the model
      /// and node pointers.
      template <typename T>
      class Fixture
      {
      private:
        std::array<T, Utilities::enum_size<Ext>()>                                   input_values_{};
        std::array<IdxT, Utilities::enum_size<Ext>()>                                input_indices_{};
        std::array<PhasorDynamics::SignalNode<T, IdxT>, Utilities::enum_size<Ext>()> input_nodes_{};

        PhasorDynamics::SignalNode<T, IdxT> iqcmd_node_;
        PhasorDynamics::SignalNode<T, IdxT> ipcmd_node_;
        bool                                commands_assigned_{true};

      public:
        explicit Fixture(const Data& data,
                         RealT       vr              = 1.0,
                         RealT       vi              = 0.0,
                         RealT       system_va_base  = kSystemBaseVa,
                         bool        assign_commands = true)
          : commands_assigned_(assign_commands),
            bus(static_cast<T>(vr), static_cast<T>(vi)),
            reecb(&bus, data)
        {
          reecb.setSystemBase(kNominalFrequency, system_va_base);
          if (commands_assigned_)
          {
            reecb.getPorts().out.template port<Data::SignalOutputs::iqcmd>().connect(&iqcmd_node_);
            reecb.getPorts().out.template port<Data::SignalOutputs::ipcmd>().connect(&ipcmd_node_);
          }
        }

        Fixture(const Fixture&)            = delete;
        Fixture& operator=(const Fixture&) = delete;

        void attachAllInputs(RealT initial_value = 0.0)
        {
          const IdxT external_index_base = reecb.size() + bus.size();
          for (auto variant : Utilities::enum_values<Ext>())
          {
            const auto port      = static_cast<size_t>(variant);
            input_values_[port]  = static_cast<T>(initial_value);
            input_indices_[port] = external_index_base + static_cast<IdxT>(port);
            input_nodes_[port].link(&input_values_[port], &input_indices_[port]);
          }

          reecb.getPorts().in.template port<Ext::pe>().connect(&input_nodes_[index(Ext::pe)]);
          reecb.getPorts().in.template port<Ext::qgen>().connect(&input_nodes_[index(Ext::qgen)]);
          reecb.getPorts().in.template port<Ext::qext>().connect(&input_nodes_[index(Ext::qext)]);
          reecb.getPorts().in.template port<Ext::pfaref>().connect(&input_nodes_[index(Ext::pfaref)]);
          reecb.getPorts().in.template port<Ext::pref>().connect(&input_nodes_[index(Ext::pref)]);
        }

        void setCommands(RealT iqcmd, RealT ipcmd)
        {
          if (commands_assigned_)
          {
            iqcmd_node_.init(static_cast<T>(iqcmd));
            ipcmd_node_.init(static_cast<T>(ipcmd));
            return;
          }
          auto* y               = reecb.y().getData();
          y[index(Vars::IQCMD)] = static_cast<T>(iqcmd);
          y[index(Vars::IPCMD)] = static_cast<T>(ipcmd);
          reecb.y().setDataUpdated();
        }

        /// Arrange the allocation, verification, bus, and command prerequisites.
        bool prepare(RealT iqcmd, RealT ipcmd)
        {
          const bool ready = (bus.allocate() == 0) && (reecb.allocate() == 0)
                             && (reecb.verify() == 0) && (bus.initialize() == 0);
          if (!ready)
          {
            std::cout << "REECB fixture preparation failed\n";
            return false;
          }
          setCommands(iqcmd, ipcmd);
          return true;
        }

        bool initialize(RealT iqcmd, RealT ipcmd)
        {
          if (!prepare(iqcmd, ipcmd))
          {
            return false;
          }
          if (reecb.initialize() != 0)
          {
            std::cout << "REECB initialization failed\n";
            return false;
          }
          return true;
        }

        int evaluate()
        {
          return reecb.evaluateResidual();
        }

        T iqcmd() const
        {
          return reecb.y().getData()[index(Vars::IQCMD)];
        }

        T ipcmd() const
        {
          return reecb.y().getData()[index(Vars::IPCMD)];
        }

        T& input(Ext port)
        {
          return input_values_[index(port)];
        }

        IdxT inputIndex(Ext port) const
        {
          return input_indices_[index(port)];
        }

        PhasorDynamics::Bus<T, IdxT>               bus;
        PhasorDynamics::Controller::Reecb<T, IdxT> reecb;
      };

      static constexpr RealT kSystemBaseVa     = static_cast<RealT>(100.0e6);
      static constexpr RealT kNominalFrequency = static_cast<RealT>(60.0);
      static constexpr RealT kStateVr          = 0.9;
      static constexpr RealT kStateVi          = 0.4;
      // The commands, current circle, and voltage give a well-conditioned
      // interior initialization point.
      static constexpr RealT kInitialIqcmd     = 0.75;
      static constexpr RealT kInitialIpcmd     = 0.75;
      static constexpr RealT kNonunitAlpha     = 0.7;

      inline static const RealT kUnitSlopeAngle = std::atan(ONE<RealT>);

      static constexpr size_t kBusVrColumn = Utilities::enum_size<Vars>();

      Data makeMinimalData() const
      {
        Data data;
        data.device_class          = "Reecb";
        data.disambiguation_string = "reecb_test";
        data.monitored_variables.insert(Mon::iqcmd);
        data.monitored_variables.insert(Mon::ipcmd);
        data.monitored_variables.insert(Mon::vmeas);
        data.monitored_variables.insert(Mon::pmeas);
        return data;
      }

      Data makeExplicitDefaultData() const
      {
        auto data = makeMinimalData();

        // These are the documented defaults; Vref0 has no fixed default and is
        // resolved from the terminal voltage.
        data.parameters[Params::mva]    = 100.0;
        data.parameters[Params::PfFlag] = false;
        data.parameters[Params::VFlag]  = false;
        data.parameters[Params::QFlag]  = false;
        data.parameters[Params::Pqflag] = false;
        data.parameters[Params::Trv]    = 0.02;
        data.parameters[Params::Tp]     = 0.0;
        data.parameters[Params::Vdip]   = 0.85;
        data.parameters[Params::Vup]    = 1.15;
        data.parameters[Params::dbd1]   = 0.0;
        data.parameters[Params::dbd2]   = 0.0;
        data.parameters[Params::kqv]    = 5.0;
        data.parameters[Params::Iql1]   = -1.1;
        data.parameters[Params::Iqh1]   = 1.1;
        data.parameters[Params::Qmax]   = 0.436;
        data.parameters[Params::Qmin]   = -0.436;
        data.parameters[Params::Kqp]    = 0.0;
        data.parameters[Params::Kqi]    = 0.1;
        data.parameters[Params::Vmax]   = 1.1;
        data.parameters[Params::Vmin]   = 0.9;
        data.parameters[Params::Kvp]    = 18.0;
        data.parameters[Params::Kvi]    = 5.0;
        data.parameters[Params::Tiq]    = 0.02;
        data.parameters[Params::Tpord]  = 0.02;
        data.parameters[Params::dPmax]  = 99.0;
        data.parameters[Params::dPmin]  = -99.0;
        data.parameters[Params::Pmax]   = 1.0;
        data.parameters[Params::Pmin]   = 0.0;
        data.parameters[Params::Imax]   = 1.3;
        return data;
      }

      /// The routine fixture: a half-size component base, wide bands, and
      /// limits far enough from the canonical commands that every smooth
      /// transition is saturated.
      Data makeData() const
      {
        auto data = makeMinimalData();

        data.parameters[Params::mva]    = 50.0;
        data.parameters[Params::PfFlag] = false;
        data.parameters[Params::VFlag]  = true;
        data.parameters[Params::QFlag]  = false;
        data.parameters[Params::Pqflag] = true;
        data.parameters[Params::Trv]    = 0.02;
        data.parameters[Params::Tp]     = 0.02;
        data.parameters[Params::Vref0]  = 1.0;
        data.parameters[Params::Vdip]   = 0.5;
        data.parameters[Params::Vup]    = 1.5;
        data.parameters[Params::dbd1]   = -0.2;
        data.parameters[Params::dbd2]   = 0.2;
        data.parameters[Params::kqv]    = 0.0;
        data.parameters[Params::Iql1]   = -1.0;
        data.parameters[Params::Iqh1]   = 1.0;
        data.parameters[Params::Qmax]   = 2.0;
        data.parameters[Params::Qmin]   = -2.0;
        data.parameters[Params::Kqp]    = 1.0;
        data.parameters[Params::Kqi]    = 0.0;
        data.parameters[Params::Vmax]   = 1.5;
        data.parameters[Params::Vmin]   = 0.5;
        data.parameters[Params::Kvp]    = 1.0;
        data.parameters[Params::Kvi]    = 0.0;
        data.parameters[Params::Tiq]    = 0.02;
        data.parameters[Params::Tpord]  = 0.02;
        data.parameters[Params::dPmax]  = 1.0;
        data.parameters[Params::dPmin]  = -1.0;
        data.parameters[Params::Pmax]   = 2.0;
        data.parameters[Params::Pmin]   = 0.0;
        data.parameters[Params::Imax]   = 2.5;
        return data;
      }

      /// Distinct nonzero values for every parameter. The bands are wide
      /// enough, and the lag reciprocals exact enough, for probe states to
      /// clear every smooth transition on an exact decimal.
      Data makeResidualData() const
      {
        auto data = makeData();

        data.parameters[Params::PfFlag] = false;
        data.parameters[Params::VFlag]  = true;
        data.parameters[Params::QFlag]  = true;
        data.parameters[Params::Pqflag] = true;
        data.parameters[Params::Trv]    = 0.2;
        data.parameters[Params::Tp]     = 0.4;
        data.parameters[Params::Vref0]  = 1.5;
        data.parameters[Params::kqv]    = 2.0;
        data.parameters[Params::Iql1]   = -0.4;
        data.parameters[Params::Iqh1]   = 0.5;
        data.parameters[Params::Qmax]   = 0.8;
        data.parameters[Params::Qmin]   = -0.7;
        data.parameters[Params::Kqp]    = 0.6;
        data.parameters[Params::Kqi]    = 0.4;
        data.parameters[Params::Vmax]   = 1.6;
        data.parameters[Params::Vmin]   = 0.4;
        data.parameters[Params::Kvp]    = 1.2;
        data.parameters[Params::Kvi]    = 0.5;
        data.parameters[Params::Tiq]    = 0.5;
        data.parameters[Params::Tpord]  = 0.25;
        data.parameters[Params::dPmax]  = 0.6;
        data.parameters[Params::dPmin]  = -0.5;
        data.parameters[Params::Pmax]   = 1.4;
        data.parameters[Params::Pmin]   = 0.1;
        data.parameters[Params::Imax]   = 1.5;
        return data;
      }

      /// The residual parameters with the injection disabled and the reactive
      /// limits widened, so the sensitivity probe sits interior everywhere.
      Data makeJacobianData() const
      {
        auto data                      = makeResidualData();
        data.parameters[Params::Vref0] = 1.0;
        data.parameters[Params::kqv]   = 0.0;
        data.parameters[Params::Qmin]  = -2.0;
        data.parameters[Params::Qmax]  = 2.0;
        data.parameters[Params::dPmin] = -0.001;
        data.parameters[Params::dPmax] = 0.1;
        data.parameters[Params::Imax]  = 2.5;
        return data;
      }

      /// The active-power reference that produces a requested pre-limit ramp
      /// rate at a given order, on the residual-parameter time constant.
      static constexpr RealT rampReference(RealT pord, RealT raw_rate)
      {
        return (pord + 0.25 * raw_rate) / 2.0;
      }

      template <typename T>
      void setAnswerKeyInputs(Fixture<T>& fixture) const
      {
        fixture.input(Ext::pe)     = static_cast<T>(0.3);
        fixture.input(Ext::qgen)   = static_cast<T>(-0.1);
        fixture.input(Ext::qext)   = static_cast<T>(0.2);
        fixture.input(Ext::pfaref) = static_cast<T>(0.15);
        fixture.input(Ext::pref)   = static_cast<T>(0.325);
      }

      /// The rich state shared by the residual answer key and the priority
      /// circle. Every smooth-transition argument keeps a saturation margin,
      /// so each row carries its ideal value.
      template <typename T>
      void setAnswerKeyState(PhasorDynamics::Controller::Reecb<T, IdxT>& reecb) const
      {
        setState(reecb,
                 {{Vars::VMEAS, 0.80},
                  {Vars::PMEAS, 0.55},
                  {Vars::XPIQ, 0.64},
                  {Vars::XPIV, -0.05},
                  {Vars::QV, 0.30},
                  {Vars::PORD, 0.60},
                  {Vars::VT, 1.00},
                  {Vars::VSAFE, 0.80},
                  {Vars::SDIP, 1.00},
                  {Vars::IQV, 0.50},
                  {Vars::QREF, 0.40},
                  {Vars::EQ, 0.60},
                  {Vars::VPIQ, 1.00},
                  {Vars::EPIV, 0.20},
                  {Vars::RPORD, 0.20},
                  {Vars::ILMAX, 1.20},
                  {Vars::ILCAP, 1.20},
                  {Vars::IQMAX, 1.20},
                  {Vars::IPMAX, 1.50},
                  {Vars::IQBASE, 0.19},
                  {Vars::IQRAW, 0.69},
                  {Vars::IQCMD, 0.25},
                  {Vars::IPCMD, 0.35}});
        setDerivative(reecb,
                      {{Vars::VMEAS, 0.01},
                       {Vars::PMEAS, -0.02},
                       {Vars::XPIQ, 0.03},
                       {Vars::XPIV, -0.03},
                       {Vars::QV, 0.05},
                       {Vars::PORD, -0.06}});
      }

      /// A neutral driven state for the control probes: unit voltage, cleared
      /// controller states, and a rested derivative.
      template <typename T>
      void setControlState(PhasorDynamics::Controller::Reecb<T, IdxT>& reecb) const
      {
        reecb.yp().setToConst(static_cast<T>(ZERO<RealT>));
        setState(reecb,
                 {{Vars::VMEAS, 1.0},
                  {Vars::PMEAS, 0.6},
                  {Vars::XPIQ, 0.0},
                  {Vars::XPIV, 0.0},
                  {Vars::QV, 0.0},
                  {Vars::PORD, 0.5},
                  {Vars::VT, 1.0},
                  {Vars::VSAFE, 1.0},
                  {Vars::SDIP, 1.0},
                  {Vars::IQV, 0.0},
                  {Vars::QREF, 0.0},
                  {Vars::EQ, 0.0},
                  {Vars::VPIQ, 0.0},
                  {Vars::EPIV, 0.0},
                  {Vars::RPORD, 0.0},
                  {Vars::ILMAX, 1.4},
                  {Vars::ILCAP, 1.4},
                  {Vars::IQMAX, 1.4},
                  {Vars::IPMAX, 1.5},
                  {Vars::IQBASE, 0.0},
                  {Vars::IQRAW, 0.0},
                  {Vars::IQCMD, 0.1},
                  {Vars::IPCMD, 0.2}});
        reecb.yp().setDataUpdated();
      }

      /// Populate every differential and algebraic variable for a Jacobian
      /// probe. The optional states place selected smooth branches explicitly.
      template <typename T>
      void setJacobianState(Fixture<T>& fixture,
                            RealT       ilmax,
                            RealT       epiv = static_cast<RealT>(0.2)) const
      {
        fixture.input(Ext::pe)     = static_cast<T>(0.25);
        fixture.input(Ext::qgen)   = static_cast<T>(0.5);
        fixture.input(Ext::qext)   = static_cast<T>(0.0);
        fixture.input(Ext::pfaref) = static_cast<T>(0.0);
        fixture.input(Ext::pref)   = static_cast<T>(0.25);

        fixture.reecb.yp().setToConst(static_cast<T>(ZERO<RealT>));
        setState(fixture.reecb,
                 {{Vars::VMEAS, 1.0},
                  {Vars::PMEAS, 0.5},
                  {Vars::XPIQ, 1.6},
                  {Vars::XPIV, 0.0},
                  {Vars::QV, 0.0},
                  {Vars::PORD, 0.5},
                  {Vars::VT, 1.0},
                  {Vars::VSAFE, 1.0},
                  {Vars::SDIP, 1.0},
                  {Vars::IQV, 0.0},
                  {Vars::QREF, 0.0},
                  {Vars::EQ, 0.5},
                  {Vars::VPIQ, 1.0},
                  {Vars::EPIV, epiv},
                  {Vars::RPORD, 0.05},
                  {Vars::ILMAX, ilmax},
                  {Vars::ILCAP, std::abs(ilmax)},
                  {Vars::IQMAX, std::abs(ilmax)},
                  {Vars::IPMAX, std::abs(ilmax)},
                  {Vars::IQBASE, 0.2},
                  {Vars::IQRAW, 0.2},
                  {Vars::IQCMD, 0.1},
                  {Vars::IPCMD, 0.2}});
        fixture.reecb.yp().setDataUpdated();
      }

      /// Omitting every parameter must give exactly the model built from the
      /// defaults the README documents, at rest and under load.
      bool defaultsMatchDocumentedValues() const
      {
        Fixture<ScalarT> implicit_defaults(makeMinimalData(), kStateVr, kStateVi);
        Fixture<ScalarT> explicit_defaults(makeExplicitDefaultData(), kStateVr, kStateVi);
        implicit_defaults.attachAllInputs();
        explicit_defaults.attachAllInputs();

        bool success = implicit_defaults.initialize(0.1, 0.2)
                       && explicit_defaults.initialize(0.1, 0.2);
        if (!success)
        {
          std::cout << "REECB documented-default comparison failed to initialize\n";
          return false;
        }

        if (implicit_defaults.evaluate() != 0 || explicit_defaults.evaluate() != 0)
        {
          success = false;
        }
        if (!vectorsMatch(implicit_defaults.reecb.y(),
                          explicit_defaults.reecb.y(),
                          "documented-default state"))
        {
          success = false;
        }
        if (!vectorsMatch(implicit_defaults.reecb.yp(),
                          explicit_defaults.reecb.yp(),
                          "documented-default derivative"))
        {
          success = false;
        }
        if (!vectorsMatch(implicit_defaults.reecb.getResidual(),
                          explicit_defaults.reecb.getResidual(),
                          "documented-default residual"))
        {
          success = false;
        }
        for (auto variable : Utilities::enum_values<Ext>())
        {
          const auto port = static_cast<size_t>(variable);
          if (!rowMatches(implicit_defaults.input(variable),
                          explicit_defaults.input(variable),
                          "documented-default signal",
                          port,
                          ""))
          {
            success = false;
          }
        }

        setAnswerKeyInputs(implicit_defaults);
        setAnswerKeyInputs(explicit_defaults);
        setAnswerKeyState(implicit_defaults.reecb);
        setAnswerKeyState(explicit_defaults.reecb);
        if (implicit_defaults.evaluate() != 0 || explicit_defaults.evaluate() != 0)
        {
          success = false;
        }
        if (!vectorsMatch(implicit_defaults.reecb.getResidual(),
                          explicit_defaults.reecb.getResidual(),
                          "documented-default dynamic residual"))
        {
          success = false;
        }
        return success;
      }

      template <typename ValueT>
      bool invalidParameterCase(Params parameter, ValueT value) const
      {
        auto data                  = makeData();
        data.parameters[parameter] = value;
        Fixture<ScalarT> fixture(data);
        return fixture.reecb.verify() > 0;
      }

      template <Ext variable>
      bool unlinkedSignalRejected() const
      {
        PhasorDynamics::SignalNode<ScalarT, IdxT> unlinked_node;
        Fixture<ScalarT>                          fixture(makeData());
        fixture.reecb.getPorts().in.template port<variable>().connect(&unlinked_node);
        return fixture.reecb.verify() > 0;
      }

      /// Fill state and derivative with a recognizable ramp, restoring the
      /// aliased commands, so any write by a rejected initialization shows.
      void poisonState(Fixture<ScalarT>& fixture, RealT iqcmd, RealT ipcmd) const
      {
        auto* y  = fixture.reecb.y().getData();
        auto* yp = fixture.reecb.yp().getData();
        for (size_t row = 0; row < Utilities::enum_size<Vars>(); ++row)
        {
          y[row]  = 0.125 + 0.01 * static_cast<RealT>(row);
          yp[row] = -0.25 - 0.01 * static_cast<RealT>(row);
        }
        fixture.setCommands(iqcmd, ipcmd);
        fixture.reecb.y().setDataUpdated();
        fixture.reecb.yp().setDataUpdated();
      }

      bool initializationRejectedAtomically(const Data& data,
                                            RealT       iqcmd,
                                            RealT       ipcmd,
                                            const char* label,
                                            RealT       pe      = 0.6,
                                            RealT       qgen    = 0.6,
                                            RealT       voltage = 1.0) const
      {
        Fixture<ScalarT> fixture(data, voltage);
        fixture.attachAllInputs(77.0);
        fixture.input(Ext::pe)   = pe;
        fixture.input(Ext::qgen) = qgen;
        if (!fixture.prepare(iqcmd, ipcmd))
        {
          return false;
        }

        poisonState(fixture, iqcmd, ipcmd);

        const auto                                     y_before   = copyVector(fixture.reecb.y());
        const auto                                     yp_before  = copyVector(fixture.reecb.yp());
        const auto                                     bus_before = copyVector(fixture.bus.y());
        std::array<RealT, Utilities::enum_size<Ext>()> inputs_before{};
        for (auto variable : Utilities::enum_values<Ext>())
        {
          const auto port     = static_cast<size_t>(variable);
          inputs_before[port] = fixture.input(variable);
        }

        bool success = true;
        if (fixture.reecb.initialize() == 0)
        {
          std::cout << "Expected REECB initialization rejection: " << label << '\n';
          success = false;
        }

        if (!scalarPreserved(fixture.iqcmd(), iqcmd, "rejected iqcmd preservation"))
        {
          success = false;
        }
        if (!scalarPreserved(fixture.ipcmd(), ipcmd, "rejected ipcmd preservation"))
        {
          success = false;
        }
        if (!vectorUnchanged(fixture.reecb.y(), y_before, "state"))
        {
          success = false;
        }
        if (!vectorUnchanged(fixture.reecb.yp(), yp_before, "derivative"))
        {
          success = false;
        }
        if (!vectorUnchanged(fixture.bus.y(), bus_before, "bus state"))
        {
          success = false;
        }
        for (auto variable : Utilities::enum_values<Ext>())
        {
          const auto port = static_cast<size_t>(variable);
          if (!valueUnchanged(fixture.input(variable),
                              inputs_before[port],
                              "external signal",
                              port))
          {
            success = false;
          }
        }
        return success;
      }

      template <typename T>
      void setState(PhasorDynamics::Controller::Reecb<T, IdxT>& reecb,
                    std::initializer_list<VariableValue>        values) const
      {
        auto* y = reecb.y().getData();
        for (const auto& [variable, value] : values)
        {
          y[index(variable)] = static_cast<T>(value);
        }
        reecb.y().setDataUpdated();
      }

      template <typename T>
      void setDerivative(PhasorDynamics::Controller::Reecb<T, IdxT>& reecb,
                         std::initializer_list<VariableValue>        values) const
      {
        auto* yp = reecb.yp().getData();
        for (const auto& [variable, value] : values)
        {
          yp[index(variable)] = static_cast<T>(value);
        }
        reecb.yp().setDataUpdated();
      }

      static const char* variableName(Vars variable)
      {
        static constexpr std::array<const char*, Utilities::enum_size<Vars>()> names{{
            "VMEAS",
            "PMEAS",
            "XPIQ",
            "XPIV",
            "QV",
            "PORD",
            "VT",
            "VSAFE",
            "SDIP",
            "IQV",
            "QREF",
            "EQ",
            "VPIQ",
            "EPIV",
            "RPORD",
            "ILMAX",
            "ILCAP",
            "IQMAX",
            "IPMAX",
            "IQBASE",
            "IQRAW",
            "IQCMD",
            "IPCMD",
        }};
        return names[index(variable)];
      }

      static bool variableMatches(RealT       actual,
                                  RealT       expected,
                                  const char* what,
                                  Vars        variable,
                                  const char* context,
                                  RealT       tolerance = kTol)
      {
        if (isEqual(actual, expected, tolerance))
        {
          return true;
        }
        std::cout << "REECB " << what << ' ' << variableName(variable);
        if (context[0] != '\0')
        {
          std::cout << ' ' << context;
        }
        std::cout << " mismatch: "
                  << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                  << actual << " != " << expected << '\n';
        return false;
      }

      static bool rowMatches(RealT       actual,
                             RealT       expected,
                             const char* what,
                             size_t      row,
                             const char* context,
                             RealT       tolerance = kTol)
      {
        if (isEqual(actual, expected, tolerance))
        {
          return true;
        }
        std::cout << "REECB " << what << " row " << row;
        if (context[0] != '\0')
        {
          std::cout << ' ' << context;
        }
        std::cout << " mismatch: "
                  << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                  << actual << " != " << expected << '\n';
        return false;
      }

      bool scalarMatches(RealT       actual,
                         RealT       expected,
                         const char* label,
                         RealT       tolerance = kTol) const
      {
        if (isEqual(actual, expected, tolerance))
        {
          return true;
        }
        std::cout << "REECB " << label << " mismatch: "
                  << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                  << actual << " != " << expected << '\n';
        return false;
      }

      /// A value retains exactly what its owner supplied, including signed
      /// infinities and NaN.
      static bool preserved(RealT actual, RealT expected)
      {
        if (std::isnan(expected))
        {
          return std::isnan(actual);
        }
        return actual == expected;
      }

      bool scalarPreserved(RealT actual, RealT expected, const char* label) const
      {
        if (preserved(actual, expected))
        {
          return true;
        }
        std::cout << "REECB " << label << " changed: "
                  << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                  << actual << " != " << expected << '\n';
        return false;
      }

      static bool valueUnchanged(RealT       actual,
                                 RealT       expected,
                                 const char* what,
                                 size_t      row)
      {
        if (preserved(actual, expected))
        {
          return true;
        }
        std::cout << "REECB " << what << ' ' << row << " changed: "
                  << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                  << actual << " != " << expected << '\n';
        return false;
      }

      template <typename VectorT, typename ValuesT>
      bool rowsMatch(const VectorT& vector,
                     const ValuesT& values,
                     const char*    what,
                     const char*    context) const
      {
        bool        success       = true;
        const auto* vector_values = vector.getData();
        for (const auto& [variable, expected] : values)
        {
          if (!variableMatches(static_cast<RealT>(vector_values[index(variable)]),
                               expected,
                               what,
                               variable,
                               context))
          {
            success = false;
          }
        }
        return success;
      }

      bool residualsMatch(const ReecbT&                        reecb,
                          std::initializer_list<VariableValue> values,
                          const char*                          context = "") const
      {
        return rowsMatch(reecb.getResidual(), values, "residual", context);
      }

      template <size_t size>
      bool residualsMatch(const ReecbT&                          reecb,
                          const std::array<VariableValue, size>& values,
                          const char*                            context = "") const
      {
        return rowsMatch(reecb.getResidual(), values, "residual", context);
      }

      bool stateMatches(const ReecbT&                        reecb,
                        std::initializer_list<VariableValue> values,
                        const char*                          context = "") const
      {
        return rowsMatch(reecb.y(), values, "state", context);
      }

      template <size_t size>
      bool stateMatches(const ReecbT&                          reecb,
                        const std::array<VariableValue, size>& values,
                        const char*                            context = "") const
      {
        return rowsMatch(reecb.y(), values, "state", context);
      }

      bool allResidualsWithinInitTolerance(const ReecbT& reecb) const
      {
        bool        success = true;
        const auto* f       = reecb.getResidual().getData();
        const auto* yp      = reecb.yp().getData();
        for (size_t row = 0; row < Utilities::enum_size<Vars>(); ++row)
        {
          const auto variable = static_cast<Vars>(row);
          if (!variableMatches(f[row],
                               0.0,
                               "residual",
                               variable,
                               "at rest",
                               ReecbT::INITIALIZATION_TOLERANCE))
          {
            success = false;
          }
          if (!valueUnchanged(yp[row], 0.0, "derivative", row))
          {
            success = false;
          }
        }
        return success;
      }

      bool allResidualsFinite(const ReecbT& reecb) const
      {
        bool        success = true;
        const auto* f       = reecb.getResidual().getData();
        for (size_t row = 0; row < Utilities::enum_size<Vars>(); ++row)
        {
          if (!std::isfinite(f[row]))
          {
            std::cout << "REECB residual " << variableName(static_cast<Vars>(row))
                      << " is not finite\n";
            success = false;
          }
        }
        return success;
      }

      bool monitorMatches(const ReecbT&               reecb,
                          const std::array<RealT, 4>& expected,
                          const char*                 context) const
      {
        RealT                                     time = 0.0;
        Model::VariableMonitorController<ScalarT> monitor(time);
        monitor.addMonitor(reecb.getMonitor());
        std::stringstream output;
        monitor.addSink({Model::VariableMonitorFormat::CSV}, output);
        monitor.start();
        monitor.print();
        monitor.stop();

        std::string header;
        std::string values_line;
        std::getline(output, header);
        std::getline(output, values_line);

        bool success = header == "t,Reecb_reecb_test_iqcmd,Reecb_reecb_test_ipcmd,"
                                 "Reecb_reecb_test_vmeas,Reecb_reecb_test_pmeas";

        const auto values = Tokenizer<RealT>(values_line, ',')();
        if (values.size() != expected.size() + 1)
        {
          std::cout << "REECB monitor emitted " << values.size()
                    << " values instead of " << expected.size() + 1 << '\n';
          return false;
        }

        for (size_t i = 0; i < expected.size(); ++i)
        {
          if (!rowMatches(values[i + 1], expected[i], "monitor", i, context))
          {
            success = false;
          }
        }
        return success;
      }

      template <typename VectorT>
      std::vector<RealT> copyVector(const VectorT& vector) const
      {
        const auto*        values = vector.getData();
        std::vector<RealT> snapshot(static_cast<size_t>(vector.getSize()));
        for (size_t row = 0; row < snapshot.size(); ++row)
        {
          snapshot[row] = static_cast<RealT>(values[row]);
        }
        return snapshot;
      }

      template <typename VectorT>
      bool vectorUnchanged(const VectorT&            vector,
                           const std::vector<RealT>& snapshot,
                           const char*               what) const
      {
        bool        success = true;
        const auto* values  = vector.getData();
        for (size_t row = 0; row < snapshot.size(); ++row)
        {
          if (!valueUnchanged(static_cast<RealT>(values[row]), snapshot[row], what, row))
          {
            success = false;
          }
        }
        return success;
      }

      template <typename LeftVectorT, typename RightVectorT>
      bool vectorsMatch(const LeftVectorT&  left,
                        const RightVectorT& right,
                        const char*         what) const
      {
        if (left.getSize() != right.getSize())
        {
          std::cout << "REECB " << what << " size mismatch\n";
          return false;
        }
        bool        success      = true;
        const auto* left_values  = left.getData();
        const auto* right_values = right.getData();
        for (size_t row = 0; row < static_cast<size_t>(left.getSize()); ++row)
        {
          if (!rowMatches(static_cast<RealT>(left_values[row]),
                          static_cast<RealT>(right_values[row]),
                          what,
                          row,
                          ""))
          {
            success = false;
          }
        }
        return success;
      }

      void numberVariables(Fixture<DependencyTracking::Variable>& fixture, RealT alpha) const
      {
        auto* y     = fixture.reecb.y().getData();
        auto* yp    = fixture.reecb.yp().getData();
        auto* bus_y = fixture.bus.y().getData();

        for (size_t row = 0; row < Utilities::enum_size<Vars>(); ++row)
        {
          y[row].setVariableNumber(row);
          yp[row].setVariableNumber(row);
          yp[row].scaleDependencies(alpha);
        }
        for (size_t row = 0; row < static_cast<size_t>(fixture.bus.size()); ++row)
        {
          bus_y[row].setVariableNumber(kBusVrColumn + row);
        }
        for (auto variable : Utilities::enum_values<Ext>())
        {
          const auto port = static_cast<size_t>(variable);
          fixture.input(variable).setVariableNumber(fixture.inputIndex(variable));
        }

        fixture.reecb.y().setDataUpdated();
        fixture.reecb.yp().setDataUpdated();
        fixture.bus.y().setDataUpdated();
      }

      std::vector<DependencyTracking::Variable::DependencyMap>
      dependencyTrackingJacobian(const Data& data,
                                 RealT       alpha,
                                 TestStatus& success,
                                 RealT       ilmax = 2.0,
                                 RealT       epiv  = static_cast<RealT>(0.2)) const
      {
        using DepVar = DependencyTracking::Variable;

        Fixture<DepVar> fixture(data, kStateVr, kStateVi);
        fixture.attachAllInputs();
        success *= fixture.prepare(0.0, 0.2);
        setJacobianState(fixture, ilmax, epiv);
        numberVariables(fixture, alpha);
        success *= (fixture.evaluate() == 0);

        std::vector<DependencyTracking::Variable::DependencyMap> rows(Utilities::enum_size<Vars>());
        const auto*                                              f = fixture.reecb.getResidual().getData();
        for (size_t row = 0; row < rows.size(); ++row)
        {
          rows[row] = f[row].getDependencies();
        }
        return rows;
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      std::vector<DependencyTracking::Variable::DependencyMap>
      enzymeJacobian(const Data& data,
                     RealT       alpha,
                     TestStatus& success,
                     RealT       ilmax = 2.0,
                     RealT       epiv  = static_cast<RealT>(0.2)) const
      {
        Fixture<ScalarT> fixture(data, kStateVr, kStateVi);
        fixture.attachAllInputs();
        success *= fixture.prepare(0.0, 0.2);

        for (IdxT row = 0; row < fixture.bus.size(); ++row)
        {
          fixture.bus.setVariableIndex(row, fixture.reecb.size() + row);
        }

        setJacobianState(fixture, ilmax, epiv);
        fixture.reecb.updateTime(0.0, alpha);
        success *= (fixture.evaluate() == 0);
        success *= (fixture.reecb.evaluateJacobian() == 0);
        success *= (fixture.reecb.constructCsr() == 0);
        return MapFromCsr(fixture.reecb.getCsrJacobian());
      }

      bool jacobiansMatch(
          const std::vector<DependencyTracking::Variable::DependencyMap>& dependency,
          const std::vector<DependencyTracking::Variable::DependencyMap>& enzyme,
          const char*                                                     probe) const
      {
        if (dependency.size() != enzyme.size())
        {
          std::cout << "REECB Jacobian row-count mismatch for " << probe << '\n';
          return false;
        }

        bool success = true;
        for (size_t row = 0; row < dependency.size(); ++row)
        {
          if (!isEqual(dependency[row], enzyme[row], kTol))
          {
            std::cout << "REECB Jacobian row " << row
                      << " mismatch between dependency tracking and Enzyme for "
                      << probe << '\n';
            success = false;
          }
        }
        return success;
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
