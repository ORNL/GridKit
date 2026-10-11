#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Load/PERC1/Perc1.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Model/VariableMonitorController.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit
{
  namespace Testing
  {
    class LoadPerc1Tests
    {
      using IdxT                  = size_t;
      using Data                  = PhasorDynamics::Perc1Data<double, IdxT>;
      using Params                = Data::Parameters;
      using Vars                  = PhasorDynamics::Perc1InternalVariables;
      using DepVar                = DependencyTracking::Variable;
      static constexpr size_t N   = Utilities::enum_size<Vars>();
      static constexpr size_t ND  = static_cast<size_t>(Vars::V);
      static constexpr double TOL = 2.0e-11;
      static constexpr double TAU = 1.0 / (Math::MU<double> * Math::MU<double>);

      template <typename T>
      struct Fixture
      {
        PhasorDynamics::Bus<T, IdxT>        bus;
        T                                   frequency{1};
        IdxT                                frequency_index{N + 2};
        PhasorDynamics::SignalNode<T, IdxT> frequency_node;
        PhasorDynamics::Perc1<T, IdxT>      model;

        explicit Fixture(const Data& data, double vr = 1.0, double vi = 0.0)
          : bus(static_cast<T>(vr), static_cast<T>(vi)), model(&bus, data)
        {
        }

        Fixture(const Fixture&)            = delete;
        Fixture& operator=(const Fixture&) = delete;

        void attachFrequency(double value)
        {
          frequency = value;
          frequency_node.link(&frequency, &frequency_index);
          model.getPorts().in.template port<Data::SignalInputs::freq>().connect(&frequency_node);
        }

        bool prepare()
        {
          return bus.allocate() == 0 && model.allocate() == 0 && bus.initialize() == 0;
        }

        bool initialize()
        {
          return prepare() && model.initialize() == 0;
        }

        bool evaluate()
        {
          return bus.evaluateResidual() == 0 && model.evaluateResidual() == 0;
        }

        T& y(Vars variable)
        {
          return model.y().getData()[static_cast<size_t>(variable)];
        }

        T& f(Vars variable)
        {
          return model.getResidual().getData()[static_cast<size_t>(variable)];
        }
      };

      static Data makeData()
      {
        Data data;
        data.device_class             = "Perc1";
        data.disambiguation_string    = "load";
        data.parameters[Params::Pnom] = 0.4;
        data.parameters[Params::Qnom] = 0.1;
        for (const auto monitor : Utilities::enum_values<Data::MonitorableVariables>())
        {
          data.monitored_variables.insert(monitor);
        }
        return data;
      }

      static bool matches(double actual, double expected, const char* label, double tolerance = TOL)
      {
        if (std::isfinite(actual) && std::abs(actual - expected) <= tolerance * std::max(1.0, std::abs(expected)))
        {
          return true;
        }
        std::cout << label << ": " << actual << " != " << expected << '\n';
        return false;
      }

      template <typename T>
      static bool zeroResidual(Fixture<T>& fixture)
      {
        bool ok = fixture.evaluate();
        for (size_t i = 0; i < N; ++i)
        {
          ok &= matches(static_cast<double>(fixture.model.getResidual().getData()[i]), 0.0, "initial residual");
        }
        return ok;
      }

    public:
      TestOutcome validation()
      {
        TestStatus success   = true;
        using Log            = Utilities::Logger;
        const auto verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::NONE);
        PhasorDynamics::Bus<double, IdxT> bus(1.0, 0.0);

        for (const auto parameter : Utilities::enum_values<Params>())
        {
          for (const auto value : {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
          {
            auto data                  = makeData();
            data.parameters[parameter] = value;
            PhasorDynamics::Perc1<double, IdxT> model(&bus, data);
            success *= model.verify() != 0;
          }
          auto data                  = makeData();
          data.parameters[parameter] = true;
          PhasorDynamics::Perc1<double, IdxT> model(&bus, data);
          success *= model.verify() != 0;
        }
        for (const auto parameter : {Params::Pnom, Params::Qnom})
        {
          auto data = makeData();
          data.parameters.erase(parameter);
          PhasorDynamics::Perc1<double, IdxT> model(&bus, data);
          success *= model.verify() != 0;
        }
        const std::vector<std::pair<Params, double>> invalid = {
            {Params::Pnom, -1}, {Params::Dbfl, 0.1}, {Params::Dbfh, -0.1}, {Params::Tvp, -1}, {Params::Tvq, -1}, {Params::Tap, -1}, {Params::Tbp, -1}, {Params::Taq, -1}, {Params::Tbq, -1}, {Params::Tt, -1}, {Params::Tv, -1}, {Params::Tf, -1}, {Params::Ipmin, 2}, {Params::Iqmin, 2}, {Params::Fcease, -0.1}, {Params::Fcease, 1.1}, {Params::Frecon, -1}, {Params::Vcease, -1}, {Params::Vrecon, 0.4}, {Params::Tcease, -1}, {Params::Tdelay, -1}, {Params::Trecon, -1}, {Params::Tramp, -1}};
        for (const auto& [parameter, value] : invalid)
        {
          auto data                  = makeData();
          data.parameters[parameter] = value;
          Fixture<double> fixture(data);
          success *= fixture.prepare();
          success *= fixture.model.verify() != 0;
          success *= fixture.model.initialize() != 0;
        }
        PhasorDynamics::Perc1<double, IdxT> busless(nullptr, makeData());
        success *= busless.verify() != 0;
        Fixture<double> unlinked(makeData());
        unlinked.model.getPorts().in.template port<Data::SignalInputs::freq>().connect(&unlinked.frequency_node);
        success *= unlinked.model.verify() != 0;
        for (const double base : {0.0, -1.0, std::numeric_limits<double>::infinity()})
        {
          Fixture<double> fixture(makeData());
          fixture.model.setSystemBase(60.0, base);
          success *= fixture.model.verify() != 0;
        }
        Log::setVerbosity(verbosity);
        return success.report(__func__);
      }

      TestOutcome initialization()
      {
        TestStatus success              = true;
        auto       data                 = makeData();
        data.parameters[Params::Kdroop] = 5.0;
        data.parameters[Params::Lfm]    = 0.0; // Documented 0.8 default.
        Fixture<double> fixture(data, 0.8, 0.6);
        fixture.attachFrequency(1.02);
        fixture.model.setSystemBase(50.0, 75.0e6);
        success *= fixture.initialize();
        success *= zeroResidual(fixture);
        success *= matches(fixture.frequency, 1.02, "known frequency preserved");
        success *= matches(fixture.bus.Vr(), 0.8, "known Vr preserved");
        success *= matches(fixture.bus.Vi(), 0.6, "known Vi preserved");
        success *= matches(fixture.y(Vars::IP), 0.8, "component active current");
        success *= matches(fixture.y(Vars::IQ), 0.2, "component reactive current");
        success *= matches(fixture.y(Vars::IQV), 0.528, "QPratio controlled path");
        success *= matches(fixture.bus.Ir(), -0.38, "real injection");
        success *= matches(fixture.bus.Ii(), -0.16, "imaginary injection");
        success *= fixture.model.tagDifferentiable() == 0;
        success *= fixture.model.setAbsoluteTolerance(1.0e-8) == 0;
        for (size_t i = 0; i < N; ++i)
        {
          success *= fixture.model.tag()[i] == (i < ND);
          success *= matches(fixture.model.absoluteTolerance().getData()[i], 1.0e-8, "absolute tolerance");
        }

        // Exercise all monitors through the same CSV controller used by systems.
        double                                   time = 0;
        Model::VariableMonitorController<double> monitor(time);
        monitor.addMonitor(fixture.model.getMonitor());
        std::ostringstream output;
        monitor.printFull(output, Model::VariableMonitorBase::Csv{});
        std::istringstream          values(output.str());
        std::string                 field;
        const std::array<double, 9> expected{0, -0.38, -0.16, std::hypot(0.38, 0.16), -0.4, -0.1, 0.8, 0.2, 1};
        for (const auto value : expected)
        {
          success *= static_cast<bool>(std::getline(values, field, ','));
          success *= matches(std::stod(field), value, "monitor");
        }

        // Reinitialization resets every history state and retains known inputs.
        fixture.y(Vars::XH)  = 1;
        fixture.y(Vars::TA)  = 0.3;
        success             *= fixture.model.initialize() == 0;
        success             *= matches(fixture.y(Vars::XH), 0, "history reset");
        success             *= zeroResidual(fixture);

        // At the cessation threshold the qualification clock has a nonzero
        // initial derivative even though no history has been accumulated.
        data.parameters[Params::Vcease] = 1.0;
        data.parameters[Params::Vrecon] = 1.1;
        Fixture<double> threshold(data);
        success *= threshold.initialize();
        success *= zeroResidual(threshold);
        success *= matches(threshold.model.yp().getData()[static_cast<size_t>(Vars::TC)], 0.5, "initial clock rate", 1.0e-9);

        auto integer_data                     = makeData();
        integer_data.parameters[Params::Pnom] = IdxT{1};
        integer_data.parameters[Params::Qnom] = IdxT{0};
        for (const auto parameter : {Params::Tv, Params::Tf, Params::Tt, Params::Tvp, Params::Tvq, Params::Tbp, Params::Tbq})
        {
          integer_data.parameters[parameter] = IdxT{0};
        }
        Fixture<double> floored(integer_data);
        success                *= floored.initialize();
        success                *= zeroResidual(floored);
        floored.y(Vars::XPLAG) -= 0.1;
        success                *= floored.evaluate();
        success                *= matches(floored.f(Vars::XPLAG), 100, "lag floor");
        return success.report(__func__);
      }

      TestOutcome initializationDomainAndFallback()
      {
        TestStatus success   = true;
        using Log            = Utilities::Logger;
        const auto verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::NONE);
        for (const auto& [parameter, value] : std::vector<std::pair<Params, double>>{
                 {Params::Ipmax, 0.8}, {Params::Iqmax, 0.5}, {Params::Iqmin, 0.6}, {Params::Lfm, 1.0e308}, {Params::Tcease, 1.0e-320}})
        {
          auto data                  = makeData();
          data.parameters[parameter] = value;
          Fixture<double> fixture(data);
          success *= fixture.prepare();
          success *= fixture.model.initialize() != 0;
        }
        Fixture<double> dead_bus(makeData(), 0.0);
        success *= dead_bus.prepare();
        success *= dead_bus.model.initialize() != 0;
        Fixture<double> bad_freq(makeData());
        bad_freq.attachFrequency(std::numeric_limits<double>::infinity());
        success *= bad_freq.prepare();
        success *= bad_freq.model.initialize() != 0;

        // A collapsed controlled-Q interval at zero is feasible. Extra Q
        // still reproduces the nonzero initial reactive consumption.
        auto data                        = makeData();
        data.parameters[Params::QPratio] = 0.0;
        data.parameters[Params::Iqmin]   = 0.0;
        data.parameters[Params::Iqmax]   = 0.0;
        // The inverse smooth clamp must recover targets close to either
        // active-current bound, where simply setting Pctrl = V * Ip fails.
        for (const double loading : {0.001, 0.999})
        {
          auto near_limit_data                    = data;
          near_limit_data.parameters[Params::Lfm] = loading;
          Fixture<double> near_limit(near_limit_data);
          success *= near_limit.initialize();
          success *= zeroResidual(near_limit);
        }
        Fixture<double> collapsed(data);
        success *= collapsed.initialize();
        success *= zeroResidual(collapsed);

        // Both a purely reactive load and a low-voltage active load use the
        // independent constant-admittance current/power relationships.
        for (const double p : {0.0, 0.4})
        {
          data                          = makeData();
          data.parameters[Params::Pnom] = p;
          Fixture<double> fixture(data, 0.3, 0.0);
          success             *= fixture.initialize();
          success             *= zeroResidual(fixture);
          fixture.bus.Vr()     = 0.6;
          fixture.y(Vars::V)   = std::sqrt(0.36 + 1.0e-12);
          fixture.y(Vars::IR) *= 2;
          fixture.y(Vars::II) *= 2;
          success             *= zeroResidual(fixture);
          success             *= matches(0.6 * fixture.y(Vars::IR), -4 * p, "fallback active power");
          success             *= matches(-0.6 * fixture.y(Vars::II), -0.4, "fallback reactive power");
        }
        data.parameters[Params::Pnom] = 0.0;
        data.parameters[Params::Qnom] = 0.0;
        Fixture<double> empty(data, 0.0);
        success *= empty.initialize();
        success *= zeroResidual(empty);
        Log::setVerbosity(verbosity);
        return success.report(__func__);
      }

      TestOutcome electricalResponse()
      {
        TestStatus success              = true;
        auto       data                 = makeData();
        data.parameters[Params::Fcease] = 0.0;
        data.parameters[Params::Kvp]    = 0.4;
        data.parameters[Params::Kvq]    = 0.2;
        data.parameters[Params::Tvp]    = 0.2;
        data.parameters[Params::Tvq]    = 0.4;
        data.parameters[Params::Tbp]    = 0.2;
        data.parameters[Params::Tbq]    = 0.4;
        data.parameters[Params::Tap]    = 0.1;
        data.parameters[Params::Taq]    = 0.1;
        data.parameters[Params::Kdroop] = 5.0;
        data.parameters[Params::Dbfl]   = -0.01;
        data.parameters[Params::Dbfh]   = 0.02;
        Fixture<double> fixture(data);
        fixture.attachFrequency(1.0);
        success          *= fixture.initialize();
        const double up0  = fixture.y(Vars::UP);
        const double uq0  = fixture.y(Vars::UQ);
        fixture.model.yp().setToZero();
        fixture.y(Vars::V)       = 0.8;
        fixture.y(Vars::XPWO)    = 0.9;
        fixture.y(Vars::XQWO)    = 0.9;
        fixture.y(Vars::VMEAS)   = 0.9;
        fixture.y(Vars::UP)      = up0 + 0.2;
        fixture.y(Vars::UQ)      = uq0 - 0.2;
        fixture.y(Vars::DFMEAS)  = 0.1;
        fixture.frequency        = 1.1;
        success                 *= fixture.evaluate();
        success                 *= matches(fixture.f(Vars::XPLAG), 1, "active lag");
        success                 *= matches(fixture.f(Vars::XQLAG), -0.5, "reactive lag");
        success                 *= matches(fixture.f(Vars::XPWO), -0.5, "active washout lag");
        success                 *= matches(fixture.f(Vars::XQWO), -0.25, "reactive washout lag");
        success                 *= matches(fixture.f(Vars::DP), -0.2, "active washout");
        success                 *= matches(fixture.f(Vars::DQ), -0.05, "reactive washout");
        success                 *= matches(fixture.f(Vars::PCTRL), 0.1, "active lead");
        success                 *= matches(fixture.f(Vars::QCTRL), -0.05, "reactive lead");
        success                 *= matches(fixture.f(Vars::VMEAS), -5, "voltage sensor");
        fixture.y(Vars::V)       = 0.6;
        success                 *= fixture.evaluate();
        success                 *= matches(fixture.f(Vars::DFMEAS), -5, "low-voltage frequency gate", 1.0e-8);
        for (size_t i = static_cast<size_t>(Vars::XD); i < ND; ++i)
        {
          success *= matches(fixture.model.getResidual().getData()[i], 0, "disabled history");
        }

        Fixture<double> fraction(makeData());
        success                  *= fraction.initialize();
        fraction.y(Vars::FRACON)  = 0.4;
        success                  *= fraction.evaluate();
        success                  *= matches(fraction.f(Vars::IP), -24, "active-current fraction and lag");
        success                  *= matches(fraction.f(Vars::IQ), -6, "extra reactive current shares fraction and lag");

        // Away from limits, voltage exponents give constant power, current,
        // and impedance. Bounds are wide enough to make smoothing negligible.
        for (double exponent : {0.0, 1.0, 2.0})
        {
          data                           = makeData();
          data.parameters[Params::Ipmax] = 3.0;
          data.parameters[Params::Iqmax] = 3.0;
          data.parameters[Params::nP]    = exponent;
          data.parameters[Params::nQ]    = exponent;
          Fixture<double> voltage(data);
          success                *= voltage.initialize();
          voltage.y(Vars::VSAFE) *= 0.8;
          voltage.y(Vars::IPV)    = 0;
          voltage.y(Vars::IQV)    = 0;
          success                *= voltage.evaluate();
          success                *= matches(voltage.f(Vars::IPV), 0.8 * std::pow(0.8, exponent - 1), "active voltage exponent");
          success                *= matches(voltage.f(Vars::IQV), 0.528 * std::pow(0.8, exponent - 1), "reactive voltage exponent");
          voltage.y(Vars::PCTRL)  = 10;
          voltage.y(Vars::QCTRL)  = -10;
          success                *= voltage.evaluate();
          success                *= matches(voltage.f(Vars::IPV), 3, "active current limit");
          success                *= matches(voltage.f(Vars::IQV), -0.66, "reactive current limit");
        }
        return success.report(__func__);
      }

      TestOutcome historyRates()
      {
        TestStatus success = true;
        for (const bool immediate : {false, true})
        {
          auto data                       = makeData();
          data.parameters[Params::Tdelay] = immediate ? 0.0 : 0.1;
          data.parameters[Params::Tcease] = immediate ? 0.0 : 0.1;
          data.parameters[Params::Trecon] = immediate ? 0.0 : 0.1;
          data.parameters[Params::Tramp]  = immediate ? 0.0 : 1.0;
          data.parameters[Params::Fcease] = 0.6;
          data.parameters[Params::Frecon] = 1.5;
          Fixture<double> fixture(data);
          success *= fixture.initialize();
          fixture.model.yp().setToZero();
          fixture.y(Vars::VMEAS)  = 0.2;
          fixture.y(Vars::TC)     = 0.2;
          success                *= fixture.evaluate();
          success                *= matches(fixture.f(Vars::XD), 1 / TAU, "cessation commit");
          success                *= matches(fixture.f(Vars::TC), immediate ? -0.2 / TAU : 1, "qualification timer");

          // A committed delay keeps running after voltage recovers.
          fixture.y(Vars::XD)     = 1;
          fixture.y(Vars::TD)     = 0.2;
          fixture.y(Vars::VMEAS)  = 1;
          success                *= fixture.evaluate();
          success                *= matches(fixture.f(Vars::XD), 0, "delay memory holds");
          success                *= matches(fixture.f(Vars::XC), 1 / TAU, "delay completes");
          success                *= matches(fixture.f(Vars::TD), immediate ? -0.2 / TAU : 1, "committed clock");

          // Reset has priority over set at each subsequent stage.
          fixture.y(Vars::XC)  = 1;
          fixture.y(Vars::XH)  = 1;
          fixture.y(Vars::XD)  = 0;
          fixture.y(Vars::TR)  = 0.2;
          success             *= fixture.evaluate();
          success             *= matches(fixture.f(Vars::XR), 1 / TAU, "reconnection release");
          success             *= matches(fixture.f(Vars::FRACON) + fixture.y(Vars::FRACON), 0.4, "ceased fraction");
          fixture.y(Vars::XD)  = 1;
          fixture.y(Vars::XR)  = 1;
          success             *= fixture.evaluate();
          success             *= matches(fixture.f(Vars::XD), -1 / TAU, "cessation resets delay");
          success             *= matches(fixture.f(Vars::XC), -1 / TAU, "release resets cessation");
          success             *= matches(fixture.f(Vars::XR), -1 / TAU, "new delay resets release");

          fixture.y(Vars::XD)  = 0;
          fixture.y(Vars::XC)  = 0;
          fixture.y(Vars::TA)  = 0.5;
          success             *= fixture.evaluate();
          success             *= matches(fixture.f(Vars::FRACON) + fixture.y(Vars::FRACON), immediate ? 1.3 : 0.85, "reconnection fraction");
          success             *= matches(fixture.f(Vars::TA), immediate ? -0.5 / TAU : 1, "ramp rate");
          fixture.y(Vars::XD)  = 1;
          success             *= fixture.evaluate();
          success             *= matches(fixture.f(Vars::TA), immediate ? -0.5 / TAU : 0, "ramp held during delay");
        }
        return success.report(__func__);
      }

      TestOutcome jsonSystem()
      {
        TestStatus         success = true;
        std::istringstream json(R"({
          "header": {"format_version": 0.2, "format_revision": 0,
            "case_name": "PERC1", "case_description": "load test", "case_comments": ""},
          "params": {"freq_base": 60, "va_base": 75000000},
          "buses": [{"number": 1, "class": "BusInfinite", "name": "terminal",
            "init": {"Vr": 1, "Vi": 0}, "params": {"kv": 138}}],
          "signals": [{"signal_id": 1, "name": "frequency"}],
          "devices": [
            {"class": "Perc1", "id": "load", "ports": {"bus": 1, "freq": 1},
             "params": {"Pnom": 0.4, "Qnom": 0.1}, "mon": ["ir", "ii", "im", "p", "q", "ip", "iq", "fracon"]},
            {"class": "ConstantSignalSource", "id": "frequency", "ports": {"sr": 1}, "params": {"Sr": 1.0}}
          ]})");
        auto               data  = PhasorDynamics::parseSystemModelData(json);
        success                 *= data.perc1.size() == 1;
        PhasorDynamics::SystemModel<double, IdxT> system(data);
        success *= system.allocate() == 0;
        success *= system.initialize() == 0;
        success *= system.evaluateResidual() == 0;
        for (size_t i = 0; i < system.size(); ++i)
        {
          success *= matches(system.getResidual().getData()[i], 0, "system initialization");
        }
        return success.report(__func__);
      }

      TestOutcome jacobian()
      {
        TestStatus success = true;
        // Compare every internal and bus row, every state/signal column, and
        // the full F_y + alpha F_yp against an independent central difference.
        for (const bool static_load : {false, true})
        {
          for (const bool immediate : {false, true})
          {
            for (const double alpha : {0.0, 3.0})
            {
              auto data                       = makeData();
              data.parameters[Params::Kdroop] = 4.0;
              data.parameters[Params::nP]     = 0.7;
              data.parameters[Params::nQ]     = 1.3;
              data.parameters[Params::Tcease] = immediate ? 0.0 : 0.1;
              data.parameters[Params::Tdelay] = immediate ? 0.0 : 0.2;
              data.parameters[Params::Trecon] = immediate ? 0.0 : 0.15;
              data.parameters[Params::Tramp]  = immediate ? 0.0 : 0.5;
              const double    vr              = static_load ? 0.3 : 1.0;
              Fixture<double> plain(data, vr);
              Fixture<DepVar> tracked(data, vr);
              plain.attachFrequency(1.03);
              tracked.attachFrequency(1.03);
              success        *= plain.initialize();
              success        *= tracked.initialize();
              auto set_point  = [](auto& fixture)
              {
                fixture.y(Vars::XD)     = 0.498;
                fixture.y(Vars::XC)     = 0.501;
                fixture.y(Vars::XR)     = 0.502;
                fixture.y(Vars::XH)     = 0.499;
                fixture.y(Vars::TC)     = 0.101;
                fixture.y(Vars::TD)     = 0.201;
                fixture.y(Vars::TR)     = 0.151;
                fixture.y(Vars::TA)     = TAU;
                fixture.y(Vars::VMEAS)  = 0.6;
                fixture.y(Vars::VSAFE)  = 0.6;
                fixture.y(Vars::DFMEAS) = 0.02;
                for (size_t i = 0; i < N; ++i)
                {
                  fixture.model.yp().getData()[i] = 0.01;
                }
              };
              set_point(plain);
              set_point(tracked);
              for (size_t i = 0; i < N; ++i)
              {
                tracked.model.y().getData()[i].setVariableNumber(2 * i);
                tracked.model.yp().getData()[i].setVariableNumber(2 * i + 1);
              }
              for (size_t i = 0; i < 2; ++i)
              {
                tracked.bus.y().getData()[i].setVariableNumber(2 * (N + i));
                plain.bus.setVariableIndex(i, N + i);
                plain.bus.setResidualIndex(i, N + i);
              }
              tracked.frequency.setVariableNumber(2 * (N + 2));
              tracked.model.updateTime(0, alpha);
              plain.model.updateTime(0, alpha);
              success      *= tracked.evaluate();
              success      *= plain.evaluate();
              using Matrix  = std::array<std::array<double, N + 3>, N + 2>;
              Matrix automatic{};
              for (size_t row = 0; row < N + 2; ++row)
              {
                const auto& f = row < N ? tracked.model.getResidual().getData()[row] : tracked.bus.getResidual().getData()[row - N];
                for (const auto& [column, value] : f.getDependencies())
                {
                  automatic[row][column / 2] += value * (column % 2 == 0 ? 1 : alpha);
                }
              }
#ifdef GRIDKIT_ENABLE_ENZYME
              success         *= plain.model.evaluateJacobian() == 0;
              const auto* coo  = plain.model.getCooJacobian();
              Matrix      enzyme{};
              for (size_t i = 0; i < coo->getNnz(); ++i)
              {
                const auto row  = coo->getRowData()[i];
                const auto col  = coo->getColData()[i];
                success        *= row < N + 2 && col < N + 3;
                if (row < N + 2 && col < N + 3)
                {
                  enzyme[row][col] += coo->getValues()[i];
                }
              }
              for (size_t row = 0; row < N + 2; ++row)
              {
                for (size_t col = 0; col < N + 3; ++col)
                {
                  success *= matches(enzyme[row][col], automatic[row][col], "Enzyme/DependencyTracking Jacobian", 1.0e-9);
                }
              }
#endif
              auto residual = [&]()
              {
                plain.evaluate();
                std::array<double, N + 2> result{};
                std::copy_n(plain.model.getResidual().getData(), N, result.data());
                std::copy_n(plain.bus.getResidual().getData(), 2, result.data() + N);
                return result;
              };
              for (size_t col = 0; col < N + 3; ++col)
              {
                double&      value      = col < N       ? plain.model.y().getData()[col]
                                          : col < N + 2 ? plain.bus.y().getData()[col - N]
                                                        : plain.frequency;
                const double initial    = value;
                const double step       = col == static_cast<size_t>(Vars::TA) ? 1.0e-10 : 1.0e-7;
                const double derivative = col < N ? plain.model.yp().getData()[col] : 0;
                value                   = initial + step;
                if (col < N)
                  plain.model.yp().getData()[col] = derivative + alpha * step;
                const auto plus = residual();
                value           = initial - step;
                if (col < N)
                  plain.model.yp().getData()[col] = derivative - alpha * step;
                const auto minus = residual();
                value            = initial;
                if (col < N)
                  plain.model.yp().getData()[col] = derivative;
                for (size_t row = 0; row < N + 2; ++row)
                {
                  // Stiff latch residuals can be O(1/tau) even when a
                  // particular partial is tiny. Bound cancellation error
                  // from the residual magnitudes and perturbation size.
                  const double roundoff = 8 * std::numeric_limits<double>::epsilon()
                                          * std::max(std::abs(plus[row]), std::abs(minus[row])) / step;
                  const double tolerance  = 1.0e-5 + roundoff / std::max(1.0, std::abs(automatic[row][col]));
                  success                *= matches((plus[row] - minus[row]) / (2 * step), automatic[row][col], "finite-difference Jacobian", tolerance);
                }
              }
            }
          }
        }
        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
