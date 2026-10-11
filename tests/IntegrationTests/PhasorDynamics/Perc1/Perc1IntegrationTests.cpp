#include <cmath>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <GridKit/Model/PhasorDynamics/Load/PERC1/Perc1.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace
{
  using namespace GridKit;
  using namespace GridKit::PhasorDynamics;
  using Params = Perc1Parameters;
  using Vars   = Perc1InternalVariables;
  using Data   = Perc1Data<double, size_t>;

  Data makeData(double frecon = 0.5)
  {
    Data data;
    data.parameters[Params::Pnom]   = 0.4;
    data.parameters[Params::Qnom]   = 0.1;
    data.parameters[Params::Fcease] = 0.6;
    data.parameters[Params::Frecon] = frecon;
    data.parameters[Params::Tv]     = 0.001;
    data.parameters[Params::Tt]     = 0.01;
    data.parameters[Params::Tcease] = 0.02;
    data.parameters[Params::Tdelay] = 0.04;
    data.parameters[Params::Trecon] = 0.03;
    data.parameters[Params::Tramp]  = 0.2;
    return data;
  }

  struct Simulation
  {
    BusInfinite<double, size_t>                    bus{1.0, 0.0};
    Perc1<double, size_t>                          load;
    SystemModel<double, size_t>                    system;
    AnalysisManager::Sundials::Ida<double, size_t> ida;
    double                                         time{0};

    explicit Simulation(const Data& data)
      : load(&bus, data), ida(&system)
    {
      system.addBus(&bus);
      system.addComponent(&load);
    }

    bool initialize()
    {
      ida.setTolerance(1.0e-9, 1.0e-11);
      ida.setMaxSteps(100000);
      return system.allocate() == 0 && system.initialize() == 0 && ida.configureSimulation() == 0;
    }

    bool advance(double duration, double voltage)
    {
      bus.Vr() = voltage;
      if (ida.initializeSimulation(time) != 0 || ida.runSimulation(time + duration, 0.005) != 0)
      {
        return false;
      }
      time              += duration;
      const auto* state  = load.y().getData();
      for (size_t i = 0; i < load.size(); ++i)
      {
        if (!std::isfinite(state[i]))
          return false;
      }
      return true;
    }

    double y(Vars variable) const
    {
      return load.y().getData()[static_cast<size_t>(variable)];
    }
  };

  bool near(double actual, double expected, const char* label, double tolerance = 0.005)
  {
    if (std::isfinite(actual) && std::abs(actual - expected) <= tolerance)
      return true;
    std::cout << label << ": " << actual << " != " << expected << '\n';
    return false;
  }

  Testing::TestOutcome repeatedFaults()
  {
    Testing::TestStatus success = true;
    Simulation          sim(makeData());
    if (!sim.initialize())
      return Testing::TestStatus(false).report(__func__);
    success *= sim.advance(0.05, 1.0);
    success *= near(sim.y(Vars::FRACON), 1, "initial fraction", 1.0e-9);
    success *= sim.advance(0.008, 0.2);
    success *= sim.advance(0.03, 1.0);
    success *= near(sim.y(Vars::XH), 0, "brief fault must not cease");
    success *= near(sim.y(Vars::TC), 0, "interrupted qualification resets", 1.0e-7);

    success *= sim.advance(0.03, 0.2);
    success *= near(sim.y(Vars::XD), 1, "committed delay");
    success *= near(sim.y(Vars::XC), 0, "delay before cessation");
    success *= sim.advance(0.02, 1.0);
    success *= near(sim.y(Vars::XD), 1, "recovery retains committed delay");
    success *= near(sim.y(Vars::XC), 0, "pending delay");
    success *= sim.advance(0.025, 1.0);
    success *= near(sim.y(Vars::XC), 1, "cessation after recovery");
    success *= near(sim.y(Vars::FRACON), 0.4, "ceased fraction");
    success *= sim.advance(0.05, 1.0);
    success *= near(sim.y(Vars::XR), 1, "qualified recovery releases");
    success *= sim.y(Vars::FRACON) > 0.4 && sim.y(Vars::FRACON) < 0.7;

    success                 *= sim.advance(0.03, 0.2);
    success                 *= near(sim.y(Vars::XD), 1, "second fault commits");
    const double ramp_clock  = sim.y(Vars::TA);
    success                 *= sim.advance(0.01, 0.2);
    success                 *= near(sim.y(Vars::TA), ramp_clock, "ramp held during committed delay", 1.0e-7);
    success                 *= sim.advance(0.06, 0.2);
    success                 *= near(sim.y(Vars::XC), 1, "second cessation");
    success                 *= near(sim.y(Vars::TA), 0, "cessation resets ramp", 1.0e-7);
    success                 *= sim.advance(0.015, 1.0);
    success                 *= sim.advance(0.01, 0.2);
    success                 *= near(sim.y(Vars::TR), 0, "interrupted recovery resets", 1.0e-7);
    success                 *= sim.advance(0.6, 1.0);
    success                 *= near(sim.y(Vars::FRACON), 0.7, "reconnected fraction");
    success                 *= near(sim.y(Vars::IR), -0.28, "recovered active consumption");
    success                 *= near(sim.y(Vars::II), 0.07, "recovered reactive consumption");
    return success.report(__func__);
  }

  Testing::TestOutcome immediateAndPartialRecovery()
  {
    Testing::TestStatus success = true;
    for (const double fraction : {0.0, 0.5, 1.5})
    {
      auto data = makeData(fraction);
      for (const auto parameter : {Params::Tcease, Params::Tdelay, Params::Trecon, Params::Tramp})
      {
        data.parameters[parameter] = 0.0;
      }
      Simulation sim(data);
      if (!sim.initialize())
        return Testing::TestStatus(false).report(__func__);
      success             *= sim.advance(0.02, 0.0);
      success             *= near(sim.y(Vars::FRACON), 0.4, "immediate cessation");
      success             *= sim.advance(0.2, 1.0);
      const double target  = 0.4 + 0.6 * fraction;
      success             *= near(sim.y(Vars::FRACON), target, "immediate recovery target");
      success             *= near(sim.y(Vars::IR), -0.4 * target, "recovery consumption");
      success             *= near(sim.y(Vars::TC), 0, "zero cessation timer", 1.0e-10);
      success             *= near(sim.y(Vars::TD), 0, "zero delay timer", 1.0e-10);
      success             *= near(sim.y(Vars::TR), 0, "zero recovery timer", 1.0e-10);
      success             *= near(sim.y(Vars::TA), 0, "zero ramp timer", 1.0e-10);
    }
    return success.report(__func__);
  }
} // namespace

int main()
{
  GridKit::Testing::TestingResults result;
  result += repeatedFaults();
  result += immediateAndPartialRecovery();
  return result.summary();
}
