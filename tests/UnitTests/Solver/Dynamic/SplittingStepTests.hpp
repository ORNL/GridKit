#pragma once

#include <cmath>
#include <deque>
#include <sstream>
#include <vector>

#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Solver/Dynamic/SplittingStep.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit::Testing
{
  class SplittingStepTests
  {
    using Ida   = AnalysisManager::Sundials::Ida<double, size_t>;
    using Split = AnalysisManager::Sundials::SplittingStep<double, size_t>;
    using Log   = Utilities::Logger;

    // x' = sign * u(t), z = u(t).
    class Region : public PhasorDynamics::Component<double, size_t>
    {
    public:
      explicit Region(double sign)
        : sign(sign)
      {
        size_ = 2;
      }

      int allocate() override
      {
        allocateVectors(size_);
        return 0;
      }

      int initialize() override
      {
        allocate();
        y_.setToConst(1.0);
        yp_.setToConst(0.0);
        yp_.getData()[0] = sign;
        yp_.setDataUpdated();
        return 0;
      }

      int tagDifferentiable() override
      {
        tag_ = {true, false};
        return 0;
      }

      int setAbsoluteTolerance(double tol) override
      {
        abs_tol_.setToConst(tol);
        return 0;
      }

      bool hasJacobian() override
      {
        return false;
      }

      int evaluateJacobian() override
      {
        return 0;
      }

      int verify() const override
      {
        return 0;
      }

      int setGridKitComponentID(size_t) override
      {
        return 0;
      }

      int evaluateResidual() override
      {
        const auto u    = input.at(time());
        f_.getData()[0] = yp_.getData()[0] - sign * u;
        f_.getData()[1] = y_.getData()[1] - u;
        f_.setDataUpdated();
        return 0;
      }

      Model::Input<double> input{1.0};
      double               sign;
    };

    /// Regions of alternating sign; the first reads the second, every other
    /// region reads the one before it. Three regions color as {1}, then {0, 2}.
    struct Chain
    {
      std::deque<Region> regions;
      std::deque<Ida>    solvers;
      Split              split;

      explicit Chain(std::size_t size, int threads = 1, double tolerance = 0.0)
      {
        for (std::size_t i = 0; i < size; ++i)
        {
          auto& solver = solvers.emplace_back(&regions.emplace_back(1.0 - 2.0 * static_cast<double>(i % 2)));
          solver.setTolerance(1e-10, 1e-12);
          solver.setMaxSteps(10000);
          solver.configureSimulation();
        }
        for (std::size_t i = 0; i < size; ++i)
        {
          std::size_t source = 1;
          if (i > 0)
          {
            source = i - 1;
          }
          split.addPartition(solvers[i], {{&regions[source], 0, &regions[i].input}}, solvers[i].context());
        }
        split.setNumThreads(threads);
        split.setFixedStep(0.125);
        split.setCouplingTolerance(tolerance);
        split.setTolerance(1e-9, 1e-11);
        split.configureSimulation();
        split.initializeSimulation(0.0);
      }

      std::vector<double> state()
      {
        std::vector<double> x;
        for (auto& region : regions)
        {
          x.push_back(region.y().getData()[0]);
        }
        return x;
      }
    };

  public:
    TestOutcome sequentialRamps()
    {
      TestStatus success = true;
      Chain      chain(2);
      chain.split.runSimulation(0.125);
      const auto y  = chain.state();
      success      *= isEqual(y[0], 1.125, 1e-7);
      success      *= isEqual(y[1], 1.0 - 0.125 - 0.5 * 0.125 * 0.125, 1e-7);
      return success.report(__func__);
    }

    TestOutcome event()
    {
      TestStatus       success = true;
      constexpr double h       = 0.125;
      Chain            chain(2);
      chain.split.runSimulation(2 * h);
      const auto start      = chain.state();
      chain.regions[0].sign = -1.0;
      chain.split.initializeSimulation(2 * h);
      success *= chain.regions[0].input.rate == 0.0 && chain.regions[1].input.rate == 0.0;
      chain.split.runSimulation(3 * h);
      const auto y  = chain.state();
      success      *= isEqual(y[0], start[0] - h * start[1], 1e-7);
      success      *= isEqual(y[1], start[1] - h * 0.5 * (start[0] + y[0]), 1e-7);
      return success.report(__func__);
    }

    /// Partitions of one color advance together; the thread count cannot change the result.
    TestOutcome rejectionAndThreads()
    {
      TestStatus success = true;
      Chain      serial(3, 1, 1e-4);
      serial.split.runSimulation(0.25);
      serial.regions[0].sign = -1.0;
      serial.split.initializeSimulation(0.25);
      serial.split.runSimulation(0.5);
      success *= serial.split.numRejectedSteps() > 0;
#ifdef GRIDKIT_ENABLE_OPENMP
      Chain parallel(3, 2, 1e-4);
      parallel.split.runSimulation(0.25);
      parallel.regions[0].sign = -1.0;
      parallel.split.initializeSimulation(0.25);
      parallel.split.runSimulation(0.5);
      success *= parallel.state() == serial.state();
      success *= parallel.split.numSteps() == serial.split.numSteps();
      success *= parallel.split.numRejectedSteps() == serial.split.numRejectedSteps();
#endif
      return success.report(__func__);
    }

    TestOutcome regionalFailure()
    {
      class FailingIda : public Ida
      {
      public:
        using Ida::Ida;

        SUNStepper createSUNStepper() override
        {
          auto stepper = Ida::createSUNStepper();
          SUNStepper_SetEvolveFn(stepper, [](SUNStepper, sunrealtype, N_Vector, sunrealtype*) -> int
                                 { return SUN_ERR_OP_FAIL; });
          return stepper;
        }
      };

      TestStatus success = true;
      Region     first(1.0), second(-1.0);
      FailingIda a(&first);
      Ida        b(&second);
      a.configureSimulation();
      b.configureSimulation();
      Split split;
      split.addPartition(a, {{&second, 0, &first.input}}, a.context());
      split.addPartition(b, {{&first, 0, &second.input}}, b.context());
      split.setFixedStep(0.125);
      split.setTolerance(1e-8, 1e-10);
      split.configureSimulation();
      split.initializeSimulation(0.0);
      int outputs = 0;
      split.setOutput([&](double)
                      { ++outputs; });
      success *= throws<std::runtime_error>([&]
                                            { split.runSimulation(0.125); });
      success *= outputs == 0 && split.numSteps() == 0;
      success *= second.y().getData()[0] == 1.0; // its color never started
      return success.report(__func__);
    }

    TestOutcome configuration()
    {
      TestStatus success = true;
      Split      split;
      success *= throws<std::invalid_argument>([&]
                                               { split.setNumThreads(0); });
#ifndef GRIDKIT_ENABLE_OPENMP
      success *= throws<std::invalid_argument>([]
                                               { Chain chain(2, 2); });
#endif
      Chain configured(2);
      success *= throws<std::logic_error>([&]
                                          { configured.split.setNumThreads(1); });

      Region first(1.0), second(-1.0);
      Ida    a(&first), b(&second);
      a.configureSimulation();
      b.configureSimulation();
      Split uncoupled;
      uncoupled.addPartition(a, {}, a.context());
      uncoupled.addPartition(b, {}, b.context());
      uncoupled.setFixedStep(0.125);
      success *= throws<std::invalid_argument>([&]
                                               { uncoupled.configureSimulation(); });
      return success.report(__func__);
    }

    TestOutcome diagnostics()
    {
      TestStatus         success = true;
      std::ostringstream output;
      const auto         verbosity = Log::verbosity();
      Log::setOutput(output);
      Log::setVerbosity(Log::WARNINGS);
#ifdef GRIDKIT_ENABLE_OPENMP
#pragma omp parallel for num_threads(2)
#endif
      for (int i = 0; i < 8; ++i)
      {
        Log::ScopedOutput capture;
        Log::warning() << "region " << i << " complete\n";
        Log::misc() << "hidden";
      }
      Log::setOutput(std::cout);
      Log::setVerbosity(verbosity);
      for (int i = 0; i < 8; ++i)
        success *= output.str().find("region " + std::to_string(i) + " complete\n") != std::string::npos;
      success *= output.str().find("hidden") == std::string::npos;
      return success.report(__func__);
    }
  };
} // namespace GridKit::Testing
