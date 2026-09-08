/**
 * @file ModulationTests.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Radial voltage limit, DC normalization, and signal gradient checks.
 */
#pragma once

#include <array>
#include <cmath>
#include <limits>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/EMT/Operators/Modulation/Modulation.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    class ModulationTests
    {
      using Model   = EMT::Modulation<double, size_t>;
      using Signal  = EMT::Signal<double, size_t>;
      using Outputs = EMT::ModulationOutputs;

      static constexpr double limit_factor = 0.95;

      template <typename Action>
      bool throws(Action action)
      {
        try
        {
          action();
        }
        catch (const std::exception&)
        {
          return true;
        }
        return false;
      }

      static Model::ModelDataT data(double limit = limit_factor)
      {
        Model::ModelDataT data;
        data.parameters[Model::ModelDataT::Parameters::Mmax] = limit;
        return data;
      }

      struct Fixture
      {
        std::array<double, 3> values;
        std::array<size_t, 3> indices{3, 5, 8};
        std::array<Signal, 3> signals;
        Model                 model;

        Fixture(std::array<double, 3> inputs, double limit = limit_factor)
          : values(inputs), model(data(limit))
        {
          for (size_t n = 0; n < 3; ++n)
            signals[n].set(&values[n], &indices[n]);
          model.attachInput({&signals[0], &signals[1]}, &signals[2]);
          model.allocate();
        }

        double read(Outputs output)
        {
          return model.outputSignal(output).read();
        }

        /// Sum gradient coefficients on one input index.
        double derivative(Outputs output, size_t n)
        {
          Signal::GradientT gradient;
          model.outputSignal(output).appendGradient(gradient);
          double total = 0;
          for (const auto& [index, coefficient] : gradient)
            if (index == indices[n])
              total += coefficient;
          return total;
        }
      };

      /// Limited voltage magnitude available from sinusoidal PWM.
      static double available(double vdc, double limit = limit_factor)
      {
        return std::sqrt(3.0 / 8) * limit * vdc;
      }

    public:
      TestOutcome limits()
      {
        TestStatus success = true;
        // Inside the limit the command is normalized and returned unchanged.
        Fixture    linear({150, -40, 400});
        success *= linear.model.verify() == 0 && linear.model.initialize() == 0;
        success *= linear.model.size() == 0 && linear.model.nnz() == 0;
        success *= std::abs(linear.read(Outputs::md) - 0.75) < 1e-12;
        success *= std::abs(linear.read(Outputs::mq) + 0.2) < 1e-12;
        success *= std::abs(linear.read(Outputs::ulimd) - 150) < 1e-9;
        success *= std::abs(linear.read(Outputs::ulimq) + 40) < 1e-9;
        // Beyond the limit the direction is preserved and the magnitude clipped.
        Fixture      clipped({300, 400, 100});
        const double ulimd  = clipped.read(Outputs::ulimd);
        const double ulimq  = clipped.read(Outputs::ulimq);
        success            *= std::abs(std::hypot(ulimd, ulimq) - available(100)) < 1e-9;
        success            *= std::abs(ulimq / ulimd - 4.0 / 3) < 1e-12;
        success            *= std::abs(std::hypot(clipped.read(Outputs::md), clipped.read(Outputs::mq)) - std::sqrt(1.5) * limit_factor) < 1e-12;
        // At zero DC voltage the limited command vanishes and the modulation sits on the boundary.
        Fixture collapsed({30, -40, 0});
        success *= collapsed.read(Outputs::ulimd) == 0 && collapsed.read(Outputs::ulimq) == 0;
        success *= std::abs(std::hypot(collapsed.read(Outputs::md), collapsed.read(Outputs::mq)) - std::sqrt(1.5) * limit_factor) < 1e-12;
        // The default limit is the full sinusoidal range.
        Fixture full({600, 0, 400}, 1.0);
        success *= std::abs(full.read(Outputs::ulimd) - available(400, 1.0)) < 1e-9;
        success *= std::abs(full.read(Outputs::md) - std::sqrt(1.5)) < 1e-12;
        // Prescribed outputs must match the computed values.
        success *= linear.model.initialize({{Outputs::md, 0.75}, {Outputs::ulimq, -40}}) == 0;
        success *= throws([&]
                          { linear.model.initialize({{Outputs::md, 0.7}}); });
        // Invalid limits and DC voltages are rejected.
        Fixture negative({150, -40, -1});
        success *= throws([&]
                          { negative.read(Outputs::md); });
        success *= throws([&]
                          { Model nonfinite(data(std::numeric_limits<double>::infinity())); });
        for (const double limit : {1.1, 0.0, -0.5})
        {
          Fixture invalid({150, -40, 400}, limit);
          success *= invalid.model.verify() != 0;
          success *= throws([&]
                            { invalid.read(Outputs::md); });
        }
        Model unattached;
        success *= unattached.verify() != 0;
        return success.report(__func__);
      }

      TestOutcome gradients()
      {
        TestStatus success       = true;
        double     maximum_error = 0;
        for (const auto& inputs : {std::array<double, 3>{150, -40, 400},
                                   std::array<double, 3>{300, 400, 100},
                                   std::array<double, 3>{200, 150, 250},
                                   std::array<double, 3>{-90, 20, 150}})
        {
          Fixture f(inputs);
          for (const auto output : {Outputs::md, Outputs::mq, Outputs::ulimd, Outputs::ulimq})
            for (size_t n = 0; n < 3; ++n)
            {
              const double original  = f.values[n];
              const double h         = 1e-6 * (1 + std::abs(original));
              f.values[n]            = original + h;
              const double plus      = f.read(output);
              f.values[n]            = original - h;
              const double minus     = f.read(output);
              f.values[n]            = original;
              const double fd        = (plus - minus) / (2 * h);
              const double error     = std::abs(f.derivative(output, n) - fd) / (1 + std::abs(fd));
              maximum_error          = std::max(maximum_error, error);
              success               *= error < 1e-7;
            }
        }
        std::cout << "Modulation gradient maximum scaled error: " << maximum_error << "\n";
        // At zero DC voltage the limited command grows linearly with the DC voltage.
        Fixture collapsed({30, -40, 0});
        success *= std::abs(collapsed.derivative(Outputs::ulimd, 2) - collapsed.read(Outputs::md) / 2) < 1e-12;
        success *= std::abs(collapsed.derivative(Outputs::ulimq, 2) - collapsed.read(Outputs::mq) / 2) < 1e-12;
        success *= std::abs(collapsed.derivative(Outputs::md, 2)) < 1e-12;
        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
