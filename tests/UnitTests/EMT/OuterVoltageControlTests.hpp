/**
 * @file OuterVoltageControlTests.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Output initialization and finite-difference Jacobian checks.
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/EMT/Component/Controller/OuterVoltageControl/OuterVoltageControl.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <typename ScalarT, typename IdxT>
    class OuterVoltageControlTests
    {
      using Outer = EMT::Controller::OuterVoltageControl<ScalarT, IdxT>;

      template <typename Action>
      bool rejects(Action action)
      {
        try
        {
          action();
        }
        catch (const std::invalid_argument&)
        {
          return true;
        }
        return false;
      }

      template <typename Model>
      struct Fixture
      {
        static constexpr size_t                    count = std::tuple_size_v<typename Model::InputSignals>;
        std::array<double, count>                  values;
        std::array<size_t, count>                  indices;
        std::array<typename Model::SignalT, count> signals;
        typename Model::SignalT                    alias;
        Model                                      model;

        Fixture(const typename Model::ModelDataT& data, std::array<double, count> inputs)
          : values(inputs), model(data)
        {
          typename Model::InputSignals attached;
          for (size_t n = 0; n < count; ++n)
          {
            indices[n] = n;
            signals[n].set(&values[n], &indices[n]);
            attached[n] = &signals[n];
          }
          model.attachInput(attached);
          model.assignOutput(static_cast<typename Model::Outputs>(0), &alias);
          model.allocate();
          model.assignGlobalIndices(count);
        }

        Fixture(const Fixture&)            = delete;
        Fixture& operator=(const Fixture&) = delete;

        bool algebraicResidualsZero()
        {
          model.evaluateResidual();
          for (size_t n = 2; n < model.size(); ++n)
            if (std::abs(model.getResidual().getData()[n]) > 1e-10)
              return false;
          return true;
        }

        bool jacobian()
        {
          bool   success       = true;
          double maximum_error = 0;
          for (const auto& [ys, yps] : {std::pair{1.0, 0.0}, std::pair{0.0, 1.0}, std::pair{2.0, 3.0}, std::pair{0.0, 0.0}})
          {
            std::map<std::pair<size_t, size_t>, double> entries;
            for (const auto& entry : model.jacobianEntries(ys, yps))
              entries[{entry.row, entry.column}] += entry.value;
            for (size_t j = 0; j < count + model.size(); ++j)
            {
              double&               y        = j < count ? values[j] : model.y().getData()[j - count];
              const double          original = y;
              const double          h        = (ys == 0 ? 1e-3 : 1e-6) * (1 + std::abs(original));
              std::array<double, 6> plus, minus;
              for (const double sign : {1.0, -1.0})
              {
                y = original + sign * h * ys;
                if (j >= count)
                  model.yp().getData()[j - count] = sign * h * yps;
                model.evaluateResidual();
                for (size_t n = 0; n < model.size(); ++n)
                  (sign > 0 ? plus : minus)[n] = model.getResidual().getData()[n];
              }
              y = original;
              if (j >= count)
                model.yp().getData()[j - count] = 0;
              for (size_t n = 0; n < model.size(); ++n)
              {
                const double fd     = (plus[n] - minus[n]) / (2 * h);
                const double error  = std::abs(entries[{count + n, j}] - fd) / (1 + std::abs(fd));
                maximum_error       = std::max(maximum_error, error);
                success            &= error < 1e-7;
              }
            }
          }
          std::cout << "Jacobian maximum scaled error: " << maximum_error << "\n";
          return success;
        }
      };

      typename Outer::ModelDataT outerData()
      {
        typename Outer::ModelDataT data;
        using P         = typename Outer::ModelDataT::Parameters;
        data.parameters = {{P::C, .0001}, {P::Kp, .04}, {P::Ki, 2.0}, {P::Kaw, 400.0}};
        return data;
      }

    public:
      Testing::TestOutcome outerControl()
      {
        Testing::TestStatus success = true;
        Fixture<Outer>      f(outerData(), {210, 5, 208, 3, 8, -2, 377, 9, 6});
        success *= f.model.verify() == 0 && f.model.size() == 4;
        success *= f.model.initialize() == 0;
        success *= std::abs(f.model.y().getData()[0]) < 1e-12;
        success *= std::abs(f.model.y().getData()[1]) < 1e-12;
        success *= f.algebraicResidualsZero();
        success *= f.model.initializeState({{"icmdd", 11}, {"icmdq", 7}}) == 0;
        success *= f.alias.read() == 11 && f.algebraicResidualsZero();
        // Independent componentwise feedforward and PI balance.
        success *= std::abs(f.model.y().getData()[0] - (11 - 8 + 377 * .0001 * 3 - .04 * 2)) < 1e-12;
        success *= std::abs(f.model.y().getData()[1] - (7 + 2 - 377 * .0001 * 208 - .04 * 2)) < 1e-12;
        success *= rejects([&]
                           { f.model.initializeState({{"etad", 0}}); });
        success *= rejects([&]
                           { f.model.initializeState({{"icmdd", std::numeric_limits<double>::infinity()}}); });
        success *= f.model.initializationPorts().inputs.size() == 7;
#ifdef GRIDKIT_ENABLE_ENZYME
        f.model.tagDifferentiable();
        success *= f.model.tag()[0] && f.model.tag()[1] && !f.model.tag()[2] && !f.model.tag()[3];
        success *= f.jacobian();
#endif
        auto invalid                                         = outerData();
        invalid.parameters[Outer::ModelDataT::Parameters::C] = -1.0;
        Fixture<Outer> bad(invalid, f.values);
        success *= bad.model.verify() != 0 && bad.model.initialize() != 0;
        return success.report("OuterVoltageControl output initialization, ownership and Jacobian");
      }
    };
  } // namespace Testing
} // namespace GridKit
