/**
 * @file InnerCurrentControlTests.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Output initialization and finite-difference Jacobian checks.
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <vector>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/EMT/Component/Controller/InnerCurrentControl/InnerCurrentControl.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <typename ScalarT, typename IdxT>
    class InnerCurrentControlTests
    {
      using Inner = EMT::Controller::InnerCurrentControl<ScalarT, IdxT>;

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
          for (size_t n = 2; n < 6; ++n)
            if (std::abs(model.getResidual().getData()[n]) > 1e-10)
              return false;
          return true;
        }

        bool jacobian()
        {
          bool                  success       = true;
          double                maximum_error = 0;
          std::array<double, 5> worst{};
          for (const auto& [ys, yps] : {std::pair{1.0, 0.0}, std::pair{0.0, 1.0}, std::pair{2.0, 3.0}, std::pair{0.0, 0.0}})
          {
            std::map<std::pair<size_t, size_t>, double> entries;
            for (const auto& entry : model.jacobianEntries(ys, yps))
              entries[{entry.row, entry.column}] += entry.value;
            for (size_t j = 0; j < count + model.size(); ++j)
            {
              double&             y        = j < count ? values[j] : model.y().getData()[j - count];
              const double        original = y;
              const double        h        = (ys == 0 ? 1e-3 : 1e-5) * (1 + std::abs(original));
              std::vector<double> plus(model.size()), minus(model.size());
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
                const double fd    = (plus[n] - minus[n]) / (2 * h);
                const double error = std::abs(entries[{count + n, j}] - fd) / (1 + std::abs(fd));
                if (error > maximum_error)
                  worst = {static_cast<double>(n), static_cast<double>(j), entries[{count + n, j}], fd, ys};
                maximum_error  = std::max(maximum_error, error);
                success       &= error < 1e-7;
              }
            }
          }
          std::cout << "Jacobian maximum scaled error: " << maximum_error
                    << " at row " << worst[0] << ", column " << worst[1]
                    << " (" << worst[2] << " vs " << worst[3] << ", y scale " << worst[4] << ")\n";
          return success;
        }
      };

      typename Inner::ModelDataT innerData()
      {
        typename Inner::ModelDataT data;
        using P         = typename Inner::ModelDataT::Parameters;
        data.parameters = {{P::L, .002}, {P::Kp, 5.0}, {P::Ki, 500.0}, {P::Kaw, 2000.0}, {P::Imax, 30.0}};
        return data;
      }

    public:
      Testing::TestOutcome innerControl()
      {
        Testing::TestStatus success = true;
        Fixture<Inner>      f(innerData(), {208, 3, 8, -2, 8, -2, 377, 211, 9});
        success *= f.model.verify() == 0 && f.model.size() == 8;
        success *= f.model.initializeState({{"ud", 211}, {"uq", 9}}) == 0;
        success *= f.algebraicResidualsZero();
        success *= std::abs(f.model.y().getData()[0] - (211 - 208 - 377 * .002 * 2)) < 1e-10;
        success *= std::abs(f.model.y().getData()[1] - (9 - 3 - 377 * .002 * 8)) < 1e-10;
        success *= rejects([&]
                           { f.model.initializeState({{"xid", 0}}); });
        success *= rejects([&]
                           { f.model.initializeState({{"ilimd", 100}}); });
        success *= rejects([&]
                           { f.model.initializeState({{"ud", std::numeric_limits<double>::infinity()}}); });
        success *= f.model.initializationPorts().inputs.size() == 7;
#ifdef GRIDKIT_ENABLE_ENZYME
        f.model.tagDifferentiable();
        for (size_t n = 0; n < 8; ++n)
          success *= f.model.tag()[n] == (n < 2 || n >= 6);
        success *= f.jacobian();
#endif
        for (const double command : {29.0, 30.0, 45.0})
        {
          f.values[4]      = command;
          f.values[5]      = 0;
          success         *= f.model.initialize() == 0 && f.algebraicResidualsZero();
          const double id  = f.model.outputSignal(Inner::Outputs::ilimd).read();
          const double iq  = f.model.outputSignal(Inner::Outputs::ilimq).read();
          const double ud  = f.model.outputSignal(Inner::Outputs::ud).read();
          const double uq  = f.model.outputSignal(Inner::Outputs::uq).read();
          success         *= std::hypot(id, iq) <= 30 + 1e-12 && iq == 0 && id > 0;
          // The default voltage command is the feedforward plus proportional action.
          success         *= std::abs(ud - (208 + 377 * .002 * 2 + 5 * (id - 8))) < 1e-9;
          success         *= std::abs(uq - (3 + 377 * .002 * 8 + 5 * (iq + 2))) < 1e-9;
          // Tracking the unclipped command leaves the pure integral rate.
          f.values[7]      = ud;
          f.values[8]      = uq;
          f.model.evaluateResidual();
          success *= std::abs(f.model.getResidual().getData()[0] - 500 * (id - 8)) < 1e-9;
#ifdef GRIDKIT_ENABLE_ENZYME
          success *= f.jacobian();
#endif
          // A clipped command from the voltage limiter restores the integral.
          f.values[7] = 100;
          f.values[8] = 0;
          f.model.evaluateResidual();
          success *= f.model.getResidual().getData()[0] < 0;
        }
        auto compensated                                         = innerData();
        compensated.parameters[Inner::ModelDataT::Parameters::C] = 1e-4;
        Fixture<Inner> c(compensated, {208, 3, 8, -2, 8 + 377 * 1e-4 * 3, -2 - 377 * 1e-4 * 208, 377, 211, 9});
        success *= c.model.initializeState({{"ud", 211}, {"uq", 9}}) == 0;
        success *= c.algebraicResidualsZero();
        c.model.tagDifferentiable();
        success *= c.model.size() == 8 && c.model.tag()[6] && c.model.tag()[7];
        success *= c.model.getResidual().getData()[6] == 0 && c.model.getResidual().getData()[7] == 0;
        success *= std::abs(c.model.outputSignal(Inner::Outputs::ilimd).read() - 8) < 1e-10;
        success *= std::abs(c.model.outputSignal(Inner::Outputs::ilimq).read() + 2) < 1e-10;
#ifdef GRIDKIT_ENABLE_ENZYME
        success *= c.jacobian();
#endif
        auto invalid                                            = innerData();
        invalid.parameters[Inner::ModelDataT::Parameters::Imax] = -1.0;
        Fixture<Inner> bad(invalid, f.values);
        success *= bad.model.verify() != 0 && bad.model.initialize() != 0;
        return success.report("InnerCurrentControl output initialization, smooth current limit and Jacobian");
      }
    };
  } // namespace Testing
} // namespace GridKit
