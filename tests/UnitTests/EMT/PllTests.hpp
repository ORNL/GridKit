/**
 * @file PllTests.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief PLL tests on a hand-assembled bus and reference operator.
 */
#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numbers>

#include <GridKit/Definitions.hpp>
#include <GridKit/LinearAlgebra/Vector/Vector.hpp>
#include <GridKit/Model/EMT/Component/Bus/Bus.hpp>
#include <GridKit/Model/EMT/Operators/Reference/PLL/Pll.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <typename ScalarT, typename IdxT>
    class PllTests
    {
      using RealT                          = ScalarT;
      using VectorT                        = LinearAlgebra::Vector<ScalarT, IdxT>;
      using BusT                           = EMT::Bus<ScalarT, IdxT>;
      using PllT                           = EMT::Pll<ScalarT, IdxT>;
      using Data                           = EMT::PllData<ScalarT, IdxT>;
      static constexpr IdxT  system_size   = 6;
      static constexpr RealT omega         = RealT{120} * std::numbers::pi_v<RealT>;
      static constexpr RealT rated_voltage = 480.0;
      static constexpr RealT kp            = 80.0;
      static constexpr RealT ki            = 2500.0;

      static Data makeData()
      {
        Data data;
        using P                = typename Data::Parameters;
        data.parameters[P::V]  = rated_voltage;
        data.parameters[P::f]  = RealT{60};
        data.parameters[P::Kp] = kp;
        data.parameters[P::Ki] = ki;
        return data;
      }

      struct Fixture
      {
        VectorT y, yp, f, abs_tol;
        BusT    bus;
        PllT    pll;

        explicit Fixture(const Data& data = makeData())
          : pll(data)
        {
          y.resize(system_size);
          yp.resize(system_size);
          f.resize(system_size);
          abs_tol.resize(system_size);
          pll.attachInput(bus.voltages());
          IdxT offset = 0;
          for (auto* component : components())
          {
            component->bind(y, yp, f, abs_tol, offset);
            component->allocate();
            component->assignGlobalIndices(offset);
            offset += component->size();
          }
          bus.initialize();
          setTerminalSamples(0.37);
          pll.initialize();
          for (auto* component : components())
            component->tagDifferentiable();
        }

        std::array<EMT::Component<ScalarT, IdxT>*, 2> components()
        {
          return {&bus, &pll};
        }

        void updateTime(RealT time, RealT alpha)
        {
          for (auto* component : components())
            component->updateTime(time, alpha);
        }

        void evaluateResidual()
        {
          for (auto* component : components())
            component->evaluateInternalResidual();
          for (auto* component : components())
            component->evaluateExternalResidual();
        }

        void setTerminalSamples(RealT angle, RealT frequency = omega)
        {
          const RealT                peak = rated_voltage * std::sqrt(RealT{2} / 3);
          const std::array<RealT, 3> offset{0, -RealT{2} * std::numbers::pi_v<RealT> / 3, RealT{2} * std::numbers::pi_v<RealT> / 3};
          for (size_t p = 0; p < 3; ++p)
          {
            y.getData()[p]  = peak * std::cos(angle + offset[p]);
            yp.getData()[p] = -frequency * peak * std::sin(angle + offset[p]);
          }
          y.setDataUpdated();
          yp.setDataUpdated();
        }

        void setProbeState()
        {
          auto* values = y.getData();
          values[0]    = 170.0;
          values[1]    = -91.0;
          values[2]    = 37.0;
          values[3]    = -0.71;
          values[4]    = 0.0037;
          values[5]    = 370.0;
          for (IdxT n = 0; n < system_size; ++n)
            yp.getData()[n] = 0.31 * static_cast<RealT>(n);
        }

        // Solve the scalar frequency row, then use the two differential rows
        // as the RHS. The prescribed stiff bus is independent of the PLL.
        std::array<RealT, 2> derivative(RealT time, const std::array<RealT, 2>& state, RealT phase, RealT frequency)
        {
          setTerminalSamples(frequency * time + phase, frequency);
          auto* values    = y.getData();
          values[3]       = state[0];
          values[4]       = state[1];
          values[5]       = omega;
          yp.getData()[3] = yp.getData()[4] = 0.0;
          evaluateResidual();
          values[5] -= f.getData()[5];
          evaluateResidual();
          return {-f.getData()[3], -f.getData()[4]};
        }

        void advance(RealT time, RealT step, RealT phase, RealT frequency)
        {
          const std::array<RealT, 2> state{y.getData()[3], y.getData()[4]};
          const auto                 shifted = [&](const auto& slope, RealT scale)
          { return std::array<RealT, 2>{state[0] + scale * slope[0], state[1] + scale * slope[1]}; };
          const auto           k1 = derivative(time, state, phase, frequency);
          const auto           k2 = derivative(time + step / 2, shifted(k1, step / 2), phase, frequency);
          const auto           k3 = derivative(time + step / 2, shifted(k2, step / 2), phase, frequency);
          const auto           k4 = derivative(time + step, shifted(k3, step), phase, frequency);
          std::array<RealT, 2> next{};
          for (size_t n = 0; n < next.size(); ++n)
            next[n] = state[n] + step / 6 * (k1[n] + 2 * k2[n] + 2 * k3[n] + k4[n]);
          derivative(time + step, next, phase, frequency);
        }
      };

    public:
      TestOutcome wiring()
      {
        TestStatus success = true;
        Fixture    fixture;
        success *= fixture.pll.size() == 3;
        success *= fixture.pll.verify() == 0;
        success *= &fixture.pll.inputSignal(EMT::PllInputs::va) == &fixture.bus.outputSignal(EMT::BusOutputs::va);
        success *= fixture.pll.tag()[0] && fixture.pll.tag()[1] && !fixture.pll.tag()[2];
        success *= isEqual(fixture.pll.outputSignal(EMT::PllOutputs::theta).read(), RealT{0.37}, 1e-12);
        success *= isEqual(fixture.pll.outputSignal(EMT::PllOutputs::omega).read(), omega, 1e-12);
        PllT invalid;
        success *= invalid.verify() == 5;
        return success.report(__func__);
      }

      TestOutcome initialState()
      {
        TestStatus success = true;
        Fixture    fixture;
        fixture.evaluateResidual();
        for (IdxT n = 3; n < system_size; ++n)
          success *= std::abs(fixture.f.getData()[n]) < 1e-12;
        using Outputs = typename PllT::Outputs;
        fixture.pll.initialize({{Outputs::theta, 0.1}});
        success               *= isEqual(fixture.y.getData()[5], omega, 1e-12);
        success               *= std::abs(fixture.y.getData()[4] + kp / ki * std::sin(0.27)) < 1e-12;
        const RealT frequency  = omega + 2 * std::numbers::pi_v<RealT>;
        fixture.pll.initializeState({{"theta", 0.1}, {"omega", frequency}});
        success *= fixture.y.getData()[5] == frequency;
        success *= std::abs(fixture.y.getData()[4] - (frequency - omega - kp * std::sin(0.27)) / ki) < 1e-12;
        fixture.evaluateResidual();
        for (IdxT n = 3; n < system_size; ++n)
          success *= std::abs(fixture.f.getData()[n]) < 1e-12;
        fixture.pll.initialize({{Outputs::omega, frequency}});
        success *= std::abs(fixture.y.getData()[3] - 0.37) < 1e-12;
        success *= std::abs(fixture.y.getData()[4] - (frequency - omega) / ki) < 1e-12;
        for (const auto& values : {std::map<std::string, RealT>{{"xi", 0.0}},
                                   std::map<std::string, RealT>{{"omega", std::numeric_limits<RealT>::infinity()}},
                                   std::map<std::string, RealT>{{"theta", std::numeric_limits<RealT>::quiet_NaN()}}})
        {
          bool rejected = false;
          try
          {
            fixture.pll.initializeState(values);
          }
          catch (const std::invalid_argument&)
          {
            rejected = true;
          }
          success *= rejected;
        }
        for (size_t n = 0; n < 3; ++n)
          fixture.y.getData()[n] = 0;
        bool rejected = false;
        try
        {
          fixture.pll.initializeState({});
        }
        catch (const std::invalid_argument&)
        {
          rejected = true;
        }
        success *= rejected;
        success *= fixture.pll.initializeState({{"theta", 0.2}, {"omega", frequency}}) == 0;
        success *= fixture.y.getData()[5] == frequency;
        success *= std::abs(fixture.y.getData()[4] - (frequency - omega) / ki) < 1e-12;
        return success.report(__func__);
      }

      TestOutcome residual()
      {
        TestStatus success = true;
        Fixture    fixture;
        fixture.setProbeState();
        fixture.evaluateResidual();
        const auto* y      = fixture.y.getData();
        const auto* yp     = fixture.yp.getData();
        const auto* f      = fixture.f.getData();
        const RealT offset = 2 * std::numbers::pi_v<RealT> / 3;
        const RealT vq     = -std::sqrt(RealT{2} / 3) / rated_voltage
                         * (y[0] * std::sin(y[3]) + y[1] * std::sin(y[3] - offset)
                            + y[2] * std::sin(y[3] + offset));
        success *= std::abs(f[3] - (yp[3] - y[5])) < 1e-12;
        success *= std::abs(f[4] - (yp[4] - vq)) < 1e-12;
        success *= std::abs(f[5] - (y[5] - omega - kp * vq - ki * y[4])) < 1e-12;
        return success.report(__func__);
      }

      TestOutcome lockIn()
      {
        TestStatus success = true;
        Fixture    fixture;
        fixture.setTerminalSamples(-0.4);
        fixture.pll.initialize();
        const RealT frequency = omega + 2 * std::numbers::pi_v<RealT>;
        const RealT step      = 1e-4;
        for (size_t n = 0; n < 10000; ++n)
          fixture.advance(step * static_cast<RealT>(n), step, 0.37, frequency);
        const RealT phase_error     = std::abs(fixture.y.getData()[3] - frequency - 0.37);
        const RealT frequency_error = std::abs(fixture.y.getData()[5] - frequency);
        const RealT integral_error  = std::abs(fixture.y.getData()[4] - (frequency - omega) / ki);
        std::cout << "PLL lock errors (rad, rad/s, s): " << phase_error << ", " << frequency_error << ", " << integral_error << "\n";
        success *= phase_error < 1e-11 && frequency_error < 1e-8 && integral_error < 1e-11;
        return success.report(__func__);
      }

      TestOutcome phaseJump()
      {
        TestStatus success = true;
        // Independent linearized phase-error solution: e'' + Kp e' + Ki e = 0,
        // e(0+) = jump, e'(0+) = -Kp jump. A small jump measures sign and gain.
        Fixture    fixture;
        fixture.setTerminalSamples(0.0);
        fixture.pll.initialize();
        const RealT jump          = 1e-4;
        const RealT decay         = kp / 2;
        const RealT damped        = std::sqrt(ki - decay * decay);
        const RealT step          = 1e-4;
        RealT       maximum_phase = 0, maximum_frequency = 0;
        for (size_t n = 0; n < 5000; ++n)
        {
          const RealT time = step * static_cast<RealT>(n);
          fixture.advance(time, step, jump, omega);
          const RealT t = time + step;
          const RealT c = std::cos(damped * t), s = std::sin(damped * t);
          const RealT envelope  = jump * std::exp(-decay * t);
          const RealT error     = envelope * (c - decay / damped * s);
          const RealT deviation = envelope * (kp * c + (ki - 2 * decay * decay) / damped * s);
          maximum_phase         = std::max(maximum_phase, std::abs(fixture.y.getData()[3] - (omega * t + jump - error)));
          maximum_frequency     = std::max(maximum_frequency, std::abs(fixture.y.getData()[5] - omega - deviation));
        }
        std::cout << "PLL small-jump maximum errors (rad, rad/s): " << maximum_phase << ", " << maximum_frequency << "\n";
        success                  *= maximum_phase < 1e-11 && maximum_frequency < 1e-9;
        // A finite jump must settle in the same direction without resetting states.
        const RealT theta_before  = fixture.y.getData()[3];
        fixture.setTerminalSamples(omega * 0.5 + 0.3);
        success *= fixture.y.getData()[3] == theta_before;
        for (size_t n = 0; n < 5000; ++n)
          fixture.advance(0.5 + step * static_cast<RealT>(n), step, 0.3, omega);
        const RealT phase_error     = std::abs(fixture.y.getData()[3] - omega - 0.3);
        const RealT frequency_error = std::abs(fixture.y.getData()[5] - omega);
        std::cout << "PLL finite-jump final errors (rad, rad/s): " << phase_error << ", " << frequency_error << "\n";
        success *= phase_error < 5e-9 && frequency_error < 2e-7;
        return success.report(__func__);
      }

      TestOutcome jacobian()
      {
        TestStatus success = true;

        const RealT alpha = 3.7;

        Fixture fixture;
        fixture.setProbeState();
        fixture.updateTime(0.0, alpha);
        fixture.evaluateResidual();

        for (auto* component : fixture.components())
        {
          component->evaluateJacobian();
        }

        std::map<std::pair<IdxT, IdxT>, RealT> enzyme_entries;
        for (auto* component : fixture.components())
        {
          auto* coo = component->getCooJacobian();
          if (coo == nullptr)
          {
            continue;
          }
          const IdxT  entry_count = coo->getNnz();
          const auto* rows        = coo->getRowData();
          const auto* cols        = coo->getColData();
          const auto* vals        = coo->getValues();
          for (IdxT i = 0; i < entry_count; ++i)
          {
            enzyme_entries[{rows[i], cols[i]}] += vals[i];
          }
        }

        success *= (!enzyme_entries.empty());

        auto* y_data  = fixture.y.getData();
        auto* yp_data = fixture.yp.getData();
        auto* f_data  = fixture.f.getData();

        RealT maximum_error = 0.0;
        for (IdxT j = 0; j < system_size; ++j)
        {
          const RealT step = 1.0e-5 * (1.0 + std::abs(y_data[j]));

          std::array<RealT, system_size> fd_column{};

          const RealT y_saved = y_data[j];
          y_data[j]           = y_saved + step;
          fixture.evaluateResidual();
          for (IdxT i = 0; i < system_size; ++i)
          {
            fd_column[static_cast<size_t>(i)] = f_data[i];
          }
          y_data[j] = y_saved - step;
          fixture.evaluateResidual();
          for (IdxT i = 0; i < system_size; ++i)
          {
            fd_column[static_cast<size_t>(i)] = (fd_column[static_cast<size_t>(i)] - f_data[i]) / (2.0 * step);
          }
          y_data[j] = y_saved;

          const RealT yp_saved = yp_data[j];
          yp_data[j]           = yp_saved + step;
          fixture.evaluateResidual();
          std::array<RealT, system_size> fp_plus{};
          for (IdxT i = 0; i < system_size; ++i)
          {
            fp_plus[static_cast<size_t>(i)] = f_data[i];
          }
          yp_data[j] = yp_saved - step;
          fixture.evaluateResidual();
          for (IdxT i = 0; i < system_size; ++i)
          {
            fd_column[static_cast<size_t>(i)] += alpha * (fp_plus[static_cast<size_t>(i)] - f_data[i]) / (2.0 * step);
          }
          yp_data[j] = yp_saved;

          for (IdxT i = 0; i < system_size; ++i)
          {
            const RealT fd_value     = fd_column[static_cast<size_t>(i)];
            const auto  it           = enzyme_entries.find({i, j});
            RealT       enzyme_value = 0.0;
            if (it != enzyme_entries.end())
            {
              enzyme_value = it->second;
            }
            const RealT error  = std::abs(enzyme_value - fd_value) / (1.0 + std::abs(fd_value));
            maximum_error      = std::max(maximum_error, error);
            success           *= (error < 5.0e-9);
          }
        }

        std::cout << "PLL Jacobian maximum scaled error: " << maximum_error << "\n";
        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
