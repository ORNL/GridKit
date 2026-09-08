/**
 * @file TransformerTests.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Transformer model tests on a hand-assembled two-bus system.
 *
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numbers>
#include <utility>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/LinearAlgebra/Vector/Vector.hpp>
#include <GridKit/Model/EMT/Component/Bus/Bus.hpp>
#include <GridKit/Model/EMT/Component/Transformer/Transformer.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    /**
     * @brief Tests for the transformer bank between two hand-assembled buses.
     *
     * The connection maps, the per-unit bases, and the two-slope magnetizing
     * characteristic exercise the terminal-current gradients and the
     * saturation derivatives.
     */
    template <typename ScalarT, typename IdxT>
    class TransformerTests
    {
    public:
      using RealT        = ScalarT;
      using VectorT      = GridKit::LinearAlgebra::Vector<ScalarT, IdxT>;
      using BusT         = GridKit::EMT::Bus<ScalarT, IdxT>;
      using TransformerT = GridKit::EMT::Transformer<ScalarT, IdxT>;

      static constexpr IdxT system_size = 33;

      TransformerTests()  = default;
      ~TransformerTests() = default;

    private:
      static constexpr RealT omega = RealT{120} * std::numbers::pi_v<RealT>;

      static GridKit::EMT::TransformerData<ScalarT, IdxT> makeResidualData()
      {
        using Data      = GridKit::EMT::TransformerData<ScalarT, IdxT>;
        using Parameter = typename Data::Parameters;
        Data data;
        data.parameters[Parameter::S]     = RealT{100.0e6};
        data.parameters[Parameter::V1]    = RealT{138.0e3};
        data.parameters[Parameter::V2]    = RealT{69.0e3};
        data.parameters[Parameter::f]     = RealT{60.0};
        data.parameters[Parameter::tap]   = RealT{1.025};
        data.parameters[Parameter::R]     = RealT{0.005};
        data.parameters[Parameter::X]     = RealT{0.1};
        data.parameters[Parameter::I0]    = RealT{0.01};
        data.parameters[Parameter::P0]    = RealT{60.0e3};
        data.parameters[Parameter::knee]  = RealT{1.2};
        data.parameters[Parameter::Lsat]  = RealT{0.25};
        data.parameters[Parameter::split] = RealT{0.4};
        return data;
      }

      /// Nominal ratio with winding 2 in delta across phases ab, bc, and ca.
      static GridKit::EMT::TransformerData<ScalarT, IdxT> makeDeltaData()
      {
        using Data                      = GridKit::EMT::TransformerData<ScalarT, IdxT>;
        using Parameter                 = typename Data::Parameters;
        Data data                       = makeResidualData();
        data.parameters[Parameter::tap] = RealT{1.0};
        data.parameters[Parameter::P2]  = GridKit::EMT::ABCMatrix<RealT>{{{{1.0, 0.0, -1.0}},
                                                                          {{-1.0, 1.0, 0.0}},
                                                                          {{0.0, -1.0, 1.0}}}};
        return data;
      }

      /**
       * @brief Two buses joined by one transformer bank.
       *
       * Variable layout: bus 1 voltage [0, 3), bus 2 voltage [3, 6),
       * transformer variables [6, 33).
       */
      struct Fixture
      {
        VectorT y;
        VectorT yp;
        VectorT f;
        VectorT abs_tol;

        BusT         bus1;
        BusT         bus2;
        TransformerT transformer;

        explicit Fixture(const GridKit::EMT::TransformerData<ScalarT, IdxT>& data = makeResidualData())
          : transformer(data)
        {
          y.resize(system_size);
          yp.resize(system_size);
          f.resize(system_size);
          abs_tol.resize(system_size);

          transformer.attachTerminal(0, bus1.voltages());
          transformer.attachTerminal(1, bus2.voltages());
          for (size_t p = 0; p < 3; ++p)
          {
            bus1.addCurrent(p, transformer.currentSignal(0, p));
            bus2.addCurrent(p, transformer.currentSignal(1, p));
          }

          IdxT offset = 0;
          for (auto* component : components())
          {
            component->bind(y, yp, f, abs_tol, offset);
            component->allocate();
            component->assignGlobalIndices(offset);
            offset += component->size();
          }

          bus1.initialize();
          bus2.initialize();
          transformer.initialize();
          for (auto* component : components())
          {
            component->tagDifferentiable();
          }
        }

        std::array<GridKit::EMT::Component<ScalarT, IdxT>*, 3> components()
        {
          return {&bus1, &bus2, &transformer};
        }

        void updateTime(RealT t, RealT alpha)
        {
          for (auto* component : components())
          {
            component->updateTime(t, alpha);
          }
        }

        void evaluateResidual()
        {
          for (auto* component : components())
          {
            component->evaluateInternalResidual();
          }
          for (auto* component : components())
          {
            component->evaluateExternalResidual();
          }
        }

        /// Balanced sinusoidal terminal samples at the initialization instant.
        void setTerminalSamples(RealT peak1, RealT angle1, RealT peak2, RealT angle2)
        {
          auto*                      y_data  = y.getData();
          auto*                      yp_data = yp.getData();
          const RealT                gamma   = RealT{2} * std::numbers::pi_v<RealT> / RealT{3};
          const std::array<RealT, 3> offset{RealT{0}, -gamma, gamma};
          for (size_t p = 0; p < 3; ++p)
          {
            y_data[p]      = peak1 * std::cos(angle1 + offset[p]);
            yp_data[p]     = -omega * peak1 * std::sin(angle1 + offset[p]);
            y_data[3 + p]  = peak2 * std::cos(angle2 + offset[p]);
            yp_data[3 + p] = -omega * peak2 * std::sin(angle2 + offset[p]);
          }
          y.setDataUpdated();
          yp.setDataUpdated();
        }

        /// Bus voltages in volts and transformer states in per unit, with
        /// flux linkages straddling the knee.
        void setProbeState()
        {
          auto* y_data  = y.getData();
          auto* yp_data = yp.getData();
          for (IdxT j = 0; j < system_size; ++j)
          {
            RealT scale = 1.0;
            if (j < 6)
            {
              scale = 1.0e4;
            }
            y_data[j]  = scale * (0.37 + 0.11 * static_cast<RealT>(j));
            yp_data[j] = scale * (-1.7 + 0.9 * static_cast<RealT>(j));
          }
          y.setDataUpdated();
          yp.setDataUpdated();
        }
      };

    public:
      /**
       * @brief Wiring smoke test: terminals, injections, and classification.
       */
      TestOutcome wiring()
      {
        TestStatus success = true;

        Fixture fixture;

        success *= (fixture.transformer.size() == 27);
        success *= &fixture.transformer.inputSignal(EMT::TransformerInputs::v1a)
                   == &fixture.bus1.outputSignal(EMT::BusOutputs::va);
        success *= &fixture.transformer.inputSignal(EMT::TransformerInputs::v2c)
                   == &fixture.bus2.outputSignal(EMT::BusOutputs::vc);
        success *= (fixture.bus1.verify() == 0);
        success *= (fixture.bus2.verify() == 0);
        success *= (fixture.transformer.verify() == 0);

        // The series current and flux linkages are differential; the node
        // voltages and currents are algebraic, as are the bus voltages
        // without a shunt.
        for (IdxT j = 0; j < 9; ++j)
        {
          success *= (fixture.transformer.tag()[j] == true);
        }
        for (IdxT j = 9; j < 27; ++j)
        {
          success *= (fixture.transformer.tag()[j] == false);
        }
        for (IdxT p = 0; p < 3; ++p)
        {
          success *= (fixture.bus1.tag()[p] == false);
          success *= (fixture.bus2.tag()[p] == false);
        }

        return success.report(__func__);
      }

      /**
       * @brief State-file keys seed the differential states and reject the rest.
       */
      TestOutcome initialState()
      {
        TestStatus success = true;

        Fixture fixture;

        success       *= (fixture.transformer.initializeState({{"i12b", 0.3}, {"psi1a", -0.7}, {"psi2c", 1.1}}) == 0);
        const auto* y  = fixture.y.getData();
        for (IdxT j = 6; j < system_size; ++j)
        {
          RealT expected = 0.0;
          if (j == 7)
          {
            expected = 0.3;
          }
          if (j == 9)
          {
            expected = -0.7;
          }
          if (j == 14)
          {
            expected = 1.1;
          }
          success *= (y[j] == expected);
        }

        bool rejected = false;
        try
        {
          fixture.transformer.validateInitialState({{"i1a", 1.0}});
        }
        catch (const std::invalid_argument&)
        {
          rejected = true;
        }
        success *= rejected;

        rejected = false;
        try
        {
          fixture.transformer.initializeState({{"psi1a", std::numeric_limits<RealT>::infinity()}});
        }
        catch (const std::invalid_argument&)
        {
          rejected = true;
        }
        success *= rejected;

        return success.report(__func__);
      }

      /**
       * @brief Assembled residual against an independent computation.
       */
      TestOutcome residual()
      {
        TestStatus success = true;

        Fixture fixture;
        fixture.setProbeState();
        fixture.updateTime(0.0, 1.0);
        fixture.evaluateResidual();

        const auto* y  = fixture.y.getData();
        const auto* yp = fixture.yp.getData();
        const auto* f  = fixture.f.getData();

        const auto data  = makeResidualData();
        using Parameter  = typename GridKit::EMT::TransformerData<ScalarT, IdxT>::Parameters;
        const auto S     = std::get<RealT>(data.parameters.at(Parameter::S));
        const auto V1    = std::get<RealT>(data.parameters.at(Parameter::V1));
        const auto V2    = std::get<RealT>(data.parameters.at(Parameter::V2));
        const auto freq  = std::get<RealT>(data.parameters.at(Parameter::f));
        const auto tap   = std::get<RealT>(data.parameters.at(Parameter::tap));
        const auto R     = std::get<RealT>(data.parameters.at(Parameter::R));
        const auto X     = std::get<RealT>(data.parameters.at(Parameter::X));
        const auto I0    = std::get<RealT>(data.parameters.at(Parameter::I0));
        const auto P0    = std::get<RealT>(data.parameters.at(Parameter::P0));
        const auto knee  = std::get<RealT>(data.parameters.at(Parameter::knee));
        const auto Lsat  = std::get<RealT>(data.parameters.at(Parameter::Lsat));
        const auto split = std::get<RealT>(data.parameters.at(Parameter::split));

        // Derived parameters for identity connection maps
        const RealT omega_base = 2.0 * std::numbers::pi_v<RealT> * freq;
        const RealT vw1        = V1 / std::sqrt(3.0);
        const RealT vw2        = V2 / std::sqrt(3.0);
        const RealT v_peak1    = std::sqrt(2.0) * vw1;
        const RealT v_peak2    = std::sqrt(2.0) * vw2;
        const RealT i_peak1    = std::sqrt(2.0) * S / (3.0 * vw1);
        const RealT i_peak2    = std::sqrt(2.0) * S / (3.0 * vw2);
        const RealT R1         = 0.5 * R;
        const RealT R2         = 0.5 * R;
        const RealT Gc         = P0 / S;
        const RealT Lm         = 1.0 / std::sqrt(I0 * I0 - Gc * Gc);
        const RealT beta1      = split;
        const RealT beta2      = 1.0 - split;

        const auto ramp = [](RealT x)
        {
          const RealT mu = GridKit::Math::MU<RealT>;
          return std::max(x, RealT{0}) + std::log1p(std::exp(-mu * std::abs(x))) / mu;
        };
        const auto g = [&](RealT psi)
        {
          return psi / Lm + (1.0 / Lsat - 1.0 / Lm) * (ramp(psi - knee) - ramp(-psi - knee));
        };

        std::array<RealT, system_size> expected{};

        for (size_t n = 0; n < 3; ++n)
        {
          const RealT i12  = y[6 + n];
          const RealT psi1 = y[9 + n];
          const RealT psi2 = y[12 + n];
          const RealT e1   = y[15 + n];
          const RealT e2   = y[18 + n];
          const RealT im1  = y[21 + n];
          const RealT im2  = y[24 + n];
          const RealT iw1  = y[27 + n];
          const RealT iw2  = y[30 + n];

          expected[6 + n]  = X / omega_base * yp[6 + n] + e2 - e1;
          expected[9 + n]  = yp[9 + n] / omega_base - e1;
          expected[12 + n] = yp[12 + n] / omega_base - e2;
          expected[15 + n] = e1 - y[n] / v_peak1 + R1 * iw1;
          expected[18 + n] = e2 - tap * y[3 + n] / v_peak2 + tap * R2 * iw2;
          expected[21 + n] = im1 - beta1 * g(psi1) - beta1 * Gc * e1;
          expected[24 + n] = im2 - beta2 * g(psi2) - beta2 * Gc * e2;
          expected[27 + n] = iw1 - im1 - i12;
          expected[30 + n] = iw2 + tap * i12 - tap * im2;

          // The bus current-balance rows hold the uncancelled injections
          expected[n]     = -i_peak1 * iw1;
          expected[3 + n] = -i_peak2 * iw2;
        }

        for (IdxT j = 0; j < system_size; ++j)
        {
          success *= isEqual(f[j], expected[static_cast<size_t>(j)], 1.0e-12);
        }

        return success.report(__func__);
      }

      /**
       * @brief Steady-state initialization satisfies the assembled residual.
       */
      TestOutcome steadyState()
      {
        TestStatus success = true;

        Fixture fixture;

        const RealT peak1 = 138.0e3 * std::sqrt(2.0 / 3.0);
        const RealT peak2 = 0.97 * 69.0e3 * std::sqrt(2.0 / 3.0);
        fixture.setTerminalSamples(peak1, 0.05, peak2, -0.12);

        success *= (fixture.transformer.initializeSteadyState(omega) == 0);
        fixture.updateTime(0.0, 1.0);
        fixture.evaluateResidual();

        const auto* y  = fixture.y.getData();
        const auto* yp = fixture.yp.getData();
        const auto* f  = fixture.f.getData();
        for (IdxT row = 6; row < system_size; ++row)
        {
          success *= (std::abs(f[row]) < 1.0e-10);
        }

        // The rated terminal voltage magnetizes the core near one per unit
        for (size_t n = 0; n < 3; ++n)
        {
          success *= (std::abs(std::hypot(y[9 + n], yp[9 + n] / omega) - 1.0) < 2.0e-2);
        }

        success *= (fixture.transformer.initializeSteadyState(0.0) != 0);
        success *= (fixture.transformer.initializeSteadyState(std::numeric_limits<RealT>::quiet_NaN()) != 0);

        return success.report(__func__);
      }

      /**
       * @brief Delta connection map: the thirty-degree shift, the no-load
       * magnetizing shares, and the folded terminal injections.
       */
      TestOutcome connection()
      {
        TestStatus success = true;

        Fixture fixture(makeDeltaData());

        const auto data  = makeDeltaData();
        using Parameter  = typename GridKit::EMT::TransformerData<ScalarT, IdxT>::Parameters;
        const auto S     = std::get<RealT>(data.parameters.at(Parameter::S));
        const auto V2    = std::get<RealT>(data.parameters.at(Parameter::V2));
        const auto I0    = std::get<RealT>(data.parameters.at(Parameter::I0));
        const auto split = std::get<RealT>(data.parameters.at(Parameter::split));

        // Winding 2 spans phases a and b, so bus 2 lags the winding voltage
        // by thirty degrees at no load
        const RealT peak1  = 138.0e3 * std::sqrt(2.0 / 3.0);
        const RealT peak2  = 69.0e3 * std::sqrt(2.0 / 3.0);
        const RealT angle1 = 0.3;
        fixture.setTerminalSamples(peak1, angle1, peak2, angle1 - std::numbers::pi_v<RealT> / 6.0);

        success *= (fixture.transformer.initializeSteadyState(omega) == 0);
        fixture.updateTime(0.0, 1.0);
        fixture.evaluateResidual();

        const auto* y  = fixture.y.getData();
        const auto* yp = fixture.yp.getData();
        const auto* f  = fixture.f.getData();
        for (IdxT row = 6; row < system_size; ++row)
        {
          success *= (std::abs(f[row]) < 1.0e-10);
        }

        // The series current carries only the resistive drop mismatch and the
        // winding currents carry the magnetizing shares
        for (size_t n = 0; n < 3; ++n)
        {
          success *= (std::hypot(y[6 + n], yp[6 + n] / omega) < 1.0e-3);
          success *= (std::abs(std::hypot(y[27 + n], yp[27 + n] / omega) - split * I0) < 2.0e-4);
          success *= (std::abs(std::hypot(y[30 + n], yp[30 + n] / omega) - (1.0 - split) * I0) < 2.0e-4);
        }

        // Terminal 2 injections fold two winding currents through the map
        const RealT i_peak2  = std::sqrt(2.0) * S / (3.0 * V2);
        success             *= isEqual(static_cast<RealT>(fixture.transformer.currentSignal(1, 0).read()), -i_peak2 * (y[30] - y[32]), 1.0e-12);
        success             *= isEqual(static_cast<RealT>(fixture.transformer.currentSignal(1, 1).read()), -i_peak2 * (y[31] - y[30]), 1.0e-12);
        success             *= isEqual(static_cast<RealT>(fixture.transformer.currentSignal(1, 2).read()), -i_peak2 * (y[32] - y[31]), 1.0e-12);

        // Without the shift the delta winding drives series current
        fixture.setTerminalSamples(peak1, angle1, peak2, angle1);
        success *= (fixture.transformer.initializeSteadyState(omega) == 0);
        for (size_t n = 0; n < 3; ++n)
        {
          success *= (std::hypot(y[6 + n], yp[6 + n] / omega) > 1.0);
        }

        return success.report(__func__);
      }

      /**
       * @brief Enzyme Jacobian against central finite differences of the
       * assembled residual, J = dF/dy + alpha * dF/dyp.
       */
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

        // The state spans per-unit flux to tens of kilovolts, so the step
        // and the comparison floor scale with each column.
        const RealT zero_floor = 1.0e-7;
        for (IdxT j = 0; j < system_size; ++j)
        {
          const RealT step = 1.0e-7 * (1.0 + std::abs(y_data[j]));

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
            if (std::abs(fd_value) > zero_floor || it != enzyme_entries.end())
            {
              success *= isEqual(enzyme_value, fd_value, 5.0e-6);
            }
          }
        }

        return success.report(__func__);
      }
    };

  } // namespace Testing
} // namespace GridKit
