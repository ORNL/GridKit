/**
 * @file LineLumpedTests.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief LineLumped model tests on a hand-assembled two-bus system.
 *
 */
#pragma once

#include <cmath>
#include <map>
#include <utility>

#include <GridKit/Definitions.hpp>
#include <GridKit/LinearAlgebra/Vector/Vector.hpp>
#include <GridKit/Model/EMT/Component/Bus/Bus.hpp>
#include <GridKit/Model/EMT/Component/Line/LineLumped/LineLumped.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    /**
     * @brief Tests for the lumped line between two hand-assembled buses.
     *
     * The full asymmetric parameter matrices exercise every Jacobian block,
     * including the alpha-scaled external-derivative entries from the shunt
     * capacitance rows.
     */
    template <typename ScalarT, typename IdxT>
    class LineLumpedTests
    {
    public:
      using RealT   = ScalarT;
      using VectorT = GridKit::LinearAlgebra::Vector<ScalarT, IdxT>;
      using BusT    = GridKit::EMT::Bus<ScalarT, IdxT>;
      using LineT   = GridKit::EMT::LineLumped<ScalarT, IdxT>;

      static constexpr IdxT system_size = 15;

      LineLumpedTests()  = default;
      ~LineLumpedTests() = default;

    private:
      static GridKit::EMT::LineLumpedData<ScalarT, IdxT> makeResidualData()
      {
        using Data      = GridKit::EMT::LineLumpedData<ScalarT, IdxT>;
        using Parameter = typename Data::Parameters;
        Data data;
        data.parameters[Parameter::conductors] = GridKit::EMT::ABCVector<IdxT>{{1, 2, 3}};
        data.parameters[Parameter::dx]         = RealT{2.5};
        data.parameters[Parameter::Rp]         = GridKit::EMT::ABCMatrix<RealT>{{{{1.7, 0.21, 0.32}},
                                                                                 {{0.13, 1.9, 0.24}},
                                                                                 {{0.35, 0.16, 2.1}}}};
        data.parameters[Parameter::Lp]         = GridKit::EMT::ABCMatrix<RealT>{{{{0.041, 0.002, 0.013}},
                                                                                 {{0.004, 0.052, 0.015}},
                                                                                 {{0.016, 0.007, 0.063}}}};
        data.parameters[Parameter::Gp]         = GridKit::EMT::ABCMatrix<RealT>{{{{0.31, 0.021, 0.012}},
                                                                                 {{0.023, 0.42, 0.034}},
                                                                                 {{0.015, 0.026, 0.53}}}};
        data.parameters[Parameter::Cp]         = GridKit::EMT::ABCMatrix<RealT>{{{{0.61, 0.041, 0.022}},
                                                                                 {{0.043, 0.72, 0.054}},
                                                                                 {{0.025, 0.046, 0.83}}}};
        return data;
      }

      /**
       * @brief Two buses joined by one lumped line.
       *
       * Variable layout: bus 1 voltage/shunt current [0, 6), bus 2
       * voltage/shunt current [6, 12), series current [12, 15).
       */
      struct Fixture
      {
        VectorT y;
        VectorT yp;
        VectorT f;
        VectorT abs_tol;

        BusT  bus1;
        BusT  bus2;
        LineT line;

        Fixture()
          : line(makeResidualData())
        {
          y.resize(system_size);
          yp.resize(system_size);
          f.resize(system_size);
          abs_tol.resize(system_size);

          const auto data = makeResidualData();
          using Parameter = typename LineT::ModelDataT::Parameters;
          typename BusT::YDataT Y;
          Y.D              = std::get<EMT::ABCMatrix<RealT>>(data.parameters.at(Parameter::Gp));
          Y.E              = std::get<EMT::ABCMatrix<RealT>>(data.parameters.at(Parameter::Cp));
          const auto scale = 0.5 * std::get<RealT>(data.parameters.at(Parameter::dx));
          using Output     = EMT::LineLumpedOutputs;
          bus1.addNorton("line_1", Y, {&line.outputSignal(Output::i21a), &line.outputSignal(Output::i21b), &line.outputSignal(Output::i21c)}, scale);
          bus2.addNorton("line_2", Y, {&line.outputSignal(Output::i12a), &line.outputSignal(Output::i12b), &line.outputSignal(Output::i12c)}, scale);
          line.attachTerminal(0, bus1.voltages());
          line.attachTerminal(1, bus2.voltages());

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
          line.initialize();
          for (auto* component : components())
          {
            component->tagDifferentiable();
          }
        }

        std::array<GridKit::EMT::Component<ScalarT, IdxT>*, 3> components()
        {
          return {&bus1, &bus2, &line};
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

        void setProbeState()
        {
          auto* y_data  = y.getData();
          auto* yp_data = yp.getData();
          for (IdxT j = 0; j < system_size; ++j)
          {
            y_data[j]  = 2.3 + 1.1 * static_cast<RealT>(j);
            yp_data[j] = -1.7 + 0.9 * static_cast<RealT>(j);
          }
          y.setDataUpdated();
          yp.setDataUpdated();
        }
      };

    public:
      /** Independent terminal states, shunt ports, and the full bus Jacobian. */
      TestOutcome nortonTerminals()
      {
        TestStatus                      success = true;
        BusT                            bus;
        EMT::VectorFitData<RealT, IdxT> Y;
        Y.poles = {{-2.0, 0.0}};
        Y.residues.resize(1);
        for (size_t p = 0; p < 3; ++p)
        {
          Y.D[p][p]           = 0.5;
          Y.E[p][p]           = 0.25;
          Y.residues[0][p][p] = 3.0;
        }
        std::array<typename BusT::SignalT, 3> incident;
        typename BusT::PhaseSignals           inputs;
        for (size_t p = 0; p < 3; ++p)
        {
          auto* voltage = &bus.outputSignal(static_cast<EMT::BusOutputs>((p + 1) % 3));
          incident[p].setComputed([voltage]
                                  { return 0.1 * voltage->read(); },
                                  [voltage](typename BusT::SignalT::GradientT& gradient, RealT scale)
                                  { voltage->appendGradient(gradient, 0.1 * scale); });
          inputs[p] = &incident[p];
        }
        auto& first     = bus.addNorton("first", Y, inputs);
        auto& second    = bus.addNorton("second", Y, {&first.outputSignal(0), &first.outputSignal(1), &first.outputSignal(2)}, 2.0, {1, 2, 0});
        success        *= &first == &bus.norton("first");
        success        *= &bus.outputSignal("first_Ish_a") == &first.outputSignal(0);
        success        *= &bus.inputSignal("first_inc_a") == inputs[0];
        success        *= bus.template component<EMT::KCL<ScalarT, IdxT>>("KCL").size() == 3;
        success        *= bus.size() == 15;
        bool duplicate  = false;
        try
        {
          bus.addShunt("first", Y);
        }
        catch (const std::invalid_argument&)
        {
          duplicate = true;
        }
        success *= duplicate && bus.size() == 15;
        success *= bus.allocate() == 0;
        success *= bus.verify() == 0;
        success *= bus.initialize({{EMT::BusOutputs::va, 2.0}, {EMT::BusOutputs::vb, -4.0}, {EMT::BusOutputs::vc, 6.0}}) == 0;
        success *= bus.initializeSteadyState(0.0) == 0;
        bus.tagDifferentiable();
        bus.evaluateResidual();
        const std::array<RealT, 3> voltage{2.0, -4.0, 6.0};
        std::array<RealT, 3>       expected_kcl{};
        for (size_t p = 0; p < 3; ++p)
        {
          const size_t q         = (p + 1) % 3;
          const RealT  incoming  = 0.1 * voltage[q];
          // At DC, Y(0) = 0.5 + 3/2 = 2. The second instance has scale 2.
          success               *= isEqual(first.shuntCurrent(p), 2.0 * voltage[p], 1e-13);
          success               *= isEqual(first.outputSignal(p).read(), 2.0 * voltage[p], 1e-13);
          success               *= isEqual(second.shuntCurrent(p), 4.0 * voltage[q], 1e-13);
          success               *= isEqual(second.outputSignal(p).read(), 4.0 * voltage[q], 1e-13);
          expected_kcl[p]       += incoming - 2.0 * voltage[p];
          expected_kcl[q]       += 2.0 * voltage[p] - 4.0 * voltage[q];
          success               *= bus.tag()[p] && !first.tag()[p] && first.tag()[3 + p];
          success               *= first.getVariableIndex(static_cast<IdxT>(3 + p))
                     != second.getVariableIndex(static_cast<IdxT>(3 + p));
        }
        for (size_t p = 0; p < 3; ++p)
          success *= isEqual(bus.getResidual().getData()[p], expected_kcl[p], 1e-13);
        for (IdxT row = 3; row < bus.size(); ++row)
          success *= std::abs(bus.getResidual().getData()[row]) < 1e-12;

        const RealT alpha = 2.7, step = 1e-6;
        bus.updateTime(0.0, alpha);
        success *= bus.evaluateJacobian() == 0;
        std::map<std::pair<IdxT, IdxT>, RealT> jacobian;
        auto*                                  coo = bus.getCooJacobian();
        for (IdxT k = 0; k < coo->getNnz(); ++k)
          jacobian[{coo->getRowData()[k], coo->getColData()[k]}] += coo->getValues()[k];
        for (IdxT col = 0; col < bus.size(); ++col)
        {
          auto*      y       = bus.y().getData();
          auto*      yp      = bus.yp().getData();
          const auto saved_y = y[col], saved_yp = yp[col];
          y[col]  = saved_y + step;
          yp[col] = saved_yp + alpha * step;
          bus.evaluateResidual();
          std::vector<RealT> plus(bus.getResidual().getData(), bus.getResidual().getData() + bus.size());
          y[col]  = saved_y - step;
          yp[col] = saved_yp - alpha * step;
          bus.evaluateResidual();
          for (IdxT row = 0; row < bus.size(); ++row)
          {
            const RealT difference  = (plus[static_cast<size_t>(row)] - bus.getResidual().getData()[row]) / (2 * step);
            success                *= std::abs(jacobian[{row, col}] - difference) < 2e-8;
          }
          y[col]  = saved_y;
          yp[col] = saved_yp;
        }
        bool frozen = false;
        try
        {
          bus.addShunt("late", Y);
        }
        catch (const std::logic_error&)
        {
          frozen = true;
        }
        success *= frozen;
        return success.report(__func__);
      }

      /**
       * @brief Wiring smoke test: shunt capacitance makes the bus voltages
       * differential.
       */
      TestOutcome wiring()
      {
        TestStatus success = true;

        Fixture fixture;

        success *= (fixture.line.size() == 3);
        success *= &fixture.line.inputSignal(EMT::LineLumpedInputs::v1a)
                   == &fixture.bus1.outputSignal(EMT::BusOutputs::va);
        success *= &fixture.line.inputSignal(EMT::LineLumpedInputs::v2c)
                   == &fixture.bus2.outputSignal(EMT::BusOutputs::vc);
        success *= &fixture.bus1.inputSignal("line_1_inc_a") == &fixture.line.outputSignal(EMT::LineLumpedOutputs::i21a);
        success *= &fixture.bus2.inputSignal("line_2_inc_c") == &fixture.line.outputSignal(EMT::LineLumpedOutputs::i12c);
        success *= (fixture.bus1.verify() == 0);
        success *= (fixture.bus2.verify() == 0);
        success *= (fixture.line.verify() == 0);

        success *= (fixture.bus1.tag()[0] == true);
        success *= (fixture.bus1.tag()[1] == true);
        success *= (fixture.bus1.tag()[2] == true);
        success *= (fixture.bus2.tag()[0] == true);
        success *= (fixture.bus2.tag()[1] == true);
        success *= (fixture.bus2.tag()[2] == true);

        success *= (fixture.line.tag()[0] == true);
        success *= (fixture.bus1.tag()[3] == false);
        success *= (fixture.bus2.tag()[3] == false);

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

        const auto data = makeResidualData();
        using Parameter = typename GridKit::EMT::LineLumpedData<ScalarT, IdxT>::Parameters;
        const auto dx   = std::get<RealT>(data.parameters.at(Parameter::dx));
        const auto Rp   = std::get<GridKit::EMT::ABCMatrix<RealT>>(data.parameters.at(Parameter::Rp));
        const auto Lp   = std::get<GridKit::EMT::ABCMatrix<RealT>>(data.parameters.at(Parameter::Lp));
        const auto Gp   = std::get<GridKit::EMT::ABCMatrix<RealT>>(data.parameters.at(Parameter::Gp));
        const auto Cp   = std::get<GridKit::EMT::ABCMatrix<RealT>>(data.parameters.at(Parameter::Cp));

        std::array<RealT, system_size> expected{};

        // The series current is positive from bus 1 to bus 2.
        for (size_t n = 0; n < 3; ++n)
        {
          expected[n]     = -y[12 + n] - y[3 + n];
          expected[6 + n] = y[12 + n] - y[9 + n];
        }

        for (size_t n = 0; n < 3; ++n)
        {
          RealT series = y[6 + n] - y[n];
          RealT shunt1 = -y[3 + n];
          RealT shunt2 = -y[9 + n];
          for (size_t k = 0; k < 3; ++k)
          {
            series += dx * Rp[n][k] * y[12 + k] + dx * Lp[n][k] * yp[12 + k];
            shunt1 += 0.5 * dx * (Gp[n][k] * y[k] + Cp[n][k] * yp[k]);
            shunt2 += 0.5 * dx * (Gp[n][k] * y[6 + k] + Cp[n][k] * yp[6 + k]);
          }
          expected[12 + n] = series;
          expected[3 + n]  = shunt1;
          expected[9 + n]  = shunt2;
        }

        for (IdxT j = 0; j < system_size; ++j)
        {
          success *= isEqual(f[j], expected[static_cast<size_t>(j)], 1.0e-14);
        }

        return success.report(__func__);
      }

      /**
       * @brief Enzyme Jacobian against central finite differences of the
       * assembled residual, J = dF/dy + alpha * dF/dyp.
       *
       * The alpha-scaled shunt capacitance entries exercise the
       * external-derivative Jacobian block.
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

        const RealT step = 1.0e-6;

        auto* y_data  = fixture.y.getData();
        auto* yp_data = fixture.yp.getData();
        auto* f_data  = fixture.f.getData();

        const RealT zero_floor = 1.0e-8;
        for (IdxT j = 0; j < system_size; ++j)
        {
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
              success *= isEqual(enzyme_value, fd_value, 1.0e-6);
            }
          }
        }

        return success.report(__func__);
      }
    };

  } // namespace Testing
} // namespace GridKit
