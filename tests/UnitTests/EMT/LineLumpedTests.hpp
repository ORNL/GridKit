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
#include <GridKit/Model/EMT/Bus/Bus.hpp>
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
       * Variable layout: bus 1 voltages [0, 3), bus 2 voltages [3, 6),
       * line variables [6, 15).
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

          line.getSignals().template attachPort<GridKit::EMT::LineLumpedExternalVariables::V1A>(&bus1.voltagePort());
          line.getSignals().template attachPort<GridKit::EMT::LineLumpedExternalVariables::V2A>(&bus2.voltagePort());

          IdxT offset = 0;
          for (auto* component : components())
          {
            component->bind(y, yp, f, abs_tol, offset);
            component->allocate();
            for (IdxT j = 0; j < component->size(); ++j)
            {
              component->setVariableIndex(j, offset + j);
              component->setResidualIndex(j, offset + j);
            }
            offset += component->size();
          }

          for (auto* component : components())
          {
            component->initialize();
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
      /**
       * @brief Wiring smoke test: shunt capacitance makes the bus voltages
       * differential.
       */
      TestOutcome wiring()
      {
        TestStatus success = true;

        Fixture fixture;

        success *= (fixture.line.size() == 9);
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
        success *= (fixture.line.tag()[3] == false);
        success *= (fixture.line.tag()[6] == false);

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

        // Bus current balances: line injections
        for (size_t n = 0; n < 3; ++n)
        {
          expected[n]     = y[9 + n] - y[6 + n];
          expected[3 + n] = y[12 + n] + y[6 + n];
        }

        // Series rows
        for (size_t n = 0; n < 3; ++n)
        {
          RealT row = y[3 + n] - y[n];
          for (size_t k = 0; k < 3; ++k)
          {
            row += dx * Rp[n][k] * y[6 + k] + dx * Lp[n][k] * yp[6 + k];
          }
          expected[6 + n] = row;
        }

        // Shunt rows
        for (size_t n = 0; n < 3; ++n)
        {
          RealT row1 = 2.0 * y[9 + n];
          RealT row2 = 2.0 * y[12 + n];
          for (size_t k = 0; k < 3; ++k)
          {
            row1 += dx * Gp[n][k] * y[k] + dx * Cp[n][k] * yp[k];
            row2 += dx * Gp[n][k] * y[3 + k] + dx * Cp[n][k] * yp[3 + k];
          }
          expected[9 + n]  = row1;
          expected[12 + n] = row2;
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
