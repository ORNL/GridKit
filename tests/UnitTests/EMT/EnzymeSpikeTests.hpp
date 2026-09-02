/**
 * @file EnzymeSpikeTests.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Hand-assembled one-bus EMT system tests for the core coupling and
 * Enzyme Jacobian machinery.
 *
 */
#pragma once

#include <cmath>
#include <map>
#include <numbers>
#include <utility>
#include <vector>

#include <GridKit/Definitions.hpp>
#include <GridKit/LinearAlgebra/Vector/Vector.hpp>
#include <GridKit/Model/EMT/Bus/Bus.hpp>
#include <GridKit/Model/EMT/Component/Load/LoadZ/LoadZ.hpp>
#include <GridKit/Model/EMT/Component/Source/VoltageSource/VoltageSource.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    /**
     * @brief Tests for a hand-assembled EMT system of one bus, one voltage
     * source, and one impedance load.
     */
    template <typename ScalarT, typename IdxT>
    class EnzymeSpikeTests
    {
    public:
      using RealT   = ScalarT;
      using VectorT = GridKit::LinearAlgebra::Vector<ScalarT, IdxT>;
      using BusT    = GridKit::EMT::Bus<ScalarT, IdxT>;
      using SourceT = GridKit::EMT::VoltageSource<ScalarT, IdxT>;
      using LoadT   = GridKit::EMT::LoadZ<ScalarT, IdxT>;

      static constexpr IdxT system_size = 12;

      EnzymeSpikeTests()  = default;
      ~EnzymeSpikeTests() = default;

    private:
      /**
       * @brief One bus with a voltage source and an impedance load attached.
       *
       * Variable layout: bus voltages [0, 3), source variables [3, 9),
       * load variables [9, 12).
       */
      struct Fixture
      {
        VectorT y;
        VectorT yp;
        VectorT f;
        VectorT abs_tol;

        BusT    bus;
        SourceT source;
        LoadT   load;

        Fixture()
          : source(makeSourceData()),
            load(makeLoadData())
        {
          y.resize(system_size);
          yp.resize(system_size);
          f.resize(system_size);
          abs_tol.resize(system_size);

          source.getSignals().template attachPort<GridKit::EMT::VoltageSourceExternalVariables::VA>(&bus.voltagePort());
          load.getSignals().template attachPort<GridKit::EMT::LoadZExternalVariables::VA>(&bus.voltagePort());

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
          return {&bus, &source, &load};
        }

        void updateTime(RealT t, RealT alpha)
        {
          for (auto* component : components())
          {
            component->updateTime(t, alpha);
          }
        }

        /// Evaluate the assembled residual: the internal phase over every
        /// component, then the external accumulation phase.
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

        /// Fill the state with distinct values away from 0 and 1.
        void setProbeState()
        {
          auto* y_data  = y.getData();
          auto* yp_data = yp.getData();
          for (IdxT j = 0; j < system_size; ++j)
          {
            y_data[j]  = 3.7 + 1.3 * static_cast<RealT>(j);
            yp_data[j] = -2.9 + 0.7 * static_cast<RealT>(j);
          }
          y.setDataUpdated();
          yp.setDataUpdated();
        }
      };

      static GridKit::EMT::VoltageSourceData<ScalarT, IdxT> makeSourceData()
      {
        using Data      = GridKit::EMT::VoltageSourceData<ScalarT, IdxT>;
        using Parameter = typename Data::Parameters;
        Data data;
        data.parameters[Parameter::E]     = GridKit::EMT::ABCVector<RealT>{{100.0, 90.0, 80.0}};
        data.parameters[Parameter::phi]   = GridKit::EMT::ABCVector<RealT>{{0.1, -2.0, 2.2}};
        data.parameters[Parameter::omega] = RealT{2.0 * std::numbers::pi_v<RealT> * 60.0};
        data.parameters[Parameter::Rs]    = GridKit::EMT::ABCMatrix<RealT>{{{{1.1, 0.11, 0.12}},
                                                                            {{0.13, 1.2, 0.14}},
                                                                            {{0.15, 0.16, 1.3}}}};
        data.parameters[Parameter::Ls]    = GridKit::EMT::ABCMatrix<RealT>{{{{0.021, 0.002, 0.003}},
                                                                            {{0.004, 0.022, 0.005}},
                                                                            {{0.006, 0.007, 0.023}}}};
        return data;
      }

      static GridKit::EMT::LoadZData<ScalarT, IdxT> makeLoadData()
      {
        using Data      = GridKit::EMT::LoadZData<ScalarT, IdxT>;
        using Parameter = typename Data::Parameters;
        Data data;
        data.parameters[Parameter::R] = GridKit::EMT::ABCMatrix<RealT>{{{{5.1, 0.21, 0.22}},
                                                                        {{0.23, 5.2, 0.24}},
                                                                        {{0.25, 0.26, 5.3}}}};
        data.parameters[Parameter::L] = GridKit::EMT::ABCMatrix<RealT>{{{{0.31, 0.011, 0.012}},
                                                                        {{0.013, 0.32, 0.014}},
                                                                        {{0.015, 0.016, 0.33}}}};
        return data;
      }

    public:
      /**
       * @brief Wiring smoke test: ports link, sizes add up, tags propagate.
       */
      TestOutcome wiring()
      {
        TestStatus success = true;

        Fixture fixture;

        success *= (fixture.bus.size() == 3);
        success *= (fixture.source.size() == 6);
        success *= (fixture.load.size() == 3);

        success *= (fixture.bus.verify() == 0);
        success *= (fixture.source.verify() == 0);
        success *= (fixture.load.verify() == 0);

        // No connected component reads the bus voltage derivative, so the
        // bus voltage is algebraic.
        success *= (fixture.bus.tag()[0] == false);
        success *= (fixture.bus.tag()[1] == false);
        success *= (fixture.bus.tag()[2] == false);

        success *= (fixture.source.tag()[0] == false);
        success *= (fixture.source.tag()[3] == true);
        success *= (fixture.load.tag()[0] == true);

        return success.report(__func__);
      }

      /**
       * @brief Assembled residual against an independent computation.
       */
      TestOutcome residual()
      {
        TestStatus success = true;

        const RealT t     = 0.017;
        const RealT alpha = 3.7;

        Fixture fixture;
        fixture.setProbeState();
        fixture.updateTime(t, alpha);
        fixture.evaluateResidual();

        const auto* y  = fixture.y.getData();
        const auto* yp = fixture.yp.getData();
        const auto* f  = fixture.f.getData();

        const auto source_data = makeSourceData();
        const auto load_data   = makeLoadData();

        using SourceParameter = typename GridKit::EMT::VoltageSourceData<ScalarT, IdxT>::Parameters;
        using LoadParameter   = typename GridKit::EMT::LoadZData<ScalarT, IdxT>::Parameters;

        const auto E     = std::get<GridKit::EMT::ABCVector<RealT>>(source_data.parameters.at(SourceParameter::E));
        const auto phi   = std::get<GridKit::EMT::ABCVector<RealT>>(source_data.parameters.at(SourceParameter::phi));
        const auto omega = std::get<RealT>(source_data.parameters.at(SourceParameter::omega));
        const auto Rs    = std::get<GridKit::EMT::ABCMatrix<RealT>>(source_data.parameters.at(SourceParameter::Rs));
        const auto Ls    = std::get<GridKit::EMT::ABCMatrix<RealT>>(source_data.parameters.at(SourceParameter::Ls));
        const auto R     = std::get<GridKit::EMT::ABCMatrix<RealT>>(load_data.parameters.at(LoadParameter::R));
        const auto L     = std::get<GridKit::EMT::ABCMatrix<RealT>>(load_data.parameters.at(LoadParameter::L));

        const RealT sqrt2 = std::numbers::sqrt2_v<RealT>;

        std::array<RealT, system_size> expected{};

        // Bus current balance: source and load injections
        for (size_t n = 0; n < 3; ++n)
        {
          expected[n] = y[6 + n] + y[9 + n];
        }

        // Source voltage rows
        for (size_t n = 0; n < 3; ++n)
        {
          expected[3 + n] = y[3 + n] - sqrt2 * E[n] * std::cos(omega * t + phi[n]);
        }

        // Source series branch rows
        for (size_t n = 0; n < 3; ++n)
        {
          RealT row = y[n] - y[3 + n];
          for (size_t k = 0; k < 3; ++k)
          {
            row += Rs[n][k] * y[6 + k] + Ls[n][k] * yp[6 + k];
          }
          expected[6 + n] = row;
        }

        // Load branch rows
        for (size_t n = 0; n < 3; ++n)
        {
          RealT row = y[n];
          for (size_t k = 0; k < 3; ++k)
          {
            row += R[n][k] * y[9 + k] + L[n][k] * yp[9 + k];
          }
          expected[9 + n] = row;
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
       */
      TestOutcome jacobian()
      {
        TestStatus success = true;

        const RealT t     = 0.017;
        const RealT alpha = 3.7;

        Fixture fixture;
        fixture.setProbeState();
        fixture.updateTime(t, alpha);
        fixture.evaluateResidual();

        for (auto* component : fixture.components())
        {
          component->evaluateJacobian();
        }

        // Collect Jacobian triplets, with duplicates accumulated
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

        // Central finite differences of the assembled residual
        const RealT                                             step = 1.0e-6;
        std::array<std::array<RealT, system_size>, system_size> fd{};

        auto* y_data  = fixture.y.getData();
        auto* yp_data = fixture.yp.getData();
        auto* f_data  = fixture.f.getData();

        for (IdxT j = 0; j < system_size; ++j)
        {
          std::array<RealT, system_size> f_plus{};
          std::array<RealT, system_size> f_minus{};

          const RealT y_saved = y_data[j];
          y_data[j]           = y_saved + step;
          fixture.evaluateResidual();
          for (IdxT i = 0; i < system_size; ++i)
          {
            f_plus[static_cast<size_t>(i)] = f_data[i];
          }
          y_data[j] = y_saved - step;
          fixture.evaluateResidual();
          for (IdxT i = 0; i < system_size; ++i)
          {
            f_minus[static_cast<size_t>(i)] = f_data[i];
          }
          y_data[j] = y_saved;

          for (IdxT i = 0; i < system_size; ++i)
          {
            fd[static_cast<size_t>(i)][static_cast<size_t>(j)] = (f_plus[static_cast<size_t>(i)] - f_minus[static_cast<size_t>(i)]) / (2.0 * step);
          }

          const RealT yp_saved = yp_data[j];
          yp_data[j]           = yp_saved + step;
          fixture.evaluateResidual();
          for (IdxT i = 0; i < system_size; ++i)
          {
            f_plus[static_cast<size_t>(i)] = f_data[i];
          }
          yp_data[j] = yp_saved - step;
          fixture.evaluateResidual();
          for (IdxT i = 0; i < system_size; ++i)
          {
            f_minus[static_cast<size_t>(i)] = f_data[i];
          }
          yp_data[j] = yp_saved;

          for (IdxT i = 0; i < system_size; ++i)
          {
            fd[static_cast<size_t>(i)][static_cast<size_t>(j)] += alpha * (f_plus[static_cast<size_t>(i)] - f_minus[static_cast<size_t>(i)]) / (2.0 * step);
          }
        }

        // Every finite-difference nonzero must be matched by an Enzyme entry
        const RealT zero_floor = 1.0e-8;
        for (IdxT i = 0; i < system_size; ++i)
        {
          for (IdxT j = 0; j < system_size; ++j)
          {
            const RealT fd_value     = fd[static_cast<size_t>(i)][static_cast<size_t>(j)];
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
