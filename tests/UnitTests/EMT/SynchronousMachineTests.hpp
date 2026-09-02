/**
 * @file SynchronousMachineTests.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief SynchronousMachine model tests on a hand-assembled one-bus system.
 *
 */
#pragma once

#include <cmath>
#include <map>
#include <utility>

#include <GridKit/Definitions.hpp>
#include <GridKit/LinearAlgebra/Vector/Vector.hpp>
#include <GridKit/Model/EMT/Bus/Bus.hpp>
#include <GridKit/Model/EMT/Component/Machine/SynchronousMachine/SynchronousMachine.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    /**
     * @brief Tests for the synchronous machine on a hand-assembled bus.
     *
     * The machine is the strongly nonlinear EMT model, so the Jacobian test
     * exercises the trigonometric, square-root, and saturation derivatives.
     */
    template <typename ScalarT, typename IdxT>
    class SynchronousMachineTests
    {
    public:
      using RealT    = ScalarT;
      using VectorT  = GridKit::LinearAlgebra::Vector<ScalarT, IdxT>;
      using BusT     = GridKit::EMT::Bus<ScalarT, IdxT>;
      using MachineT = GridKit::EMT::SynchronousMachine<ScalarT, IdxT>;

      static constexpr IdxT system_size = 27;

      SynchronousMachineTests()  = default;
      ~SynchronousMachineTests() = default;

    private:
      static GridKit::EMT::SynchronousMachineData<ScalarT, IdxT> makeResidualData()
      {
        using Data      = GridKit::EMT::SynchronousMachineData<ScalarT, IdxT>;
        using Parameter = typename Data::Parameters;
        Data data;
        data.parameters[Parameter::S]    = RealT{100.0e6};
        data.parameters[Parameter::V]    = RealT{13800.0};
        data.parameters[Parameter::f]    = RealT{60.0};
        data.parameters[Parameter::H]    = RealT{3.7};
        data.parameters[Parameter::F]    = RealT{0.013};
        data.parameters[Parameter::Rs]   = RealT{0.003};
        data.parameters[Parameter::Ll]   = RealT{0.15};
        data.parameters[Parameter::Lmd]  = RealT{1.66};
        data.parameters[Parameter::Lmq]  = RealT{1.61};
        data.parameters[Parameter::L0]   = RealT{0.17};
        data.parameters[Parameter::Rfd]  = RealT{0.0006};
        data.parameters[Parameter::Llfd] = RealT{0.165};
        data.parameters[Parameter::R1d]  = RealT{0.0284};
        data.parameters[Parameter::Ll1d] = RealT{0.1713};
        data.parameters[Parameter::R1q]  = RealT{0.0062};
        data.parameters[Parameter::Ll1q] = RealT{0.7252};
        data.parameters[Parameter::R2q]  = RealT{0.0237};
        data.parameters[Parameter::Ll2q] = RealT{0.125};
        data.parameters[Parameter::S10]  = RealT{0.1};
        data.parameters[Parameter::S12]  = RealT{0.5};
        data.parameters[Parameter::p0]   = RealT{50.0e6};
        data.parameters[Parameter::q0]   = RealT{10.0e6};
        return data;
      }

      /**
       * @brief One bus with a synchronous machine attached.
       *
       * Variable layout: bus voltages [0, 3), machine variables [3, 27).
       */
      struct Fixture
      {
        VectorT y;
        VectorT yp;
        VectorT f;
        VectorT abs_tol;

        BusT     bus;
        MachineT machine;

        Fixture()
          : bus(ScalarT{11267.65281680262},
                ScalarT{-5633.82640840131},
                ScalarT{-5633.82640840131}),
            machine(makeResidualData())
        {
          y.resize(system_size);
          yp.resize(system_size);
          f.resize(system_size);
          abs_tol.resize(system_size);

          machine.getSignals().template attachPort<GridKit::EMT::SynchronousMachineExternalVariables::VA>(&bus.voltagePort());

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

        std::array<GridKit::EMT::Component<ScalarT, IdxT>*, 2> components()
        {
          return {&bus, &machine};
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

        /// Perturb the initialized operating point so no residual row or
        /// Jacobian entry sits on a special value.
        void setProbeState()
        {
          auto* y_data  = y.getData();
          auto* yp_data = yp.getData();
          for (IdxT j = 0; j < system_size; ++j)
          {
            y_data[j]  *= 1.0 + 0.01 * static_cast<RealT>(j + 1) / static_cast<RealT>(system_size);
            yp_data[j] += 0.03 * static_cast<RealT>(j) - 0.2;
          }
          y.setDataUpdated();
          yp.setDataUpdated();
        }
      };

    public:
      /**
       * @brief Initialization satisfies the assembled residual exactly.
       */
      TestOutcome initialization()
      {
        TestStatus success = true;

        Fixture fixture;

        success *= (fixture.machine.size() == 24);
        success *= (fixture.bus.verify() == 0);
        success *= (fixture.machine.verify() == 0);

        fixture.updateTime(0.0, 1.0);
        fixture.evaluateResidual();

        // The machine rows are exactly satisfied; the bus current-balance
        // rows hold the uncancelled injections because the machine is the
        // only device.
        const auto* f             = fixture.f.getData();
        RealT       residual_norm = 0.0;
        for (IdxT j = 3; j < system_size; ++j)
        {
          residual_norm += f[j] * f[j];
        }
        residual_norm  = std::sqrt(residual_norm);
        success       *= (residual_norm < 1.0e-12);

        success *= (std::abs(f[0]) > 1.0e2);

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

        // The state spans radians to hundreds of kilovolts, so the step and
        // the comparison floor scale with each column.
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
