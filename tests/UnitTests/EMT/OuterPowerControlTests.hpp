/**
 * @file OuterPowerControlTests.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Hand-assembled power and current controller fixture.
 */
#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numbers>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/CommonMath.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/LinearAlgebra/Vector/Vector.hpp>
#include <GridKit/Model/EMT/Component/Controller/InnerCurrentControl/InnerCurrentControl.hpp>
#include <GridKit/Model/EMT/Component/Controller/OuterPowerControl/OuterPowerControl.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <typename ScalarT, typename IdxT>
    class OuterPowerControlTests
    {
      using RealT                        = ScalarT;
      using VectorT                      = LinearAlgebra::Vector<ScalarT, IdxT>;
      using OuterT                       = EMT::Controller::OuterPowerControl<ScalarT, IdxT>;
      using InnerT                       = EMT::Controller::InnerCurrentControl<ScalarT, IdxT>;
      using SignalT                      = EMT::Signal<ScalarT, IdxT>;
      using Data                         = typename OuterT::ModelDataT;
      static constexpr IdxT  system_size = 18;
      static constexpr RealT omega       = RealT{120} * std::numbers::pi_v<RealT>;

      static Data makeData()
      {
        Data data;
        using P                  = typename Data::Parameters;
        data.parameters[P::V]    = RealT{208};
        data.parameters[P::Pref] = RealT{2496};
        data.parameters[P::Qref] = RealT{1040};
        data.parameters[P::Kp]   = RealT{0.3};
        data.parameters[P::Ki]   = RealT{40};
        data.parameters[P::Kaw]  = RealT{200};
        return data;
      }

      static typename InnerT::ModelDataT makeInnerData()
      {
        typename InnerT::ModelDataT data;
        using P                  = typename InnerT::ModelDataT::Parameters;
        data.parameters[P::L]    = RealT{0.002};
        data.parameters[P::Kp]   = RealT{5};
        data.parameters[P::Ki]   = RealT{500};
        data.parameters[P::Kaw]  = RealT{2000};
        data.parameters[P::Imax] = RealT{30};
        data.parameters[P::Mmax] = RealT{0.95};
        return data;
      }

      struct Fixture
      {
        VectorT                y, yp, f, abs_tol;
        std::array<SignalT, 8> input;
        std::array<IdxT, 8>    indices{};
        OuterT                 outer;
        InnerT                 inner;

        Fixture()
          : outer(makeData()), inner(makeInnerData())
        {
          y.resize(system_size);
          yp.resize(system_size);
          f.resize(system_size);
          abs_tol.resize(system_size);
          y.setToConst(0.0);
          yp.setToConst(0.0);
          f.setToConst(0.0);
          for (size_t n = 0; n < input.size(); ++n)
          {
            indices[n] = static_cast<IdxT>(n);
            input[n].set(&y.getData()[n], &yp.getData()[n], &f.getData()[n], &indices[n], &indices[n]);
          }
          using O = typename OuterT::Outputs;
          using I = typename InnerT::Outputs;
          outer.attachInput({&input[4], &input[5], &inner.outputSignal(I::ilimd), &inner.outputSignal(I::ilimq)});
          inner.attachInput({&input[2], &input[3], &input[4], &input[5], &outer.outputSignal(O::icmdd), &outer.outputSignal(O::icmdq), &input[6], &input[7]});
          IdxT offset = 8;
          for (auto* component : components())
          {
            component->bind(y, yp, f, abs_tol, offset);
            component->allocate();
            component->assignGlobalIndices(offset);
            offset += component->size();
          }
          setProbeState();
          outer.initialize();
          inner.initialize();
          for (auto* component : components())
            component->tagDifferentiable();
        }

        std::array<EMT::Component<ScalarT, IdxT>*, 2> components()
        {
          return {&outer, &inner};
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

        void setProbeState()
        {
          const std::array<RealT, system_size> state{12, -5, 208, 13, 8, -3, omega, 400, 0.3, -0.2, 29, 8, 1.3, -0.8, 28, 7, 200, 20};
          for (IdxT n = 0; n < system_size; ++n)
          {
            y.getData()[n]  = state[n];
            yp.getData()[n] = RealT{0.31} * static_cast<RealT>(n);
          }
          y.setDataUpdated();
          yp.setDataUpdated();
        }
      };

    public:
      TestOutcome wiring()
      {
        TestStatus success = true;
        Fixture    fixture;
        success *= fixture.outer.verify() == 0;
        success *= fixture.outer.size() == 4;
        success *= fixture.outer.tag()[0] && fixture.outer.tag()[1]
                   && !fixture.outer.tag()[2] && !fixture.outer.tag()[3];
        success *= &fixture.outer.inputSignal(EMT::Controller::OuterPowerControlInputs::ilimd)
                   == &fixture.inner.outputSignal(InnerT::Outputs::ilimd);
        OuterT invalid;
        success *= invalid.verify() == 5;
        return success.report(__func__);
      }

      TestOutcome initialState()
      {
        TestStatus success = true;
        Fixture    fixture;
        using Outputs  = typename OuterT::Outputs;
        const RealT ed = RealT{12} - 8;
        const RealT eq = RealT{-5} + 3;
        auto*       y  = fixture.y.getData();
        std::cout << "OuterPowerControl default integral (A): " << y[8] << ", " << y[9] << "\n";
        success *= std::abs(y[8]) < 1e-12 && std::abs(y[9]) < 1e-12;
        success *= std::abs(y[10] - RealT{0.3} * ed) < 1e-12;
        success *= std::abs(y[11] - RealT{0.3} * eq) < 1e-12;
        fixture.outer.initialize({{Outputs::icmdd, 2.0}, {Outputs::icmdq, -3.0}});
        success *= y[10] == 2 && y[11] == -3;
        success *= std::abs(y[8] - (2 - RealT{0.3} * ed)) < 1e-12;
        success *= std::abs(y[9] - (-3 - RealT{0.3} * eq)) < 1e-12;
        fixture.evaluateResidual();
        success *= std::abs(fixture.f.getData()[10]) < 1e-12;
        success *= std::abs(fixture.f.getData()[11]) < 1e-12;
        fixture.outer.initializeState({{"icmdd", 4.0}});
        success *= y[10] == 4 && std::abs(y[9]) < 1e-12;
        success *= std::abs(y[11] - RealT{0.3} * eq) < 1e-12;
        success *= y[0] == 12 && y[1] == -5;
        for (size_t n = 8; n < 12; ++n)
          success *= fixture.yp.getData()[n] == 0;
        for (const auto& values : {std::map<std::string, RealT>{{"etad", 0.0}},
                                   std::map<std::string, RealT>{{"etaq", 0.0}},
                                   std::map<std::string, RealT>{{"icmdd", std::numeric_limits<RealT>::infinity()}},
                                   std::map<std::string, RealT>{{"icmdq", std::numeric_limits<RealT>::quiet_NaN()}}})
        {
          bool rejected = false;
          try
          {
            fixture.outer.initializeState(values);
          }
          catch (const std::invalid_argument&)
          {
            rejected = true;
          }
          success *= rejected;
        }
        y[4]          = std::numeric_limits<RealT>::quiet_NaN();
        bool rejected = false;
        try
        {
          fixture.outer.initialize();
        }
        catch (const std::invalid_argument&)
        {
          rejected = true;
        }
        success *= rejected;
        return success.report(__func__);
      }

      TestOutcome residual()
      {
        TestStatus success = true;
        Fixture    fixture;
        fixture.setProbeState();
        RealT maximum_error = 0;
        // Exercise tracking below, near, and above the shared current circle.
        for (const RealT reference : {RealT{8}, RealT{29}, RealT{45}})
        {
          auto* y = fixture.y.getData();
          y[10]   = reference;
          fixture.inner.initialize();
          fixture.evaluateResidual();
          const auto* yp   = fixture.yp.getData();
          const auto* f    = fixture.f.getData();
          const RealT mu   = Math::MU<RealT>;
          const auto  ramp = [mu](RealT x)
          {
            return std::max(x, RealT{0}) + std::log1p(std::exp(-mu * std::abs(x))) / mu;
          };
          const std::array<RealT, 2> e{RealT{12} - y[4], RealT{-5} - y[5]};
          const RealT                limiter = std::sqrt(1 + ramp((y[10] * y[10] + y[11] * y[11]) / 900 - 1));
          for (size_t n = 0; n < 2; ++n)
          {
            const RealT expected = -yp[8 + n] + 40 * e[n] + 200 * (y[10 + n] / limiter - y[10 + n]);
            maximum_error        = std::max(maximum_error, std::abs(f[8 + n] - expected));
            maximum_error        = std::max(maximum_error, std::abs(f[10 + n] - (y[10 + n] - 0.3 * e[n] - y[8 + n])));
          }
        }
        std::cout << "OuterPowerControl residual maximum absolute error: " << maximum_error << "\n";
        success *= maximum_error < 1e-12;
        return success.report(__func__);
      }

      TestOutcome steadyState()
      {
        TestStatus success = true;
        Fixture    fixture;
        auto*      y = fixture.y.getData();
        fixture.yp.setToConst(0);
        const RealT vd = 208, igd = 12, igq = -5;
        y[0]  = igd;
        y[1]  = igq;
        y[2]  = vd;
        y[3]  = 0;
        y[4]  = igd;
        y[5]  = igq;
        y[8]  = igd;
        y[9]  = igq + omega * 0.0001 * vd;
        y[10] = igd;
        y[11] = igq + omega * 0.0001 * vd;
        fixture.inner.initialize();
        fixture.evaluateResidual();
        for (size_t n = 8; n < 12; ++n)
          success *= std::abs(fixture.f.getData()[n]) < 1e-10;
        // At identical current error, saturation supplies the restoring tracking term.
        y[10] = 45;
        y[11] = 0;
        fixture.inner.initialize();
        fixture.evaluateResidual();
        success *= fixture.f.getData()[8] < -2900;
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
          const RealT step = 1.0e-6 * (1.0 + std::abs(y_data[j]));

          std::array<RealT, system_size> fd_column{};

          const RealT y_saved = y_data[j];
          y_data[j]           = y_saved + step;
          fixture.evaluateResidual();
          for (IdxT i = 8; i < system_size; ++i)
          {
            fd_column[static_cast<size_t>(i)] = f_data[i];
          }
          y_data[j] = y_saved - step;
          fixture.evaluateResidual();
          for (IdxT i = 8; i < system_size; ++i)
          {
            fd_column[static_cast<size_t>(i)] = (fd_column[static_cast<size_t>(i)] - f_data[i]) / (2.0 * step);
          }
          y_data[j] = y_saved;

          const RealT yp_saved = yp_data[j];
          yp_data[j]           = yp_saved + step;
          fixture.evaluateResidual();
          std::array<RealT, system_size> fp_plus{};
          for (IdxT i = 8; i < system_size; ++i)
          {
            fp_plus[static_cast<size_t>(i)] = f_data[i];
          }
          yp_data[j] = yp_saved - step;
          fixture.evaluateResidual();
          for (IdxT i = 8; i < system_size; ++i)
          {
            fd_column[static_cast<size_t>(i)] += alpha * (fp_plus[static_cast<size_t>(i)] - f_data[i]) / (2.0 * step);
          }
          yp_data[j] = yp_saved;

          for (IdxT i = 8; i < system_size; ++i)
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
            success           *= (error < 3.0e-7);
          }
        }

        std::cout << "OuterPowerControl Jacobian maximum scaled error: " << maximum_error << "\n";
        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
