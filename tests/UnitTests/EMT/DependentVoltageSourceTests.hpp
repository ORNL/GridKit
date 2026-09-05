/**
 * @file DependentVoltageSourceTests.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief EMT dependent-voltage-source tests.
 */
#pragma once

#include <array>
#include <cmath>
#include <complex>
#include <map>
#include <sstream>
#include <utility>

#include <GridKit/Definitions.hpp>
#include <GridKit/LinearAlgebra/Vector/Vector.hpp>
#include <GridKit/Model/EMT/Component/Bus/Bus.hpp>
#include <GridKit/Model/EMT/Component/Source/DependentVoltageSource/DependentVoltageSource.hpp>
#include <GridKit/Model/EMT/SystemModel.hpp>
#include <GridKit/Model/EMT/SystemModelData.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <typename ScalarT, typename IdxT>
    class DependentVoltageSourceTests
    {
    public:
      using RealT     = ScalarT;
      using VectorT   = LinearAlgebra::Vector<ScalarT, IdxT>;
      using BusT      = EMT::Bus<ScalarT, IdxT>;
      using SourceT   = EMT::DependentVoltageSource<ScalarT, IdxT>;
      using DataT     = EMT::DependentVoltageSourceData<RealT, IdxT>;
      using Component = EMT::Component<ScalarT, IdxT>;

    private:
      static DataT seriesData()
      {
        using Parameter = typename DataT::Parameters;

        DataT data;
        data.parameters[Parameter::N]  = IdxT{3};
        data.parameters[Parameter::Rs] = EMT::ABCMatrix<RealT>{{{{1.1, 0.12, 0.13}},
                                                                {{0.21, 1.2, 0.23}},
                                                                {{0.31, 0.32, 1.3}}}};
        data.parameters[Parameter::Ls] = EMT::ABCMatrix<RealT>{{{{0.041, 0.002, 0.003}},
                                                                {{0.004, 0.052, 0.006}},
                                                                {{0.007, 0.008, 0.063}}}};
        return data;
      }

      static DataT rationalData()
      {
        using ComplexT  = std::complex<RealT>;
        using Parameter = typename DataT::Parameters;

        DataT data;
        data.parameters[Parameter::N] = IdxT{3};

        EMT::VectorFitData<RealT, IdxT> fit;
        fit.D = {{{{1.7, 0.21, 0.32}}, {{0.13, 1.9, 0.24}}, {{0.35, 0.16, 2.1}}}};
        fit.poles.emplace_back(-3.1, 0.0);
        fit.residues.push_back(EMT::ABCMatrix<ComplexT>{{{{ComplexT{0.61}, ComplexT{0.041}, ComplexT{0.022}}},
                                                         {{ComplexT{0.043}, ComplexT{0.72}, ComplexT{0.054}}},
                                                         {{ComplexT{0.025}, ComplexT{0.046}, ComplexT{0.83}}}}});
        data.Y = fit;
        return data;
      }

      struct Fixture
      {
        VectorT y;
        VectorT yp;
        VectorT f;
        VectorT abs_tol;

        BusT    drive;
        BusT    terminal;
        SourceT source;

        explicit Fixture(const DataT& data)
          : source(data)
        {
          const IdxT system_size = drive.size() + terminal.size() + source.size();
          y.resize(system_size);
          yp.resize(system_size);
          f.resize(system_size);
          abs_tol.resize(system_size);

          source.getSignals().template attachPort<EMT::DependentVoltageSourceExternalVariables::VA>(
              &terminal.voltagePort());
          source.getSignals().template attachPort<EMT::DependentVoltageSourceExternalVariables::EA>(
              &drive.voltagePort());

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

        std::array<Component*, 3> components()
        {
          return {&drive, &terminal, &source};
        }

        void setProbeState()
        {
          auto* y_data  = y.getData();
          auto* yp_data = yp.getData();
          for (IdxT j = 0; j < y.getSize(); ++j)
          {
            y_data[j]  = 0.7 + 0.31 * static_cast<RealT>(j);
            yp_data[j] = -1.3 + 0.23 * static_cast<RealT>(j);
          }
          y.setDataUpdated();
          yp.setDataUpdated();
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
      };

    public:
      TestOutcome seriesResidual()
      {
        TestStatus success = true;

        const auto data = seriesData();
        Fixture    fixture(data);
        fixture.setProbeState();
        fixture.evaluateResidual();

        success *= (fixture.source.verify() == 0);
        success *= (fixture.source.size() == 3);

        const auto* y  = fixture.y.getData();
        const auto* yp = fixture.yp.getData();
        const auto* f  = fixture.f.getData();
        const auto  Rs = std::get<EMT::ABCMatrix<RealT>>(
            data.parameters.at(DataT::Parameters::Rs));
        const auto Ls = std::get<EMT::ABCMatrix<RealT>>(
            data.parameters.at(DataT::Parameters::Ls));

        for (size_t n = 0; n < 3; ++n)
        {
          success *= (fixture.source.tag()[n] == true);
          success *= isEqual(f[n], 0.0, 1.0e-14);
          success *= isEqual(f[3 + n], y[6 + n], 1.0e-14);

          RealT expected = y[3 + n] - y[n];
          for (size_t k = 0; k < 3; ++k)
          {
            expected += Rs[n][k] * y[6 + k] + Ls[n][k] * yp[6 + k];
          }
          success *= isEqual(f[6 + n], expected, 1.0e-14);
        }

        return success.report(__func__);
      }

      TestOutcome rationalResidual()
      {
        TestStatus success = true;

        const auto data = rationalData();
        Fixture    fixture(data);
        fixture.setProbeState();
        fixture.evaluateResidual();

        success *= (fixture.source.verify() == 0);
        success *= (fixture.source.size() == 6);

        const auto* y   = fixture.y.getData();
        const auto* yp  = fixture.yp.getData();
        const auto* f   = fixture.f.getData();
        const auto& fit = *data.Y;

        for (size_t n = 0; n < 3; ++n)
        {
          success *= (fixture.source.tag()[n] == false);
          success *= (fixture.source.tag()[3 + n] == true);
          success *= isEqual(f[n], 0.0, 1.0e-14);
          success *= isEqual(f[6 + n], y[6 + n] + y[3 + n] - y[n], 1.0e-14);
          success *= isEqual(f[9 + n],
                             -yp[9 + n] + fit.poles[0].real() * y[9 + n] + y[6 + n],
                             1.0e-14);

          RealT current = 0.0;
          for (size_t k = 0; k < 3; ++k)
          {
            current += fit.D[n][k] * y[6 + k]
                       + fit.residues[0][n][k].real() * y[9 + k];
          }
          success *= isEqual(f[3 + n], current, 1.0e-14);
        }

        constexpr RealT tolerance = 1.0e-6;
        fixture.source.setAbsoluteTolerance(tolerance);
        const auto* abs_tol = fixture.abs_tol.getData();
        for (size_t n = 0; n < 3; ++n)
        {
          success *= isEqual(abs_tol[6 + n], tolerance, 1.0e-18);
          success *= isEqual(abs_tol[9 + n],
                             tolerance / std::abs(fit.poles[0].real()),
                             1.0e-18);
        }

        return success.report(__func__);
      }

      TestOutcome validation()
      {
        TestStatus success = true;

        auto wrong_n                             = seriesData();
        wrong_n.parameters[DataT::Parameters::N] = IdxT{2};
        Fixture wrong_n_fixture(wrong_n);

        using Log                     = Utilities::Logger;
        const auto previous_verbosity = Log::verbosity();
        Log::setVerbosity(Log::Verbosity::NONE);
        success *= (wrong_n_fixture.source.verify() > 0);
        Log::setVerbosity(previous_verbosity);

        auto resistive                              = seriesData();
        resistive.parameters[DataT::Parameters::Ls] = EMT::ABCMatrix<RealT>{};
        Fixture resistive_fixture(resistive);
        success *= (resistive_fixture.source.verify() == 0);
        for (size_t n = 0; n < 3; ++n)
        {
          success *= (resistive_fixture.source.tag()[n] == false);
        }

        auto partial_ls                                       = EMT::ABCMatrix<RealT>{};
        partial_ls[0][1]                                      = 0.1;
        auto partially_inductive                              = seriesData();
        partially_inductive.parameters[DataT::Parameters::Ls] = partial_ls;
        Fixture partially_inductive_fixture(partially_inductive);
        success *= (partially_inductive_fixture.source.tag()[0] == false);
        success *= (partially_inductive_fixture.source.tag()[1] == true);
        success *= (partially_inductive_fixture.source.tag()[2] == false);

        return success.report(__func__);
      }

      TestOutcome parseAndBuild()
      {
        TestStatus success = true;

        std::istringstream stream(R"({
          "header": {
            "case_name": "Dependent source smoke case",
            "case_description": "Parser and wiring coverage",
            "case_comments": ""
          },
          "signals": [
            { "id": "ea" },
            { "id": "eb" },
            { "id": "ec" }
          ],
          "devices": [
            { "class": "Bus", "id": "bus_1" },
            {
              "class": "DependentVoltageSource",
              "id": "source_1",
              "params": {
                "N": 3,
                "Rs": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                "Ls": [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.1]]
              },
              "inputs": {
                "bus": "bus_1",
                "ea": "ea",
                "eb": "eb",
                "ec": "ec"
              }
            }
          ]
        })");

        const auto data  = EMT::parseSystemModelData(stream);
        success         *= (data.dependent_voltage_source.size() == 1);

        using Input         = EMT::DependentVoltageSourceInputs;
        const auto& inputs  = data.dependent_voltage_source[0].inputs;
        success            *= (data.bus.size() == 1);
        success            *= (data.bus[0].id == "bus_1");
        success            *= (data.signal.size() == 3);
        success            *= (data.signal[0].id == "ea");
        success            *= (data.signal[1].id == "eb");
        success            *= (data.signal[2].id == "ec");
        success            *= (inputs.at(Input::bus) == "bus_1");
        success            *= (inputs.at(Input::ea) == "ea");
        success            *= (inputs.at(Input::eb) == "eb");
        success            *= (inputs.at(Input::ec) == "ec");

        EMT::SystemModel<ScalarT, IdxT> system(data);
        auto*                           source  = dynamic_cast<SourceT*>(system.getComponent(1));
        success                                *= (source != nullptr);
        if (source != nullptr)
        {
          const auto& signals  = source->getSignals();
          success             *= signals.template isAttached<EMT::DependentVoltageSourceExternalVariables::VA>();
          success             *= signals.template isAttached<EMT::DependentVoltageSourceExternalVariables::EA>();
          success             *= signals.template isAttached<EMT::DependentVoltageSourceExternalVariables::EB>();
          success             *= signals.template isAttached<EMT::DependentVoltageSourceExternalVariables::EC>();
        }

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      TestOutcome jacobian()
      {
        TestStatus success = true;

        constexpr RealT alpha = 3.7;
        const auto      data  = seriesData();
        Fixture         fixture(data);
        fixture.setProbeState();
        fixture.source.updateTime(0.0, alpha);
        fixture.source.evaluateJacobian();

        std::map<std::pair<IdxT, IdxT>, RealT> entries;
        auto*                                  coo  = fixture.source.getCooJacobian();
        success                                    *= (coo != nullptr);
        if (coo != nullptr)
        {
          for (IdxT j = 0; j < coo->getNnz(); ++j)
          {
            entries[{coo->getRowData()[j], coo->getColData()[j]}] += coo->getValues()[j];
          }
        }

        const auto Rs = std::get<EMT::ABCMatrix<RealT>>(
            data.parameters.at(DataT::Parameters::Rs));
        const auto Ls = std::get<EMT::ABCMatrix<RealT>>(
            data.parameters.at(DataT::Parameters::Ls));
        for (IdxT n = 0; n < 3; ++n)
        {
          success *= isEqual(entries[{3 + n, 6 + n}], 1.0, 1.0e-14);
          success *= isEqual(entries[{6 + n, n}], -1.0, 1.0e-14);
          success *= isEqual(entries[{6 + n, 3 + n}], 1.0, 1.0e-14);
          for (IdxT k = 0; k < 3; ++k)
          {
            success *= isEqual(entries[{6 + n, 6 + k}], Rs[n][k] + alpha * Ls[n][k], 1.0e-14);
          }
        }

        return success.report(__func__);
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
