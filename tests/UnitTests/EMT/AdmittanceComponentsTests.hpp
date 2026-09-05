#pragma once

#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <map>
#include <type_traits>
#include <vector>

#include <GridKit/Model/EMT/Component/Bus/Bus.hpp>
#include <GridKit/Model/EMT/Component/Load/LoadZ/LoadZ.hpp>
#include <GridKit/Model/EMT/Component/Source/DependentVoltageSource/DependentVoltageSource.hpp>
#include <GridKit/Model/EMT/Component/Source/VoltageSource/VoltageSource.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit::Testing
{
  template <typename ScalarT, typename IdxT>
  class AdmittanceComponentsTests
  {
    using RealT      = ScalarT;
    using ComplexT   = std::complex<RealT>;
    using ComponentT = EMT::Component<ScalarT, IdxT>;
    using BusT       = EMT::Bus<ScalarT, IdxT>;
    using SourceT    = EMT::VoltageSource<ScalarT, IdxT>;
    using DependentT = EMT::DependentVoltageSource<ScalarT, IdxT>;
    using LoadT      = EMT::LoadZ<ScalarT, IdxT>;
    using FitDataT   = EMT::VectorFitData<RealT, IdxT>;
    using Phasors    = EMT::ABCVector<ComplexT>;

    static FitDataT fitData()
    {
      FitDataT fit;
      fit.poles = {{-2.7, 0.0}, {-1.3, 4.1}, {-1.3, -4.1}};
      fit.residues.resize(3);
      for (size_t n = 0; n < 3; ++n)
      {
        for (size_t k = 0; k < 3; ++k)
        {
          fit.D[n][k]           = n == k ? 1.7 + 0.1 * static_cast<RealT>(n) : 0.05 * static_cast<RealT>(1 + n + k);
          fit.residues[0][n][k] = n == k ? 0.9 + 0.2 * static_cast<RealT>(n) : 0.08 * static_cast<RealT>(1 + n + k);
          fit.residues[1][n][k] = {0.07 * static_cast<RealT>(1 + n + k), 0.04 * static_cast<RealT>(1 + 2 * n + k)};
          fit.residues[2][n][k] = std::conj(fit.residues[1][n][k]);
        }
      }
      return fit;
    }

    // Independent frequency-domain oracle evaluated directly from the input data.
    static Phasors applyTransfer(const FitDataT& fit, RealT omega, const Phasors& input)
    {
      Phasors        output{};
      const ComplexT s{0.0, omega};
      for (size_t n = 0; n < 3; ++n)
      {
        for (size_t k = 0; k < 3; ++k)
        {
          ComplexT value = fit.D[n][k] + s * fit.E[n][k];
          for (size_t q = 0; q < fit.poles.size(); ++q)
          {
            value += fit.residues[q][n][k] / (s - fit.poles[q]);
          }
          output[n] += value * input[k];
        }
      }
      return output;
    }

    template <typename ModelT>
    struct Fixture
    {
      typename ComponentT::VectorT y, yp, f, tolerance;
      BusT                         drive, terminal;
      ModelT                       model;

      explicit Fixture(const typename ModelT::ModelDataT& data)
        : model(data)
      {
        const auto size = drive.size() + terminal.size() + model.size();
        y.resize(size);
        yp.resize(size);
        f.resize(size);
        tolerance.resize(size);
        if constexpr (std::is_same_v<ModelT, DependentT>)
        {
          model.getSignals().template attachPort<EMT::DependentVoltageSourceExternalVariables::VA>(&terminal.voltagePort());
          model.getSignals().template attachPort<EMT::DependentVoltageSourceExternalVariables::EA>(&drive.voltagePort());
        }
        else if constexpr (std::is_same_v<ModelT, SourceT>)
        {
          model.getSignals().template attachPort<EMT::VoltageSourceExternalVariables::VA>(&terminal.voltagePort());
        }
        else
        {
          model.getSignals().template attachPort<EMT::LoadZExternalVariables::VA>(&terminal.voltagePort());
        }
        IdxT offset = 0;
        for (auto* component : components())
        {
          component->bind(y, yp, f, tolerance, offset);
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

      std::array<ComponentT*, 3> components()
      {
        return {&drive, &terminal, &model};
      }

      void setPhasors(IdxT offset, RealT omega, const Phasors& values)
      {
        for (IdxT n = 0; n < 3; ++n)
        {
          y.getData()[offset + n]  = values[n].real();
          yp.getData()[offset + n] = -omega * values[n].imag();
        }
        y.setDataUpdated();
        yp.setDataUpdated();
      }

      void evaluate()
      {
        for (auto* component : components())
          component->evaluateInternalResidual();
        for (auto* component : components())
          component->evaluateExternalResidual();
      }

      bool zeroInternalResidual()
      {
        evaluate();
        for (IdxT n = 6; n < y.getSize(); ++n)
        {
          if (std::abs(f.getData()[n]) > 2.0e-11)
            return false;
        }
        return true;
      }
    };

  public:
    TestOutcome dependentSteadyState()
    {
      TestStatus                      success = true;
      typename DependentT::ModelDataT data;
      data.Y = fitData();
      const Phasors voltage{{{1.2, -0.3}, {-0.7, 1.1}, {-0.5, -0.8}}};
      const Phasors drop{{{0.3, 0.4}, {-0.2, 0.5}, {0.7, -0.6}}};
      for (RealT omega : {0.9, 7.3, 21.0})
      {
        Fixture<DependentT> fixture(data);
        Phasors             drive = voltage;
        for (size_t n = 0; n < 3; ++n)
          drive[n] += drop[n];
        fixture.setPhasors(0, omega, drive);
        fixture.setPhasors(3, omega, voltage);
        success             *= (fixture.model.initializeSteadyState(omega) == 0);
        success             *= fixture.zeroInternalResidual();
        const auto expected  = applyTransfer(*data.Y, omega, drop);
        for (size_t n = 0; n < 3; ++n)
        {
          success *= isEqual(fixture.f.getData()[3 + n], expected[n].real(), 2.0e-12);
        }
      }
      return success.report(__func__);
    }

    TestOutcome loadSteadyState()
    {
      TestStatus                 success = true;
      typename LoadT::ModelDataT data;
      data.Z          = fitData();
      // Only current column 1 has derivative terms; the other currents are algebraic.
      data.Z->E[0][1] = 0.03;
      data.Z->E[1][1] = 0.2;
      data.Z->E[2][1] = 0.05;
      const Phasors current{{{0.3, 0.4}, {-0.2, 0.5}, {0.7, -0.6}}};
      for (RealT omega : {0.9, 7.3, 21.0})
      {
        Fixture<LoadT> fixture(data);
        auto           voltage = applyTransfer(*data.Z, omega, current);
        for (auto& phase : voltage)
          phase = -phase;
        fixture.setPhasors(3, omega, voltage);
        success *= (fixture.model.verify() == 0);
        success *= (fixture.model.tag()[0] == false);
        success *= (fixture.model.tag()[1] == true);
        success *= (fixture.model.tag()[2] == false);
        success *= (fixture.model.initializeSteadyState(omega) == 0);
        success *= fixture.zeroInternalResidual();
        for (size_t n = 0; n < 3; ++n)
        {
          success *= isEqual(fixture.y.getData()[6 + n], current[n].real(), 2.0e-12);
          success *= isEqual(fixture.yp.getData()[6 + n], -omega * current[n].imag(), 2.0e-12);
          success *= isEqual(fixture.f.getData()[3 + n], current[n].real(), 2.0e-12);
        }
        success *= (fixture.model.initializeSteadyState(0.0) != 0);
        success *= (fixture.model.initializeSteadyState(std::numeric_limits<RealT>::quiet_NaN()) != 0);
      }
      typename LoadT::ModelDataT singular;
      singular.Z = FitDataT{};
      Fixture<LoadT> fixture(singular);
      success *= (fixture.model.initializeSteadyState(7.3) != 0);
      return success.report(__func__);
    }

    TestOutcome voltageSteadyState()
    {
      TestStatus                   success = true;
      typename SourceT::ModelDataT data;
      using Parameter = typename SourceT::ModelDataT::Parameters;
      const EMT::ABCVector<RealT> magnitude{{1.0, 1.2, 0.8}};
      const EMT::ABCVector<RealT> phase{{0.1, -2.0, 2.1}};
      constexpr RealT             omega = 7.3;
      constexpr RealT             time  = 0.23;
      data.parameters[Parameter::E]     = magnitude;
      data.parameters[Parameter::phi]   = phase;
      data.parameters[Parameter::omega] = omega;
      data.Y                            = fitData();
      Fixture<SourceT> fixture(data);
      fixture.model.updateTime(time, 0.0);
      const Phasors voltage{{{1.2, -0.3}, {-0.7, 1.1}, {-0.5, -0.8}}};
      fixture.setPhasors(3, omega, voltage);
      Phasors drop{};
      for (size_t n = 0; n < 3; ++n)
      {
        drop[n] = std::sqrt(2.0) * magnitude[n] * std::exp(ComplexT{0.0, omega * time + phase[n]}) - voltage[n];
      }
      success             *= (fixture.model.initializeSteadyState(omega) == 0);
      success             *= fixture.zeroInternalResidual();
      const auto expected  = applyTransfer(*data.Y, omega, drop);
      for (size_t n = 0; n < 3; ++n)
      {
        success *= isEqual(fixture.f.getData()[3 + n], expected[n].real(), 2.0e-12);
        success *= (fixture.model.tag()[n] == false);
        success *= (fixture.model.tag()[3 + n] == false);
      }
      return success.report(__func__);
    }

    TestOutcome legacySteadyState()
    {
      TestStatus            success = true;
      FitDataT              impedance;
      EMT::ABCMatrix<RealT> R{};
      EMT::ABCMatrix<RealT> L{};
      for (size_t n = 0; n < 3; ++n)
      {
        for (size_t k = 0; k < 3; ++k)
        {
          R[n][k]           = n == k ? 1.7 + 0.1 * static_cast<RealT>(n) : 0.05 * static_cast<RealT>(1 + n + k);
          L[n][k]           = n == k ? 0.2 + 0.02 * static_cast<RealT>(n) : 0.003 * static_cast<RealT>(1 + n + k);
          impedance.D[n][k] = R[n][k];
          impedance.E[n][k] = L[n][k];
        }
      }
      typename DependentT::ModelDataT source_data;
      source_data.parameters[DependentT::ModelDataT::Parameters::Rs] = R;
      source_data.parameters[DependentT::ModelDataT::Parameters::Ls] = L;
      typename LoadT::ModelDataT load_data;
      load_data.parameters[LoadT::ModelDataT::Parameters::R] = R;
      load_data.parameters[LoadT::ModelDataT::Parameters::L] = L;
      const Phasors       current{{{0.3, 0.4}, {-0.2, 0.5}, {0.7, -0.6}}};
      constexpr RealT     omega   = 7.3;
      auto                voltage = applyTransfer(impedance, omega, current);
      Fixture<DependentT> source(source_data);
      source.setPhasors(0, omega, voltage);
      success *= (source.model.size() == 3);
      success *= (source.model.verify() == 0);
      success *= (source.model.initializeSteadyState(omega) == 0);
      success *= source.zeroInternalResidual();
      for (auto& phase : voltage)
        phase = -phase;
      Fixture<LoadT> load(load_data);
      load.setPhasors(3, omega, voltage);
      success *= (load.model.size() == 3);
      success *= (load.model.verify() == 0);
      success *= (load.model.initializeSteadyState(omega) == 0);
      success *= load.zeroInternalResidual();
      for (size_t n = 0; n < 3; ++n)
      {
        success *= isEqual(source.y.getData()[6 + n], current[n].real(), 2.0e-12);
        success *= isEqual(source.yp.getData()[6 + n], -omega * current[n].imag(), 2.0e-12);
        success *= isEqual(load.y.getData()[6 + n], current[n].real(), 2.0e-12);
        success *= isEqual(load.yp.getData()[6 + n], -omega * current[n].imag(), 2.0e-12);
      }
      return success.report(__func__);
    }

    TestOutcome validation()
    {
      TestStatus                 success = true;
      typename LoadT::ModelDataT load;
      load.Z          = fitData();
      load.Z->E[0][0] = 0.2;
      load.Z->E[0][1] = 0.4;
      Fixture<LoadT> dependent_columns(load);
      const auto     previous = EMT::Log::verbosity();
      EMT::Log::setVerbosity(EMT::Log::Verbosity::NONE);
      success   *= (dependent_columns.model.verify() > 0);
      load.Z->E  = EMT::RationalMatrix<RealT>(2, 3);
      LoadT malformed_coefficient(load);
      BusT  bus;
      bus.allocate();
      malformed_coefficient.getSignals().template attachPort<EMT::LoadZExternalVariables::VA>(&bus.voltagePort());
      success *= (malformed_coefficient.allocate() != 0);
      success *= (malformed_coefficient.verify() > 0);
      typename SourceT::ModelDataT source;
      source.Y                                                  = fitData();
      source.parameters[SourceT::ModelDataT::Parameters::omega] = std::numeric_limits<RealT>::infinity();
      Fixture<SourceT> nonfinite_frequency(source);
      success                                                   *= (nonfinite_frequency.model.verify() > 0);
      source.parameters[SourceT::ModelDataT::Parameters::omega]  = 7.3;
      source.parameters[SourceT::ModelDataT::Parameters::E]      = EMT::ABCVector<RealT>{{1.0, -0.1, 1.0}};
      Fixture<SourceT> negative_magnitude(source);
      success                                               *= (negative_magnitude.model.verify() > 0);
      source.parameters[SourceT::ModelDataT::Parameters::E]  = EMT::ABCVector<RealT>{{1.0, 1.0, 1.0}};
      source.parameters[SourceT::ModelDataT::Parameters::N]  = IdxT{2};
      Fixture<SourceT> wrong_phases(source);
      success *= (wrong_phases.model.verify() > 0);
      EMT::Log::setVerbosity(previous);
      return success.report(__func__);
    }

    TestOutcome idealShortInitialization()
    {
      TestStatus                 success = true;
      typename LoadT::ModelDataT legacy_data;
      typename LoadT::ModelDataT rational_data;
      rational_data.Z = FitDataT{};
      Fixture<LoadT> legacy(legacy_data);
      Fixture<LoadT> rational(rational_data);
      const Phasors  voltage{{{1.2, -0.3}, {-0.7, 1.1}, {-0.5, -0.8}}};
      legacy.setPhasors(3, 7.3, voltage);
      rational.setPhasors(3, 7.3, voltage);
      success *= (legacy.model.verify() == 0);
      success *= (rational.model.verify() == 0);
      success *= (legacy.model.initialize() == 0);
      success *= (rational.model.initialize() == 0);
      legacy.evaluate();
      rational.evaluate();
      for (size_t n = 0; n < 3; ++n)
      {
        success *= (legacy.model.tag()[n] == false);
        success *= (rational.model.tag()[n] == false);
        success *= isEqual(legacy.y.getData()[6 + n], 0.0, 1.0e-14);
        success *= isEqual(rational.y.getData()[6 + n], 0.0, 1.0e-14);
        success *= isEqual(legacy.f.getData()[6 + n], voltage[n].real(), 1.0e-14);
        success *= isEqual(rational.f.getData()[6 + n], voltage[n].real(), 1.0e-14);
      }
      success *= (legacy.model.initializeSteadyState(7.3) != 0);
      success *= (rational.model.initializeSteadyState(7.3) != 0);
      return success.report(__func__);
    }

    template <typename ModelT>
    static bool standaloneStorage(const typename ModelT::ModelDataT& data, size_t branch_offset, size_t fit_offset)
    {
      BusT bus;
      bus.allocate();
      bus.initialize();
      ModelT model(data);
      if constexpr (std::is_same_v<ModelT, SourceT>)
      {
        model.getSignals().template attachPort<EMT::VoltageSourceExternalVariables::VA>(&bus.voltagePort());
      }
      else
      {
        model.getSignals().template attachPort<EMT::LoadZExternalVariables::VA>(&bus.voltagePort());
      }
      if (model.allocate() != 0 || model.initialize() != 0 || model.y().getSize() != model.size())
        return false;
      for (size_t n = 0; n < 3; ++n)
      {
        bus.y().getData()[n]                   = 7.0 + static_cast<RealT>(n);
        model.y().getData()[branch_offset + n] = 0.5 + static_cast<RealT>(n);
        model.y().getData()[fit_offset + n]    = 0.2 + static_cast<RealT>(n);
        model.yp().getData()[fit_offset + n]   = -0.1 * static_cast<RealT>(n + 1);
      }
      bus.y().setDataUpdated();
      model.y().setDataUpdated();
      model.yp().setDataUpdated();
      bus.evaluateInternalResidual();
      model.evaluateResidual();
      for (size_t n = 0; n < 3; ++n)
      {
        const RealT input      = model.y().getData()[branch_offset + n];
        const RealT state      = model.y().getData()[fit_offset + n];
        const RealT derivative = model.yp().getData()[fit_offset + n];
        if (std::abs(model.getResidual().getData()[fit_offset + n] + derivative + 3.0 * state - input) > 1.0e-13)
          return false;
        const RealT output = 2.0 * input + 4.0 * state;
        if constexpr (std::is_same_v<ModelT, SourceT>)
        {
          if (std::abs(bus.getResidual().getData()[n] - output) > 1.0e-13)
            return false;
        }
        else
        {
          if (std::abs(model.getResidual().getData()[n] - bus.y().getData()[n] - output) > 1.0e-13)
            return false;
        }
      }
      return true;
    }

    TestOutcome standaloneStorage()
    {
      TestStatus success = true;
      FitDataT   fit;
      fit.poles.emplace_back(-3.0, 0.0);
      fit.residues.resize(1);
      for (size_t n = 0; n < 3; ++n)
      {
        fit.D[n][n]           = 2.0;
        fit.residues[0][n][n] = 4.0;
      }
      typename LoadT::ModelDataT load;
      load.Z   = fit;
      success *= standaloneStorage<LoadT>(load, 0, 3);
      typename SourceT::ModelDataT source;
      source.Y                                                   = fit;
      source.parameters[SourceT::ModelDataT::Parameters::omega]  = 7.3;
      success                                                   *= standaloneStorage<SourceT>(source, 3, 6);
      return success.report(__func__);
    }

    template <typename ModelT>
    static bool finiteDifferenceJacobian(const typename ModelT::ModelDataT& data)
    {
      Fixture<ModelT> fixture(data);
      const IdxT      size = fixture.y.getSize();
      for (IdxT n = 0; n < size; ++n)
      {
        fixture.y.getData()[n]  = 0.7 + 0.13 * static_cast<RealT>(n);
        fixture.yp.getData()[n] = -0.4 + 0.21 * static_cast<RealT>(n);
      }
      fixture.y.setDataUpdated();
      fixture.yp.setDataUpdated();
      for (RealT alpha : {1.7, 5.3})
      {
        fixture.model.updateTime(0.17, alpha);
        fixture.model.evaluateJacobian();
        std::map<std::pair<IdxT, IdxT>, RealT> entries;
        auto*                                  coo = fixture.model.getCooJacobian();
        if (coo == nullptr)
          return false;
        for (IdxT n = 0; n < coo->getNnz(); ++n)
          entries[{coo->getRowData()[n], coo->getColData()[n]}] += coo->getValues()[n];
        constexpr RealT step = 1.0e-5;
        for (IdxT k = 0; k < size; ++k)
        {
          const RealT original    = fixture.y.getData()[k];
          const RealT derivative  = fixture.yp.getData()[k];
          fixture.y.getData()[k]  = original + step;
          fixture.yp.getData()[k] = derivative + alpha * step;
          fixture.y.setDataUpdated();
          fixture.yp.setDataUpdated();
          fixture.evaluate();
          std::vector<RealT> plus(fixture.f.getData(), fixture.f.getData() + size);
          fixture.y.getData()[k]  = original - step;
          fixture.yp.getData()[k] = derivative - alpha * step;
          fixture.y.setDataUpdated();
          fixture.yp.setDataUpdated();
          fixture.evaluate();
          for (IdxT n = 0; n < size; ++n)
          {
            const RealT expected = (plus[n] - fixture.f.getData()[n]) / (2.0 * step);
            if (std::abs(entries[{n, k}] - expected) > 2.0e-8)
              return false;
          }
          fixture.y.getData()[k]  = original;
          fixture.yp.getData()[k] = derivative;
        }
      }
      return true;
    }

    TestOutcome jacobians()
    {
      TestStatus                      success = true;
      typename DependentT::ModelDataT dependent;
      dependent.Y  = fitData();
      success     *= finiteDifferenceJacobian<DependentT>(dependent);
      typename LoadT::ModelDataT load;
      load.Z           = fitData();
      load.Z->E[0][1]  = 0.03;
      load.Z->E[1][1]  = 0.2;
      success         *= finiteDifferenceJacobian<LoadT>(load);
      typename SourceT::ModelDataT source;
      source.Y                                                   = fitData();
      source.parameters[SourceT::ModelDataT::Parameters::omega]  = 7.3;
      success                                                   *= finiteDifferenceJacobian<SourceT>(source);
      return success.report(__func__);
    }
  };
} // namespace GridKit::Testing
