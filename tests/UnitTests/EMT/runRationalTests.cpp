#include <cmath>
#include <complex>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>

#include <GridKit/Model/EMT/Operators/Rational/StateSpace/StateSpace.hpp>
#include <GridKit/Model/EMT/Operators/Rational/StateSpace/StateSpaceDataJSONParser.hpp>
#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFit.hpp>
#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFitDataJSONParser.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace
{
  using namespace GridKit::EMT;
  using Complex = std::complex<double>;
  using VFData  = VectorFitData<double, size_t>;
  using SSData  = StateSpaceData<double, size_t>;
  using VF      = VectorFit<double, size_t>;
  using SS      = StateSpace<double, size_t>;
  using GridKit::Testing::isEqual;
  using GridKit::Testing::TestStatus;

  template <typename Model>
  struct Fixture
  {
    using Vector = GridKit::LinearAlgebra::Vector<double, size_t>;
    Model                               model;
    Vector                              y, yp, f, tolerance;
    std::vector<double>                 u, udot, out;
    std::vector<size_t>                 input_indices, output_indices;
    std::vector<Signal<double, size_t>> input_signals, output_signals;

    explicit Fixture(const typename Model::ModelDataT& data)
      : model(data), u(static_cast<size_t>(data.cols)), udot(u.size()), out(static_cast<size_t>(data.rows)), input_indices(u.size()), output_indices(out.size()), input_signals(u.size()), output_signals(out.size())
    {
      y.resize(model.size());
      yp.resize(model.size());
      f.resize(model.size());
      tolerance.resize(model.size());
      std::vector<Signal<double, size_t>*> inputs, outputs;
      for (size_t k = 0; k < u.size(); ++k)
      {
        input_indices[k] = model.size() + k;
        input_signals[k].set(&u[k], &udot[k], nullptr, &input_indices[k], nullptr);
        inputs.push_back(&input_signals[k]);
      }
      for (size_t n = 0; n < out.size(); ++n)
      {
        output_indices[n] = model.size() + n;
        output_signals[n].set(nullptr, nullptr, &out[n], nullptr, &output_indices[n]);
        outputs.push_back(&output_signals[n]);
      }
      model.attachInput(inputs);
      model.attachOutput(outputs);
      model.bind(y, yp, f, tolerance, 0);
      if (model.allocate() || model.initialize())
      {
        throw std::runtime_error("Rational fixture allocation failed");
      }
      model.tagDifferentiable();
    }

    std::vector<double> residual()
    {
      std::fill(out.begin(), out.end(), 0.0);
      model.evaluateResidual();
      std::vector<double> result;
      for (size_t j = 0; j < model.size(); ++j)
      {
        result.push_back(f.getData()[j]);
      }
      result.insert(result.end(), out.begin(), out.end());
      return result;
    }

    bool jacobian(double alpha)
    {
      model.updateTime(0.0, alpha);
      model.evaluateJacobian();
      auto*                                       coo = model.getCooJacobian();
      std::map<std::pair<size_t, size_t>, double> entries;
      for (size_t j = 0; j < coo->getNnz(); ++j)
      {
        entries[{coo->getRowData()[j], coo->getColData()[j]}] += coo->getValues()[j];
      }
      const double h       = 1.e-5;
      bool         success = true;
      for (size_t k = 0; k < model.size() + u.size(); ++k)
      {
        double&      value      = k < model.size() ? y.getData()[k] : u[k - model.size()];
        double&      derivative = k < model.size() ? yp.getData()[k] : udot[k - model.size()];
        const double original = value, original_derivative = derivative;
        value            = original + h;
        derivative       = original_derivative + alpha * h;
        const auto plus  = residual();
        value            = original - h;
        derivative       = original_derivative - alpha * h;
        const auto minus = residual();
        value            = original;
        derivative       = original_derivative;
        for (size_t n = 0; n < plus.size(); ++n)
        {
          success &= isEqual(entries[{n, k}], (plus[n] - minus[n]) / (2 * h), 1.e-8);
        }
      }
      return success;
    }
  };

  VFData vectorData()
  {
    VFData data;
    data.rows  = 2;
    data.cols  = 4;
    data.D     = RationalMatrix<double>(2, 4);
    data.E     = RationalMatrix<double>(2, 4);
    data.poles = {{-3, 0}, {-2, 7}, {-2, -7}};
    data.residues.assign(3, RationalMatrix<Complex>(2, 4));
    for (size_t n = 0; n < 2; ++n)
    {
      for (size_t k = 0; k < 4; ++k)
      {
        data.D[n][k]           = 0.2 + 0.1 * static_cast<double>(n + k);
        data.E[n][k]           = k == 1 ? 0.03 * static_cast<double>(n + 1) : 0;
        data.residues[0][n][k] = 0.3 + 0.05 * static_cast<double>(n + k);
        data.residues[1][n][k] = {0.2 * static_cast<double>(n + 1), 0.04 * static_cast<double>(k + 1)};
        data.residues[2][n][k] = std::conj(data.residues[1][n][k]);
      }
    }
    return data;
  }

  SSData stateData()
  {
    SSData data;
    data.rows  = 2;
    data.cols  = 4;
    data.D     = RationalMatrix<double>(2, 4);
    data.E     = RationalMatrix<double>(2, 4);
    // Repeated real poles give a rank-two combined residue; one conjugate pair.
    data.poles = {{-3, 0}, {-3, 0}, {-2, 7}, {-2, -7}};
    data.C     = RationalMatrix<Complex>(2, 4);
    data.B     = RationalMatrix<Complex>(4, 4);
    for (size_t n = 0; n < 2; ++n)
    {
      data.C[n][n] = 1.0;
      data.C[n][2] = {0.5 + static_cast<double>(n), -0.2 * static_cast<double>(n + 1)};
      data.C[n][3] = std::conj(data.C[n][2]);
    }
    for (size_t k = 0; k < 4; ++k)
    {
      data.B[0][k] = 0.2 + 0.1 * static_cast<double>(k);
      data.B[1][k] = (k % 2) ? 0.4 : -0.3;
      data.B[2][k] = {0.4 + 0.1 * static_cast<double>(k), 0.3 - 0.05 * static_cast<double>(k)};
      data.B[3][k] = std::conj(data.B[2][k]);
      data.D[0][k] = 0.2;
      data.D[1][k] = -0.3;
      data.E[0][k] = k == 1 ? 0.03 : 0.0;
    }
    return data;
  }

  template <typename Data>
  Complex expected(const Data& data, size_t n, size_t k, double omega)
  {
    Complex H{data.D[n][k], omega * data.E[n][k]};
    for (size_t q = 0; q < data.poles.size(); ++q)
    {
      Complex residue;
      if constexpr (std::is_same_v<Data, VFData>)
      {
        residue = data.residues[q][n][k];
      }
      else
      {
        residue = data.C[n][q] * data.B[q][k];
      }
      H += residue / (Complex{0, omega} - data.poles[q]);
    }
    return H;
  }

  template <typename Model>
  bool response(const typename Model::ModelDataT& data)
  {
    Fixture<Model> fixture(data);
    bool           success = true;
    for (double omega : {0.0, 0.3, 7.0, 50.0, -11.0})
    {
      for (size_t k = 0; k < fixture.u.size(); ++k)
      {
        fixture.u[k]    = 0.7 + 0.2 * static_cast<double>(k);
        fixture.udot[k] = -omega * (0.4 - 0.1 * static_cast<double>(k));
      }
      success             &= fixture.model.initializeSteadyState(omega, fixture.u, fixture.udot) == 0;
      const auto residual  = fixture.residual();
      for (size_t j = 0; j < fixture.model.size(); ++j)
      {
        success &= isEqual(residual[j], 0.0, 2.e-13);
      }
      RationalMatrix<double> re, im;
      fixture.model.transfer(omega, re, im);
      for (size_t n = 0; n < fixture.out.size(); ++n)
      {
        double output = 0.0;
        for (size_t k = 0; k < fixture.u.size(); ++k)
        {
          const auto H  = expected(data, n, k, omega);
          success      &= isEqual(re[n][k], H.real(), 2.e-13) && isEqual(im[n][k], H.imag(), 2.e-13);
          const Complex U{fixture.u[k], omega == 0 ? 0 : -fixture.udot[k] / omega};
          output += (H * U).real();
        }
        success &= isEqual(fixture.out[n], output, 2.e-13);
      }
      success &= fixture.jacobian(3.7) && fixture.jacobian(-0.8);
    }
    return success;
  }

  bool validation()
  {
    bool success            = true;
    auto data               = vectorData();
    data.residues[0][0][0]  = {std::numeric_limits<double>::infinity(), 0};
    success                &= data.validate() != 0;
    data                    = vectorData();
    data.residues[1][1].pop_back();
    success           &= data.validate() != 0;
    auto ss            = stateData();
    ss.B[3][0]        += Complex{0.1, 0};
    success           &= ss.validate() != 0;
    ss                 = stateData();
    ss.C[0][0]        += Complex{0, 0.1};
    success           &= ss.validate() != 0;
    const auto j       = nlohmann::json::parse(R"({"rows":1,"cols":2,"D":[[1,2]],"poles":[[-3,0]],"residues":[[[[0.2,0],[0.4,0]]]]})");
    const auto parsed  = j.get<VFData>();
    success           &= parsed.rows == 1 && parsed.cols == 2 && parsed.validate() == 0;
    auto reused        = parsed;
    from_json(j, reused);
    success &= reused.poles.size() == 1;
    for (const auto& dimension : {nlohmann::json(-1), nlohmann::json(0), nlohmann::json(1.5), nlohmann::json("3")})
    {
      auto bad      = j;
      bad["rows"]   = dimension;
      bool rejected = false;
      try
      {
        (void) bad.get<VFData>();
      }
      catch (const std::exception&)
      {
        rejected = true;
      }
      success  &= rejected;
      rejected  = false;
      try
      {
        (void) bad.get<SSData>();
      }
      catch (const std::exception&)
      {
        rejected = true;
      }
      success &= rejected;
    }
    auto bad = j;
    bad["D"][0].push_back(3);
    bool rejected = false;
    try
    {
      (void) bad.get<VFData>();
    }
    catch (const std::exception&)
    {
      rejected = true;
    }
    success            &= rejected;
    const auto ss_json  = nlohmann::json::parse(R"({"rows":1,"cols":2,"poles":[[-3,0]],"C":[[[2,0]]],"B":[[[3,0],[4,0]]]})");
    success            &= ss_json.get<SSData>().validate() == 0;
    VFData scaled;
    scaled.D[0][0] = std::numeric_limits<double>::max();
    VF                    invalid(scaled, 2.0);
    Port3<double, size_t> unused_port;
    invalid.attachInput(&unused_port);
    invalid.attachOutput(&unused_port);
    success &= invalid.verify() != 0;
    return success;
  }

  bool boundaries()
  {
    bool        success = true;
    auto        data    = vectorData();
    Fixture<VF> fixture(data);
    success &= fixture.model.size() == 12;
    for (size_t k = 0; k < 4; ++k)
    {
      success &= fixture.input_signals[k].hasDerivativeCoupling() == (k == 1);
    }
    std::fill(fixture.u.begin(), fixture.u.end(), 2.0);
    success            &= fixture.model.initialize() == 0;
    const auto initial  = fixture.residual();
    for (size_t j = 0; j < fixture.model.size(); ++j)
    {
      success &= isEqual(initial[j], 0.0, 1.e-14);
    }
    data.poles.clear();
    data.residues.clear();
    Fixture<VF> feedthrough(data);
    success       &= feedthrough.model.size() == 0 && feedthrough.jacobian(0.0) && feedthrough.jacobian(2.0);
    data           = vectorData();
    data.poles[0]  = {0, 0};
    Fixture<VF> integrator(data);
    success       &= integrator.model.initializeSteadyState(0, integrator.u, integrator.udot) != 0;
    data           = vectorData();
    data.poles[1]  = {0, 7};
    data.poles[2]  = {0, -7};
    Fixture<VF> resonance(data);
    success &= resonance.model.initializeSteadyState(7, resonance.u, resonance.udot) != 0;
    fixture.input_signals[0].setComputed(
        [&]()
        { return fixture.u[0] * fixture.u[0] + 2.0 * fixture.u[2]; },
        [&](auto& gradient, double scale)
        {
          gradient.emplace_back(fixture.input_indices[0], scale * 2.0 * fixture.u[0]);
          gradient.emplace_back(fixture.input_indices[2], scale * 2.0);
        });
    success      &= fixture.jacobian(1.7);
    fixture.u[0]  = 0.0;
    success      &= fixture.jacobian(-2.0);
    fixture.input_signals[1].setComputed([]()
                                         { return 0.0; },
                                         [](auto&, double) {});
    success &= fixture.model.verify() != 0;
    Fixture<SS> states(stateData());
    success &= states.model.size() == 4;
    return success;
  }
} // namespace

int main()
{
  GridKit::Testing::TestingResults results;
  TestStatus                       vf  = response<VF>(vectorData());
  results                             += vf.report("rectangular VectorFit response and Jacobian");
  TestStatus ss                        = response<SS>(stateData());
  results                             += ss.report("factorized StateSpace response and Jacobian");
  TestStatus valid                     = validation();
  results                             += valid.report("rational coefficient and JSON validation");
  TestStatus edge                      = boundaries();
  results                             += edge.report("rational zero-state and singular E cases");
  return results.summary();
}
