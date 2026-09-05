/**
 * @file VectorFitTests.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief VectorFit rational operator tests on a hand-bound fixture.
 *
 */
#pragma once

#include <array>
#include <cmath>
#include <map>
#include <numbers>
#include <utility>

#include <GridKit/Definitions.hpp>
#include <GridKit/LinearAlgebra/Vector/Vector.hpp>
#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFit.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace Testing
  {
    /**
     * @brief Tests for the vector-fitted rational operator.
     *
     * The fixture binds the input and output ports over raw arrays, so the
     * operator is exercised standalone with no auxiliary components: the
     * input triple plays the consumer variables and the output rows play the
     * consumer residual rows.
     */
    template <typename ScalarT, typename IdxT>
    class VectorFitTests
    {
    public:
      using RealT      = ScalarT;
      using VectorT    = GridKit::LinearAlgebra::Vector<ScalarT, IdxT>;
      using VectorFitT = GridKit::EMT::VectorFit<ScalarT, IdxT>;
      using DataT      = GridKit::EMT::VectorFitData<RealT, IdxT>;
      using ComplexT   = std::complex<RealT>;

      VectorFitTests()  = default;
      ~VectorFitTests() = default;

    private:
      static constexpr RealT test_scale = 2.5;

      /// One real pole, one conjugate pair, and full asymmetric D and E.
      static DataT makeResidualData()
      {
        DataT data;
        data.D = {{{{1.7, 0.21, 0.32}}, {{0.13, 1.9, 0.24}}, {{0.35, 0.16, 2.1}}}};
        data.E = {{{{0.041, 0.002, 0.013}}, {{0.004, 0.052, 0.015}}, {{0.016, 0.007, 0.063}}}};

        data.poles.emplace_back(-3.1, 0.0);
        auto& real_residue = data.residues.emplace_back();
        real_residue       = {{{{ComplexT{0.61, 0.0}, ComplexT{0.041, 0.0}, ComplexT{0.022, 0.0}}},
                               {{ComplexT{0.043, 0.0}, ComplexT{0.72, 0.0}, ComplexT{0.054, 0.0}}},
                               {{ComplexT{0.025, 0.0}, ComplexT{0.046, 0.0}, ComplexT{0.83, 0.0}}}}};

        const ComplexT pair_pole{-2.2, 7.3};
        data.poles.push_back(pair_pole);
        data.poles.push_back(std::conj(pair_pole));
        auto& pair_residue           = data.residues.emplace_back();
        pair_residue                 = {{{{ComplexT{0.31, 0.15}, ComplexT{0.011, 0.021}, ComplexT{0.012, 0.032}}},
                                         {{ComplexT{0.013, 0.043}, ComplexT{0.32, 0.25}, ComplexT{0.014, 0.054}}},
                                         {{ComplexT{0.015, 0.065}, ComplexT{0.016, 0.076}, ComplexT{0.33, 0.35}}}}};
        const auto pair_copy         = pair_residue;
        auto&      conjugate_residue = data.residues.emplace_back();
        for (size_t n = 0; n < 3; ++n)
        {
          for (size_t k = 0; k < 3; ++k)
          {
            conjugate_residue[n][k] = std::conj(pair_copy[n][k]);
          }
        }

        return data;
      }

      /**
       * @brief The operator alone, with ports over raw arrays.
       *
       * Variable layout: operator states [0, 9), input triple [9, 12); the
       * output rows share the input triple's global row identifiers, playing
       * a consumer whose rows are its own variables.
       */
      struct Fixture
      {
        static constexpr IdxT state_size  = 9;
        static constexpr IdxT system_size = 12;

        VectorT y;
        VectorT yp;
        VectorT f;
        VectorT abs_tol;

        std::array<ScalarT, 3> u{};
        std::array<ScalarT, 3> u_dot{};
        std::array<ScalarT, 3> out_row{};
        std::array<IdxT, 3>    io_index{9, 10, 11};

        GridKit::EMT::Port3<ScalarT, IdxT> io_port;

        VectorFitT vf;

        Fixture()
          : vf(makeResidualData(), test_scale)
        {
          y.resize(state_size);
          yp.resize(state_size);
          f.resize(state_size);
          abs_tol.resize(state_size);

          for (size_t n = 0; n < 3; ++n)
          {
            io_port.signals[n].set(&u[n], &u_dot[n], &out_row[n], &io_index[n], &io_index[n]);
          }

          vf.attachInput(&io_port);
          vf.attachOutput(&io_port);

          vf.bind(y, yp, f, abs_tol, 0);
          vf.allocate();
          for (IdxT j = 0; j < vf.size(); ++j)
          {
            vf.setVariableIndex(j, j);
            vf.setResidualIndex(j, j);
          }
          vf.initialize();
          vf.tagDifferentiable();
        }

        /// Assembled residual: state rows in f, output rows accumulated from
        /// zero into the raw output array.
        void evaluateResidual()
        {
          out_row[0] = 0.0;
          out_row[1] = 0.0;
          out_row[2] = 0.0;
          vf.evaluateInternalResidual();
          vf.evaluateExternalResidual();
        }

        /// Fill states and inputs with distinct values away from 0 and 1.
        void setProbeState()
        {
          auto* y_data  = y.getData();
          auto* yp_data = yp.getData();
          for (IdxT j = 0; j < state_size; ++j)
          {
            y_data[j]  = 0.7 + 0.31 * static_cast<RealT>(j);
            yp_data[j] = -1.3 + 0.23 * static_cast<RealT>(j);
          }
          for (size_t n = 0; n < 3; ++n)
          {
            u[n]     = 2.1 + 0.7 * static_cast<RealT>(n);
            u_dot[n] = -0.9 + 0.4 * static_cast<RealT>(n);
          }
          y.setDataUpdated();
          yp.setDataUpdated();
        }
      };

    public:
      /**
       * @brief Specification violations are counted and clean data passes.
       */
      TestOutcome validation()
      {
        TestStatus success = true;

        using Log = GridKit::Utilities::Logger;
        Log::setVerbosity(Log::Verbosity::NONE);

        success *= (makeResidualData().validate() == 0);

        // A nonreal pole without its conjugate
        auto broken_pair = makeResidualData();
        broken_pair.poles.pop_back();
        broken_pair.residues.pop_back();
        success *= (broken_pair.validate() > 0);

        // A complex residue on a real pole
        auto broken_real               = makeResidualData();
        broken_real.residues[0][1][2]  = {0.022, 0.5};
        success                       *= (broken_real.validate() > 0);

        // A conjugate residue mismatch
        auto broken_residue               = makeResidualData();
        broken_residue.residues[2][0][0] += ComplexT{0.0, 0.1};
        success                          *= (broken_residue.validate() > 0);

        Log::setVerbosity(Log::Verbosity::EVERYTHING);

        return success.report(__func__);
      }

      /**
       * @brief Assembled residual against an independent computation from
       * the original complex data.
       */
      TestOutcome residual()
      {
        TestStatus success = true;

        Fixture fixture;
        fixture.setProbeState();
        fixture.vf.updateTime(0.0, 1.0);
        fixture.evaluateResidual();

        const auto* y  = fixture.y.getData();
        const auto* yp = fixture.yp.getData();
        const auto* f  = fixture.f.getData();

        const auto data = makeResidualData();

        // State rows: real section [0, 3), pair section [3, 9)
        for (size_t n = 0; n < 3; ++n)
        {
          const RealT a_r  = data.poles[0].real();
          success         *= isEqual<RealT>(f[n], -yp[n] + a_r * y[n] + fixture.u[n], 1.0e-14);

          const RealT a_c  = data.poles[1].real();
          const RealT w_c  = data.poles[1].imag();
          success         *= isEqual<RealT>(f[3 + n],
                                    -yp[3 + n] + a_c * y[3 + n] - w_c * y[6 + n] + fixture.u[n],
                                    1.0e-14);
          success         *= isEqual<RealT>(f[6 + n],
                                    -yp[6 + n] + w_c * y[3 + n] + a_c * y[6 + n],
                                    1.0e-14);
        }

        // Output rows: scale (D u + E du) + scale A_r w + 2 scale (A_c w - B_c v)
        for (size_t n = 0; n < 3; ++n)
        {
          RealT expected = 0.0;
          for (size_t k = 0; k < 3; ++k)
          {
            expected += test_scale * (data.D[n][k] * fixture.u[k] + data.E[n][k] * fixture.u_dot[k]);
            expected += test_scale * data.residues[0][n][k].real() * y[k];
            expected += TWO<RealT> * test_scale
                        * (data.residues[1][n][k].real() * y[3 + k]
                           - data.residues[1][n][k].imag() * y[6 + k]);
          }
          success *= isEqual<RealT>(fixture.out_row[n], expected, 1.0e-14);
        }

        return success.report(__func__);
      }

      /**
       * @brief Sinusoidal steady state: initialized states satisfy the state
       * rows exactly, and the output matches the frequency response.
       */
      TestOutcome steadyState()
      {
        TestStatus success = true;

        const RealT omega0 = 2.0 * std::numbers::pi_v<RealT> * 60.0;
        const RealT phase  = 0.4;

        Fixture fixture;

        // Instantaneous input pair of u_n(t) = U_n cos(omega0 t + phase_n)
        GridKit::EMT::ABCVector<RealT> amplitude{3.7, 2.9, 4.3};
        GridKit::EMT::ABCVector<RealT> u{};
        GridKit::EMT::ABCVector<RealT> u_dot{};
        for (size_t n = 0; n < 3; ++n)
        {
          const RealT phi  = phase + 0.5 * static_cast<RealT>(n);
          u[n]             = amplitude[n] * std::cos(phi);
          u_dot[n]         = -omega0 * amplitude[n] * std::sin(phi);
          fixture.u[n]     = u[n];
          fixture.u_dot[n] = u_dot[n];
        }

        fixture.vf.initializeSteadyState(omega0, u, u_dot);
        fixture.vf.updateTime(0.0, 1.0);
        fixture.evaluateResidual();

        // The state rows are exactly consistent
        const auto* f = fixture.f.getData();
        for (IdxT j = 0; j < Fixture::state_size; ++j)
        {
          success *= isEqual<RealT>(f[j], 0.0, 1.0e-13);
        }

        // The output matches Re[H(j omega0) U] with U the rotating phasor
        // sampled at the initialization instant
        GridKit::EMT::ABCMatrix<RealT> H_re{};
        GridKit::EMT::ABCMatrix<RealT> H_im{};
        fixture.vf.transfer(omega0, H_re, H_im);

        for (size_t n = 0; n < 3; ++n)
        {
          RealT expected = 0.0;
          for (size_t k = 0; k < 3; ++k)
          {
            const RealT u_im  = -u_dot[k] / omega0;
            expected         += H_re[n][k] * u[k] - H_im[n][k] * u_im;
          }
          success *= isEqual<RealT>(fixture.out_row[n], expected, 1.0e-13);
        }

        return success.report(__func__);
      }

      /**
       * @brief Enzyme Jacobian against central finite differences of the
       * assembled residual, J = dF/dy + alpha * dF/dyp, covering the
       * discovery pass, the cached alpha refresh, and the reset path.
       */
      TestOutcome jacobian()
      {
        TestStatus success = true;

        Fixture fixture;
        fixture.setProbeState();

        success *= checkJacobian(fixture, 3.7);  // discovery at alpha = 3.7
        success *= checkJacobian(fixture, -1.3); // cached refresh at a new alpha

        fixture.vf.resetJacobianStructure();
        success *= checkJacobian(fixture, 0.9); // rediscovery after reset

        return success.report(__func__);
      }

    private:
      bool checkJacobian(Fixture& fixture, RealT alpha)
      {
        TestStatus success = true;

        fixture.vf.updateTime(0.0, alpha);
        fixture.vf.evaluateJacobian();

        std::map<std::pair<IdxT, IdxT>, RealT> enzyme_entries;
        auto*                                  coo  = fixture.vf.getCooJacobian();
        success                                    *= (coo != nullptr);
        if (coo != nullptr)
        {
          const IdxT  entry_count = coo->getNnz();
          const auto* rows        = coo->getRowData();
          const auto* cols        = coo->getColData();
          const auto* vals        = coo->getValues();
          for (IdxT i = 0; i < entry_count; ++i)
          {
            enzyme_entries[{rows[i], cols[i]}] += vals[i];
          }
        }

        // Central finite differences over states and the input triple
        const RealT step       = 1.0e-6;
        const RealT zero_floor = 1.0e-9;

        auto* y_data  = fixture.y.getData();
        auto* yp_data = fixture.yp.getData();

        for (IdxT j = 0; j < Fixture::system_size; ++j)
        {
          std::array<RealT, Fixture::system_size> fd_column{};

          auto value_at = [&](IdxT i) -> RealT
          {
            if (i < Fixture::state_size)
            {
              return fixture.f.getData()[i];
            }
            return fixture.out_row[static_cast<size_t>(i - Fixture::state_size)];
          };

          auto perturb = [&](RealT delta, bool derivative)
          {
            if (j < Fixture::state_size)
            {
              if (derivative)
              {
                yp_data[j] += delta;
              }
              else
              {
                y_data[j] += delta;
              }
              return;
            }
            const auto n = static_cast<size_t>(j - Fixture::state_size);
            if (derivative)
            {
              fixture.u_dot[n] += delta;
            }
            else
            {
              fixture.u[n] += delta;
            }
          };

          for (const bool derivative : {false, true})
          {
            perturb(step, derivative);
            fixture.evaluateResidual();
            std::array<RealT, Fixture::system_size> f_plus{};
            for (IdxT i = 0; i < Fixture::system_size; ++i)
            {
              f_plus[static_cast<size_t>(i)] = value_at(i);
            }
            perturb(-2.0 * step, derivative);
            fixture.evaluateResidual();
            RealT weight = 1.0;
            if (derivative)
            {
              weight = alpha;
            }
            for (IdxT i = 0; i < Fixture::system_size; ++i)
            {
              fd_column[static_cast<size_t>(i)] +=
                  weight * (f_plus[static_cast<size_t>(i)] - value_at(i)) / (2.0 * step);
            }
            perturb(step, derivative);
          }

          for (IdxT i = 0; i < Fixture::system_size; ++i)
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

        return static_cast<bool>(success);
      }
    };

  } // namespace Testing
} // namespace GridKit
