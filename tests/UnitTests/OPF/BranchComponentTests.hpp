#pragma once

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>

#include <GridKit/Model/OPF/Branch/Branch.hpp>
#include <GridKit/Model/StateData.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    class OPFBranchComponentTests
    {
    public:
      using IdxT    = std::size_t;
      using RealT   = double;
      using BranchT = OPF::Branch<RealT, IdxT>;

      TestOutcome terminalPowerAndDerivativeStructures()
      {
        auto                      data  = branchData();
        auto                      state = branchState(false);
        const std::array<IdxT, 6> rated_rows{{10, 11, 12, 13, 20, 21}};
        BranchT                   rated(data, state, {0, 1, 2, 3}, rated_rows);

        const auto point    = branchPoint();
        TestStatus success  = rated.gatherVariables(point.data()) == 0;
        success            *= rated.evaluateConstraints() == 0;

        const auto expected  = complexReference(data, state, point);
        success             *= rated.constraintValues().size() == 6;
        for (std::size_t row = 0; row < expected.size(); ++row)
        {
          success *= near(rated.constraintValues()[row], expected[row]);
        }

        success *= rated.hasJacobian() && rated.hasHessian();
        success *= rated.jacobianEntries().size() == 24;
        for (IdxT variable = 0; variable < 4; ++variable)
        {
          for (IdxT constraint = 0; constraint < 6; ++constraint)
          {
            const auto entry  = rated.jacobianEntries()[6 * variable + constraint];
            success          *= entry.variable == variable;
            success          *= entry.constraint == constraint;
          }
        }
        success *= validFullLowerTriangle(rated.hessianEntries());

        data.smax.reset();
        const std::array<IdxT, 4> unrated_rows{{10, 11, 12, 13}};
        BranchT                   unrated(data, state, {0, 1, 2, 3}, unrated_rows);
        success *= unrated.gatherVariables(point.data()) == 0;
        success *= unrated.evaluateConstraints() == 0;
        success *= unrated.constraintIndices().size() == 4;
        success *= unrated.constraintValues().size() == 4;
        success *= unrated.jacobianEntries().size() == 16;
        for (IdxT variable = 0; variable < 4; ++variable)
        {
          for (IdxT constraint = 0; constraint < 4; ++constraint)
          {
            const auto entry =
                unrated.jacobianEntries()[4 * variable + constraint];
            success *= entry.variable == variable;
            success *= entry.constraint == constraint;
          }
        }
        success *= validFullLowerTriangle(unrated.hessianEntries());

        return success.report(__func__);
      }

      TestOutcome exactJacobianAndHessian()
      {
        auto                      data  = branchData();
        auto                      state = branchState(false);
        const std::array<IdxT, 6> rows{{0, 1, 2, 3, 4, 5}};
        BranchT                   rated(data, state, {0, 1, 2, 3}, rows);
        const auto                point = branchPoint();

        TestStatus success  = rated.gatherVariables(point.data()) == 0;
        success            *= rated.evaluateJacobian() == 0;
        success            *= same(rated.jacobianValues(), ratedJacobianGolden());

        const std::array<RealT, 6> thermal_multipliers{{0.3, -0.2, 0.1, 0.4, 0.5, -0.7}};
        success *= rated.evaluateHessian(1.3, thermal_multipliers) == 0;
        success *= same(rated.hessianValues(), thermalHessianGolden());

        const std::array<RealT, 6> flow_multipliers{{-0.6, 0.25, 0.45, -0.15, 0.0, 0.0}};
        success *= rated.evaluateHessian(0.0, flow_multipliers) == 0;
        success *= same(rated.hessianValues(), flowHessianGolden());

        data.smax.reset();
        const std::array<IdxT, 4> unrated_rows{{0, 1, 2, 3}};
        BranchT                   unrated(data, state, {0, 1, 2, 3}, unrated_rows);
        success *= unrated.gatherVariables(point.data()) == 0;
        success *= unrated.evaluateJacobian() == 0;
        success *= same(unrated.jacobianValues(), unratedJacobianGolden());
        const std::array<RealT, 4> unrated_multipliers{{-0.6, 0.25, 0.45, -0.15}};
        success *= unrated.evaluateHessian(2.0, unrated_multipliers) == 0;
        success *= same(unrated.hessianValues(), flowHessianGolden());

        state.open = true;
        BranchT open(data, state, {0, 1, 2, 3}, unrated_rows);
        success *= open.gatherVariables(point.data()) == 0;
        success *= open.evaluateConstraints() == 0;
        success *= open.evaluateJacobian() == 0;
        success *= open.evaluateHessian(1.0, unrated_multipliers) == 0;
        success *= allNear(open.constraintValues(), 0.0);
        success *= allNear(open.jacobianValues(), 0.0);
        success *= allNear(open.hessianValues(), 0.0);
        success *= open.jacobianEntries().size()
                   == unrated.jacobianEntries().size();
        success *= open.hessianEntries().size()
                   == unrated.hessianEntries().size();
        success *= open.hasJacobian() && open.hasHessian();

        return success.report(__func__);
      }

    private:
      static OPF::BranchData<RealT, IdxT> branchData()
      {
        OPF::BranchData<RealT, IdxT> data;
        data.id   = "BR";
        data.from = 1;
        data.to   = 2;
        data.R    = 0.02;
        data.X    = 0.17;
        data.G    = 0.006;
        data.B    = 0.04;
        data.smax = 1.7;
        return data;
      }

      static Model::DeviceState branchState(bool open)
      {
        Model::DeviceState state;
        state.open  = open;
        state.tap   = 1.08;
        state.phase = 0.13;
        return state;
      }

      static std::array<RealT, 4> branchPoint()
      {
        return {1.04, 0.21, 0.97, -0.08};
      }

      static std::array<RealT, 6> complexReference(
          const OPF::BranchData<RealT, IdxT>& data,
          const Model::DeviceState&           state,
          const std::array<RealT, 4>&         point)
      {
        using ComplexT          = std::complex<RealT>;
        const ComplexT series_y = RealT{1} / ComplexT(data.R, data.X);
        const ComplexT shunt_y(data.G, data.B);
        const ComplexT tap = std::polar(state.tap.value_or(1.0),
                                        state.phase.value_or(0.0));

        const ComplexT yff = (series_y + RealT{0.5} * shunt_y)
                             / std::norm(tap);
        const ComplexT yft = -series_y / std::conj(tap);
        const ComplexT ytf = -series_y / tap;
        const ComplexT ytt = series_y + RealT{0.5} * shunt_y;
        const ComplexT vf  = std::polar(point[0], point[1]);
        const ComplexT vt  = std::polar(point[2], point[3]);
        const ComplexT sf  = vf * std::conj(yff * vf + yft * vt);
        const ComplexT st  = vt * std::conj(ytf * vf + ytt * vt);

        return {-sf.real(),
                -sf.imag(),
                -st.real(),
                -st.imag(),
                std::norm(sf),
                std::norm(st)};
      }

      static std::array<RealT, 24> ratedJacobianGolden()
      {
        return {
            -1.44757838391509935,
            -5.06871010548153076,
            1.43546269650411838,
            5.04686530844155329,
            1.60340297367896322,
            0.500141505932413902,
            -5.45190033897217627,
            -0.233982350275884098,
            5.24873992077921542,
            -1.49288120436428312,
            9.44116345671534697,
            9.47186456053162201,
            -0.241218917810189792,
            5.62051581337337760,
            0.209000706093152789,
            -5.80610061296819212,
            1.43376774699602875,
            2.57915220884384983,
            5.45190033897217627,
            0.233982350275884098,
            -5.24873992077921542,
            1.49288120436428312,
            -9.44116345671534697,
            -9.47186456053162201};
      }

      static std::array<RealT, 16> unratedJacobianGolden()
      {
        const auto            rated = ratedJacobianGolden();
        std::array<RealT, 16> unrated{};
        for (std::size_t variable = 0; variable < 4; ++variable)
        {
          for (std::size_t constraint = 0; constraint < 4; ++constraint)
          {
            unrated[4 * variable + constraint] =
                rated[6 * variable + constraint];
          }
        }
        return unrated;
      }

      static std::array<RealT, 10> thermalHessianGolden()
      {
        return {
            -8.99850909063832667,
            5.64460372249674739,
            13.8711522178072628,
            -5.64460372249674739,
            -13.3306696753754134,
            -17.3519442990340211,
            13.3306696753754134,
            -21.8449230280343184,
            17.3519442990340211,
            -13.3306696753754134};
      }

      static std::array<RealT, 10> flowHessianGolden()
      {
        return {
            -1.77324353806467319,
            5.57549015463425398,
            1.37574350590366734,
            -5.57549015463425398,
            -1.38785004875561962,
            5.97784511424703520,
            1.38785004875561962,
            1.11757986348122867,
            -5.97784511424703520,
            -1.38785004875561962};
      }

      template <typename EntrySpanT>
      static bool validFullLowerTriangle(const EntrySpanT& entries)
      {
        if (entries.size() != 10)
        {
          return false;
        }
        std::size_t entry = 0;
        for (IdxT column = 0; column < 4; ++column)
        {
          for (IdxT row = column; row < 4; ++row)
          {
            if (entries[entry].row != row
                || entries[entry].column != column)
            {
              return false;
            }
            ++entry;
          }
        }
        return true;
      }

      static bool near(RealT value, RealT expected, RealT tolerance = 2.0e-10)
      {
        return std::abs(value - expected)
               <= tolerance * (RealT{1} + std::abs(expected));
      }

      template <typename SpanT, typename ArrayT>
      static bool same(const SpanT&  values,
                       const ArrayT& expected,
                       RealT         tolerance = 2.0e-10)
      {
        if (values.size() != expected.size())
        {
          return false;
        }
        for (std::size_t entry = 0; entry < values.size(); ++entry)
        {
          if (!near(values[entry], expected[entry], tolerance))
          {
            return false;
          }
        }
        return true;
      }

      template <typename SpanT>
      static bool allNear(const SpanT& values,
                          RealT        expected,
                          RealT        tolerance = 2.0e-10)
      {
        for (const RealT value : values)
        {
          if (!near(value, expected, tolerance))
          {
            return false;
          }
        }
        return true;
      }
    };
  } // namespace Testing
} // namespace GridKit
