#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <type_traits>
#include <utility>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/OptimalPowerFlow/Branch/Branch.hpp>
#include <GridKit/Model/OptimalPowerFlow/Bus/Bus.hpp>
#include <GridKit/Model/OptimalPowerFlow/Generator/Generator.hpp>
#include <GridKit/Model/OptimalPowerFlow/Load/Load.hpp>
#include <GridKit/Model/OptimalPowerFlow/Shunt/Shunt.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    namespace OPF = GridKit::OptimalPowerFlow;

    /**
     * @brief Component derivatives against dependency tracking and closed forms
     */
    class ComponentTests
    {
      using RealT         = double;
      using IdxT          = size_t;
      using VariableT     = DependencyTracking::Variable;
      using DependencyMap = VariableT::DependencyMap;
      using EntryMap      = std::map<std::pair<IdxT, IdxT>, RealT>;

    public:
      /// Power into each bus at V1 = 1 and V2 = j, from the admittances by hand
      TestOutcome branchPower()
      {
        TestStatus success = true;

        OPF::Branch<RealT, IdxT> branch(branchData());
        const auto               x = variables(branch, POWER_X);

        std::vector<RealT> g(branch.sizeConstraints());
        branch.evaluateConstraints(x.data(), g.data());

        // S1 = -3.9 - j3.45 and S2 = -7.6 - j7.3
        success *= isEqual(g[0], 27.1125, POWER_TOL);
        success *= isEqual(g[1], 111.05, POWER_TOL);
        success *= isEqual(g[2], -3.9, POWER_TOL);
        success *= isEqual(g[3], -3.45, POWER_TOL);
        success *= isEqual(g[4], -7.6, POWER_TOL);
        success *= isEqual(g[5], -7.3, POWER_TOL);

        return success.report(__func__);
      }

      /// Enzyme Jacobian entries equal the dependency-tracking derivatives
      TestOutcome branchJacobian()
      {
        TestStatus success = true;

        OPF::Branch<RealT, IdxT> branch(branchData());
        const auto               x = variables(branch, BRANCH_X);
        branch.evaluateJacobian(x.data());

        OPF::Branch<VariableT, IdxT> tracked(branchData());
        const auto                   x_tracked = variables(tracked, BRANCH_X);
        std::vector<VariableT>       g(tracked.sizeConstraints());
        tracked.evaluateConstraints(x_tracked.data(), g.data());

        const auto enzyme = rowsOf(branch.jacobian(), branch.sizeConstraints());
        for (IdxT i = 0; i < g.size(); ++i)
        {
          success *= isEqual(enzyme[i], g[i].getDependencies(), JACOBIAN_TOL);
        }
        success *= branch.jacobian().values.size() == BRANCH_JACOBIAN_SIZE;

        return success.report(__func__);
      }

      /// Enzyme Hessian entries equal the closed-form Cartesian second derivatives
      TestOutcome branchHessian()
      {
        TestStatus success = true;

        OPF::Branch<RealT, IdxT> branch(branchData());
        const auto               x = variables(branch, BRANCH_X);
        branch.evaluateHessian(x.data(), ONE<RealT>, BRANCH_LAMBDA.data());

        // Rows 2 to 5: P1, Q1, P2, Q2
        EntryMap expected;
        addPowerHessian(expected, 0, BUS1_ADMITTANCE, BRANCH_LAMBDA[2], BRANCH_LAMBDA[3]);
        addPowerHessian(expected, 1, BUS2_ADMITTANCE, BRANCH_LAMBDA[4], BRANCH_LAMBDA[5]);

        // Rows 0 and 1: |S|^2 = P^2 + Q^2 by the chain rule over tracked gradients
        OPF::Branch<VariableT, IdxT> tracked(branchData());
        const auto                   x_tracked = variables(tracked, BRANCH_X);
        std::vector<VariableT>       g(tracked.sizeConstraints());
        tracked.evaluateConstraints(x_tracked.data(), g.data());

        for (IdxT terminal = 0; terminal < 2; ++terminal)
        {
          const RealT     weight = TWO<RealT> * BRANCH_LAMBDA[terminal];
          const VariableT p      = g[2 + 2 * terminal];
          const VariableT q      = g[3 + 2 * terminal];
          addOuterProduct(expected, p.getDependencies(), weight);
          addOuterProduct(expected, q.getDependencies(), weight);
        }
        const RealT weight1 = TWO<RealT> * BRANCH_LAMBDA[0];
        const RealT weight2 = TWO<RealT> * BRANCH_LAMBDA[1];
        addPowerHessian(expected, 0, BUS1_ADMITTANCE, weight1 * g[2].getValue(), weight1 * g[3].getValue());
        addPowerHessian(expected, 1, BUS2_ADMITTANCE, weight2 * g[4].getValue(), weight2 * g[5].getValue());

        success *= isEqual(sum(branch.hessian()), expected, HESSIAN_TOL);
        success *= branch.hessian().values.size() == BRANCH_HESSIAN_SIZE;

        return success.report(__func__);
      }

      /// Entry counts are structural, so an out-of-service branch at a flat
      /// start with zero multipliers has the same entries, all zero
      TestOutcome branchPattern()
      {
        TestStatus success = true;

        auto data       = branchData();
        data.in_service = false;

        OPF::Branch<RealT, IdxT> branch(data);
        const auto               x = variables(branch, FLAT_START);

        const std::vector<RealT> lambda(branch.sizeConstraints(), ZERO<RealT>);
        branch.evaluateJacobian(x.data());
        branch.evaluateHessian(x.data(), ONE<RealT>, lambda.data());

        success *= branch.jacobian().values.size() == BRANCH_JACOBIAN_SIZE;
        success *= branch.hessian().values.size() == BRANCH_HESSIAN_SIZE;
        for (const RealT value : branch.jacobian().values)
        {
          success *= isEqual(value, ZERO<RealT>);
        }
        for (const RealT value : branch.hessian().values)
        {
          success *= isEqual(value, ZERO<RealT>);
        }

        return success.report(__func__);
      }

      /// Generator cost, bounds, and derivatives in closed form
      TestOutcome generator()
      {
        TestStatus success = true;

        OPF::Generator<RealT, IdxT> generator(generatorData());
        std::vector<RealT>          x(generator.size());
        numberLocally(generator);
        success *= generator.initialize(x.data()) == 0;
        success *= isEqual(x[0], ZERO<RealT>) && isEqual(x[1], ZERO<RealT>);

        const auto& lower = generator.variableLower();
        const auto& upper = generator.variableUpper();

        success *= isEqual(lower[0], PMIN) && isEqual(upper[0], PMAX);
        success *= isEqual(lower[1], QMIN) && isEqual(upper[1], QMAX);

        x             = GENERATOR_X;
        const RealT p = x[0];

        RealT f = ZERO<RealT>;
        generator.evaluateObjective(x.data(), f);
        success *= isEqual(f, C0 + C1 * p + C2 * p * p);

        std::vector<RealT> gradient(generator.size(), ZERO<RealT>);
        generator.evaluateGradient(x.data(), gradient.data());
        success *= isEqual(gradient[0], C1 + TWO<RealT> * C2 * p);

        generator.evaluateJacobian(x.data());
        success *= isEqual(sum(generator.jacobian()), EntryMap{{{0, 0}, ONE<RealT>}, {{1, 1}, ONE<RealT>}});

        const std::vector<RealT> lambda(generator.sizeConstraints(), ONE<RealT>);
        generator.evaluateHessian(x.data(), SIGMA, lambda.data());
        success *= isEqual(sum(generator.hessian()), EntryMap{{{0, 0}, TWO<RealT> * C2 * SIGMA}});

        return success.report(__func__);
      }

      /// An out-of-service generator is fixed at zero with no cost and the same entries
      TestOutcome outOfServiceGenerator()
      {
        TestStatus success = true;

        auto data       = generatorData();
        data.in_service = false;

        OPF::Generator<RealT, IdxT> generator(data);
        numberLocally(generator);

        const auto& lower = generator.variableLower();
        const auto& upper = generator.variableUpper();

        success *= isEqual(lower[0], ZERO<RealT>) && isEqual(upper[0], ZERO<RealT>);
        success *= isEqual(lower[1], ZERO<RealT>) && isEqual(upper[1], ZERO<RealT>);

        std::vector<RealT> x = GENERATOR_X;
        RealT              f = ZERO<RealT>;
        generator.evaluateObjective(x.data(), f);
        success *= isEqual(f, ZERO<RealT>);

        const std::vector<RealT> lambda(generator.sizeConstraints(), ZERO<RealT>);
        generator.evaluateJacobian(x.data());
        generator.evaluateHessian(x.data(), ONE<RealT>, lambda.data());
        success *= generator.jacobian().values.size() == 2;
        success *= generator.hessian().values.size() == 1;

        return success.report(__func__);
      }

      /// Shunt derivatives in closed form, and no power out of service
      TestOutcome shunt()
      {
        TestStatus success = true;

        OPF::ShuntData<RealT, IdxT> data;
        data.id                                  = "shunt_1";
        data.buses[OPF::ShuntBuses::bus]         = 1;
        data.parameters[OPF::ShuntParameters::G] = SHUNT_G;
        data.parameters[OPF::ShuntParameters::B] = SHUNT_B;

        OPF::Shunt<RealT, IdxT> shunt(data);
        const auto              x  = variables(shunt, SHUNT_X);
        const RealT             vr = x[0];
        const RealT             vi = x[1];

        std::vector<RealT> g(shunt.sizeConstraints(), ZERO<RealT>);
        shunt.evaluateConstraints(x.data(), g.data());
        success *= isEqual(g[0], -SHUNT_G * (vr * vr + vi * vi));
        success *= isEqual(g[1], SHUNT_B * (vr * vr + vi * vi));

        shunt.evaluateJacobian(x.data());
        success *= isEqual(sum(shunt.jacobian()),
                           EntryMap{{{0, 0}, -TWO<RealT> * SHUNT_G * vr},
                                    {{0, 1}, -TWO<RealT> * SHUNT_G * vi},
                                    {{1, 0}, TWO<RealT> * SHUNT_B * vr},
                                    {{1, 1}, TWO<RealT> * SHUNT_B * vi}});

        shunt.evaluateHessian(x.data(), ONE<RealT>, SHUNT_LAMBDA.data());
        const RealT hessian = TWO<RealT> * (-SHUNT_G * SHUNT_LAMBDA[0] + SHUNT_B * SHUNT_LAMBDA[1]);

        success *= isEqual(sum(shunt.hessian()), EntryMap{{{0, 0}, hessian}, {{1, 1}, hessian}});

        data.in_service = false;
        OPF::Shunt<RealT, IdxT> out_of_service(data);
        const auto              x_out = variables(out_of_service, SHUNT_X);
        std::vector<RealT>      g_out(out_of_service.sizeConstraints(), ZERO<RealT>);
        out_of_service.evaluateConstraints(x_out.data(), g_out.data());
        success *= isEqual(g_out[0], ZERO<RealT>) && isEqual(g_out[1], ZERO<RealT>);

        return success.report(__func__);
      }

      /// Bus rows follow the limits, and a reference bus fixes its imaginary voltage at zero
      TestOutcome bus()
      {
        TestStatus success = true;

        OPF::BusData<RealT, IdxT> data;
        data.number                               = 3;
        data.parameters[OPF::BusParameters::Vmin] = VMIN;
        data.parameters[OPF::BusParameters::Vmax] = VMAX;

        OPF::Bus<RealT, IdxT> bus(data);
        std::vector<RealT>    x(bus.size());
        numberLocally(bus);

        bus.initialize(x.data());
        success *= isEqual(x[0], ONE<RealT>) && isEqual(x[1], ZERO<RealT>);
        for (IdxT j = 0; j < bus.sizeInternalVariables(); ++j)
        {
          success *= std::isinf(bus.variableLower()[j]) && std::isinf(bus.variableUpper()[j]);
        }

        const auto& lower = bus.constraintLower();
        const auto& upper = bus.constraintUpper();

        success *= isEqual(lower[0], ZERO<RealT>) && isEqual(upper[0], ZERO<RealT>);
        success *= isEqual(lower[1], ZERO<RealT>) && isEqual(upper[1], ZERO<RealT>);
        success *= isEqual(lower[2], VMIN * VMIN) && isEqual(upper[2], VMAX * VMAX);

        x = BUS_X;
        std::vector<RealT> g(bus.sizeConstraints());
        bus.evaluateConstraints(x.data(), g.data());
        success *= isEqual(g[2], x[0] * x[0] + x[1] * x[1]);

        bus.evaluateJacobian(x.data());
        success *= isEqual(sum(bus.jacobian()), EntryMap{{{2, 0}, TWO<RealT> * x[0]}, {{2, 1}, TWO<RealT> * x[1]}});

        bus.evaluateHessian(x.data(), ONE<RealT>, BUS_LAMBDA.data());
        const RealT hessian = TWO<RealT> * BUS_LAMBDA[2];

        success *= isEqual(sum(bus.hessian()), EntryMap{{{0, 0}, hessian}, {{1, 1}, hessian}});

        data.reference = true;
        OPF::Bus<RealT, IdxT> reference(data);
        success *= std::isinf(reference.variableLower()[0]) && std::isinf(reference.variableUpper()[0]);
        success *= isEqual(reference.variableLower()[1], ZERO<RealT>) && isEqual(reference.variableUpper()[1], ZERO<RealT>);

        return success.report(__func__);
      }

      /// Load demand comes from its parameters, with no entries, and is zero
      /// out of service
      TestOutcome load()
      {
        TestStatus success = true;

        OPF::LoadData<RealT, IdxT> data;
        data.id                                 = "load_2";
        data.buses[OPF::LoadBuses::bus]         = 2;
        data.parameters[OPF::LoadParameters::P] = LOAD_P;
        data.parameters[OPF::LoadParameters::Q] = LOAD_Q;

        OPF::Load<RealT, IdxT> load(data);
        const auto             x = variables(load, LOAD_X);

        std::vector<RealT> g(load.sizeConstraints(), ZERO<RealT>);
        load.evaluateConstraints(x.data(), g.data());
        success *= isEqual(g[0], -LOAD_P) && isEqual(g[1], -LOAD_Q);

        const std::vector<RealT> lambda(load.sizeConstraints(), ONE<RealT>);
        load.evaluateJacobian(x.data());
        load.evaluateHessian(x.data(), ONE<RealT>, lambda.data());
        success *= load.jacobian().values.empty() && load.hessian().values.empty();

        data.in_service = false;
        OPF::Load<RealT, IdxT> out_of_service(data);
        const auto             x_out = variables(out_of_service, LOAD_X);
        std::vector<RealT>     g_out(out_of_service.sizeConstraints(), ZERO<RealT>);
        out_of_service.evaluateConstraints(x_out.data(), g_out.data());
        success *= isEqual(g_out[0], ZERO<RealT>) && isEqual(g_out[1], ZERO<RealT>);

        return success.report(__func__);
      }

    private:
      static constexpr const char* BRANCH_ID = "branch_1_2";

      static constexpr RealT R     = 0.1;
      static constexpr RealT X     = 0.1;
      static constexpr RealT G     = 0.2;
      static constexpr RealT B     = 0.4;
      static constexpr RealT GMAG  = 0.05;
      static constexpr RealT BMAG  = 0.1;
      static constexpr RealT TAP   = 2.0;
      static constexpr RealT PHASE = std::numbers::pi / 2.0;
      static constexpr RealT SMAX  = 1.2;

      /// Admittances g_k1, b_k1, g_k2, b_k2 of the current into bus k: series
      /// admittance 5 - j5, tap 2 at 90 degrees, and the shunts
      static constexpr std::array<RealT, 4> BUS1_ADMITTANCE = {-1.4, 0.95, 2.5, 2.5};
      static constexpr std::array<RealT, 4> BUS2_ADMITTANCE = {-2.5, -2.5, -5.1, 4.8};

      static constexpr IdxT BRANCH_JACOBIAN_SIZE = 24;
      static constexpr IdxT BRANCH_HESSIAN_SIZE  = 10;

      static constexpr RealT PMIN  = 0.1;
      static constexpr RealT PMAX  = 0.9;
      static constexpr RealT QMIN  = -0.4;
      static constexpr RealT QMAX  = 0.4;
      static constexpr RealT C0    = 3.0;
      static constexpr RealT C1    = 11.0;
      static constexpr RealT C2    = 0.7;
      static constexpr RealT SIGMA = 1.7;

      static constexpr RealT SHUNT_G = 0.05;
      static constexpr RealT SHUNT_B = 0.3;

      static constexpr RealT VMIN = 0.9;
      static constexpr RealT VMAX = 1.1;

      static constexpr RealT LOAD_P = 0.6;
      static constexpr RealT LOAD_Q = 0.2;

      static inline const std::vector<RealT> POWER_X       = {1.0, 0.0, 0.0, 1.0};
      static inline const std::vector<RealT> BRANCH_X      = {1.02, 0.12, 0.96, -0.08};
      static inline const std::vector<RealT> BRANCH_LAMBDA = {0.3, -0.2, 1.5, -0.7, 2.1, 0.4};
      static inline const std::vector<RealT> FLAT_START    = {1.0, 0.0, 1.0, 0.0};
      static inline const std::vector<RealT> GENERATOR_X   = {0.3, 0.1, 1.01, 0.05};
      static inline const std::vector<RealT> SHUNT_X       = {0.98, 0.1};
      static inline const std::vector<RealT> SHUNT_LAMBDA  = {1.3, -0.6};
      static inline const std::vector<RealT> BUS_X         = {1.02, -0.05};
      static inline const std::vector<RealT> BUS_LAMBDA    = {0.7, -0.4, 1.3};
      static inline const std::vector<RealT> LOAD_X        = {1.01, -0.05};

      // Allow roundoff in the independent power and derivative references.
      static constexpr RealT POWER_TOL    = 2 * std::numeric_limits<RealT>::epsilon();
      static constexpr RealT JACOBIAN_TOL = 2 * std::numeric_limits<RealT>::epsilon();
      static constexpr RealT HESSIAN_TOL  = 9 * std::numeric_limits<RealT>::epsilon();

      static OPF::BranchData<RealT, IdxT> branchData()
      {
        OPF::BranchData<RealT, IdxT> data;
        data.id                                       = BRANCH_ID;
        data.buses[OPF::BranchBuses::bus1]            = 1;
        data.buses[OPF::BranchBuses::bus2]            = 2;
        data.parameters[OPF::BranchParameters::R]     = R;
        data.parameters[OPF::BranchParameters::X]     = X;
        data.parameters[OPF::BranchParameters::G]     = G;
        data.parameters[OPF::BranchParameters::B]     = B;
        data.parameters[OPF::BranchParameters::Gmag]  = GMAG;
        data.parameters[OPF::BranchParameters::Bmag]  = BMAG;
        data.parameters[OPF::BranchParameters::tap]   = TAP;
        data.parameters[OPF::BranchParameters::phase] = PHASE;
        data.parameters[OPF::BranchParameters::Smax]  = SMAX;
        return data;
      }

      static OPF::GeneratorData<RealT, IdxT> generatorData()
      {
        OPF::GeneratorData<RealT, IdxT> data;
        data.id                                         = "gen_1";
        data.buses[OPF::GeneratorBuses::bus]            = 1;
        data.parameters[OPF::GeneratorParameters::Pmin] = PMIN;
        data.parameters[OPF::GeneratorParameters::Pmax] = PMAX;
        data.parameters[OPF::GeneratorParameters::Qmin] = QMIN;
        data.parameters[OPF::GeneratorParameters::Qmax] = QMAX;
        data.parameters[OPF::GeneratorParameters::c0]   = C0;
        data.parameters[OPF::GeneratorParameters::c1]   = C1;
        data.parameters[OPF::GeneratorParameters::c2]   = C2;
        return data;
      }

      /// Local indices are global indices
      template <typename ComponentT>
      static void numberLocally(ComponentT& component)
      {
        for (IdxT j = 0; j < component.size(); ++j)
        {
          component.variableIndices()[j] = j;
        }
        for (IdxT i = 0; i < component.sizeConstraints(); ++i)
        {
          component.constraintIndices()[i] = i;
        }
      }

      /// Initialized component variables set to `values`, tracked by index
      template <typename ComponentT>
      static std::vector<typename ComponentT::ScalarT> variables(ComponentT& component, const std::vector<RealT>& values)
      {
        using ScalarT = typename ComponentT::ScalarT;

        numberLocally(component);
        std::vector<ScalarT> x(values.size());
        component.initialize(x.data());

        for (IdxT j = 0; j < values.size(); ++j)
        {
          x[j] = values[j];
          if constexpr (std::is_same_v<ScalarT, VariableT>)
          {
            x[j].setVariableNumber(j);
          }
        }
        return x;
      }

      /// Entries summed by coordinate
      static EntryMap sum(const Enzyme::Sparse::CooEntries& entries)
      {
        EntryMap result;
        for (IdxT k = 0; k < entries.values.size(); ++k)
        {
          result[{entries.rows[k], entries.cols[k]}] += entries.values[k];
        }
        return result;
      }

      /// Rows of the summed entries
      static std::vector<DependencyMap> rowsOf(const Enzyme::Sparse::CooEntries& entries, IdxT count)
      {
        std::vector<DependencyMap> rows(count);
        for (const auto& [coordinate, value] : sum(entries))
        {
          rows[coordinate.first][coordinate.second] = value;
        }
        return rows;
      }

      /// Coordinate maps with the same keys and values within `tol`
      static bool isEqual(const EntryMap& actual, const EntryMap& expected, RealT tol = std::numeric_limits<RealT>::epsilon())
      {
        bool equal = actual.size() == expected.size();
        for (const auto& [coordinate, value] : expected)
        {
          const auto entry = actual.find(coordinate);
          if (entry == actual.end())
          {
            std::cerr << "Entry (" << coordinate.first << ", " << coordinate.second << ") is missing\n";
            equal = false;
          }
          else
          {
            equal = isEqual(entry->second, value, tol) && equal;
          }
        }
        return equal;
      }

      /// Values within `tol` relative to one plus the reference magnitude
      static bool isEqual(RealT value, RealT reference, RealT tol = std::numeric_limits<RealT>::epsilon())
      {
        if (Testing::isEqual(value, reference, tol))
        {
          return true;
        }
        std::cerr << std::setprecision(std::numeric_limits<RealT>::max_digits10) << value << " differs from "
                  << reference << " by " << std::abs(value - reference) / (ONE<RealT> + std::abs(reference)) / std::numeric_limits<RealT>::epsilon()
                  << " eps\n";
        return false;
      }

      /// Dependency maps with the same keys and values within `tol`
      static bool isEqual(const DependencyMap& actual, const DependencyMap& expected, RealT tol)
      {
        bool equal = actual.size() == expected.size();
        for (const auto& [variable, value] : expected)
        {
          const auto entry = actual.find(variable);
          if (entry == actual.end())
          {
            std::cerr << "Derivative with respect to " << variable << " is missing\n";
            equal = false;
          }
          else
          {
            equal = isEqual(entry->second, value, tol) && equal;
          }
        }
        return equal;
      }

      /// Add `value` at the lower-triangle coordinate of (i, j)
      static void addEntry(EntryMap& hessian, IdxT i, IdxT j, RealT value)
      {
        hessian[{std::max(i, j), std::min(i, j)}] += value;
      }

      /// Add `weight` times the lower triangle of the outer product of `gradient`
      static void addOuterProduct(EntryMap& hessian, const DependencyMap& gradient, RealT weight)
      {
        for (const auto& [i, di] : gradient)
        {
          for (const auto& [j, dj] : gradient)
          {
            if (i >= j)
            {
              hessian[{i, j}] += weight * di * dj;
            }
          }
        }
      }

      /**
       * @brief Add the Hessian of \f$w_P P_k + w_Q Q_k\f$ for the current
       * \f$I_k = (g_{k1} + j b_{k1}) V_1 + (g_{k2} + j b_{k2}) V_2\f$
       *
       * \f$P_k + j Q_k = V_k I_k^*\f$ is bilinear in the voltages, so each
       * admittance adds a constant block. Local indices \f$2m\f$ and
       * \f$2m + 1\f$ hold the real and imaginary voltage of bus \f$m + 1\f$.
       *
       * @param[in] admittance - \f$g_{k1}, b_{k1}, g_{k2}, b_{k2}\f$
       */
      static void addPowerHessian(EntryMap&                   hessian,
                                  IdxT                        k,
                                  const std::array<RealT, 4>& admittance,
                                  RealT                       weight_p,
                                  RealT                       weight_q)
      {
        for (IdxT m = 0; m < 2; ++m)
        {
          const RealT g = admittance[2 * m];
          const RealT b = admittance[2 * m + 1];

          // Real and imaginary parts of the weighted admittance (w_P + j w_Q)(g + j b)
          const RealT re = weight_p * g - weight_q * b;
          const RealT im = weight_p * b + weight_q * g;

          addProduct(hessian, 2 * k, 2 * m, re);
          addProduct(hessian, 2 * k, 2 * m + 1, -im);
          addProduct(hessian, 2 * k + 1, 2 * m, im);
          addProduct(hessian, 2 * k + 1, 2 * m + 1, re);
        }
      }

      /// Add the Hessian of \f$c x_i x_j\f$
      static void addProduct(EntryMap& hessian, IdxT i, IdxT j, RealT c)
      {
        RealT scale = ONE<RealT>;
        if (i == j)
        {
          scale = TWO<RealT>;
        }
        addEntry(hessian, i, j, scale * c);
      }
    };
  } // namespace Testing
} // namespace GridKit
