#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <type_traits>
#include <utility>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/OptimalPowerFlow/Branch/Branch.hpp>
#include <GridKit/Model/OptimalPowerFlow/Bus/Bus.hpp>
#include <GridKit/Model/OptimalPowerFlow/Generator/Generator.hpp>
#include <GridKit/Model/OptimalPowerFlow/Load/Load.hpp>
#include <GridKit/Model/OptimalPowerFlow/Shunt/Shunt.hpp>
#include <GridKit/Model/PhasorDynamics/Branch/Branch.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/StateData.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    namespace OPF = GridKit::OptimalPowerFlow;
    namespace PD  = GridKit::PhasorDynamics;

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
      /// Power into each bus is V conj(I) with I from `PhasorDynamics::Branch`
      TestOutcome branchPower()
      {
        TestStatus success = true;

        OPF::Branch<RealT, IdxT> branch(opfBranchData());
        const auto               x = variables(branch, Model::StateData{}, BRANCH_X);

        std::vector<RealT> g(branch.sizeConstraints());
        branch.evaluateConstraints(x.data(), g.data());

        const RealT vr1 = BRANCH_X[0];
        const RealT vi1 = BRANCH_X[1];
        const RealT vr2 = BRANCH_X[2];
        const RealT vi2 = BRANCH_X[3];

        RealT ir1 = 0.0;
        RealT ii1 = 0.0;
        RealT ir2 = 0.0;
        RealT ii2 = 0.0;
        pdCurrents(vr1, vi1, vr2, vi2, ir1, ii1, ir2, ii2);

        const RealT p1 = vr1 * ir1 + vi1 * ii1;
        const RealT q1 = vi1 * ir1 - vr1 * ii1;
        const RealT p2 = vr2 * ir2 + vi2 * ii2;
        const RealT q2 = vi2 * ir2 - vr2 * ii2;

        success *= isEqual(g[0], p1 * p1 + q1 * q1, POWER_TOL);
        success *= isEqual(g[1], p2 * p2 + q2 * q2, POWER_TOL);
        success *= isEqual(g[2], p1, POWER_TOL);
        success *= isEqual(g[3], q1, POWER_TOL);
        success *= isEqual(g[4], p2, POWER_TOL);
        success *= isEqual(g[5], q2, POWER_TOL);

        return success.report(__func__);
      }

      /// Enzyme Jacobian entries equal the dependency-tracking derivatives
      TestOutcome branchJacobian()
      {
        TestStatus success = true;

        OPF::Branch<RealT, IdxT> branch(opfBranchData());
        const auto               x = variables(branch, Model::StateData{}, BRANCH_X);
        branch.evaluateJacobian(x.data());

        OPF::Branch<VariableT, IdxT> tracked(opfBranchData());
        const auto                   x_tracked = variables(tracked, Model::StateData{}, BRANCH_X);
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

        OPF::Branch<RealT, IdxT> branch(opfBranchData());
        const auto               x = variables(branch, Model::StateData{}, BRANCH_X);
        branch.evaluateHessian(x.data(), ONE<RealT>, BRANCH_LAMBDA.data());

        // Admittances of PhasorDynamics::Branch from unit voltages
        RealT g11 = 0.0;
        RealT b11 = 0.0;
        RealT g12 = 0.0;
        RealT b12 = 0.0;
        RealT g21 = 0.0;
        RealT b21 = 0.0;
        RealT g22 = 0.0;
        RealT b22 = 0.0;
        pdCurrents(ONE<RealT>, ZERO<RealT>, ZERO<RealT>, ZERO<RealT>, g11, b11, g21, b21);
        pdCurrents(ZERO<RealT>, ZERO<RealT>, ONE<RealT>, ZERO<RealT>, g12, b12, g22, b22);

        // Rows 2 to 5: P1, Q1, P2, Q2
        EntryMap expected;
        addPowerHessian(expected, 0, {g11, b11, g12, b12}, BRANCH_LAMBDA[2], BRANCH_LAMBDA[3]);
        addPowerHessian(expected, 1, {g21, b21, g22, b22}, BRANCH_LAMBDA[4], BRANCH_LAMBDA[5]);

        // Rows 0 and 1: |S|^2 = P^2 + Q^2 by the chain rule over tracked gradients
        OPF::Branch<VariableT, IdxT> tracked(opfBranchData());
        const auto                   x_tracked = variables(tracked, Model::StateData{}, BRANCH_X);
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
        addPowerHessian(expected, 0, {g11, b11, g12, b12}, weight1 * g[2].getValue(), weight1 * g[3].getValue());
        addPowerHessian(expected, 1, {g21, b21, g22, b22}, weight2 * g[4].getValue(), weight2 * g[5].getValue());

        success *= isEqual(sum(branch.hessian()), expected, HESSIAN_TOL);
        success *= branch.hessian().values.size() == BRANCH_HESSIAN_SIZE;

        return success.report(__func__);
      }

      /// Entry counts are structural, so an open branch at a flat start with
      /// zero multipliers has the same entries
      TestOutcome branchPattern()
      {
        TestStatus success = true;

        Model::StateData state;
        state.devices[BRANCH_ID].flags["open"] = true;

        OPF::Branch<RealT, IdxT> branch(opfBranchData());
        const auto               x = variables(branch, state, FLAT_START);

        const std::vector<RealT> lambda(branch.sizeConstraints(), ZERO<RealT>);
        branch.evaluateJacobian(x.data());
        branch.evaluateHessian(x.data(), ONE<RealT>, lambda.data());

        success *= branch.jacobian().values.size() == BRANCH_JACOBIAN_SIZE;
        success *= branch.hessian().values.size() == BRANCH_HESSIAN_SIZE;
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

        Model::StateData state;
        state.buses[Model::busKey(1)].values = {{"vr", 1.0}, {"vi", 0.0}};

        OPF::Generator<RealT, IdxT> generator(generatorData());
        std::vector<RealT>          x(generator.size());
        std::vector<RealT>          x_lower(generator.size());
        std::vector<RealT>          x_upper(generator.size());
        numberLocally(generator);
        success *= generator.initialize(state, x.data(), x_lower.data(), x_upper.data()) == 0;

        success *= isEqual(x_lower[0], PMIN) && isEqual(x_upper[0], PMAX);
        success *= isEqual(x_lower[1], QMIN) && isEqual(x_upper[1], QMAX);

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

      /// An offline generator is fixed at zero with no cost and the same entries
      TestOutcome offlineGenerator()
      {
        TestStatus success = true;

        Model::StateData state;
        state.devices["gen_1"].flags["online"] = false;

        OPF::Generator<RealT, IdxT> generator(generatorData());
        std::vector<RealT>          x(generator.size());
        std::vector<RealT>          x_lower(generator.size());
        std::vector<RealT>          x_upper(generator.size());
        numberLocally(generator);
        generator.initialize(state, x.data(), x_lower.data(), x_upper.data());

        success *= isEqual(x_lower[0], ZERO<RealT>) && isEqual(x_upper[0], ZERO<RealT>);
        success *= isEqual(x_lower[1], ZERO<RealT>) && isEqual(x_upper[1], ZERO<RealT>);

        x       = GENERATOR_X;
        RealT f = ZERO<RealT>;
        generator.evaluateObjective(x.data(), f);
        success *= isEqual(f, ZERO<RealT>);

        const std::vector<RealT> lambda(generator.sizeConstraints(), ZERO<RealT>);
        generator.evaluateJacobian(x.data());
        generator.evaluateHessian(x.data(), ONE<RealT>, lambda.data());
        success *= generator.jacobian().values.size() == 2;
        success *= generator.hessian().values.size() == 1;

        return success.report(__func__);
      }

      /// Shunt derivatives in closed form
      TestOutcome shunt()
      {
        TestStatus success = true;

        OPF::ShuntData<RealT, IdxT> data;
        data.id                                  = "shunt_1";
        data.buses[OPF::ShuntBuses::bus]         = 1;
        data.parameters[OPF::ShuntParameters::G] = SHUNT_G;
        data.parameters[OPF::ShuntParameters::B] = SHUNT_B;

        OPF::Shunt<RealT, IdxT> shunt(data);
        const auto              x  = variables(shunt, Model::StateData{}, SHUNT_X);
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

        return success.report(__func__);
      }

      /// Bus rows follow the limits and the reference and infinite settings
      TestOutcome bus()
      {
        TestStatus success = true;

        Model::StateData state;
        state.buses[Model::busKey(3)].values = {{"vr", BUS_VR}, {"vi", BUS_VI}};

        OPF::BusData<RealT, IdxT> data;
        data.number                               = 3;
        data.parameters[OPF::BusParameters::Vmin] = VMIN;
        data.parameters[OPF::BusParameters::Vmax] = VMAX;

        OPF::Bus<RealT, IdxT> bus(data);
        std::vector<RealT>    x(bus.size());
        std::vector<RealT>    x_lower(bus.size());
        std::vector<RealT>    x_upper(bus.size());
        numberLocally(bus);

        bus.initialize(state, x.data(), x_lower.data(), x_upper.data());
        success *= isEqual(x[0], BUS_VR) && isEqual(x[1], BUS_VI);
        for (IdxT j = 0; j < bus.size(); ++j)
        {
          success *= std::isinf(x_lower[j]) && std::isinf(x_upper[j]);
        }

        const auto& lower = bus.constraintLower();
        const auto& upper = bus.constraintUpper();

        success *= isEqual(lower[0], ZERO<RealT>) && isEqual(upper[0], ZERO<RealT>);
        success *= isEqual(lower[1], ZERO<RealT>) && isEqual(upper[1], ZERO<RealT>);
        success *= isEqual(lower[2], VMIN * VMIN) && isEqual(upper[2], VMAX * VMAX);
        success *= std::isinf(lower[3]) && std::isinf(upper[3]);

        bus.setReference();
        success *= isEqual(lower[3], ZERO<RealT>) && isEqual(upper[3], ZERO<RealT>);

        x = BUS_X;
        std::vector<RealT> g(bus.sizeConstraints());
        bus.evaluateConstraints(x.data(), g.data());
        success *= isEqual(g[2], x[0] * x[0] + x[1] * x[1]);
        success *= isEqual(g[3], BUS_VR * x[1] - BUS_VI * x[0]);

        bus.evaluateJacobian(x.data());
        success *= isEqual(sum(bus.jacobian()),
                           EntryMap{{{2, 0}, TWO<RealT> * x[0]}, {{2, 1}, TWO<RealT> * x[1]}, {{3, 0}, -BUS_VI}, {{3, 1}, BUS_VR}});

        bus.evaluateHessian(x.data(), ONE<RealT>, BUS_LAMBDA.data());
        const RealT hessian = TWO<RealT> * BUS_LAMBDA[2];

        success *= isEqual(sum(bus.hessian()), EntryMap{{{0, 0}, hessian}, {{1, 1}, hessian}});

        data.infinite = true;
        OPF::Bus<RealT, IdxT> infinite(data);
        numberLocally(infinite);
        infinite.initialize(state, x.data(), x_lower.data(), x_upper.data());
        success *= isEqual(x_lower[0], BUS_VR) && isEqual(x_upper[0], BUS_VR);
        success *= isEqual(x_lower[1], BUS_VI) && isEqual(x_upper[1], BUS_VI);
        for (IdxT i = 0; i < infinite.sizeConstraints(); ++i)
        {
          success *= std::isinf(infinite.constraintLower()[i]) && std::isinf(infinite.constraintUpper()[i]);
        }

        return success.report(__func__);
      }

      /// Load demand comes from the state, with no entries
      TestOutcome load()
      {
        TestStatus success = true;

        Model::StateData state;
        state.buses[Model::busKey(2)].values = {{"vr", 1.01}, {"vi", -0.05}};
        Model::setTerminalCurrent(state, "load_2", 2, 0, 1, -0.6, -0.2);

        OPF::LoadData<RealT, IdxT> data;
        data.id                         = "load_2";
        data.buses[OPF::LoadBuses::bus] = 2;

        OPF::Load<RealT, IdxT> load(data);
        std::vector<RealT>     x = {1.0, 0.0};
        numberLocally(load);
        success *= load.initialize(state, x.data(), nullptr, nullptr) == 0;

        std::vector<RealT> g(load.sizeConstraints(), ZERO<RealT>);
        load.evaluateConstraints(x.data(), g.data());
        success *= isEqual(g[0], -0.6) && isEqual(g[1], -0.2);

        const std::vector<RealT> lambda(load.sizeConstraints(), ONE<RealT>);
        load.evaluateJacobian(x.data());
        load.evaluateHessian(x.data(), ONE<RealT>, lambda.data());
        success *= load.jacobian().values.empty() && load.hessian().values.empty();

        OPF::Load<RealT, IdxT> missing(data);
        numberLocally(missing);
        success *= missing.initialize(Model::StateData{}, x.data(), nullptr, nullptr) != 0;

        return success.report(__func__);
      }

    private:
      static constexpr const char* BRANCH_ID = "branch_1_2";

      static constexpr RealT R     = 0.02;
      static constexpr RealT X     = 0.15;
      static constexpr RealT G     = 0.01;
      static constexpr RealT B     = 0.3;
      static constexpr RealT GMAG  = 0.004;
      static constexpr RealT BMAG  = -0.02;
      static constexpr RealT TAP   = 1.05;
      static constexpr RealT PHASE = 0.1;
      static constexpr RealT SMAX  = 1.2;

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

      static constexpr RealT VMIN   = 0.9;
      static constexpr RealT VMAX   = 1.1;
      static constexpr RealT BUS_VR = 0.99;
      static constexpr RealT BUS_VI = 0.1;

      static inline const std::vector<RealT> BRANCH_X      = {1.02, 0.12, 0.96, -0.08};
      static inline const std::vector<RealT> BRANCH_LAMBDA = {0.3, -0.2, 1.5, -0.7, 2.1, 0.4};
      static inline const std::vector<RealT> FLAT_START    = {1.0, 0.0, 1.0, 0.0};
      static inline const std::vector<RealT> GENERATOR_X   = {0.3, 0.1, 1.01, 0.05};
      static inline const std::vector<RealT> SHUNT_X       = {0.98, 0.1};
      static inline const std::vector<RealT> SHUNT_LAMBDA  = {1.3, -0.6};
      static inline const std::vector<RealT> BUS_X         = {1.02, -0.05};
      static inline const std::vector<RealT> BUS_LAMBDA    = {0.7, -0.4, 1.3, 0.5};

      // Measured worst errors: power 1.5 eps, Jacobian 54.4 eps, Hessian 11.3 eps
      static constexpr RealT POWER_TOL    = 2 * std::numeric_limits<RealT>::epsilon();
      static constexpr RealT JACOBIAN_TOL = 64 * std::numeric_limits<RealT>::epsilon();
      static constexpr RealT HESSIAN_TOL  = 16 * std::numeric_limits<RealT>::epsilon();

      static OPF::BranchData<RealT, IdxT> opfBranchData()
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

      /**
       * @brief Currents into bus 1 and bus 2 from `PhasorDynamics::Branch`
       */
      static void pdCurrents(RealT  vr1,
                             RealT  vi1,
                             RealT  vr2,
                             RealT  vi2,
                             RealT& ir1,
                             RealT& ii1,
                             RealT& ir2,
                             RealT& ii2)
      {
        PD::BranchData<RealT, IdxT> data;
        data.parameters[PD::BranchParameters::R]     = R;
        data.parameters[PD::BranchParameters::X]     = X;
        data.parameters[PD::BranchParameters::G]     = G;
        data.parameters[PD::BranchParameters::B]     = B;
        data.parameters[PD::BranchParameters::Gmag]  = GMAG;
        data.parameters[PD::BranchParameters::Bmag]  = BMAG;
        data.parameters[PD::BranchParameters::tap]   = TAP;
        data.parameters[PD::BranchParameters::phase] = PHASE;

        PD::Bus<RealT, IdxT> bus1(vr1, vi1);
        PD::Bus<RealT, IdxT> bus2(vr2, vi2);
        bus1.allocate();
        bus1.initialize();
        bus1.evaluateResidual();
        bus2.allocate();
        bus2.initialize();
        bus2.evaluateResidual();

        PD::Branch<RealT, IdxT> branch(&bus1, &bus2, data);
        branch.allocate();
        branch.initialize();
        branch.evaluateResidual();

        ir1 = bus1.Ir();
        ii1 = bus1.Ii();
        ir2 = bus2.Ir();
        ii2 = bus2.Ii();
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
      static std::vector<typename ComponentT::ScalarT> variables(ComponentT&               component,
                                                                 const Model::StateData&   state,
                                                                 const std::vector<RealT>& values)
      {
        using ScalarT = typename ComponentT::ScalarT;

        numberLocally(component);
        std::vector<ScalarT> x(values.size());
        std::vector<RealT>   x_lower(values.size());
        std::vector<RealT>   x_upper(values.size());
        component.initialize(state, x.data(), x_lower.data(), x_upper.data());

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
