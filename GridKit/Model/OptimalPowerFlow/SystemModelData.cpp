/**
 * @file SystemModelData.cpp
 * @brief Optimal power flow data from a MATPOWER case.
 */

#include <array>
#include <format>
#include <fstream>
#include <numbers>
#include <ratio>
#include <stdexcept>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/OptimalPowerFlow/SystemModelData.hpp>
#include <GridKit/Model/PowerFlow/MatpowerParser.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace
    {
      using SystemModelDataT = SystemModelData<double, size_t>;
      using MatpowerDataT    = PowerFlowData::SystemModelData<double, size_t>;

      /// MATPOWER bus types
      constexpr size_t REFERENCE_BUS = 3;
      constexpr size_t ISOLATED_BUS  = 4;

      /// MATPOWER polynomial cost model
      constexpr size_t POLYNOMIAL = 2;

      /// MATPOWER angle-difference limit of no effect [deg]
      constexpr double NO_ANGLE_LIMIT = 360.0;

      constexpr double RADIANS_PER_DEGREE = std::numbers::pi / 180.0;

      bool isZero(double value)
      {
        return !(value < ZERO<double>) && !(value > ZERO<double>);
      }

      /// Buses, and the loads and shunts of their demand and shunt columns
      void addBuses(SystemModelDataT& data, const MatpowerDataT& matpower)
      {
        const double mva_base = matpower.baseMVA;
        for (size_t i = 0; i < matpower.bus.size(); ++i)
        {
          const auto& row = matpower.bus[i];
          if (row.type == ISOLATED_BUS)
          {
            throw std::invalid_argument(std::format("MATPOWER bus {}: isolated buses are not supported", row.bus_i));
          }

          auto& bus                           = data.bus.emplace_back();
          bus.number                          = row.bus_i;
          bus.reference                       = row.type == REFERENCE_BUS;
          bus.parameters[BusParameters::Vmin] = row.Vmin;
          bus.parameters[BusParameters::Vmax] = row.Vmax;

          const auto& demand = matpower.load[i];
          if (!isZero(demand.Pd) || !isZero(demand.Qd))
          {
            auto& load                         = data.load.emplace_back();
            load.id                            = std::format("load_{}", row.bus_i);
            load.buses[LoadBuses::bus]         = row.bus_i;
            load.parameters[LoadParameters::P] = demand.Pd / mva_base;
            load.parameters[LoadParameters::Q] = demand.Qd / mva_base;
          }

          if (!isZero(row.Gs) || !isZero(row.Bs))
          {
            auto& shunt                          = data.shunt.emplace_back();
            shunt.id                             = std::format("shunt_{}", row.bus_i);
            shunt.buses[ShuntBuses::bus]         = row.bus_i;
            shunt.parameters[ShuntParameters::G] = row.Gs / mva_base;
            shunt.parameters[ShuntParameters::B] = row.Bs / mva_base;
          }
        }
      }

      /// Branches, numbered by MATPOWER row
      void addBranches(SystemModelDataT& data, const MatpowerDataT& matpower)
      {
        const double mva_base = matpower.baseMVA;
        for (size_t k = 0; k < matpower.branch.size(); ++k)
        {
          const auto& row = matpower.branch[k];

          const bool no_angle_limit = (row.angmin <= -NO_ANGLE_LIMIT && row.angmax >= NO_ANGLE_LIMIT)
                                      || (isZero(row.angmin) && isZero(row.angmax));
          if (!no_angle_limit)
          {
            throw std::invalid_argument(std::format("MATPOWER branch {}: angle difference limits are not supported", k + 1));
          }

          // A zero ratio is a line
          double tap = ONE<double>;
          if (!isZero(row.ratio))
          {
            tap = row.ratio;
          }

          auto& branch                               = data.branch.emplace_back();
          branch.id                                  = std::format("branch_{}", k + 1);
          branch.in_service                          = row.status > 0;
          branch.buses[BranchBuses::bus1]            = row.fbus;
          branch.buses[BranchBuses::bus2]            = row.tbus;
          branch.parameters[BranchParameters::R]     = row.r;
          branch.parameters[BranchParameters::X]     = row.x;
          branch.parameters[BranchParameters::B]     = row.b;
          branch.parameters[BranchParameters::tap]   = tap;
          branch.parameters[BranchParameters::phase] = row.angle * RADIANS_PER_DEGREE;

          // MATPOWER scales the line charging at the tapped bus by 1 / tap^2
          branch.parameters[BranchParameters::Bmag] = HALF<double> * row.b * (ONE<double> / (tap * tap) - ONE<double>);

          // A zero rating is unlimited
          if (row.rateA > ZERO<double>)
          {
            branch.parameters[BranchParameters::Smax] = row.rateA / mva_base;
          }
        }
      }

      /**
       * @brief Costs from a MATPOWER polynomial in P [MW], whose coefficients
       * run from the highest order down
       */
      void setCost(SystemModelDataT::GeneratorDataT& generator, const MatpowerDataT::GenCostDataT& cost, double mva_base)
      {
        if (cost.kind != POLYNOMIAL || cost.rest.size() != cost.n)
        {
          throw std::invalid_argument(generator.id + ": only polynomial MATPOWER costs are supported");
        }

        std::array<double, 3> coefficients{};
        double                scale = ONE<double>;
        for (size_t k = 0; k < cost.n; ++k)
        {
          const double coefficient = cost.rest[cost.n - 1 - k] * scale;
          if (k < coefficients.size())
          {
            coefficients[k] = coefficient;
          }
          else if (!isZero(coefficient))
          {
            throw std::invalid_argument(generator.id + ": MATPOWER cost is above quadratic");
          }
          scale *= mva_base;
        }

        generator.parameters[GeneratorParameters::c0] = coefficients[0];
        generator.parameters[GeneratorParameters::c1] = coefficients[1];
        generator.parameters[GeneratorParameters::c2] = coefficients[2];
      }

      /// Generators, numbered by MATPOWER row
      void addGenerators(SystemModelDataT& data, const MatpowerDataT& matpower)
      {
        if (matpower.gencost.size() != matpower.gen.size())
        {
          throw std::invalid_argument("MATPOWER case needs one active power cost per generator");
        }

        const double mva_base = matpower.baseMVA;
        for (size_t i = 0; i < matpower.gen.size(); ++i)
        {
          const auto& row = matpower.gen[i];

          auto& generator                                 = data.generator.emplace_back();
          generator.id                                    = std::format("gen_{}", i + 1);
          generator.in_service                            = row.status > 0;
          generator.buses[GeneratorBuses::bus]            = row.bus;
          generator.parameters[GeneratorParameters::Pmin] = row.Pmin / mva_base;
          generator.parameters[GeneratorParameters::Pmax] = row.Pmax / mva_base;
          generator.parameters[GeneratorParameters::Qmin] = row.Qmin / mva_base;
          generator.parameters[GeneratorParameters::Qmax] = row.Qmax / mva_base;
          setCost(generator, matpower.gencost[i], mva_base);
        }
      }
    } // namespace

    SystemModelData<double, size_t> parseMatpowerCase(std::istream& stream)
    {
      MatpowerDataT matpower{};
      readMatPower(matpower, stream);
      if (matpower.baseMVA <= ZERO<double>)
      {
        throw std::invalid_argument("MATPOWER case needs a positive baseMVA");
      }

      SystemModelDataT data;
      data.va_base = matpower.baseMVA * std::mega::num;

      addBuses(data, matpower);
      addBranches(data, matpower);
      addGenerators(data, matpower);

      return data;
    }

    SystemModelData<double, size_t> parseMatpowerCase(const std::filesystem::path& file)
    {
      std::ifstream stream(file);
      if (!stream)
      {
        throw std::runtime_error("Could not open MATPOWER case " + file.string());
      }
      return parseMatpowerCase(stream);
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
