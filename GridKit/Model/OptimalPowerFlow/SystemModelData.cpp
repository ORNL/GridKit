/**
 * @file SystemModelData.cpp
 * @brief Optimal power flow limits and costs from a MATPOWER case.
 */

#include <algorithm>
#include <array>
#include <map>
#include <ratio>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <GridKit/Model/OptimalPowerFlow/SystemModelData.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    namespace
    {
      namespace Column = MatpowerColumns;

      using SystemModelDataT = SystemModelData<double, size_t>;
      using BusPair          = std::pair<size_t, size_t>;

      /// MATPOWER polynomial cost model
      constexpr double POLYNOMIAL = 2.0;

      /// Buses of a branch in increasing order
      BusPair busPair(size_t bus1, size_t bus2)
      {
        return {std::min(bus1, bus2), std::max(bus1, bus2)};
      }

      std::string name(const BusPair& buses)
      {
        return "Branch " + std::to_string(buses.first) + "-" + std::to_string(buses.second);
      }

      /// Voltage limits of the buses, matched by number
      void applyBuses(SystemModelDataT& data, const MatpowerMatrix& rows)
      {
        std::map<size_t, SystemModelDataT::BusDataT*> buses;
        for (auto& bus : data.bus)
        {
          buses[bus.number] = &bus;
        }

        for (const auto& row : rows)
        {
          const auto number = static_cast<size_t>(row.at(Column::BUS_I));
          const auto bus    = buses.find(number);
          if (bus == buses.end())
          {
            throw std::invalid_argument("MATPOWER bus " + std::to_string(number) + " is not in the case");
          }
          bus->second->parameters[BusParameters::Vmin] = row.at(Column::VMIN);
          bus->second->parameters[BusParameters::Vmax] = row.at(Column::VMAX);
        }

        if (rows.size() != data.bus.size())
        {
          throw std::invalid_argument("MATPOWER case has " + std::to_string(rows.size()) + " buses, but the case has "
                                      + std::to_string(data.bus.size()));
        }
      }

      /// Ratings of the branches, matched in order by their buses
      void applyBranches(SystemModelDataT& data, const MatpowerMatrix& rows, double mva_base)
      {
        std::map<BusPair, std::vector<SystemModelDataT::BranchDataT*>> branches;
        for (auto& branch : data.branch)
        {
          branches[busPair(branch.buses.at(BranchBuses::bus1), branch.buses.at(BranchBuses::bus2))].push_back(&branch);
        }

        std::map<BusPair, size_t> matched;
        for (const auto& row : rows)
        {
          if (row.at(Column::BR_STATUS) <= 0.0)
          {
            continue;
          }

          const BusPair buses = busPair(static_cast<size_t>(row.at(Column::F_BUS)), static_cast<size_t>(row.at(Column::T_BUS)));
          size_t&       count = matched[buses];
          if (count == branches[buses].size())
          {
            throw std::invalid_argument(name(buses) + " has more in-service MATPOWER branches than the case");
          }

          // A zero rating is unlimited
          auto& branch = *branches[buses][count++];
          if (row.at(Column::RATE_A) > 0.0)
          {
            branch.parameters[BranchParameters::Smax] = row.at(Column::RATE_A) / mva_base;
          }
        }

        for (const auto& [buses, list] : branches)
        {
          if (matched[buses] != list.size())
          {
            throw std::invalid_argument(name(buses) + " has fewer in-service MATPOWER branches than the case");
          }
        }
      }

      /**
       * @brief Costs from a MATPOWER polynomial in P [MW], whose coefficients
       * run from the highest order down
       */
      void setCost(SystemModelDataT::GeneratorDataT& generator, const std::vector<double>& cost, double mva_base)
      {
        if (cost.at(Column::MODEL) != POLYNOMIAL)
        {
          throw std::invalid_argument(generator.id + ": only polynomial MATPOWER costs are supported");
        }

        const auto            n = static_cast<size_t>(cost.at(Column::NCOST));
        std::array<double, 3> coefficients{};
        double                scale = 1.0;
        for (size_t k = 0; k < n; ++k)
        {
          const double coefficient = cost.at(Column::COST + n - 1 - k) * scale;
          if (k < coefficients.size())
          {
            coefficients[k] = coefficient;
          }
          else if (coefficient != 0.0)
          {
            throw std::invalid_argument(generator.id + ": MATPOWER cost is above quadratic");
          }
          scale *= mva_base;
        }

        generator.parameters[GeneratorParameters::c0] = coefficients[0];
        generator.parameters[GeneratorParameters::c1] = coefficients[1];
        generator.parameters[GeneratorParameters::c2] = coefficients[2];
      }

      /**
       * @brief Limits and costs of the generators, matched in order by bus
       *
       * MATPOWER generators at a bus without generators are skipped.
       */
      void applyGenerators(SystemModelDataT& data, const MatpowerMatrix& rows, const MatpowerMatrix& costs, double mva_base)
      {
        std::map<size_t, std::vector<SystemModelDataT::GeneratorDataT*>> generators;
        for (auto& generator : data.generator)
        {
          generators[generator.buses.at(GeneratorBuses::bus)].push_back(&generator);
        }

        std::map<size_t, size_t> matched;
        for (size_t i = 0; i < rows.size(); ++i)
        {
          const auto& row    = rows[i];
          const auto  number = static_cast<size_t>(row.at(Column::GEN_BUS));
          const auto  entry  = generators.find(number);
          if (row.at(Column::GEN_STATUS) <= 0.0 || entry == generators.end())
          {
            continue;
          }

          size_t& count = matched[number];
          if (count == entry->second.size())
          {
            throw std::invalid_argument("Bus " + std::to_string(number) + " has more in-service MATPOWER generators than the case");
          }

          auto& generator                                 = *entry->second[count++];
          generator.parameters[GeneratorParameters::Pmin] = row.at(Column::PMIN) / mva_base;
          generator.parameters[GeneratorParameters::Pmax] = row.at(Column::PMAX) / mva_base;
          generator.parameters[GeneratorParameters::Qmin] = row.at(Column::QMIN) / mva_base;
          generator.parameters[GeneratorParameters::Qmax] = row.at(Column::QMAX) / mva_base;
          setCost(generator, costs.at(i), mva_base);
        }

        for (const auto& [number, list] : generators)
        {
          if (matched[number] != list.size())
          {
            throw std::invalid_argument("Bus " + std::to_string(number) + " has fewer in-service MATPOWER generators than the case");
          }
        }
      }
    } // namespace

    void applyMatpowerData(SystemModelData<double, size_t>& data, const MatpowerData& matpower)
    {
      const double mva_base = data.va_base / std::mega::num;

      applyBuses(data, matpower.matrix("bus"));
      applyBranches(data, matpower.matrix("branch"), mva_base);
      applyGenerators(data, matpower.matrix("gen"), matpower.matrix("gencost"), mva_base);
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
