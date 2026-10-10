/**
 * @file MatpowerData.hpp
 * @brief Numeric data of a MATPOWER case file.
 */

#pragma once

#include <cstddef>
#include <filesystem>
#include <istream>
#include <map>
#include <string>
#include <vector>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /// Rows of a MATPOWER matrix
    using MatpowerMatrix = std::vector<std::vector<double>>;

    /**
     * @brief Numeric matrices of a MATPOWER case file by field name
     *
     * Holds every `mpc.<name> = [ ... ];` matrix. Rows end at `;` or at the
     * end of a line and may have any number of columns. Scalars and cell
     * arrays such as `mpc.bus_name` are skipped.
     */
    struct MatpowerData
    {
      std::map<std::string, MatpowerMatrix> matrices;

      /// Matrix `mpc.<name>`
      const MatpowerMatrix& matrix(const std::string& name) const;
    };

    /// Zero-based columns of MATPOWER `idx_bus`, `idx_brch`, `idx_gen`, and `idx_cost`
    namespace MatpowerColumns
    {
      inline constexpr size_t BUS_I      = 0;
      inline constexpr size_t VMAX       = 11;
      inline constexpr size_t VMIN       = 12;
      inline constexpr size_t F_BUS      = 0;
      inline constexpr size_t T_BUS      = 1;
      inline constexpr size_t RATE_A     = 5;
      inline constexpr size_t BR_STATUS  = 10;
      inline constexpr size_t GEN_BUS    = 0;
      inline constexpr size_t QMAX       = 3;
      inline constexpr size_t QMIN       = 4;
      inline constexpr size_t GEN_STATUS = 7;
      inline constexpr size_t PMAX       = 8;
      inline constexpr size_t PMIN       = 9;
      inline constexpr size_t MODEL      = 0;
      inline constexpr size_t NCOST      = 3;
      inline constexpr size_t COST       = 4;
    } // namespace MatpowerColumns

    MatpowerData parseMatpowerData(std::istream& stream);
    MatpowerData parseMatpowerData(const std::filesystem::path& file);
  } // namespace OptimalPowerFlow
} // namespace GridKit
