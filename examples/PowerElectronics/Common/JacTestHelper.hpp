#pragma once

#include <cmath>
#include <iostream>
#include <optional>
#include <unordered_map>
#include <vector>

#include <GridKit/LinearAlgebra/SparseMatrix/CsrMatrix.hpp>
#include <GridKit/MemoryUtilities/MemoryUtils.hpp>
#include <GridKit/Model/PowerElectronics/SubsystemModel.hpp>
#include <GridKit/Testing/TestHelpers.hpp>

#include "GridKit/Testing/Testing.hpp"

namespace GridKit
{
  namespace Testing
  {
    /**
     * @brief Verify that the subsystem Jacobian is an exact subset of the
     * full-system Jacobian to ensure accuracy of subsystem Jacobians.
     *
     * Each subsystem Jacobian entry should match the corresponding entry in the
     * full-system Jacobian after converting local subsystem indices back to global
     * indices. The check verifies that:
     *
     * - every subsystem entry exists in the full-system Jacobian,
     * - every full-system entry whose row and column belong to the subsystem
     * exists in the subsystem Jacobian, and
     * - matching entries agree within the specified tolerance.
     *
     * @note When the components in a subsystem are evaluated in a different order
     *       than in the monolithic reference system, a larger tolerance may be
     *       required to account for floating-point roundoff introduced by the
     *       different summation order.
     *
     * @param full_jac Full-system Jacobian.
     * @param sub_jac Subsystem Jacobian.
     * @param subsystem Subsystem associated with sub_jac.
     * @param tolerance Maximum allowed difference between matching entries.
     *
     * @return true if the subsystem Jacobian is an exact subset of the
     * full-system Jacobian, false otherwise.
     */

    template <typename RealT, typename IdxT>
    bool verifySubsystemJacobian(
        GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>& full_jac,
        GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>& sub_jac,
        GridKit::SubsystemModel<RealT, IdxT>&           subsystem,
        std::optional<RealT>                            tolerance = std::nullopt)
    {
      constexpr auto host      = memory::HOST;
      // Full-system CSR data.
      const IdxT*    full_rows = full_jac.getRowData(host);
      const IdxT*    full_cols = full_jac.getColData(host);
      const RealT*   full_vals = full_jac.getValues(host);

      // Subsystem CSR data.
      const IdxT*  sub_rows = sub_jac.getRowData(host);
      const IdxT*  sub_cols = sub_jac.getColData(host);
      const RealT* sub_vals = sub_jac.getValues(host);

      const IdxT num_sub_rows = sub_jac.getNumRows();
      const IdxT num_sub_cols = sub_jac.getNumColumns();

      // The subsystem internal map stores global-to-local indices.
      const auto& global_to_local = subsystem.getInternalMap();

      // The internal map should contain one entry for every subsystem row.
      if (global_to_local.size() != num_sub_rows)
      {
        std::cout << "Internal map size mismatch: map has "
                  << global_to_local.size()
                  << " entries, but subsystem Jacobian has "
                  << num_sub_rows
                  << " rows\n";

        return false;
      }

      // The subsystem Jacobian is expected to be square.
      if (num_sub_rows != num_sub_cols)
      {
        std::cout << "Subsystem Jacobian is not square: "
                  << num_sub_rows
                  << " rows and "
                  << num_sub_cols
                  << " columns\n";

        return false;
      }

      // Build the reverse map from local subsystem indices to global indices.
      // This is needed to compare subsystem rows and columns with the
      // corresponding entries in the full-system Jacobian.
      std::vector<IdxT> local_to_global(num_sub_rows);
      std::vector<bool> local_index_found(num_sub_rows, false);

      for (const auto& [global_index, local_index] : global_to_local)
      {

        // Make sure the global index is valid for the full-system Jacobian.
        if (global_index >= full_jac.getNumRows())
        {
          std::cout << "Invalid global index "
                    << global_index
                    << " in subsystem internal map\n";

          return false;
        }

        // Make sure the local index is valid for the subsystem Jacobian.
        if (local_index >= num_sub_rows)
        {
          std::cout << "Invalid local index "
                    << local_index
                    << " mapped from global index "
                    << global_index << '\n';

          return false;
        }

        // Each local index should correspond to only one global index.
        if (local_index_found[local_index])
        {
          std::cout << "Duplicate local index "
                    << local_index
                    << " in subsystem internal map\n";

          return false;
        }

        local_to_global[local_index]   = global_index;
        local_index_found[local_index] = true;
      }

      // Make sure every subsystem local index has a global index.
      for (IdxT local_index = 0; local_index < num_sub_rows; ++local_index)
      {
        if (!local_index_found[local_index])
        {
          std::cout << "No global index maps to local index "
                    << local_index << '\n';

          return false;
        }
      }

      bool matches = true;

      // Compare each subsystem row with the corresponding full-system row.
      for (IdxT local_row = 0; local_row < num_sub_rows; ++local_row)
      {
        const IdxT global_row = local_to_global[local_row];

        const IdxT sub_begin  = sub_rows[local_row];
        const IdxT sub_end    = sub_rows[local_row + 1];
        const IdxT full_begin = full_rows[global_row];
        const IdxT full_end   = full_rows[global_row + 1];

        /*
         * Store the current subsystem row using global column indices.
         *
         * The subsystem Jacobian uses local column indices, while the full-system
         * Jacobian uses global column indices. Convert each local column to its
         * corresponding global column so the two rows can be compared directly.
         */
        std::unordered_map<IdxT, RealT> sub_row_entries;

        for (IdxT sub_index = sub_begin; sub_index < sub_end; ++sub_index)
        {
          const IdxT local_column = sub_cols[sub_index];

          // Fail if the subsystem column index is outside the valid local range.
          if (local_column >= num_sub_cols)
          {
            std::cout << "Invalid subsystem column index "
                      << local_column
                      << " in local row "
                      << local_row << '\n';

            matches = false;
            continue;
          }

          const IdxT global_column = local_to_global[local_column];

          // Fail if the subsystem row contains more than one entry for the same column.
          if (sub_row_entries.find(global_column) != sub_row_entries.end())
          {
            std::cout << "Duplicate subsystem entry at ("
                      << global_row << ", "
                      << global_column << ")\n";

            matches = false;
            continue;
          }

          sub_row_entries[global_column] = sub_vals[sub_index];
        }

        // Compare full-system entries whose columns belong to the subsystem.
        for (IdxT full_index = full_begin; full_index < full_end; ++full_index)
        {
          const IdxT global_column = full_cols[full_index];

          // No need to check if column belongs outside the subsystem.
          if (global_to_local.find(global_column) == global_to_local.end())
          {
            continue;
          }

          auto sub_entry = sub_row_entries.find(global_column);

          // Then it must exist in the subsystem Jacobian; fail otherwise.
          if (sub_entry == sub_row_entries.end())
          {
            std::cout << "Entry exists only in full Jacobian at ("
                      << global_row << ", "
                      << global_column << ")\n";

            matches = false;
            continue;
          }

          const RealT full_value = full_vals[full_index];
          const RealT sub_value  = sub_entry->second;
          const RealT difference = std::abs(full_value - sub_value);

          // Different component evaluation orders can change the order of floating-point
          // summation and introduce small roundoff differences. Use an appropriate
          // tolerance when comparing against the monolithic reference.
          // if tolerance is not supplied we simply use default machine precision provided
          // by GridKit's Test::isEqual
          auto isEqual = [&tolerance](RealT value, RealT reference)
          {
            if (tolerance)
            {
              return GridKit::Testing::isEqual(value, reference, *tolerance);
            }

            return GridKit::Testing::isEqual(value, reference);
          };

          // Then the values must agree, fail otherwise
          if (!isEqual(sub_value, full_value))
          {
            std::cout << "Jacobian value mismatch at ("
                      << global_row << ", "
                      << global_column << "): "
                      << "full = " << full_value
                      << ", subsystem = " << sub_value
                      << ", difference = " << difference << '\n';

            matches = false;
          }

          // Remove the matched entry
          sub_row_entries.erase(sub_entry);
        }

        // Any entries remaining must be missing from the
        // full-system Jacobian, so this subsystem row contains incorrect entries, fail!
        for (const auto& entry : sub_row_entries)
        {
          const IdxT global_column = entry.first;

          std::cout << "Entry exists only in subsystem Jacobian at ("
                    << global_row << ", "
                    << global_column << ")\n";

          matches = false;
        }
      }

      return matches;
    }

  } // namespace Testing
} // namespace GridKit
