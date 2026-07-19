#pragma once

#include <cmath>
#include <iostream>
#include <unordered_map>
#include <vector>

#include <GridKit/LinearAlgebra/SparseMatrix/CsrMatrix.hpp>
#include <GridKit/MemoryUtilities/MemoryUtils.hpp>
#include <GridKit/Model/PowerElectronics/SubsystemModel.hpp>

namespace GridKit
{
  namespace Testing
  {

    template <typename RealT, typename IdxT>
    bool verifySubsystemJacobian(
        GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>& full_jac,
        GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>& sub_jac,
        GridKit::SubsystemModel<RealT, IdxT>&           subsystem,
        RealT                                           tolerance = static_cast<RealT>(1e-12))
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

      /*
       * The subsystem internal map stores:
       *
       *     global index -> local index
       */
      const auto& global_to_local = subsystem.getInternalMap();

      if (global_to_local.size() != num_sub_rows)
      {
        std::cout << "Internal map size mismatch: map has "
                  << global_to_local.size()
                  << " entries, but subsystem Jacobian has "
                  << num_sub_rows
                  << " rows\n";

        return false;
      }

      if (num_sub_rows != num_sub_cols)
      {
        std::cout << "Subsystem Jacobian is not square: "
                  << num_sub_rows
                  << " rows and "
                  << num_sub_cols
                  << " columns\n";

        return false;
      }

      /*
       * Build the reverse mapping:
       *
       *     local index -> global index
       *
       * This is needed because the subsystem Jacobian stores local row
       * and column indices.
       */
      std::vector<IdxT> local_to_global(num_sub_rows);
      std::vector<bool> local_index_found(num_sub_rows, false);

      for (const auto& [global_index, local_index] : global_to_local)
      {

        if (global_index >= full_jac.getNumRows())
        {
          std::cout << "Invalid global index "
                    << global_index
                    << " in subsystem internal map\n";

          return false;
        }

        if (local_index >= num_sub_rows)
        {
          std::cout << "Invalid local index "
                    << local_index
                    << " mapped from global index "
                    << global_index << '\n';

          return false;
        }

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

      // Verify that every local index has a corresponding global index.
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

      /*
       * Compare each subsystem row against the corresponding full-system row.
       */
      for (IdxT local_row = 0; local_row < num_sub_rows; ++local_row)
      {
        const IdxT global_row = local_to_global[local_row];

        const IdxT sub_begin  = sub_rows[local_row];
        const IdxT sub_end    = sub_rows[local_row + 1];
        const IdxT full_begin = full_rows[global_row];
        const IdxT full_end   = full_rows[global_row + 1];

        /*
         * Store this subsystem row using global column indices.
         *
         * The subsystem CSR row is sorted by local column index. After
         * converting local columns to global columns, the global columns
         * may no longer be sorted. Therefore, we use a lookup map instead
         * of comparing both CSR rows sequentially.
         */
        std::unordered_map<IdxT, RealT> sub_row_entries;

        for (IdxT sub_index = sub_begin; sub_index < sub_end; ++sub_index)
        {
          const IdxT local_column = sub_cols[sub_index];

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

          // Check for duplicate columns in this subsystem row.
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

        /*
         * Compare every full-system entry whose column belongs to this
         * subsystem.
         */
        for (IdxT full_index = full_begin;
             full_index < full_end;
             ++full_index)
        {
          const IdxT global_column = full_cols[full_index];

          // Ignore columns owned by another subsystem.
          if (global_to_local.find(global_column) == global_to_local.end())
          {
            continue;
          }

          auto sub_entry = sub_row_entries.find(global_column);

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

          if (difference > tolerance)
          {
            std::cout << "Jacobian value mismatch at ("
                      << global_row << ", "
                      << global_column << "): "
                      << "full = " << full_value
                      << ", subsystem = " << sub_value
                      << ", difference = " << difference << '\n';

            matches = false;
          }

          /*
           * Remove the matched subsystem entry.
           *
           * After the full row has been checked, any remaining entries
           * exist only in the subsystem Jacobian.
           */
          sub_row_entries.erase(sub_entry);
        }

        // Report subsystem entries that were not found in the full Jacobian.
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
