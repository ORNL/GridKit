#pragma once

#include <algorithm>
#include <cstddef>
#include <memory>
#include <numeric>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

#include <GridKit/LinearAlgebra/SparseMatrix/CsrMatrix.hpp>
#include <GridKit/Optimization/DerivativeStructure.hpp>

namespace GridKit
{
  namespace Optimization
  {
    /**
     * @brief Assemble fixed local sparse contributions into a deduplicated CSR matrix.
     *
     * Entry groups describe contributor-owned derivative buffers. allocate()
     * sorts and deduplicates their coordinates once and binds each group to a
     * stable set of CSR value slots. Runtime evaluations clear the values and
     * scatter-add through those bindings without changing the structure.
     */
    template <typename real_type, typename index_type>
    class SparseMatrixAssembler
    {
    public:
      using RealT      = real_type;
      using IdxT       = index_type;
      using CsrMatrixT = LinearAlgebra::CsrMatrix<RealT, IdxT>;
      using EntrySpan  = std::span<const SparseEntry<IdxT>>;

      SparseMatrixAssembler() = default;

      // Contribution bindings use the assembler address as part of their identity.
      SparseMatrixAssembler(const SparseMatrixAssembler&)            = delete;
      SparseMatrixAssembler& operator=(const SparseMatrixAssembler&) = delete;
      SparseMatrixAssembler(SparseMatrixAssembler&&)                 = delete;
      SparseMatrixAssembler& operator=(SparseMatrixAssembler&&)      = delete;

      /**
       * @brief Opaque binding between one contributor and its CSR value slots.
       */
      class Contribution
      {
      public:
        std::size_t size() const
        {
          return slots_.size();
        }

      private:
        friend class SparseMatrixAssembler;

        const SparseMatrixAssembler* assembler_{nullptr};
        std::size_t                  generation_{0};
        std::vector<IdxT>            slots_;
      };

      int allocate(IdxT                       row_count,
                   IdxT                       column_count,
                   std::span<const EntrySpan> entry_groups,
                   std::span<Contribution>    bindings)
      {
        if (!validDimension(row_count)
            || !validDimension(column_count)
            || entry_groups.size() != bindings.size())
        {
          return 1;
        }

        std::vector<SparseEntry<IdxT>> entries;
        std::vector<std::size_t>       group_offsets;
        group_offsets.reserve(entry_groups.size() + 1);
        group_offsets.push_back(0);

        for (const auto group : entry_groups)
        {
          for (const auto& entry : group)
          {
            if (!validIndex(entry.row, row_count)
                || !validIndex(entry.column, column_count))
            {
              return 1;
            }
            entries.push_back(entry);
          }
          group_offsets.push_back(entries.size());
        }

        std::vector<std::size_t> order(entries.size());
        std::iota(order.begin(), order.end(), std::size_t{0});
        std::stable_sort(order.begin(), order.end(), [&entries](std::size_t lhs, std::size_t rhs)
                         {
                           const auto& left  = entries[lhs];
                           const auto& right = entries[rhs];
                           if (left.row != right.row)
                           {
                             return left.row < right.row;
                           }
                           return left.column < right.column; });

        std::vector<IdxT> contribution_to_csr(entries.size(), IdxT{0});
        std::vector<IdxT> row_data(static_cast<std::size_t>(row_count) + 1, IdxT{0});
        std::vector<IdxT> column_data;
        column_data.reserve(entries.size());

        bool              have_previous = false;
        SparseEntry<IdxT> previous{};
        IdxT              slot{0};

        for (const auto original : order)
        {
          const auto entry = entries[original];
          if (!have_previous
              || entry.row != previous.row
              || entry.column != previous.column)
          {
            previous      = entry;
            have_previous = true;
            slot          = static_cast<IdxT>(column_data.size());
            column_data.push_back(entry.column);
            ++row_data[static_cast<std::size_t>(entry.row) + 1];
          }
          contribution_to_csr[original] = slot;
        }

        for (std::size_t row = 1; row < row_data.size(); ++row)
        {
          row_data[row] += row_data[row - 1];
        }

        auto rows = std::make_unique<IdxT[]>(row_data.size());
        auto cols = std::make_unique<IdxT[]>(column_data.size());
        auto vals = std::make_unique<RealT[]>(column_data.size());

        std::copy(row_data.begin(), row_data.end(), rows.get());
        std::copy(column_data.begin(), column_data.end(), cols.get());

        auto  nonzero_count = static_cast<IdxT>(column_data.size());
        auto* row_pointer   = rows.get();
        auto* col_pointer   = cols.get();
        auto* val_pointer   = vals.get();
        auto  matrix        = std::make_unique<CsrMatrixT>(row_count,
                                                   column_count,
                                                   nonzero_count,
                                                   &row_pointer,
                                                   &col_pointer,
                                                   &val_pointer);
        rows.release();
        cols.release();
        vals.release();

        const std::size_t         new_generation = generation_ + 1;
        std::vector<Contribution> new_bindings(entry_groups.size());
        for (std::size_t group = 0; group < entry_groups.size(); ++group)
        {
          auto& binding       = new_bindings[group];
          binding.assembler_  = this;
          binding.generation_ = new_generation;
          binding.slots_.assign(
              contribution_to_csr.begin()
                  + static_cast<std::ptrdiff_t>(group_offsets[group]),
              contribution_to_csr.begin()
                  + static_cast<std::ptrdiff_t>(group_offsets[group + 1]));
        }

        matrix_     = std::move(matrix);
        generation_ = new_generation;
        std::move(new_bindings.begin(), new_bindings.end(), bindings.begin());
        return 0;
      }

      int clearValues()
      {
        if (matrix_ == nullptr)
        {
          return 1;
        }

        auto* values = matrix_->getValues(memory::HOST);
        std::fill_n(values, static_cast<std::size_t>(matrix_->getNnz()), RealT{0});
        return matrix_->setUpdated(memory::HOST);
      }

      int addValues(const Contribution& binding, std::span<const RealT> values)
      {
        if (matrix_ == nullptr
            || binding.assembler_ != this
            || binding.generation_ != generation_
            || binding.size() != values.size())
        {
          return 1;
        }

        auto* matrix_values = matrix_->getValues(memory::HOST);
        for (std::size_t entry = 0; entry < values.size(); ++entry)
        {
          matrix_values[static_cast<std::size_t>(binding.slots_[entry])] += values[entry];
        }
        return matrix_->setUpdated(memory::HOST);
      }

      CsrMatrixT* matrix() const
      {
        return matrix_.get();
      }

    private:
      static constexpr bool validDimension(IdxT dimension)
      {
        if constexpr (std::is_signed_v<IdxT>)
        {
          return dimension >= 0;
        }
        return true;
      }

      static constexpr bool validIndex(IdxT index, IdxT extent)
      {
        if constexpr (std::is_signed_v<IdxT>)
        {
          if (index < 0)
          {
            return false;
          }
        }
        return index < extent;
      }

      std::unique_ptr<CsrMatrixT> matrix_;
      std::size_t                 generation_{0};
    };
  } // namespace Optimization
} // namespace GridKit
