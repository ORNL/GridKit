#include <LinearAlgebra/SparseMatrix/CsrMatrix.hpp>

namespace GridKit
{
  namespace LinearAlgebra
  {
    /**
     * @brief Represents the sparse matrix operation \f$\sum_i P_i A_i P_i^\top\f$ AKA a permuted sum.
     *        Each summand \f$A_i\f$ has its rows and columns permuted before being summed. This class
     *        has two steps: 1) an expensive `analyzeSparsity()` step which only needs to be computed once
     *        as long as the summands, their sparsity patterns, and their permutations are unchanged;
     *        and 2) an inexpensive `apply()` step which computes the resultant sum matrix using a
     *        given `CsrBuilder`. Due to the use of the `CsrBuilder` class, if a previously constructed
     *        sum is used, then its sparsity pattern can be re-used and the subsequent `apply()` steps
     *        will be even faster.
     *
     *        Useful for e.g. computing the system jacobian from component jacobians, where the permutations
     *        \f$P_i\f$ are simply mappings from component variables to system variables.
     */
    template <class ScalarT, typename IdxT>
    class CsrPermutedSum
    {
    private:
      using CsrMatrix = LinearAlgebra::CsrMatrix<ScalarT, IdxT>;

      /**
       * @brief Represents the index of an element in a summand
       */
      struct LocalSummandIdx
      {
        /// @brief The index of the summand matrix that this refers to
        size_t summand_idx_;
        /// @brief The index of the row/column inside that summand that this refers to
        IdxT   local_idx_;
      };

      /**
       * @brief A plan for constructing a row of the final sum which is calculated
       *        exclusively by one of the summands. This happens, for instance,
       *        during system jacobian construction when a system variable is an
       *        internal variable of a component.
       */
      struct ExclusiveRowPlan
      {
        /**
         * @brief A single element in a summand's row.
         */
        struct Element
        {
          /**
           * @brief The index of the element in the summand.
           *        Used to copy the value of the element into the final sum matrix.
           */
          size_t elem_idx_;
          /**
           * @brief The column of the element in the final sum. Used to construct the
           *        sparsity pattern of this row in the final sum matrix.
           */
          IdxT   sum_col_;
        };

        /**
         * @brief The index of the summand that this row "belongs" to.
         */
        LocalSummandIdx summand_idx_;

        /**
         * @brief The expected number of nonzero elements in the row of the summand.
         *        Currently only used to verify the sparsity of the summand has not changed since being
         *        analyzed, but could enable a quick memcpy in the future (currently limited by the fact
         *        that some columns in component jacobians do not map to columns in the system jacobian).
         */
        IdxT                 row_nnz_;
        /**
         * @brief The elements which belong to this row. Iterating this will tell you when you need to call
         *        `CsrBuilder::elem()` to construct this row.
         */
        std::vector<Element> elements_;

        /**
         * @brief Apply this plan by using a builder to construct this row and all of its elements.
         *        Assumes `CsrBuilder::row()` has been called already and only calls `CsrBuilder::elem()`.
         * @tparam CsrBuilder Expected to be some sort of `CsrBuilder`.
         * @param summands A list of the matrices being summed.
         * @param builder The builder which is constructing the final sum matrix.
         */
        template <class CsrBuilder>
        void apply(const std::vector<std::reference_wrapper<CsrMatrix>>& summands, CsrBuilder& builder) const
        {
          const CsrMatrix& summand = summands[summand_idx_.summand_idx_];

          IdxT row_start = summand.rowIndices()[summand_idx_.local_idx_];
          IdxT row_end   = summand.rowIndices()[summand_idx_.local_idx_ + 1];

          // Verify the sparsity pattern is as we expect
          assert(summand.rowIndices().size() >= summand_idx_.local_idx_ + 2);
          assert(row_end - row_start == row_nnz_);

          // Construct the elements - simply copy the values from the summand, since this row is exclusively
          // constructed by this summand.
          for (auto element : elements_)
          {
            builder.elem(element.sum_col_, summand.values()[element.elem_idx_]);
          }
        }
      };

      /**
       * @brief A plan for constructing a row of the final sum which is calculated
       *        by multiple summands. This happens, for instance, during system jacobian
       *        construction when a system variable is an external variable, such as a bus voltage.
       */
      struct SharedRowPlan
      {
        /**
         * @brief A single element in the row. An element is the sum of the elements of (potentially) several summands
         */
        struct SumElement
        {
          /**
           * @brief A list of all elements in summands which sum to this element in the final sum.
           */
          std::vector<LocalSummandIdx> summand_elements_;
          /**
           * @brief The column of this element in the final sum.
           */
          IdxT                         column_;
        };

        /**
         * @brief A list of all elements on this row of the final sum.
         */
        std::vector<SumElement> elements_;

        /**
         * @brief Apply this plan by using a builder to construct this row and all of its elements.
         *        Assumes `CsrBuilder::row()` has been called already and only calls `CsrBuilder::elem()`.
         * @tparam CsrBuilder Expected to be some sort of `CsrBuilder`.
         * @param summands A list of the matrices being summed.
         * @param builder The builder which is constructing the final sum matrix.
         */
        template <class CsrBuilder>
        void apply(const std::vector<std::reference_wrapper<CsrMatrix>>& summands, CsrBuilder& builder) const
        {
          for (auto& element : elements_)
          {
            ScalarT sum = 0;

            // Sum over the elements of all summands which contribute to this element in the final sum
            for (auto& element : element.summand_elements_)
            {
              const CsrMatrix& summand = summands[element.summand_idx_];

              sum += summand.values()[element.local_idx_];
            }

            builder.elem(element.column_, sum);
          }
        }
      };

      /**
       * @brief The plan for constructing a single row in the final sum matrix.
       */
      using RowPlan = std::variant<ExclusiveRowPlan, SharedRowPlan>;

      /**
       * @brief A collection of plans for constructing each row in the final sum matrix.
       */
      std::vector<RowPlan> row_plans_;

      /**
       * @brief Number of summands which are part of the sum.
       */
      size_t num_summands_;

      /**
       * @brief The dimension of the final sum matrix. Due to the nature of the operation being performed
       *        (rows and columns of summands are permuted the same way) - the size of the final sum must be
       *        square.
       */
      size_t summand_size_;

      /**
       * @brief Pass-through constructor. Row plans are not expected to be created by the user, so make sure to call
       *        `analyzeSparsity()` to properly construct a `CsrPermutedSum`.
       */
      CsrPermutedSum(std::vector<RowPlan> row_plans, size_t num_summands, size_t summand_size) : row_plans_(row_plans), num_summands_(num_summands), summand_size_(summand_size)
      {
      }

      /**
       * @brief Create a plan for constructing a shared row of the final sum matrix, based on the given sparsities and permutations.
       * @param summands A list of sparsity patterns for summands.
       * @param summand_contributions A list of all summands and their rows that map to this row.
       * @param permutations A list of the permutations \f$P_i\f$.
       * @param size The size of the resulting sum matrix.
       */
      static SharedRowPlan createSharedRowPlan(
          const std::vector<std::reference_wrapper<CsrMatrix>>& summands,
          const std::vector<LocalSummandIdx>&                   summand_contributions,
          const std::vector<std::vector<IdxT>>&                 permutations,
          size_t                                                size)
      {
        // Keep track of the number of nonzero elements in this row
        size_t nnz = 0;

        // Sort the columns of all of the permuted rows of the summands.
        // Sort used is a radix sort, where the radix is equal to the number of columns of the final sum.
        // Due to this choice, only a single step of radix sort needs to be performed,
        // but we allocate a vector of size equal to the number of columns in the matrix.
        std::vector<std::vector<LocalSummandIdx>> summand_elements(size);

        // Loop over all summands which have rows that map to this row in the final sum
        for (const auto& contribution : summand_contributions)
        {
          auto [summand_idx, local_idx] = contribution;
          const CsrMatrix& comp_jac     = summands[summand_idx];

          // Loop over all elements in the row that maps to this row in the summand under inspection
          for (size_t elem_idx = comp_jac.rowIndices()[local_idx]; elem_idx < comp_jac.rowIndices()[local_idx + 1]; elem_idx++)
          {
            IdxT local_col_idx  = comp_jac.colIndices()[elem_idx];
            IdxT global_col_idx = permutations[summand_idx][local_col_idx];

            // Verify that this element should be added to the final sum.
            // In system jacobian construction, not all component variables map to system variables.
            if (global_col_idx != static_cast<IdxT>(-1))
            {
              // If we haven't encountered an element in one of the summands that maps to this column in the final sum
              // yet, then this is an additional non-zero element.
              if (summand_elements[global_col_idx].empty())
              {
                nnz++;
              }

              // Insert into the radix sort bucket
              summand_elements[global_col_idx].push_back({.summand_idx_ = summand_idx,
                                                          .local_idx_   = elem_idx});
            }
          }
        }

        // Now that columns have been sorted, construct the final row plan
        SharedRowPlan rowPlan;
        rowPlan.elements_.reserve(nnz);

        for (size_t col_idx = 0; col_idx < summand_elements.size(); col_idx++)
        {
          if (!summand_elements[col_idx].empty())
          {
            rowPlan.elements_.push_back({
                .summand_elements_ = std::move(summand_elements[col_idx]),
                .column_           = col_idx,
            });
          }
        }
        return rowPlan;
      }

      /**
       * @brief Create a plan for constructing an exclusive row of the final sum matrix, based on the given sparsity and permutation.
       * @param summand The sole summand which has a row which maps to this row.
       * @param summand_contribution How the summand maps to this row.
       * @param summand_permutation The permutation \f$P_i\f$ of the summand.
       */
      static ExclusiveRowPlan createExclusiveRowPlan(
          const CsrMatrix&         summand,
          LocalSummandIdx          summand_contribution,
          const std::vector<IdxT>& summand_permutation)
      {
        IdxT row_idx      = summand.rowIndices()[summand_contribution.local_idx_];
        IdxT next_row_idx = summand.rowIndices()[summand_contribution.local_idx_ + 1];

        using Element = ExclusiveRowPlan::Element;
        std::vector<Element> elements;
        elements.reserve(next_row_idx - row_idx);

        // Remove all columns which do not map to a column in the final sum
        for (size_t i = row_idx; i < next_row_idx; i++)
        {
          IdxT local_col  = summand.colIndices()[i];
          IdxT global_col = summand_permutation[local_col];

          if (global_col != static_cast<IdxT>(-1))
          {
            elements.push_back({
                .elem_idx_ = i,
                .sum_col_  = global_col,
            });
          }
        }

        return {
            .summand_idx_ = summand_contribution,
            .row_nnz_     = next_row_idx - row_idx,
            .elements_    = elements,
        };
      }

      /**
       * @brief "Inverts" `permutations` and returns a mapping from rows/columns in the final sum to rows/columns in summands.
       * @param size The size of the final sum matrix.
       */
      static std::vector<std::vector<LocalSummandIdx>> invertPermutations(const std::vector<std::vector<IdxT>>& permutations, size_t size)
      {
        std::vector<std::vector<LocalSummandIdx>> inverse(size);

        // Loop over all summands
        for (size_t summand_idx = 0; summand_idx < permutations.size(); summand_idx++)
        {
          auto summand_perm = permutations[summand_idx];
          // Loop over all local variables
          for (IdxT local_idx = 0; local_idx < summand_perm.size(); local_idx++)
          {
            // The index in the final sum
            IdxT sum_idx = summand_perm[local_idx];

            // Remove rows/columns which are not present in the final sum
            if (sum_idx != static_cast<IdxT>(-1))
            {
              inverse[sum_idx].push_back({.summand_idx_ = summand_idx, .local_idx_ = local_idx});
            }
          }
        }

        return inverse;
      }

    public:
      /**
       * @brief Construct a `CsrPermutedSum` by analyzing the sparsity and permutations of some summands.
       * @param summands A list of sparsity patterns of eventual summands.
       * @param permutations The mappings of local rows/columns in summands to global rows/columns in the final sum matrix.
       *                     Contains a separate permutation for each summand. Each permutation is a mapping from local index
       *                     to global index. A special value of -1 indicates that this row/column does not appear in the final
       *                     sum.
       * @param size The size of the final sum matrix.
       * @return A `CsrPermutedSum` object, which can be applied to compute the permuted sum of summands for as long as the sparsity patterns
       *         and permutations do of those summands do not change.
       */
      static CsrPermutedSum analyzeSparsity(const std::vector<std::reference_wrapper<CsrMatrix>>& summands, const std::vector<std::vector<IdxT>>& permutations, size_t size)
      {
        std::vector<typename CsrPermutedSum::RowPlan> row_plans;
        row_plans.reserve(size);

        assert(summands.size() == permutations.size());

        std::vector<std::vector<LocalSummandIdx>> inverse_permutation = invertPermutations(permutations, size);

        // Construct a plan for constructing each row of the final sum
        for (size_t row = 0; row < size; row++)
        {
          auto component_contributions = inverse_permutation[row];
          if (component_contributions.size() > 1)
          {
            row_plans.push_back(createSharedRowPlan(summands, component_contributions, permutations, size));
          }
          else
          {
            auto contribution = component_contributions.front();
            row_plans.push_back(createExclusiveRowPlan(summands[contribution.summand_idx_], contribution, permutations[contribution.summand_idx_]));
          }
        }

        return CsrPermutedSum(
            row_plans,
            summands.size(),
            size);
      }

      /**
       * @brief Apply the `CsrPermutedSum` by constructing the final sum matrix, based on the values contained in the `summands`.
       * @tparam CsrBuilder Some type of `CsrBuilder`
       * @param summands A list of matrices which contain the values to be summed. The sparsity patterns should match the ones given
       *                 when `analyzeSparsity()` was used to construct the `CsrPermutedSum`.
       * @param builder The builder which is used to construct the matrix.
       * @return The same `builder`.
       */
      template <class CsrBuilder>
      CsrBuilder apply(const std::vector<std::reference_wrapper<CsrMatrix>>& summands, CsrBuilder&& builder)
      {
        assert(summands.size() == num_summands_);

        // Apply all of the row plans
        for (size_t row = 0; row < row_plans_.size(); row++)
        {
          builder.row(row);

          auto row_plan = row_plans_[row];
          std::visit(
              Utility::OverloadVisitor{
                  [&](const ExclusiveRowPlan& row_plan)
                  { row_plan.apply(summands, builder); },
                  [&](const SharedRowPlan& row_plan)
                  { row_plan.apply(summands, builder); }},
              row_plan);
        }

        return std::move(builder);
      }
    };
  } // namespace LinearAlgebra
} // namespace GridKit