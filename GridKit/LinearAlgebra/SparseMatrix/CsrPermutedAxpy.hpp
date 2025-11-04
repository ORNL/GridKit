#include <GridKit/LinearAlgebra/SparseMatrix/CsrMatrix.hpp>

namespace GridKit
{
  namespace LinearAlgebra
  {
    template <class ScalarT, typename IdxT>
    class CsrPermutedAxpy
    {
    private:
      using CsrMatrix = LinearAlgebra::CsrMatrix<ScalarT, IdxT>;

      struct LocalSummandIdx
      {
        // The index of the summand matrix that this refers to
        size_t summand_idx_;
        // The index of the row/column inside that summand that this refers to
        IdxT   local_idx_;
      };

      /**
       * @brief A plan for constructing a row of the final sum which is calculated exclusively by one of the summands
       */
      struct ExclusiveRowPlan
      {
        struct Element
        {
          /**
           * @brief The index of the element in the component jacobian
           */
          size_t elem_idx_;
          /**
           * @brief The column of the element in the final sum
           */
          IdxT   sum_col_;
        };

        /**
         * @brief The index of the summand that this equation belongs to
         */
        size_t               summand_idx_;
        /**
         * @brief The index of the row of the summand that this row corresponds to.
         */
        size_t               row_idx_;
        /**
         * @brief The expected number of nonzero elements in the row of the summand.
         */
        IdxT                 row_nnz_;
        /**
         * @brief The elements which belong to this row
         */
        std::vector<Element> elements_;

        template <class CsrBuilder>
        void apply(const std::vector<std::reference_wrapper<CsrMatrix>>& summands, CsrBuilder& builder) const
        {
          const CsrMatrix& summand = summands[summand_idx_];

          IdxT row_start = summand.rowIndices()[row_idx_];
          IdxT row_end   = summand.rowIndices()[row_idx_ + 1];

          assert(summand.rowIndices().size() >= row_idx_ + 2);
          assert(row_end - row_start == row_nnz_);

          for (auto element : elements_)
          {
            builder.elem(element.sum_col_, summand.values()[element.elem_idx_]);
          }
        }
      };

      /**
       * @brief A plan for constructing a row of the final sum which is calculated by multiple summands
       */
      struct SharedRowPlan
      {
        /**
         * @brief A single element in the row. An element is the sum of the elements of (potentially) several summands
         */
        struct SumElement
        {
          /**
           * @brief A single element in a summand's row.
           */
          struct SummandElement
          {
            /**
             * @brief The index of the summand that this element belongs to.
             */
            size_t summand_idx_;
            /**
             * @brief The index of the element inside the summand's `CsrMatrix::values()` and `CsrMatrix::column_indices()` lists.
             */
            size_t element_idx_;
            /**
             * @brief The expected column of the element inside the summand.
             */
            IdxT   column_idx_;
          };

          /**
           * @brief A list of all elements in summands which sum to this element in the final sum.
           */
          std::vector<SummandElement> summand_elements_;
          /**
           * @brief The column of this element in the final sum.
           */
          IdxT                        column_;
        };

        /**
         * @brief A list of all elements on this row of the final sum.
         */
        std::vector<SumElement> elements_;

        template <class CsrBuilder>
        void apply(const std::vector<std::reference_wrapper<CsrMatrix>>& summands, CsrBuilder& builder) const
        {
          for (auto& element : elements_)
          {
            ScalarT sum = 0;

            for (auto& element : element.summand_elements_)
            {
              const CsrMatrix& summand = summands[element.summand_idx_];

              assert(summand.colIndices()[element.element_idx_] == element.column_idx_);

              sum += summand.values()[element.element_idx_];
            }

            builder.elem(element.column_, sum);
          }
        }
      };

      using RowPlan = std::variant<ExclusiveRowPlan, SharedRowPlan>;

      /**
       * @brief A plan for assembling each row
       */
      std::vector<RowPlan> row_plans_;

      /**
       * @brief Number of summands which are part of the sum.
       */
      size_t num_summands_;

      /**
       * @brief The dimension of the sum
       */
      size_t summand_size_;

      CsrPermutedAxpy(std::vector<RowPlan> row_plans, size_t num_summands, size_t summand_size) : row_plans_(row_plans), num_summands_(num_summands), summand_size_(summand_size)
      {
      }

      static SharedRowPlan createSharedRowPlan(
          const std::vector<std::reference_wrapper<CsrMatrix>>& summands,
          const std::vector<LocalSummandIdx>&                   summand_contributions,
          const std::vector<std::vector<IdxT>>&                 permutations,
          size_t                                                size)
      {
        using SummandElement = SharedRowPlan::SumElement::SummandElement;
        std::vector<std::vector<SummandElement>>
               summand_elements(size);
        size_t nnz = 0;

        for (const auto& contribution : summand_contributions)
        {
          auto [summand_idx, local_idx] = contribution;
          const CsrMatrix& comp_jac     = summands[summand_idx];

          for (size_t elem_idx = comp_jac.rowIndices()[local_idx]; elem_idx < comp_jac.rowIndices()[local_idx + 1]; elem_idx++)
          {
            IdxT local_col_idx  = comp_jac.colIndices()[elem_idx];
            IdxT global_col_idx = permutations[summand_idx][local_col_idx];

            if (global_col_idx != static_cast<IdxT>(-1))
            {
              if (summand_elements[global_col_idx].empty())
              {
                nnz++;
              }

              summand_elements[global_col_idx].push_back({.summand_idx_ = summand_idx,
                                                          .element_idx_ = elem_idx,
                                                          .column_idx_  = local_col_idx});
            }
          }
        }

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

      static ExclusiveRowPlan createExclusiveRowPlan(
          const CsrMatrix&         summand,
          size_t                   summand_idx,
          IdxT                     local_idx,
          const std::vector<IdxT>& summand_permutation)
      {
        IdxT row_idx      = summand.rowIndices()[local_idx];
        IdxT next_row_idx = summand.rowIndices()[local_idx + 1];

        using Element = ExclusiveRowPlan::Element;
        std::vector<Element> elements;
        elements.reserve(next_row_idx - row_idx);

        for (size_t i = row_idx; i < next_row_idx; i++)
        {
          IdxT local_col  = summand.colIndices()[i];
          IdxT global_col = summand_permutation[local_col];

          if (global_col != static_cast<IdxT>(-1))
          {
            elements.push_back(Element{
                .elem_idx_ = i,
                .sum_col_  = global_col,
            });
          }
        }

        return ExclusiveRowPlan{
            .summand_idx_ = summand_idx,
            .row_idx_     = local_idx,
            .row_nnz_     = next_row_idx - row_idx,
            .elements_    = elements,
        };
      }

      /**
       * @brief Returns a vector which indicates which system equations are external (true)
       */
      static std::vector<std::vector<LocalSummandIdx>> invertPermutations(const std::vector<std::vector<IdxT>>& permutations, size_t size)
      {
        std::vector<std::vector<LocalSummandIdx>> inverse(size);

        for (size_t summand_idx = 0; summand_idx < permutations.size(); summand_idx++)
        {
          auto summand_perm = permutations[summand_idx];
          for (IdxT local_idx = 0; local_idx < summand_perm.size(); local_idx++)
          {
            // The index in the final sum
            IdxT sum_idx = summand_perm[local_idx];

            if (sum_idx != static_cast<IdxT>(-1))
            {
              inverse[sum_idx].push_back({.summand_idx_ = summand_idx, .local_idx_ = local_idx});
            }
          }
        }

        return inverse;
      }

    public:
      static CsrPermutedAxpy analyzeSparsity(const std::vector<std::reference_wrapper<CsrMatrix>>& summands, const std::vector<std::vector<IdxT>>& permutations, size_t size)
      {
        std::vector<typename CsrPermutedAxpy::RowPlan> row_plans;
        row_plans.reserve(size);

        std::vector<std::vector<LocalSummandIdx>> inverse_permutation = invertPermutations(permutations, size);

        for (size_t row = 0; row < inverse_permutation.size(); row++)
        {
          auto component_contributions = inverse_permutation[row];
          if (component_contributions.size() > 1)
          {
            row_plans.push_back(createSharedRowPlan(summands, component_contributions, permutations, size));
          }
          else
          {
            auto [comp_idx, local_idx] = component_contributions.front();
            row_plans.push_back(createExclusiveRowPlan(summands[comp_idx], comp_idx, local_idx, permutations[comp_idx]));
          }
        }

        return CsrPermutedAxpy(
            row_plans,
            permutations.size(),
            size);
      }

      template <class CsrBuilder>
      CsrBuilder apply(const std::vector<std::reference_wrapper<CsrMatrix>>& summands, CsrBuilder&& builder)
      {
        assert(summands.size() == num_summands_);

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