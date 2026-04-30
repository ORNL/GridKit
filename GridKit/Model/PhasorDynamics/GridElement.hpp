#pragma once

#include <vector>

#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Utilities/Errors.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename ScalarP, typename IdxP>
    class GridElement : public Model::Evaluator<ScalarP, IdxP>
    {
    public:
      using ScalarT = ScalarP;
      using IdxT    = IdxP;
      using RealT   = typename Model::Evaluator<ScalarT, IdxT>::RealT;
      using MatrixT = typename Model::Evaluator<ScalarT, IdxT>::MatrixT;

      GridElement() = default;

      virtual ~GridElement();

      virtual int verify() const = 0;

      IdxT size() override final
      {
        return size_;
      }

      IdxT nnz() override final
      {
        return nnz_;
      }

      void setTolerances(RealT& rtol, RealT& atol) const override;

      void setMaxSteps(IdxT& msa) const override;

      std::vector<ScalarT>& y() override
      {
        return y_;
      }

      const std::vector<ScalarT>& y() const override
      {
        return y_;
      }

      std::vector<ScalarT>& yp() override
      {
        return yp_;
      }

      const std::vector<ScalarT>& yp() const override
      {
        return yp_;
      }

      std::vector<bool>& tag() override
      {
        return tag_;
      }

      const std::vector<bool>& tag() const override
      {
        return tag_;
      }

      std::vector<ScalarT>& getResidual() override
      {
        return f_;
      }

      const std::vector<ScalarT>& getResidual() const override
      {
        return f_;
      }

      MatrixT& getJacobian() override
      {
        return J_;
      }

      const MatrixT& getJacobian() const override
      {
        return J_;
      }

      void setVariableIndex(IdxT local_index, IdxT global_index);

      IdxT& getVariableIndex(IdxT local_index);

      const std::vector<IdxT>& getVariableIndices() const
      {
        return variable_indices_;
      }

      void setResidualIndex(IdxT local_index, IdxT global_index);

      IdxT& getResidualIndex(IdxT local_index);

      const std::vector<IdxT>& getResidualIndices() const
      {
        return residual_indices_;
      }

    protected:
      IdxT              size_{0};
      IdxT              nnz_{0};
      /// Global (system-level) variable indices
      std::vector<IdxT> variable_indices_;
      /// Global (system-level) residual indices
      std::vector<IdxT> residual_indices_;

      std::vector<ScalarT> y_;
      std::vector<ScalarT> yp_;
      std::vector<bool>    tag_;
      std::vector<ScalarT> f_;
      std::vector<ScalarT> g_;

      MatrixT J_;
      IdxT*   J_rows_buffer_{nullptr};
      IdxT*   J_cols_buffer_{nullptr};
      RealT*  J_vals_buffer_{nullptr};

      RealT rel_tol_;
      RealT abs_tol_;

      IdxT max_steps_;

      //
      // Adjoint sensitivity members
      //

      std::vector<ScalarT> yB_{};
      std::vector<ScalarT> ypB_{};
      std::vector<ScalarT> fB_{};
      std::vector<ScalarT> gB_{};

      std::vector<ScalarT> param_{};
      std::vector<ScalarT> param_up_{};
      std::vector<ScalarT> param_lo_{};
    };

  } // namespace PhasorDynamics
} // namespace GridKit
