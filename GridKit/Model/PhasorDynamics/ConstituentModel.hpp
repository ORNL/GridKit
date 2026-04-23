#pragma once

#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Utilities/Errors.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using Log = ::GridKit::Utilities::Logger;

    /**
     * @brief Model base class for all system constituents
     */
    template <class ScalarT, typename IdxT>
    class ConstituentModel : public Model::Evaluator<ScalarT, IdxT>
    {
    public:
      using RealT   = typename Model::Evaluator<ScalarT, IdxT>::RealT;
      using MatrixT = typename Model::Evaluator<ScalarT, IdxT>::MatrixT;

      ConstituentModel()
        : size_{0}
      {
      }

      virtual ~ConstituentModel()
      {
        if (J_rows_buffer_ != nullptr)
        {
          delete[] J_rows_buffer_;
          delete[] J_cols_buffer_;
          delete[] J_vals_buffer_;
          J_rows_buffer_ = nullptr;
          J_cols_buffer_ = nullptr;
          J_vals_buffer_ = nullptr;
        }
      }

      virtual int verify() const = 0;

      IdxT size() override final
      {
        return size_;
      }

      IdxT nnz() override final
      {
        return nnz_;
      }

      void setTolerances(RealT& rtol, RealT& atol) const override
      {
        rtol = rel_tol_;
        atol = abs_tol_;
      }

      void setMaxSteps(IdxT& msa) const override
      {
        msa = max_steps_;
      }

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

      int setVariableIndex(IdxT local_index, IdxT global_index)
      {
        variable_indices_[static_cast<size_t>(local_index)] = global_index;
        return 0;
      }

      IdxT& getVariableIndex(IdxT local_index)
      {
        return variable_indices_[static_cast<size_t>(local_index)];
      }

      const std::vector<IdxT>& getVariableIndices() const
      {
        return variable_indices_;
      }

      int setResidualIndex(IdxT local_index, IdxT global_index)
      {
        residual_indices_[static_cast<size_t>(local_index)] = global_index;
        return 0;
      }

      IdxT& getResidualIndex(IdxT local_index)
      {
        return residual_indices_[static_cast<size_t>(local_index)];
      }

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

      using NotImplementedError = GridKit::Utilities::NotImplementedError;

    public:
      [[noreturn]] IdxT sizeQuadrature() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] IdxT sizeParams() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] std::vector<ScalarT>& yB() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const std::vector<ScalarT>& yB() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] std::vector<ScalarT>& ypB() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const std::vector<ScalarT>& ypB() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] std::vector<ScalarT>& param() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const std::vector<ScalarT>& param() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] std::vector<ScalarT>& param_up() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const std::vector<ScalarT>& param_up() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] std::vector<ScalarT>& param_lo() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const std::vector<ScalarT>& param_lo() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] int evaluateIntegrand() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] int initializeAdjoint() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] int evaluateAdjointResidual() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] int evaluateAdjointIntegrand() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] std::vector<ScalarT>& getIntegrand() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const std::vector<ScalarT>& getIntegrand() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] std::vector<ScalarT>& getAdjointResidual() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const std::vector<ScalarT>& getAdjointResidual() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] std::vector<ScalarT>& getAdjointIntegrand() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const std::vector<ScalarT>& getAdjointIntegrand() const override
      {
        throw NotImplementedError(__func__);
      }
    };

  } // namespace PhasorDynamics
} // namespace GridKit
