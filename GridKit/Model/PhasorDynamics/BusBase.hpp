#pragma once

#include <memory>
#include <set>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/CommonMath.hpp>
#include <GridKit/Constants.hpp>
#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>
#include <GridKit/Utilities/Errors.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using Log = ::GridKit::Utilities::Logger;

    /*!
     * @brief BusBase model implementation base class.
     *
     */
    template <typename ScalarT, typename IdxT>
    class BusBase : public Model::Evaluator<ScalarT, IdxT>
    {
    public:
      using RealT    = typename Model::Evaluator<ScalarT, IdxT>::RealT;
      using MatrixT  = typename Model::Evaluator<ScalarT, IdxT>::MatrixT;
      using BusTypeT = typename BusData<RealT, IdxT>::BusType;
      using MonitorT = Model::VariableMonitor<BusBase, BusData>;

      BusBase() = default;

      virtual ~BusBase();

      virtual int verify() const
      {
        return 0;
      }

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

      /// Pure virtual function, returns bus type (DEFAULT or SLACK).
      virtual BusTypeT BusType() const
      {
        return BusTypeT::DEFAULT;
      }

      bool hasJacobian() override
      {
        return false;
      }

      void updateTime(RealT /* t */, RealT /* a */) override
      {
        // No time to update in bus models
      }

      virtual ScalarT&       Vr()       = 0;
      virtual const ScalarT& Vr() const = 0;
      virtual ScalarT&       Vi()       = 0;
      virtual const ScalarT& Vi() const = 0;
      virtual ScalarT&       Ir()       = 0;
      virtual const ScalarT& Ir() const = 0;
      virtual ScalarT&       Ii()       = 0;
      virtual const ScalarT& Ii() const = 0;

      virtual int setBusID(IdxT) = 0;

      virtual const IdxT busID() const
      {
        return bus_id_;
      }

      const Model::VariableMonitorBase* getMonitor() const override;

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

      IdxT bus_id_{INVALID_INDEX<IdxT>};

      /// Variable monitor
      std::unique_ptr<MonitorT> monitor_;

      using NotImplementedError = GridKit::Utilities::NotImplementedError;

    public:
      // TODO: evaluate how this complies with xSDK guidelines

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
