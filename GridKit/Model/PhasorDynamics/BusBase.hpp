#pragma once

#include <cassert>
#include <cstddef>
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
    template <typename scalar_type, typename index_type>
    class BusBase : public Model::Evaluator<scalar_type, index_type>
    {
    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename Model::Evaluator<ScalarT, IdxT>::RealT;
      using CsrMatrixT = typename Model::Evaluator<ScalarT, IdxT>::CsrMatrixT;
      using CooMatrixT = typename Model::Evaluator<ScalarT, IdxT>::CooMatrixT;
      using BusTypeT   = typename BusData<RealT, IdxT>::BusType;
      using MonitorT   = Model::VariableMonitor<BusBase, BusData>;
      using VectorT  = typename Model::Evaluator<ScalarT, IdxT>::VectorT;

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

      VectorT& tag() override
      {
        return tag_;
      }

      const VectorT& tag() const override
      {
        return tag_;
      }

      VectorT& absoluteTolerance() override
      {
        return abs_tol_;
      }

      const VectorT& absoluteTolerance() const override
      {
        return abs_tol_;
      }

      VectorT& getResidual() override
      {
        return f_;
      }

      const VectorT& getResidual() const override
      {
        return f_;
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
      using Model::Evaluator<ScalarT, IdxT>::y_;
      using Model::Evaluator<ScalarT, IdxT>::yp_;
      using Model::Evaluator<ScalarT, IdxT>::f_;
      using Model::Evaluator<ScalarT, IdxT>::tag_;
      using Model::Evaluator<ScalarT, IdxT>::abs_tol_;
      using Model::Evaluator<ScalarT, IdxT>::offset_;
      using Model::Evaluator<ScalarT, IdxT>::allocated_;
      using Model::Evaluator<ScalarT, IdxT>::allocateVectors;

    public:
      /// Pure virtual function, returns bus type (DEFAULT or SLACK).
      virtual BusTypeT BusType() const
      {
        return BusTypeT::DEFAULT;
      }

      CsrMatrixT* getCsrJacobian() const override
      {
        return csr_jac_;
      }

      CooMatrixT* getCooJacobian() const
      {
        return coo_jac_;
      }

      bool hasJacobian() override
      {
        return true;
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

      std::vector<ScalarT> g_;

      CsrMatrixT* csr_jac_{nullptr};
      CooMatrixT* coo_jac_{nullptr};

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

      [[noreturn]] VectorT& yB() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& yB() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& ypB() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& ypB() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& param() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& param() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& param_up() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& param_up() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& param_lo() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& param_lo() const override
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

      [[noreturn]] VectorT& getIntegrand() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& getIntegrand() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& getAdjointResidual() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& getAdjointResidual() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& getAdjointIntegrand() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& getAdjointIntegrand() const override
      {
        throw NotImplementedError(__func__);
      }
    };

  } // namespace PhasorDynamics
} // namespace GridKit
