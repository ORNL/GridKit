#pragma once

#include <memory>
#include <set>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusData.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
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

      explicit BusBase(const BusData<RealT, IdxT>& data);

      virtual ~BusBase();

      /// Pure virtual function, returns bus type (DEFAULT or SLACK).
      virtual BusTypeT BusType() const
      {
        return BusTypeT::DEFAULT;
      }

      virtual IdxT size() override
      {
        return size_;
      }

      virtual IdxT nnz() override
      {
        return nnz_;
      }

      virtual bool hasJacobian() override
      {
        return false;
      }

      virtual void updateTime(RealT /* t */, RealT /* a */) override
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

      std::vector<ScalarT>& absoluteTolerance() override
      {
        return absTol_;
      }

      const std::vector<ScalarT>& absoluteTolerance() const override
      {
        return absTol_;
      }

      MatrixT& getJacobian() override
      {
        return J_;
      }

      const MatrixT& getJacobian() const override
      {
        return J_;
      }

      virtual int setBusID(IdxT) = 0;

      virtual const IdxT busID() const
      {
        return bus_id_;
      }

      int setVariableIndex(IdxT local_index, IdxT global_index)
      {
        variable_indices_[local_index] = global_index;
        return 0;
      }

      IdxT getVariableIndex(IdxT local_index) const
      {
        return variable_indices_.at(local_index);
      }

      const std::map<IdxT, IdxT>& getVariableIndices() const
      {
        return variable_indices_;
      }

      int setResidualIndex(IdxT local_index, IdxT global_index)
      {
        residual_indices_[local_index] = global_index;
        return 0;
      }

      IdxT getResidualIndex(IdxT local_index) const
      {
        return residual_indices_.at(local_index);
      }

      const std::map<IdxT, IdxT>& getResidualIndices() const
      {
        return residual_indices_;
      }

      const Model::VariableMonitorBase* getMonitor() const override;

    protected:
      IdxT bus_id_{static_cast<IdxT>(-1)};

      IdxT                 size_{0};
      IdxT                 nnz_{0};
      std::map<IdxT, IdxT> variable_indices_; ///< Map between local and global (system-level)
                                              /// variable indices
      std::map<IdxT, IdxT> residual_indices_; ///< Map between local and global (system-level)
                                              /// residual indices

      /// Variable monitor
      std::unique_ptr<MonitorT> monitor_;

      std::vector<ScalarT> y_;
      std::vector<ScalarT> yp_;
      std::vector<bool>    tag_;
      std::vector<ScalarT> absTol_;
      std::vector<ScalarT> f_;

      MatrixT J_;

      RealT time_;
      RealT alpha_;

      RealT rtol_;
      RealT atol_;

      IdxT max_steps_;

      //
      // Adjoint sensitivity members
      //

      std::vector<ScalarT> g_{};
      std::vector<ScalarT> yB_{};
      std::vector<ScalarT> ypB_{};
      std::vector<ScalarT> fB_{};
      std::vector<ScalarT> gB_{};

      std::vector<ScalarT> param_{};
      std::vector<ScalarT> param_up_{};
      std::vector<ScalarT> param_lo_{};

      //
      // Public adjoint sensitivity methods (not yet implemented in components)
      //

    public:
      virtual IdxT sizeQuadrature() override
      {
        throw "ERROR: Method not implemented!\n";
        return 0;
      }

      virtual IdxT sizeParams() override
      {
        throw "ERROR: Method not implemented!\n";
        return 0;
      }

      std::vector<ScalarT>& yB() override
      {
        throw "ERROR: Method not implemented!\n";
        return yB_;
      }

      const std::vector<ScalarT>& yB() const override
      {
        throw "ERROR: Method not implemented!\n";
        return yB_;
      }

      std::vector<ScalarT>& ypB() override
      {
        throw "ERROR: Method not implemented!\n";
        return ypB_;
      }

      const std::vector<ScalarT>& ypB() const override
      {
        throw "ERROR: Method not implemented!\n";
        return ypB_;
      }

      std::vector<ScalarT>& param() override
      {
        throw "ERROR: Method not implemented!\n";
        return param_;
      }

      const std::vector<ScalarT>& param() const override
      {
        throw "ERROR: Method not implemented!\n";
        return param_;
      }

      std::vector<ScalarT>& param_up() override
      {
        throw "ERROR: Method not implemented!\n";
        return param_up_;
      }

      const std::vector<ScalarT>& param_up() const override
      {
        throw "ERROR: Method not implemented!\n";
        return param_up_;
      }

      std::vector<ScalarT>& param_lo() override
      {
        throw "ERROR: Method not implemented!\n";
        return param_lo_;
      }

      const std::vector<ScalarT>& param_lo() const override
      {
        throw "ERROR: Method not implemented!\n";
        return param_lo_;
      }

      int evaluateIntegrand() override
      {
        throw "ERROR: Method not implemented!\n";
        return 1;
      }

      int initializeAdjoint() override
      {
        throw "ERROR: Method not implemented!\n";
        return 1;
      }

      int evaluateAdjointResidual() override
      {
        throw "ERROR: Method not implemented!\n";
        return 1;
      }

      int evaluateAdjointIntegrand() override
      {
        throw "ERROR: Method not implemented!\n";
        return 1;
      }

      std::vector<ScalarT>& getResidual() override
      {
        return f_;
      }

      const std::vector<ScalarT>& getResidual() const override
      {
        return f_;
      }

      std::vector<ScalarT>& getIntegrand() override
      {
        throw "ERROR: Method not implemented!\n";
        return g_;
      }

      const std::vector<ScalarT>& getIntegrand() const override
      {
        throw "ERROR: Method not implemented!\n";
        return g_;
      }

      std::vector<ScalarT>& getAdjointResidual() override
      {
        throw "ERROR: Method not implemented!\n";
        return fB_;
      }

      const std::vector<ScalarT>& getAdjointResidual() const override
      {
        throw "ERROR: Method not implemented!\n";
        return fB_;
      }

      std::vector<ScalarT>& getAdjointIntegrand() override
      {
        throw "ERROR: Method not implemented!\n";
        return gB_;
      }

      const std::vector<ScalarT>& getAdjointIntegrand() const override
      {
        throw "ERROR: Method not implemented!\n";
        return gB_;
      }
    };

  } // namespace PhasorDynamics
} // namespace GridKit
