#pragma once

#include <memory>

#include <GridKit/Model/Evaluator.hpp>

namespace GridKit
{
  namespace PowerElectronics
  {
    template <typename ScalarT, typename IdxT>
    class Bus : public Model::Evaluator<ScalarT, IdxT>
    {
    public:
      using RealT   = typename Model::Evaluator<ScalarT, IdxT>::RealT;
      using MatrixT = typename Model::Evaluator<ScalarT, IdxT>::MatrixT;

      IdxT size() final
      {
        return size_;
      }

      IdxT nnz() final
      {
        return nnz_;
      }

      bool hasJacobian() final
      {
        return false;
      }

      void updateTime(RealT, RealT) final
      {
      }

      void setTolerances(RealT& rtol, RealT& atol) const final
      {
        atol = atol_;
        rtol = rtol_;
      }

      void setMaxSteps(IdxT& msa) const final
      {
        msa = max_steps_;
      }

      std::vector<ScalarT>& y() final
      {
        return y_;
      }

      const std::vector<ScalarT>& y() const final
      {
        return y_;
      }

      std::vector<ScalarT>& yp() final
      {
        return yp_;
      }

      const std::vector<ScalarT>& yp() const final
      {
        return yp_;
      }

      std::vector<bool>& tag() final
      {
        return tag_;
      }

      const std::vector<bool>& tag() const final
      {
        return tag_;
      }

      MatrixT& getJacobian() final
      {
        return J_;
      }

      const MatrixT& getJacobian() const final
      {
        return J_;
      }

      int setBusID(IdxT bus_id)
      {
        bus_id_ = bus_id;
      }

      IdxT busID() const
      {
        return bus_id_;
      }

      /**
       * @brief Create the mappings from local to global indices
       *
       * @param local_index
       * @param global_index
       * @return int
       */
      int setExternalConnectionNodes(IdxT local_index, IdxT global_index)
      {
        if (!connection_nodes_)
        {
          connection_nodes_ = std::make_unique<IdxT[]>(size_);
        }

        connection_nodes_[local_index] = global_index;
        return 0;
      }

      /**
       * @brief Given the location of value in the local vector map to global index
       *
       * f(local_index) = global_index
       *
       * @param local_index index of local value in vector
       * @return IdxT Index of the same value in the global vector
       */
      IdxT getNodeConnection(IdxT local_index) const
      {
        return connection_nodes_[local_index];
      }

      int allocate() final
      {
        // Temporary while we use std::vector in the code
        size_t size = static_cast<size_t>(size_);

        // Resize component model data
        f_.resize(size);
        y_.resize(size);
        yp_.resize(size);
        tag_.resize(size);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        return 0;
      }

      int tagDifferentiable() final
      {
        return 0;
      }

      int initialize() final
      {
        // TODO: fill this in
        return 0;
      }

      int evaluateResidual() final
      {
        return 0;
      }

      int evaluateJacobian() final
      {
        return 0;
      }

    private:
      IdxT bus_id_{INVALID_INDEX<IdxT>};

      IdxT              size_{0};
      IdxT              nnz_{0};
      std::vector<IdxT> variable_indices_; ///< Global (system-level) variable indices
      std::vector<IdxT> residual_indices_; ///< Global (system-level) residual indices

      std::vector<ScalarT> y_;
      std::vector<ScalarT> yp_;
      std::vector<bool>    tag_;
      std::vector<ScalarT> f_;

      MatrixT J_;
      IdxT*   J_rows_buffer_{nullptr};
      IdxT*   J_cols_buffer_{nullptr};
      RealT*  J_vals_buffer_{nullptr};

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

      std::unique_ptr<IdxT[]> connection_nodes_;

    public:
      virtual IdxT sizeQuadrature() final
      {
        throw "ERROR: Method not implemented!\n";
        return 0;
      }

      virtual IdxT sizeParams() final
      {
        throw "ERROR: Method not implemented!\n";
        return 0;
      }

      std::vector<ScalarT>& yB() final
      {
        throw "ERROR: Method not implemented!\n";
        return yB_;
      }

      const std::vector<ScalarT>& yB() const final
      {
        throw "ERROR: Method not implemented!\n";
        return yB_;
      }

      std::vector<ScalarT>& ypB() final
      {
        throw "ERROR: Method not implemented!\n";
        return ypB_;
      }

      const std::vector<ScalarT>& ypB() const final
      {
        throw "ERROR: Method not implemented!\n";
        return ypB_;
      }

      std::vector<ScalarT>& param() final
      {
        throw "ERROR: Method not implemented!\n";
        return param_;
      }

      const std::vector<ScalarT>& param() const final
      {
        throw "ERROR: Method not implemented!\n";
        return param_;
      }

      std::vector<ScalarT>& param_up() final
      {
        throw "ERROR: Method not implemented!\n";
        return param_up_;
      }

      const std::vector<ScalarT>& param_up() const final
      {
        throw "ERROR: Method not implemented!\n";
        return param_up_;
      }

      std::vector<ScalarT>& param_lo() final
      {
        throw "ERROR: Method not implemented!\n";
        return param_lo_;
      }

      const std::vector<ScalarT>& param_lo() const final
      {
        throw "ERROR: Method not implemented!\n";
        return param_lo_;
      }

      int evaluateIntegrand() final
      {
        throw "ERROR: Method not implemented!\n";
        return 1;
      }

      int initializeAdjoint() final
      {
        throw "ERROR: Method not implemented!\n";
        return 1;
      }

      int evaluateAdjointResidual() final
      {
        throw "ERROR: Method not implemented!\n";
        return 1;
      }

      int evaluateAdjointIntegrand() final
      {
        throw "ERROR: Method not implemented!\n";
        return 1;
      }

      std::vector<ScalarT>& getResidual() final
      {
        return f_;
      }

      const std::vector<ScalarT>& getResidual() const final
      {
        return f_;
      }

      std::vector<ScalarT>& getIntegrand() final
      {
        throw "ERROR: Method not implemented!\n";
        return g_;
      }

      const std::vector<ScalarT>& getIntegrand() const final
      {
        throw "ERROR: Method not implemented!\n";
        return g_;
      }

      std::vector<ScalarT>& getAdjointResidual() final
      {
        throw "ERROR: Method not implemented!\n";
        return fB_;
      }

      const std::vector<ScalarT>& getAdjointResidual() const final
      {
        throw "ERROR: Method not implemented!\n";
        return fB_;
      }

      std::vector<ScalarT>& getAdjointIntegrand() final
      {
        throw "ERROR: Method not implemented!\n";
        return gB_;
      }

      const std::vector<ScalarT>& getAdjointIntegrand() const final
      {
        throw "ERROR: Method not implemented!\n";
        return gB_;
      }
    };
  } // namespace PowerElectronics
} // namespace GridKit