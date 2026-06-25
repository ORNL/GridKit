#pragma once

#include <memory>

#include <GridKit/Model/Evaluator.hpp>

namespace GridKit
{
  namespace PowerElectronics
  {
    template <typename ScalarT, typename IdxT>
    class NodeBase : public Model::Evaluator<ScalarT, IdxT>
    {
    public:
      using RealT = typename Model::Evaluator<ScalarT, IdxT>::RealT;

      NodeBase(size_t n_intern, size_t n_extern)
        : n_intern_(n_intern), n_extern_(n_extern)
      {
      }

      IdxT size() final
      {
        return static_cast<IdxT>(n_extern_ + n_intern_);
      }

      size_t getExternSize()
      {
        return n_extern_;
      }

      size_t getInternalSize()
      {
        return n_intern_;
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

      std::vector<ScalarT>& absoluteTolerance() final
      {
        return abs_tol_;
      }

      const std::vector<ScalarT>& absoluteTolerance() const final
      {
        return abs_tol_;
      }

      int setBusID(IdxT bus_id)
      {
        bus_id_ = bus_id;
        return 0;
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
        return connection_nodes_[static_cast<size_t>(local_index)];
      }

      int allocate() override
      {
        // Temporary while we use std::vector in the code
        size_t size = static_cast<size_t>(n_intern_ + n_extern_);

        // Resize component model data
        f_.resize(size);
        y_.resize(size);
        yp_.resize(size);
        tag_.resize(size);
        abs_tol_.resize(size);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        connection_nodes_ = std::make_unique<IdxT[]>(size);

        return 0;
      }

      int tagDifferentiable() final
      {
        return 0;
      }

      /**
       * @brief Compute the absolute tolerance for each variable in the model
       *
       * @param rel_tol The relative tolerance which can be used to pick the
       *        absolute tolerance.
       * @tparam ScalarT Scalar data type
       * @tparam IdxT Index data type
       * @return int 0 if successful, non-zero otherwise.
       *
       * This represents a "noise" level close to zero for which pure relative
       * error cannot be used.
       */
      int setAbsoluteTolerance(RealT rel_tol) final
      {
        std::fill(abs_tol_.begin(), abs_tol_.end(), rel_tol);
        return 0;
      }

      int initialize() override
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

      size_t n_intern_;
      size_t n_extern_;

      IdxT              nnz_{0};
      std::vector<IdxT> variable_indices_; ///< Global (system-level) variable indices
      std::vector<IdxT> residual_indices_; ///< Global (system-level) residual indices

      std::vector<ScalarT> y_;
      std::vector<ScalarT> yp_;
      std::vector<bool>    tag_;
      std::vector<ScalarT> abs_tol_;
      std::vector<ScalarT> f_;

      IdxT*  J_rows_buffer_{nullptr};
      IdxT*  J_cols_buffer_{nullptr};
      RealT* J_vals_buffer_{nullptr};

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
