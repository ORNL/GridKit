#pragma once

#include <memory>
#include <vector>

#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Model/PowerElectronics/ExternalConnection.hpp>

namespace GridKit
{
  namespace PowerElectronics
  {
    template <typename ScalarT, typename IdxT>
    class NodeBase : public Model::Evaluator<ScalarT, IdxT>
    {
    public:
      using RealT   = typename Model::Evaluator<ScalarT, IdxT>::RealT;
      using VectorT = typename Model::Evaluator<ScalarT, IdxT>::VectorT;

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

      VectorT& y() final
      {
        return y_;
      }

      const VectorT& y() const final
      {
        return y_;
      }

      VectorT& yp() final
      {
        return yp_;
      }

      const VectorT& yp() const final
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

      VectorT& absoluteTolerance() final
      {
        return abs_tol_;
      }

      const VectorT& absoluteTolerance() const final
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
       * @brief Create the mappings from local to global indices for a node variable (either internal or external),
       * to be used from an attached component. \see CircuitComponent::setExternalConnectionNodes()
       *
       * @param local_index The index of the local variable
       * @param connection The necessary connection information for the variable
       *
       * @pre `local_index` *must* be the index of an external variable. As of now, using this method
       * to set information for a local variable will silently discard the unnecessary information, but
       * this may change in the future.
       */
      int setExternalConnectionNodes(IdxT local_index, ExternalConnection<ScalarT, IdxT> connection)
      {
        y_ext_[local_index]            = connection.y_;
        yp_ext_[local_index]           = connection.yp_;
        f_ext_[local_index]            = connection.f_;
        connection_nodes_[local_index] = connection.idx_;
        return 0;
      }

      /**
       * @brief Get connection information for a particular variable, to be consumed by an attached
       * component so they can properly access their externals.
       */
      ExternalConnection<ScalarT, IdxT> getNodeConnection(size_t local_index) const
      {
        return ExternalConnection<ScalarT, IdxT>{
            .y_   = y_ext_[local_index],
            .yp_  = yp_ext_[local_index],
            .f_   = f_ext_[local_index],
            .idx_ = connection_nodes_[local_index],
        };
      }

      int allocate() override
      {
        size_t size = static_cast<size_t>(n_intern_ + n_extern_);

        if (!allocated_)
        {
          allocateVectors(static_cast<IdxT>(size));
        }

        tag_.resize(size, false);
        variable_indices_.resize(size);
        residual_indices_.resize(size);

        y_ext_            = std::make_unique<const ScalarT*[]>(size);
        yp_ext_           = std::make_unique<const ScalarT*[]>(size);
        f_ext_            = std::make_unique<ScalarT*[]>(size);
        connection_nodes_ = std::make_unique<IdxT[]>(size);

        allocated_ = true;
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
        abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
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

    protected:
      IdxT bus_id_{INVALID_INDEX<IdxT>};

      size_t n_intern_;
      size_t n_extern_;

      IdxT              nnz_{0};
      std::vector<IdxT> variable_indices_; ///< Global (system-level) variable indices
      std::vector<IdxT> residual_indices_; ///< Global (system-level) residual indices

      std::vector<bool> tag_;
      VectorT           abs_tol_;

      /// @brief A pointer to the internal variables of this component.
      const ScalarT* y_int_;
      /// @brief A pointer to the internal derivatives of this component.
      const ScalarT* yp_int_;
      /// @brief A pointer to the internal residuals of this component
      ScalarT*       f_int_;

      std::unique_ptr<const ScalarT*[]> y_ext_;
      std::unique_ptr<const ScalarT*[]> yp_ext_;
      std::unique_ptr<ScalarT*[]>       f_ext_;

      IdxT*  J_rows_buffer_{nullptr};
      IdxT*  J_cols_buffer_{nullptr};
      RealT* J_vals_buffer_{nullptr};

      //
      // Adjoint sensitivity members
      //

      VectorT g_{};
      VectorT yB_{};
      VectorT ypB_{};
      VectorT fB_{};
      VectorT gB_{};

      VectorT param_{};
      VectorT param_up_{};
      VectorT param_lo_{};

      std::unique_ptr<IdxT[]> connection_nodes_;

      bool allocated_{false};

    private:
      VectorT y_;
      VectorT yp_;
      VectorT f_;

    public:
      virtual IdxT
      sizeQuadrature() final
      {
        throw "ERROR: Method not implemented!\n";
        return 0;
      }

      virtual IdxT sizeParams() final
      {
        throw "ERROR: Method not implemented!\n";
        return 0;
      }

      VectorT& yB() final
      {
        throw "ERROR: Method not implemented!\n";
        return yB_;
      }

      const VectorT& yB() const final
      {
        throw "ERROR: Method not implemented!\n";
        return yB_;
      }

      VectorT& ypB() final
      {
        throw "ERROR: Method not implemented!\n";
        return ypB_;
      }

      const VectorT& ypB() const final
      {
        throw "ERROR: Method not implemented!\n";
        return ypB_;
      }

      VectorT& param() final
      {
        throw "ERROR: Method not implemented!\n";
        return param_;
      }

      const VectorT& param() const final
      {
        throw "ERROR: Method not implemented!\n";
        return param_;
      }

      VectorT& param_up() final
      {
        throw "ERROR: Method not implemented!\n";
        return param_up_;
      }

      const VectorT& param_up() const final
      {
        throw "ERROR: Method not implemented!\n";
        return param_up_;
      }

      VectorT& param_lo() final
      {
        throw "ERROR: Method not implemented!\n";
        return param_lo_;
      }

      const VectorT& param_lo() const final
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

      VectorT& getResidual() final
      {
        return f_;
      }

      const VectorT& getResidual() const final
      {
        return f_;
      }

      VectorT& getIntegrand() final
      {
        throw "ERROR: Method not implemented!\n";
        return g_;
      }

      const VectorT& getIntegrand() const final
      {
        throw "ERROR: Method not implemented!\n";
        return g_;
      }

      VectorT& getAdjointResidual() final
      {
        throw "ERROR: Method not implemented!\n";
        return fB_;
      }

      const VectorT& getAdjointResidual() const final
      {
        throw "ERROR: Method not implemented!\n";
        return fB_;
      }

      VectorT& getAdjointIntegrand() final
      {
        throw "ERROR: Method not implemented!\n";
        return gB_;
      }

      const VectorT& getAdjointIntegrand() const final
      {
        throw "ERROR: Method not implemented!\n";
        return gB_;
      }

    private:
      void allocateVectors(IdxT n)
      {
        abs_tol_.resize(n);
      }
    };
  } // namespace PowerElectronics
} // namespace GridKit
