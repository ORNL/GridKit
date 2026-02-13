#pragma once

#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using Log = ::GridKit::Utilities::Logger;

    /**
     * @brief Component model implementation base class.
     */
    template <class ScalarT, typename IdxT>
    class Component : public Model::Evaluator<ScalarT, IdxT>
    {
    public:
      using RealT   = typename Model::Evaluator<ScalarT, IdxT>::RealT;
      using MatrixT = typename Model::Evaluator<ScalarT, IdxT>::MatrixT;

      Component()
        : size_(0)
      {
      }

      virtual int verify() const = 0;

      virtual IdxT size() override
      {
        return size_;
      }

      virtual IdxT nnz() override
      {
        return nnz_;
      }

      /// @todo Remove this method. It should be part of DynamicSolver class.
      virtual bool hasJacobian() override
      {
        return true;
      }

      virtual void updateTime(RealT t, RealT a) override
      {
        time_  = t;
        alpha_ = a;
        // std::cout << "updateTime: t = " << time_ << ", alpha = " << alpha_ << "\n";
      }

      std::vector<ScalarT>& absoluteTolerance() override
      {
        return abs_tol_;
      }

      const std::vector<ScalarT>& absoluteTolerance() const override
      {
        return abs_tol_;
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

      virtual int setGridKitComponentID(IdxT) = 0;

      IdxT getGridKitComponentID() const
      {
        return gridkit_component_id_;
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
      IdxT              gridkit_component_id_{0};
      IdxT              size_{0};
      IdxT              nnz_{0};
      std::vector<IdxT> variable_indices_; ///< Global (system-level) variable indices
      std::vector<IdxT> residual_indices_; ///< Global (system-level) residual indices

      std::vector<ScalarT> y_;
      std::vector<ScalarT> yp_;
      std::vector<bool>    tag_;
      std::vector<ScalarT> abs_tol_;
      std::vector<ScalarT> f_;
      std::vector<ScalarT> g_;
      std::vector<ScalarT> wb_;
      std::vector<ScalarT> h_;

      MatrixT J_;

      RealT time_;
      RealT alpha_;

      /*

      ------ WARNING: Temporary ------

      The protected variable mva_system_base_ is temporarily
      hard coded. This eventually needs to be configured
      from the input JSON format, which specifies the system MVA base.

      */

      RealT mva_system_base_{100.0};

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
