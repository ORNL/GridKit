#pragma once

#include <optional>
#include <vector>

#include <AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <Model/Evaluator.hpp>
#include <Model/PhasorDynamics/SignalNode/SignalNode.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Component model implementation base class.
     */
    template <class ScalarT, typename IdxT>
    class Component : public Model::Evaluator<ScalarT, IdxT>
    {
    public:
      using real_type = typename Model::Evaluator<ScalarT, IdxT>::real_type;

      Component()
        : size_(0)
      {
      }

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
        return false;
      }

      /// Attaches a signal node to an external variable on this component
      ///
      /// @pre The object is initialized
      /// @post The signal node index corresponding to the variable specified is set to the
      ///       provided signal
      virtual void attachSignalNode([[maybe_unused]] size_t                           variable,
                                    [[maybe_unused]] const SignalNode<ScalarT, IdxT>* signal)
      {
        /// TODO: replace this interface with something significantly better. this interface has
        ///       many, many issues, among them being that there is no type safety in the way
        ///       the variable is specified; the callee must cast the variable enumeration to
        ///       the variable itself (which requires more work) or specify the integer directly
        ///       (also bad, as it makes the code harder to read and doesn't automatically update
        ///       the index when changes are made to the enumeration. additionally, since errors
        ///       resulting from this will only be caught at runtime, it is harder to debug).
        ///
        ///       another problem with this interface is that like many other methods in gridkit,
        ///       we assume that it semantically makes sense for all components to provide them
        ///       and narrow them down with runtime errors.
        ///
        ///       finally, this implementation also suffers from the problem that it requires the
        ///       non-stub implementations to write lots of boilerplate code to actually implement
        ///       it.
        ///
        ///       a prototype of a better implementation for this can be found here:
        ///       https://github.com/ORNL/GridKit/pull/193/commits/d9158691c8e4de3d5bd0269c415a76f3b7ca76c0
        ///
        ///       this implementation suffers from none of the issues described above.

        throw "No signals exist for this component";
      }

      // virtual void updateTime(real_type t, real_type a)
      // {
      //     time_ = t;
      //     alpha_ = a;
      //     std::cout << "updateTime: t = " << time_ << "\n";
      // }

      virtual void setTolerances(real_type& rtol, real_type& atol) const override
      {
        rtol = rel_tol_;
        atol = abs_tol_;
      }

      virtual void setMaxSteps(IdxT& msa) const override
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

      GridKit::LinearAlgebra::COO_Matrix<ScalarT, IdxT>& getJacobian() override
      {
        return J_;
      }

      const GridKit::LinearAlgebra::COO_Matrix<ScalarT, IdxT>& getJacobian() const override
      {
        return J_;
      }

    protected:
      IdxT size_;
      IdxT nnz_;

      std::vector<ScalarT> y_;
      std::vector<ScalarT> yp_;
      std::vector<bool>    tag_;
      std::vector<ScalarT> f_;
      std::vector<ScalarT> g_;

      GridKit::LinearAlgebra::COO_Matrix<ScalarT, IdxT> J_;

      real_type time_;
      real_type alpha_;

      real_type rel_tol_;
      real_type abs_tol_;

      IdxT max_steps_;

      IdxT component_id_;

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

      /// Vector containing signals attached to external variables
      std::vector<std::optional<const SignalNode<ScalarT, IdxT>*>> external_variable_signals_;

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

      IdxT getComponentID() const
      {
        return component_id_;
      }
    };

  } // namespace PhasorDynamics
} // namespace GridKit
