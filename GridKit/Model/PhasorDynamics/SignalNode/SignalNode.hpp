#pragma once

#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/Evaluator.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename RealT, typename IdxT>
    struct SignalNodeData;
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    /*!
     * @brief SignalNode model implementation base class.
     *
     */
    template <class ScalarT, typename IdxT>
    class SignalNode : public Model::Evaluator<ScalarT, IdxT>
    {
    public:
      using RealT   = typename Model::Evaluator<ScalarT, IdxT>::RealT;
      using MatrixT = typename Model::Evaluator<ScalarT, IdxT>::MatrixT;

      SignalNode();
      SignalNode(const SignalNodeData<RealT, IdxT>& data);

      virtual ~SignalNode() = default;

      virtual int allocate() override final;
      virtual int initialize() override final;
      virtual int tagDifferentiable() override final;
      virtual int setAbsoluteTolerance(RealT) override final;
      virtual int evaluateResidual() override final;
      virtual int evaluateJacobian() override final;

      void    set(ScalarT* signal_in, IdxT* global_index);
      bool    linked() const;
      ScalarT read() const;
      void    init(ScalarT signal_in);

      const IdxT signalId() const
      {
        return signal_id_;
      }

      virtual IdxT size() override final
      {
        return size_;
      }

      virtual IdxT nnz() override final
      {
        return nnz_;
      }

      virtual bool hasJacobian() override final
      {
        return false;
      }

      virtual void updateTime(RealT /* t */, RealT /* a */) override final
      {
        // No time to update in signal nodes
      }

      std::vector<ScalarT>& y() override final
      {
        return y_;
      }

      const std::vector<ScalarT>& y() const override final
      {
        return y_;
      }

      std::vector<ScalarT>& yp() override final
      {
        return yp_;
      }

      const std::vector<ScalarT>& yp() const override final
      {
        return yp_;
      }

      std::vector<bool>& tag() override final
      {
        return tag_;
      }

      const std::vector<bool>& tag() const override final
      {
        return tag_;
      }

      std::vector<ScalarT>& absoluteTolerance() override
      {
        return abs_tol_;
      }

      const std::vector<ScalarT>& absoluteTolerance() const override
      {
        return abs_tol_;
      }

      MatrixT& getJacobian() override final
      {
        return J_;
      }

      const MatrixT& getJacobian() const override final
      {
        return J_;
      }

      IdxT getVariableIndex() const
      {
        return *variable_index_;
      }

    private:
      ScalarT* signal_{nullptr};
      IdxT     signal_id_{0};

    protected:
      const IdxT bus_id_{INVALID_INDEX<IdxT>};

      IdxT  size_{0};
      IdxT  nnz_{0};
      IdxT* variable_index_{nullptr};

      std::vector<ScalarT> y_;
      std::vector<ScalarT> yp_;
      std::vector<bool>    tag_;
      std::vector<ScalarT> abs_tol_;
      std::vector<ScalarT> f_;

      MatrixT J_;

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
      virtual IdxT sizeQuadrature() override final
      {
        throw "ERROR: Method not implemented!\n";
        return 0;
      }

      virtual IdxT sizeParams() override final
      {
        throw "ERROR: Method not implemented!\n";
        return 0;
      }

      std::vector<ScalarT>& yB() override final
      {
        throw "ERROR: Method not implemented!\n";
        return yB_;
      }

      const std::vector<ScalarT>& yB() const override final
      {
        throw "ERROR: Method not implemented!\n";
        return yB_;
      }

      std::vector<ScalarT>& ypB() override final
      {
        throw "ERROR: Method not implemented!\n";
        return ypB_;
      }

      const std::vector<ScalarT>& ypB() const override final
      {
        throw "ERROR: Method not implemented!\n";
        return ypB_;
      }

      std::vector<ScalarT>& param() override final
      {
        throw "ERROR: Method not implemented!\n";
        return param_;
      }

      const std::vector<ScalarT>& param() const override final
      {
        throw "ERROR: Method not implemented!\n";
        return param_;
      }

      std::vector<ScalarT>& param_up() override final
      {
        throw "ERROR: Method not implemented!\n";
        return param_up_;
      }

      const std::vector<ScalarT>& param_up() const override final
      {
        throw "ERROR: Method not implemented!\n";
        return param_up_;
      }

      std::vector<ScalarT>& param_lo() override final
      {
        throw "ERROR: Method not implemented!\n";
        return param_lo_;
      }

      const std::vector<ScalarT>& param_lo() const override final
      {
        throw "ERROR: Method not implemented!\n";
        return param_lo_;
      }

      int evaluateIntegrand() override final
      {
        throw "ERROR: Method not implemented!\n";
        return 1;
      }

      int initializeAdjoint() override final
      {
        throw "ERROR: Method not implemented!\n";
        return 1;
      }

      int evaluateAdjointResidual() override final
      {
        throw "ERROR: Method not implemented!\n";
        return 1;
      }

      int evaluateAdjointIntegrand() override final
      {
        throw "ERROR: Method not implemented!\n";
        return 1;
      }

      std::vector<ScalarT>& getResidual() override final
      {
        return f_;
      }

      const std::vector<ScalarT>& getResidual() const override final
      {
        return f_;
      }

      std::vector<ScalarT>& getIntegrand() override final
      {
        throw "ERROR: Method not implemented!\n";
        return g_;
      }

      const std::vector<ScalarT>& getIntegrand() const override final
      {
        throw "ERROR: Method not implemented!\n";
        return g_;
      }

      std::vector<ScalarT>& getAdjointResidual() override final
      {
        throw "ERROR: Method not implemented!\n";
        return fB_;
      }

      const std::vector<ScalarT>& getAdjointResidual() const override final
      {
        throw "ERROR: Method not implemented!\n";
        return fB_;
      }

      std::vector<ScalarT>& getAdjointIntegrand() override final
      {
        throw "ERROR: Method not implemented!\n";
        return gB_;
      }

      const std::vector<ScalarT>& getAdjointIntegrand() const override final
      {
        throw "ERROR: Method not implemented!\n";
        return gB_;
      }

      virtual const IdxT busID() const
      {
        return bus_id_;
      }
    };

  } // namespace PhasorDynamics
} // namespace GridKit
