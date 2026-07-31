#pragma once

#include <limits>
#include <vector>

#include <GridKit/Constants.hpp>
#include <GridKit/LinearAlgebra/SparseMatrix/CooMatrix.hpp>
#include <GridKit/LinearAlgebra/SparseMatrix/CsrMatrix.hpp>
#include <GridKit/LinearAlgebra/Vector/Vector.hpp>
#include <GridKit/Model/VariableMonitor.hpp>
#include <GridKit/ScalarTraits.hpp>
#include <GridKit/Utilities/Errors.hpp>

namespace GridKit
{
  namespace Model
  {
    /*!
     * @brief Abstract class describing a model.
     *
     */
    template <class scalar_type, typename index_type>
    class Evaluator
    {
    public:
      using ScalarT     = scalar_type;
      using IdxT        = index_type;
      using RealT       = typename GridKit::ScalarTraits<ScalarT>::RealT;
      using CsrMatrixT  = GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>;
      using CooMatrixT  = GridKit::LinearAlgebra::CooMatrix<RealT, IdxT>;
      using VectorT     = GridKit::LinearAlgebra::Vector<ScalarT, IdxT>;
      using RealVectorT = GridKit::LinearAlgebra::Vector<RealT, IdxT>;

      Evaluator()
      {
      }

      virtual ~Evaluator()
      {
      }

      virtual int allocate()   = 0;
      virtual int initialize() = 0;

      virtual int tagDifferentiable()
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      /**
       * @brief Compute the absolute tolerance for each variable in the model
       *
       * @param rel_tol The relative tolerance which can be used to pick the
       *        absolute tolerance.
       * @return int 0 if successful, non-zero otherwise.
       *
       * This represents a "noise" level close to zero for which pure relative
       * error cannot be used.
       */
      virtual int setAbsoluteTolerance([[maybe_unused]] RealT rel_tol)
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      virtual int evaluateResidual() = 0;

      /**
       * @brief Evaluate the model's first derivatives.
       *
       * Optimization evaluators populate both the constraint Jacobian and
       * objective gradient in this pass.
       */
      virtual int evaluateJacobian() = 0;

      virtual int evaluateIntegrand()
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      virtual int initializeAdjoint()
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      virtual int evaluateAdjointResidual()
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      // virtual int evaluateAdjointJacobian() = 0;
      virtual int evaluateAdjointIntegrand()
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      virtual IdxT size() = 0;

      virtual IdxT sizeResidual()
      {
        return size();
      }

      virtual IdxT nnz() = 0;

      /**
       * @brief Return bounds on the variables exposed through y().
       *
       * Models without explicit variable bounds are unbounded by default.
       */
      virtual int getVariableBounds(RealVectorT& lower, RealVectorT& upper)
      {
        int status  = lower.resize(size());
        status     |= upper.resize(size());

        const RealT infinity  = std::numeric_limits<RealT>::infinity();
        status               |= lower.setToConst(-infinity);
        status               |= upper.setToConst(infinity);
        return status;
      }

      /**
       * @brief Return lower and upper bounds on the residual equations.
       *
       * Residuals are equality constraints by default.
       */
      virtual int getResidualBounds(RealVectorT& lower, RealVectorT& upper)
      {
        int status  = lower.resize(sizeResidual());
        status     |= upper.resize(sizeResidual());
        status     |= lower.setToZero();
        status     |= upper.setToZero();
        return status;
      }

      /**
       * @brief Whether this evaluator defines a scalar optimization objective.
       */
      virtual bool hasObjective() const
      {
        return false;
      }

      /**
       * @brief Evaluate the scalar optimization objective.
       */
      virtual int evaluateObjective()
      {
        return 0;
      }

      virtual RealT objective() const
      {
        return RealT{0};
      }

      virtual VectorT& getObjectiveGradient()
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      virtual const VectorT& getObjectiveGradient() const
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      /**
       * @brief Is there something to monitor? Defaults to `false`
       */
      virtual bool monitoring() const
      {
        return false;
      }

      /**
       * @brief Print variables at current state
       */
      virtual void printMonitoredVariables() const
      {
      }

      /**
       * @brief Get non-owning reference to monitor
       */
      virtual const VariableMonitorBase* getMonitor() const
      {
        return nullptr;
      }

      /**
       * @brief Get monitor ready for output
       */
      virtual void startMonitor()
      {
      }

      /**
       * @brief Tell monitor to wrap up
       */
      virtual void stopMonitor()
      {
      }

      /**
       * @brief Return a pointer to the CSR Jacobian
       *
       * @todo Remove this and use CsrMatrix for jac_
       */
      virtual CsrMatrixT* getCsrJacobian() const
      {
        return nullptr;
      }

      /**
       * @brief Is the Jacobian defined. Used in IDA to determine wether DQ is used or not
       *
       * @return true
       * @return false
       */
      virtual bool hasJacobian() = 0;

      virtual IdxT sizeQuadrature()
      {
        return 0;
      }

      virtual IdxT sizeParams()
      {
        return 0;
      }

      virtual void updateTime([[maybe_unused]] RealT t, [[maybe_unused]] RealT a)
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      /**
       * @brief Get the absolute tolerance for each variable in the model
       *
       * @return a reference to the absolute tolerance vector.
       *
       * @pre `setAbsoluteTolerance` must have been called first.
       */
      virtual VectorT& absoluteTolerance()
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      /**
       * @brief Get the absolute tolerance for each variable in the model
       *
       * @return a const reference to the absolute tolerance vector.
       *
       * @pre `setAbsoluteTolerance` must have been called first.
       */
      virtual const VectorT& absoluteTolerance() const
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      virtual VectorT&       y()       = 0;
      virtual const VectorT& y() const = 0;

      virtual VectorT& yp()
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      virtual const VectorT& yp() const
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      virtual std::vector<bool>& tag()
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      virtual const std::vector<bool>& tag() const
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      virtual VectorT& yB()
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      virtual const VectorT& yB() const
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      virtual VectorT& ypB()
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      virtual const VectorT& ypB() const
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      virtual VectorT& param()
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      virtual const VectorT& param() const
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      virtual VectorT& param_up()
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      virtual const VectorT& param_up() const
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      virtual VectorT& param_lo()
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      virtual const VectorT& param_lo() const
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      virtual VectorT&       getResidual()       = 0;
      virtual const VectorT& getResidual() const = 0;

      virtual VectorT& getIntegrand()
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      virtual const VectorT& getIntegrand() const
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      virtual VectorT& getAdjointResidual()
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      virtual const VectorT& getAdjointResidual() const
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      virtual VectorT& getAdjointIntegrand()
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }

      virtual const VectorT& getAdjointIntegrand() const
      {
        throw GridKit::Utilities::NotImplementedError(__func__);
      }
    };

  } // namespace Model
} // namespace GridKit
