#pragma once

#include <vector>

#include <GridKit/Constants.hpp>
#include <GridKit/LinearAlgebra/SparseMatrix/CooMatrix.hpp>
#include <GridKit/LinearAlgebra/SparseMatrix/CsrMatrix.hpp>
#include <GridKit/LinearAlgebra/Vector/Vector.hpp>
#include <GridKit/Model/VariableMonitor.hpp>
#include <GridKit/ScalarTraits.hpp>

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
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename GridKit::ScalarTraits<ScalarT>::RealT;
      using CsrMatrixT = GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>;
      using CooMatrixT = GridKit::LinearAlgebra::CooMatrix<RealT, IdxT>;
      using VectorT    = GridKit::LinearAlgebra::Vector<ScalarT, IdxT>;

      Evaluator()
      {
      }

      virtual ~Evaluator()
      {
      }

      /**
       * @brief Allocate all internal memory used by the model.
       *
       * @post Changing the model by changing its internal topology (i.e. \ref size(), \ref nnz(), sparsity pattern
       * in \ref getCsrJacobian(), etc.) may invalidate the model until `allocate()` is called again. Calling any
       * function that requires `allocate()` to be called in between modifying the model and reallocating may invoke
       * undefined behaviour.
       */
      virtual int allocate()                          = 0;
      virtual int initialize()                        = 0;
      /**
       * @brief Fill the \ref tag() vector with proper indicators.
       *
       * @note For some models, this may be a no-op, and the vector returned by \ref tag() will still be empty.
       *
       * @pre \ref allocate() must have been called.
       */
      virtual int tagDifferentiable()                 = 0;
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
      virtual int setAbsoluteTolerance(RealT rel_tol) = 0;

      /**
       * @brief Evaluate the model's residual function.
       *
       * For models which can be appropriately modeled in semi-explicit form, this will evaluate
       * \f[F(t, y, \dot{y}) = f(t,y) - M \dot{y},\f]
       * where \f(M\f) is the mass matrix. For models which can be modeled in Hessenberg form,
       * \f(M\f) should be a simple indicator function for if the variable is differential
       * (if \f(\dot{y}\f) appears in its equation). This indicator should be given by \ref tag().
       */
      virtual int evaluateResidual()  = 0;
      virtual int evaluateJacobian()  = 0;
      virtual int evaluateIntegrand() = 0;

      virtual int initializeAdjoint()        = 0;
      virtual int evaluateAdjointResidual()  = 0;
      // virtual int evaluateAdjointJacobian() = 0;
      virtual int evaluateAdjointIntegrand() = 0;

      /// How many variables exist in the model
      virtual IdxT size() = 0;
      /**
       * @brief How many non-zero elements are in this model's Jacobian.
       *
       * This knowledge may be needed to allocate the buffers for the model's Jacobian, so it
       * it should not depend on \ref allocate() having been called.
       *
       * @see getCsrJacobian()
       */
      virtual IdxT nnz()  = 0;

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
       * @pre \ref allocate() must have been called.
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

      virtual IdxT sizeQuadrature()             = 0;
      virtual IdxT sizeParams()                 = 0;
      virtual void updateTime(RealT t, RealT a) = 0;

      /**
       * @brief Get the absolute tolerance for each variable in the model
       *
       * @return a reference to the absolute tolerance vector.
       *
       * @pre `setAbsoluteTolerance` must have been called first.
       */
      virtual VectorT&       absoluteTolerance()       = 0;
      /**
       * @brief Get the absolute tolerance for each variable in the model
       *
       * @return a const reference to the absolute tolerance vector.
       *
       * @pre `setAbsoluteTolerance` must have been called first.
       */
      virtual const VectorT& absoluteTolerance() const = 0;

      virtual VectorT&       y()       = 0;
      virtual const VectorT& y() const = 0;

      virtual VectorT&       yp()       = 0;
      virtual const VectorT& yp() const = 0;

      /**
       * @brief An indicator for which variables in the model are differential (true) or algebraic (false).
       *
       * This information may be required for certain integrator features (such as the \ref AnalysisManager::NativeDynamicSolver::Rosenbrock
       * integrator), but is not required to be set. A reference to an empty vector may be returned.
       *
       * @pre \ref tagDifferentiable() must have been called.
       */
      virtual std::vector<bool>&       tag()       = 0;
      virtual const std::vector<bool>& tag() const = 0;

      virtual VectorT&       yB()       = 0;
      virtual const VectorT& yB() const = 0;

      virtual VectorT&       ypB()       = 0;
      virtual const VectorT& ypB() const = 0;

      virtual VectorT&       param()       = 0;
      virtual const VectorT& param() const = 0;

      virtual VectorT&       param_up()       = 0;
      virtual const VectorT& param_up() const = 0;

      virtual VectorT&       param_lo()       = 0;
      virtual const VectorT& param_lo() const = 0;

      virtual VectorT&       getResidual()       = 0;
      virtual const VectorT& getResidual() const = 0;

      virtual VectorT&       getIntegrand()       = 0;
      virtual const VectorT& getIntegrand() const = 0;

      virtual VectorT&       getAdjointResidual()       = 0;
      virtual const VectorT& getAdjointResidual() const = 0;

      virtual VectorT&       getAdjointIntegrand()       = 0;
      virtual const VectorT& getAdjointIntegrand() const = 0;
    };

  } // namespace Model
} // namespace GridKit
