#pragma once

#include <cassert>

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
    namespace memory = GridKit::LinearAlgebra::memory;

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

      virtual int allocate()                          = 0;
      virtual int initialize()                        = 0;
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
      virtual int evaluateResidual()                  = 0;
      virtual int evaluateJacobian()                  = 0;
      virtual int evaluateIntegrand()                 = 0;

      virtual int initializeAdjoint()        = 0;
      virtual int evaluateAdjointResidual()  = 0;
      // virtual int evaluateAdjointJacobian() = 0;
      virtual int evaluateAdjointIntegrand() = 0;

      virtual IdxT size() = 0;
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

      virtual IdxT sizeQuadrature()             = 0;
      virtual IdxT sizeParams()                 = 0;
      virtual void updateTime(RealT t, RealT a) = 0;

      /**
       * @brief Get the absolute tolerance for each variable in the model
       *
       * @return a view of the absolute tolerance vector.
       *
       * @pre `setAbsoluteTolerance` must have been called first.
       */
      virtual VectorT& absoluteTolerance()
      {
        return abs_tol_;
      }

      /**
       * @brief Get the absolute tolerance for each variable in the model
       *
       * @return a const view of the absolute tolerance vector.
       *
       * @pre `setAbsoluteTolerance` must have been called first.
       */
      virtual const VectorT& absoluteTolerance() const
      {
        return abs_tol_;
      }

      virtual VectorT& y()
      {
        return y_;
      }

      virtual const VectorT& y() const
      {
        return y_;
      }

      virtual VectorT& yp()
      {
        return yp_;
      }

      virtual const VectorT& yp() const
      {
        return yp_;
      }

      /**
       * @brief Get the differential/algebraic tag for each variable.
       *
       * Entries are scalar values for direct transfer to solver vectors.
       */
      virtual VectorT& tag()
      {
        return tag_;
      }

      /**
       * @brief Get the differential/algebraic tag for each variable.
       */
      virtual const VectorT& tag() const
      {
        return tag_;
      }

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

      virtual VectorT& getResidual()
      {
        return f_;
      }

      virtual const VectorT& getResidual() const
      {
        return f_;
      }

      virtual VectorT&       getIntegrand()       = 0;
      virtual const VectorT& getIntegrand() const = 0;

      virtual VectorT&       getAdjointResidual()       = 0;
      virtual const VectorT& getAdjointResidual() const = 0;

      virtual VectorT&       getAdjointIntegrand()       = 0;
      virtual const VectorT& getAdjointIntegrand() const = 0;

      /**
       * @brief Bind this evaluator's state and residual vectors.
       */
      virtual void bind(VectorT& y, VectorT& yp, VectorT& f, VectorT& tag, VectorT& abs_tol, IdxT offset)
      {
        const IdxT n = size();
        assert(static_cast<std::size_t>(offset) + static_cast<std::size_t>(n) <= y.size());
        assert(static_cast<std::size_t>(offset) + static_cast<std::size_t>(n) <= yp.size());
        assert(static_cast<std::size_t>(offset) + static_cast<std::size_t>(n) <= f.size());
        assert(static_cast<std::size_t>(offset) + static_cast<std::size_t>(n) <= tag.size());
        assert(static_cast<std::size_t>(offset) + static_cast<std::size_t>(n) <= abs_tol.size());

        y_ = VectorT(n);
        y_.setData(y.data() + offset);

        yp_ = VectorT(n);
        yp_.setData(yp.data() + offset);

        f_ = VectorT(n);
        f_.setData(f.data() + offset);

        tag_ = VectorT(n);
        tag_.setData(tag.data() + offset);

        abs_tol_ = VectorT(n);
        abs_tol_.setData(abs_tol.data() + offset);

        offset_    = offset;
        allocated_ = true;
      }

    protected:
      /**
       * @brief Allocate this evaluator's state and residual vectors.
       */
      void allocateVectors(IdxT n)
      {
        y_ = VectorT(n);
        y_.allocate(memory::HOST);
        y_.setDataUpdated(memory::HOST);

        yp_ = VectorT(n);
        yp_.allocate(memory::HOST);
        yp_.setDataUpdated(memory::HOST);

        f_ = VectorT(n);
        f_.allocate(memory::HOST);
        f_.setDataUpdated(memory::HOST);

        tag_ = VectorT(n);
        tag_.allocate(memory::HOST);
        tag_.setDataUpdated(memory::HOST);

        abs_tol_ = VectorT(n);
        abs_tol_.allocate(memory::HOST);
        abs_tol_.setDataUpdated(memory::HOST);

        offset_    = 0;
        allocated_ = true;
      }

      VectorT y_;
      VectorT yp_;
      VectorT f_;
      VectorT tag_;
      VectorT abs_tol_;
      IdxT    offset_{0};
      bool    allocated_{false};
    };

  } // namespace Model
} // namespace GridKit
