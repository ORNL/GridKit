#pragma once

#include <vector>

#include <GridKit/LinearAlgebra/SparseMatrix/COO_Matrix.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  namespace Model
  {
    /*!
     * @brief Abstract class describing a model.
     *
     */
    template <class ScalarT, typename IdxT>
    class Evaluator
    {
    public:
      typedef typename GridKit::ScalarTraits<ScalarT>::real_type real_type;

      Evaluator()
      {
      }

      virtual ~Evaluator()
      {
      }

      virtual int allocate()          = 0;
      virtual int initialize()        = 0;
      virtual int tagDifferentiable() = 0;
      virtual int evaluateResidual()  = 0;
      virtual int evaluateJacobian()  = 0;
      virtual int evaluateIntegrand() = 0;

      virtual int initializeAdjoint()        = 0;
      virtual int evaluateAdjointResidual()  = 0;
      // virtual int evaluateAdjointJacobian() = 0;
      virtual int evaluateAdjointIntegrand() = 0;

      virtual IdxT size() = 0;
      virtual IdxT nnz()  = 0;

      /**
       * @brief Is the Jacobian defined. Used in IDA to determine wether DQ is used or not
       *
       * @return true
       * @return false
       */
      virtual bool hasJacobian() = 0;

      virtual IdxT sizeQuadrature()                                      = 0;
      virtual IdxT sizeParams()                                          = 0;
      virtual void updateTime(real_type t, real_type a)                  = 0;
      virtual void setTolerances(real_type& rtol, real_type& atol) const = 0;
      virtual void setMaxSteps(IdxT& msa) const                          = 0;

      virtual std::vector<ScalarT>&       y()       = 0;
      virtual const std::vector<ScalarT>& y() const = 0;

      virtual std::vector<ScalarT>&       yp()       = 0;
      virtual const std::vector<ScalarT>& yp() const = 0;

      virtual std::vector<bool>&       tag()       = 0;
      virtual const std::vector<bool>& tag() const = 0;

      virtual std::vector<ScalarT>&       yB()       = 0;
      virtual const std::vector<ScalarT>& yB() const = 0;

      virtual std::vector<ScalarT>&       ypB()       = 0;
      virtual const std::vector<ScalarT>& ypB() const = 0;

      virtual std::vector<ScalarT>&       param()       = 0;
      virtual const std::vector<ScalarT>& param() const = 0;

      virtual std::vector<ScalarT>&       param_up()       = 0;
      virtual const std::vector<ScalarT>& param_up() const = 0;

      virtual std::vector<ScalarT>&       param_lo()       = 0;
      virtual const std::vector<ScalarT>& param_lo() const = 0;

      virtual std::vector<ScalarT>&       getResidual()       = 0;
      virtual const std::vector<ScalarT>& getResidual() const = 0;

      /// \todo Use a different approach to store and set Jacobians
      virtual GridKit::LinearAlgebra::COO_Matrix<ScalarT, IdxT>&       getJacobian()       = 0;
      virtual const GridKit::LinearAlgebra::COO_Matrix<ScalarT, IdxT>& getJacobian() const = 0;

      virtual std::vector<ScalarT>&       getIntegrand()       = 0;
      virtual const std::vector<ScalarT>& getIntegrand() const = 0;

      virtual std::vector<ScalarT>&       getAdjointResidual()       = 0;
      virtual const std::vector<ScalarT>& getAdjointResidual() const = 0;

      virtual std::vector<ScalarT>&       getAdjointIntegrand()       = 0;
      virtual const std::vector<ScalarT>& getAdjointIntegrand() const = 0;
    };

  } // namespace Model
} // namespace GridKit
