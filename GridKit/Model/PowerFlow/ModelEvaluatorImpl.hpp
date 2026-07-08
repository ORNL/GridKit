
#pragma once

#include <vector>

#include <GridKit/Model/Evaluator.hpp>

namespace GridKit
{

  /*!
   * @brief Model implementation base class.
   *
   */
  template <class ScalarT, typename IdxT>
  class ModelEvaluatorImpl : public Model::Evaluator<ScalarT, IdxT>
  {
  protected:
    using Model::Evaluator<ScalarT, IdxT>::y_;
    using Model::Evaluator<ScalarT, IdxT>::yp_;
    using Model::Evaluator<ScalarT, IdxT>::f_;
    using Model::Evaluator<ScalarT, IdxT>::tag_;
    using Model::Evaluator<ScalarT, IdxT>::abs_tol_;
    using Model::Evaluator<ScalarT, IdxT>::allocateVectors;

  public:
    using RealT   = typename Model::Evaluator<ScalarT, IdxT>::RealT;
    using VectorT = typename Model::Evaluator<ScalarT, IdxT>::VectorT;

    ModelEvaluatorImpl()
      : size_(0),
        size_quad_(0),
        size_opt_(0)
    {
    }

    ModelEvaluatorImpl(IdxT size, IdxT size_quad, IdxT size_opt)
      : size_(size),
        nnz_(0),
        size_quad_(size_quad),
        size_opt_(size_opt),
        g_(size_quad_),
        yB_(size_),
        ypB_(size_),
        fB_(size_),
        gB_(size_opt_),
        param_(size_opt_),
        param_up_(size_opt_),
        param_lo_(size_opt_)
    {
      allocateVectors(size_);

      g_.allocate(memory::HOST);
      g_.setDataUpdated(memory::HOST);

      yB_.allocate(memory::HOST);
      yB_.setDataUpdated(memory::HOST);

      ypB_.allocate(memory::HOST);
      ypB_.setDataUpdated(memory::HOST);

      fB_.allocate(memory::HOST);
      fB_.setDataUpdated(memory::HOST);

      gB_.allocate(memory::HOST);
      gB_.setDataUpdated(memory::HOST);

      param_.allocate(memory::HOST);
      param_.setDataUpdated(memory::HOST);

      param_up_.allocate(memory::HOST);
      param_up_.setDataUpdated(memory::HOST);

      param_lo_.allocate(memory::HOST);
      param_lo_.setDataUpdated(memory::HOST);
    }

    virtual IdxT size()
    {
      return size_;
    }

    virtual IdxT nnz()
    {
      return nnz_;
    }

    virtual bool hasJacobian()
    {
      return false;
    }

    virtual IdxT sizeQuadrature()
    {
      return size_quad_;
    }

    virtual IdxT sizeParams()
    {
      return size_opt_;
    }

    // virtual void updateTime(RealT t, RealT a)
    // {
    //     time_ = t;
    //     alpha_ = a;
    //     std::cout << "updateTime: t = " << time_ << "\n";
    // }

    virtual void setMaxSteps(IdxT& msa) const
    {
      msa = max_steps_;
    }

    VectorT& y()
    {
      return y_;
    }

    const VectorT& y() const
    {
      return y_;
    }

    VectorT& yp()
    {
      return yp_;
    }

    const VectorT& yp() const
    {
      return yp_;
    }

    VectorT& tag()
    {
      return tag_;
    }

    const VectorT& tag() const
    {
      return tag_;
    }

    VectorT& absoluteTolerance()
    {
      return abs_tol_;
    }

    const VectorT& absoluteTolerance() const
    {
      return abs_tol_;
    }

    VectorT& yB()
    {
      return yB_;
    }

    const VectorT& yB() const
    {
      return yB_;
    }

    VectorT& ypB()
    {
      return ypB_;
    }

    const VectorT& ypB() const
    {
      return ypB_;
    }

    VectorT& param()
    {
      return param_;
    }

    const VectorT& param() const
    {
      return param_;
    }

    VectorT& param_up()
    {
      return param_up_;
    }

    const VectorT& param_up() const
    {
      return param_up_;
    }

    VectorT& param_lo()
    {
      return param_lo_;
    }

    const VectorT& param_lo() const
    {
      return param_lo_;
    }

    VectorT& getResidual()
    {
      return f_;
    }

    const VectorT& getResidual() const
    {
      return f_;
    }

    VectorT& getIntegrand()
    {
      return g_;
    }

    const VectorT& getIntegrand() const
    {
      return g_;
    }

    VectorT& getAdjointResidual()
    {
      return fB_;
    }

    const VectorT& getAdjointResidual() const
    {
      return fB_;
    }

    VectorT& getAdjointIntegrand()
    {
      return gB_;
    }

    const VectorT& getAdjointIntegrand() const
    {
      return gB_;
    }

    //@todo Fix ID naming
    IdxT getIDcomponent()
    {
      return idc_;
    }

  protected:
    IdxT size_;
    IdxT nnz_;
    IdxT size_quad_;
    IdxT size_opt_;

    VectorT g_;

    VectorT yB_;
    VectorT ypB_;
    VectorT fB_;
    VectorT gB_;

    VectorT param_;
    VectorT param_up_;
    VectorT param_lo_;

    RealT time_;
    RealT alpha_;

    IdxT max_steps_;

    IdxT idc_;
  };

} // namespace GridKit
