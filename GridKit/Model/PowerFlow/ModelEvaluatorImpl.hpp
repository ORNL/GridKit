
#pragma once

#include <vector>

#include <GridKit/Model/Evaluator.hpp>

namespace GridKit
{

  /*!
   * @brief Model implementation base class.
   *
   */
  template <typename scalar_type, typename index_type>
  class ModelEvaluatorImpl : public Model::Evaluator<scalar_type, index_type>
  {
  public:
    using ScalarT = scalar_type;
    using IdxT    = index_type;
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
        size_quad_(size_quad),
        size_opt_(size_opt),
        y_(size_),
        yp_(size_),
        abs_tol_(size_),
        f_(size_),
        g_(size_quad_),
        yB_(size_),
        ypB_(size_),
        fB_(size_),
        gB_(size_opt_),
        param_(size_opt_),
        param_up_(size_opt_),
        param_lo_(size_opt_)
    {
      y_.resize(size_);
      yp_.resize(size_);
      abs_tol_.resize(size_);
      f_.resize(size_);
      g_.resize(size_quad_);
      yB_.resize(size_);
      ypB_.resize(size_);
      fB_.resize(size_);
      gB_.resize(size_opt_);
      param_.resize(size_opt_);
      param_up_.resize(size_opt_);
      param_lo_.resize(size_opt_);
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

    std::vector<bool>& tag()
    {
      return tag_;
    }

    const std::vector<bool>& tag() const
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
    void allocateVectors(IdxT size)
    {
      y_.resize(size);
      yp_.resize(size);
      abs_tol_.resize(size);
      f_.resize(size);
    }

    IdxT size_;
    IdxT nnz_;
    IdxT size_quad_;
    IdxT size_opt_;

    VectorT           y_;
    VectorT           yp_;
    std::vector<bool> tag_;
    VectorT           abs_tol_;
    VectorT           f_;
    VectorT           g_;

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
