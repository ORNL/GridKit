
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
      allocateHost(y_);
      allocateHost(yp_);
      allocateHost(abs_tol_);
      allocateHost(f_);
      allocateHost(g_);
      allocateHost(yB_);
      allocateHost(ypB_);
      allocateHost(fB_);
      allocateHost(gB_);
      allocateHost(param_);
      allocateHost(param_up_);
      allocateHost(param_lo_);
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
      resizeAndAllocateHost(y_, size);
      resizeAndAllocateHost(yp_, size);
      resizeAndAllocateHost(abs_tol_, size);
      resizeAndAllocateHost(f_, size);
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

  private:
    static void allocateHost(VectorT& vector)
    {
      vector.allocate(memory::HOST);
      vector.setToZero(memory::HOST);
    }

    static void resizeAndAllocateHost(VectorT& vector, IdxT size)
    {
      vector.resize(size);
      allocateHost(vector);
    }
  };

} // namespace GridKit
