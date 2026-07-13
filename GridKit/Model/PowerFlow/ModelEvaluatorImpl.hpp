
#pragma once

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
    using RealT      = typename Model::Evaluator<ScalarT, IdxT>::RealT;
    using VectorT    = typename Model::Evaluator<ScalarT, IdxT>::VectorT;
    using TagVectorT = typename Model::Evaluator<ScalarT, IdxT>::TagVectorT;

    ModelEvaluatorImpl()
      : size_(0),
        nnz_(0),
        size_quad_(0),
        size_opt_(0)
    {
    }

    ModelEvaluatorImpl(IdxT size, IdxT size_quad, IdxT size_opt)
      : size_(size),
        nnz_(0),
        size_quad_(size_quad),
        size_opt_(size_opt),
        y_(size_),
        yp_(size_),
        tag_(size_),
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
      allocateHost(tag_);
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

    IdxT size() override
    {
      return size_;
    }

    IdxT nnz() override
    {
      return nnz_;
    }

    bool hasJacobian() override
    {
      return false;
    }

    IdxT sizeQuadrature() override
    {
      return size_quad_;
    }

    IdxT sizeParams() override
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

    VectorT& y() override
    {
      return y_;
    }

    const VectorT& y() const override
    {
      return y_;
    }

    VectorT& yp() override
    {
      return yp_;
    }

    const VectorT& yp() const override
    {
      return yp_;
    }

    TagVectorT& tag() override
    {
      return tag_;
    }

    const TagVectorT& tag() const override
    {
      return tag_;
    }

    VectorT& absoluteTolerance() override
    {
      return abs_tol_;
    }

    const VectorT& absoluteTolerance() const override
    {
      return abs_tol_;
    }

    VectorT& yB() override
    {
      return yB_;
    }

    const VectorT& yB() const override
    {
      return yB_;
    }

    VectorT& ypB() override
    {
      return ypB_;
    }

    const VectorT& ypB() const override
    {
      return ypB_;
    }

    VectorT& param() override
    {
      return param_;
    }

    const VectorT& param() const override
    {
      return param_;
    }

    VectorT& param_up() override
    {
      return param_up_;
    }

    const VectorT& param_up() const override
    {
      return param_up_;
    }

    VectorT& param_lo() override
    {
      return param_lo_;
    }

    const VectorT& param_lo() const override
    {
      return param_lo_;
    }

    VectorT& getResidual() override
    {
      return f_;
    }

    const VectorT& getResidual() const override
    {
      return f_;
    }

    VectorT& getIntegrand() override
    {
      return g_;
    }

    const VectorT& getIntegrand() const override
    {
      return g_;
    }

    VectorT& getAdjointResidual() override
    {
      return fB_;
    }

    const VectorT& getAdjointResidual() const override
    {
      return fB_;
    }

    VectorT& getAdjointIntegrand() override
    {
      return gB_;
    }

    const VectorT& getAdjointIntegrand() const override
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
      resizeAndAllocateHost(tag_, size);
      resizeAndAllocateHost(abs_tol_, size);
      resizeAndAllocateHost(f_, size);
    }

    IdxT size_;
    IdxT nnz_;
    IdxT size_quad_;
    IdxT size_opt_;

    VectorT    y_;
    VectorT    yp_;
    TagVectorT tag_;
    VectorT    abs_tol_;
    VectorT    f_;
    VectorT    g_;

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
    template <typename ValueT>
    static void allocateHost(GridKit::LinearAlgebra::Vector<ValueT, IdxT>& vector)
    {
      vector.allocate(memory::HOST);
      vector.setToZero(memory::HOST);
    }

    template <typename ValueT>
    static void resizeAndAllocateHost(GridKit::LinearAlgebra::Vector<ValueT, IdxT>& vector, IdxT size)
    {
      vector.resize(size);
      allocateHost(vector);
    }
  };

} // namespace GridKit
