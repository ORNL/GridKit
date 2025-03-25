
#ifndef _LOAD_HPP_
#define _LOAD_HPP_

#include <Model/PowerFlow/ModelEvaluatorImpl.hpp>
#include <PowerSystemData.hpp>

namespace GridKit
{
  template <class ScalarT, typename IdxT>
  class BaseBus;
}

namespace GridKit
{
  /*!
   * @brief Declaration of a passive load class.
   *
   */
  template <class ScalarT, typename IdxT>
  class Load : public ModelEvaluatorImpl<ScalarT, IdxT>
  {
    using ModelEvaluatorImpl<ScalarT, IdxT>::size_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::nnz_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::time_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::alpha_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::y_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::yp_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::tag_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::f_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::g_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::yB_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::ypB_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::fB_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::gB_;
    using ModelEvaluatorImpl<ScalarT, IdxT>::param_;

    // typedef typename ModelEvaluatorImpl<ScalarT, IdxT>::real_type real_type;
    // typedef BaseBus<ScalarT, IdxT> bus_type;
    using bus_type  = BaseBus<ScalarT, IdxT>;
    using real_type = typename ModelEvaluatorImpl<ScalarT, IdxT>::real_type;
    using LoadData  = GridKit::PowerSystemData::LoadData<real_type, IdxT>;

  public:
    Load(bus_type* bus, ScalarT P, ScalarT Q);
    Load(bus_type* bus, LoadData& data);
    virtual ~Load();

    int allocate();
    int initialize();
    int tagDifferentiable();
    int evaluateResidual();
    int evaluateJacobian();
    int evaluateIntegrand();

    int initializeAdjoint();
    int evaluateAdjointResidual();
    // int evaluateAdjointJacobian();
    int evaluateAdjointIntegrand();

    void updateTime(real_type t, real_type a)
    {
      time_  = t;
      alpha_ = a;
    }

  private:
    ScalarT    P_;
    ScalarT    Q_;
    const IdxT busID_;
    bus_type*  bus_;
  };
} // namespace GridKit

#endif // _LOAD_HPP_
