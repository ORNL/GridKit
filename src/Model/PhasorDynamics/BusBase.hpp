#pragma once

#include <vector>

#include <Model/Evaluator.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /*!
     * @brief BusBase model implementation base class.
     *
     */
    template <class ScalarT, typename IdxT>
    class BusBase : public Model::Evaluator<ScalarT, IdxT>
    {
    public:
      using real_type = typename Model::Evaluator<ScalarT, IdxT>::real_type;

      enum BusType
      {
        DEFAULT = 1,
        SLACK
      };

      BusBase()
        : size_(0),
          size_quad_(0),
          size_param_(0)
      {
      }

      BusBase(IdxT bus_id)
        : bus_id_(bus_id)
      {
      }

      BusBase(IdxT size, IdxT size_quad, IdxT size_param)
        : size_(size),
          size_quad_(size_quad),
          size_param_(size_param),
          y_(size_),
          yp_(size_),
          f_(size_),
          g_(size_quad_),
          yB_(size_),
          ypB_(size_),
          fB_(size_),
          gB_(size_param_),
          J_(COO_Matrix<ScalarT, IdxT>()),
          param_(size_param_),
          param_up_(size_param_),
          param_lo_(size_param_)
      {
      }

      virtual ~BusBase()
      {
      }

      /// Pure virtual function, returns bus type (DEFAULT or SLACK).
      virtual int BusType() const = 0;

      virtual IdxT size() override
      {
        return size_;
      }

      virtual IdxT nnz() override
      {
        return nnz_;
      }

      virtual bool hasJacobian() override
      {
        return false;
      }

      virtual IdxT sizeQuadrature() override
      {
        return size_quad_;
      }

      virtual IdxT sizeParams() override
      {
        return size_param_;
      }

      virtual void updateTime(real_type /* t */, real_type /* a */) override
      {
        // No time to update in bus models
      }

      virtual void setTolerances(real_type& rtol, real_type& atol) const override
      {
        rtol = rtol_;
        atol = atol_;
      }

      virtual void setMaxSteps(IdxT& msa) const override
      {
        msa = max_steps_;
      }

      virtual ScalarT&       Vr()       = 0;
      virtual const ScalarT& Vr() const = 0;
      virtual ScalarT&       Vi()       = 0;
      virtual const ScalarT& Vi() const = 0;
      virtual ScalarT&       Ir()       = 0;
      virtual const ScalarT& Ir() const = 0;
      virtual ScalarT&       Ii()       = 0;
      virtual const ScalarT& Ii() const = 0;

      std::vector<ScalarT>& y() override
      {
        return y_;
      }

      const std::vector<ScalarT>& y() const override
      {
        return y_;
      }

      std::vector<ScalarT>& yp() override
      {
        return yp_;
      }

      const std::vector<ScalarT>& yp() const override
      {
        return yp_;
      }

      std::vector<bool>& tag() override
      {
        return tag_;
      }

      const std::vector<bool>& tag() const override
      {
        return tag_;
      }

      std::vector<ScalarT>& yB() override
      {
        return yB_;
      }

      const std::vector<ScalarT>& yB() const override
      {
        return yB_;
      }

      std::vector<ScalarT>& ypB() override
      {
        return ypB_;
      }

      const std::vector<ScalarT>& ypB() const override
      {
        return ypB_;
      }

      std::vector<ScalarT>& param() override
      {
        return param_;
      }

      const std::vector<ScalarT>& param() const override
      {
        return param_;
      }

      std::vector<ScalarT>& param_up() override
      {
        return param_up_;
      }

      const std::vector<ScalarT>& param_up() const override
      {
        return param_up_;
      }

      std::vector<ScalarT>& param_lo() override
      {
        return param_lo_;
      }

      const std::vector<ScalarT>& param_lo() const override
      {
        return param_lo_;
      }

      std::vector<ScalarT>& getResidual() override
      {
        return f_;
      }

      const std::vector<ScalarT>& getResidual() const override
      {
        return f_;
      }

      COO_Matrix<ScalarT, IdxT>& getJacobian() override
      {
        return J_;
      }

      const COO_Matrix<ScalarT, IdxT>& getJacobian() const override
      {
        return J_;
      }

      std::vector<ScalarT>& getIntegrand() override
      {
        return g_;
      }

      const std::vector<ScalarT>& getIntegrand() const override
      {
        return g_;
      }

      std::vector<ScalarT>& getAdjointResidual() override
      {
        return fB_;
      }

      const std::vector<ScalarT>& getAdjointResidual() const override
      {
        return fB_;
      }

      std::vector<ScalarT>& getAdjointIntegrand() override
      {
        return gB_;
      }

      const std::vector<ScalarT>& getAdjointIntegrand() const override
      {
        return gB_;
      }

      virtual const IdxT busID() const
      {
        return bus_id_;
      }

    protected:
      const IdxT bus_id_{static_cast<IdxT>(-1)};

      IdxT size_{0};
      IdxT nnz_{0};
      IdxT size_quad_{0};
      IdxT size_param_{0};

      std::vector<ScalarT> y_;
      std::vector<ScalarT> yp_;
      std::vector<bool>    tag_;
      std::vector<ScalarT> f_;
      std::vector<ScalarT> g_;

      std::vector<ScalarT> yB_;
      std::vector<ScalarT> ypB_;
      std::vector<ScalarT> fB_;
      std::vector<ScalarT> gB_;

      COO_Matrix<ScalarT, IdxT> J_;

      std::vector<ScalarT> param_;
      std::vector<ScalarT> param_up_;
      std::vector<ScalarT> param_lo_;

      real_type time_;
      real_type alpha_;

      real_type rtol_;
      real_type atol_;

      IdxT max_steps_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
