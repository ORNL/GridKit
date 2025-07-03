#pragma once

#include <vector>

#include <AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <Model/Evaluator.hpp>
#include <Model/PhasorDynamics/Bus/BusData.hpp>

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
      using BusTypeT  = typename BusData<real_type, IdxT>::BusType;

      BusBase()
        : size_(0)
      {
      }

      BusBase(IdxT bus_id)
        : bus_id_(bus_id)
      {
      }

      virtual ~BusBase()
      {
      }

      /// Pure virtual function, returns bus type (DEFAULT or SLACK).
      virtual BusTypeT BusType() const
      {
        return BusTypeT::DEFAULT;
      }

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

      GridKit::LinearAlgebra::COO_Matrix<ScalarT, IdxT>& getJacobian() override
      {
        return J_;
      }

      const GridKit::LinearAlgebra::COO_Matrix<ScalarT, IdxT>& getJacobian() const override
      {
        return J_;
      }

    protected:
      const IdxT bus_id_{static_cast<IdxT>(-1)};

      IdxT size_{0};
      IdxT nnz_{0};

      std::vector<ScalarT> y_;
      std::vector<ScalarT> yp_;
      std::vector<bool>    tag_;
      std::vector<ScalarT> f_;

      GridKit::LinearAlgebra::COO_Matrix<ScalarT, IdxT> J_;

      real_type time_;
      real_type alpha_;

      real_type rtol_;
      real_type atol_;

      IdxT max_steps_;

      //
      // Adjoint sensitivity members
      //

      std::vector<ScalarT> g_{};
      std::vector<ScalarT> yB_{};
      std::vector<ScalarT> ypB_{};
      std::vector<ScalarT> fB_{};
      std::vector<ScalarT> gB_{};

      std::vector<ScalarT> param_{};
      std::vector<ScalarT> param_up_{};
      std::vector<ScalarT> param_lo_{};

      //
      // Public adjoint sensitivity methods (not yet implemented in components)
      //

    public:
      virtual IdxT sizeQuadrature() override
      {
        throw "ERROR: Method not implemented!\n";
        return 0;
      }

      virtual IdxT sizeParams() override
      {
        throw "ERROR: Method not implemented!\n";
        return 0;
      }

      std::vector<ScalarT>& yB() override
      {
        throw "ERROR: Method not implemented!\n";
        return yB_;
      }

      const std::vector<ScalarT>& yB() const override
      {
        throw "ERROR: Method not implemented!\n";
        return yB_;
      }

      std::vector<ScalarT>& ypB() override
      {
        throw "ERROR: Method not implemented!\n";
        return ypB_;
      }

      const std::vector<ScalarT>& ypB() const override
      {
        throw "ERROR: Method not implemented!\n";
        return ypB_;
      }

      std::vector<ScalarT>& param() override
      {
        throw "ERROR: Method not implemented!\n";
        return param_;
      }

      const std::vector<ScalarT>& param() const override
      {
        throw "ERROR: Method not implemented!\n";
        return param_;
      }

      std::vector<ScalarT>& param_up() override
      {
        throw "ERROR: Method not implemented!\n";
        return param_up_;
      }

      const std::vector<ScalarT>& param_up() const override
      {
        throw "ERROR: Method not implemented!\n";
        return param_up_;
      }

      std::vector<ScalarT>& param_lo() override
      {
        throw "ERROR: Method not implemented!\n";
        return param_lo_;
      }

      const std::vector<ScalarT>& param_lo() const override
      {
        throw "ERROR: Method not implemented!\n";
        return param_lo_;
      }

      int evaluateIntegrand() override
      {
        throw "ERROR: Method not implemented!\n";
        return 1;
      }

      int initializeAdjoint() override
      {
        throw "ERROR: Method not implemented!\n";
        return 1;
      }

      int evaluateAdjointResidual() override
      {
        throw "ERROR: Method not implemented!\n";
        return 1;
      }

      int evaluateAdjointIntegrand() override
      {
        throw "ERROR: Method not implemented!\n";
        return 1;
      }

      std::vector<ScalarT>& getResidual() override
      {
        return f_;
      }

      const std::vector<ScalarT>& getResidual() const override
      {
        return f_;
      }

      std::vector<ScalarT>& getIntegrand() override
      {
        throw "ERROR: Method not implemented!\n";
        return g_;
      }

      const std::vector<ScalarT>& getIntegrand() const override
      {
        throw "ERROR: Method not implemented!\n";
        return g_;
      }

      std::vector<ScalarT>& getAdjointResidual() override
      {
        throw "ERROR: Method not implemented!\n";
        return fB_;
      }

      const std::vector<ScalarT>& getAdjointResidual() const override
      {
        throw "ERROR: Method not implemented!\n";
        return fB_;
      }

      std::vector<ScalarT>& getAdjointIntegrand() override
      {
        throw "ERROR: Method not implemented!\n";
        return gB_;
      }

      const std::vector<ScalarT>& getAdjointIntegrand() const override
      {
        throw "ERROR: Method not implemented!\n";
        return gB_;
      }

      virtual const IdxT busID() const
      {
        return bus_id_;
      }
    };

  } // namespace PhasorDynamics
} // namespace GridKit
