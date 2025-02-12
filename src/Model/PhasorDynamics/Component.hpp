#pragma once

#include <vector>
#include <ModelEvaluator.hpp>

namespace GridKit
{
namespace PhasorDynamics
{

    /*!
     * @brief Component model implementation base class.
     *
     */
    template <class ScalarT, typename IdxT>
    class Component : public ModelEvaluator<ScalarT, IdxT>
    {
    public:
        using real_type = typename ModelEvaluator<ScalarT, IdxT>::real_type;

        Component()
          : size_(0),
            size_quad_(0),
            size_param_(0)
        {}

        Component(IdxT size, IdxT size_quad, IdxT size_opt)
          : size_(size),
            size_quad_(size_quad),
            size_param_(size_opt),
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

        virtual IdxT size_quad()
        {
            return size_quad_;
        }

        virtual IdxT size_opt()
        {
            return size_param_;
        }

        // virtual void updateTime(real_type t, real_type a)
        // {
        //     time_ = t;
        //     alpha_ = a;
        //     std::cout << "updateTime: t = " << time_ << "\n";
        // }

        virtual void setTolerances(real_type& rtol, real_type& atol) const
        {
            rtol = rtol_;
            atol = atol_;
        }

        virtual void setMaxSteps(IdxT& msa) const
        {
            msa = max_steps_;
        }

        std::vector<ScalarT>& y()
        {
            return y_;
        }

        const std::vector<ScalarT>& y() const
        {
            return y_;
        }

        std::vector<ScalarT>& yp()
        {
            return yp_;
        }

        const std::vector<ScalarT>& yp() const
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

        std::vector<ScalarT>& yB()
        {
            return yB_;
        }

        const std::vector<ScalarT>& yB() const
        {
            return yB_;
        }

        std::vector<ScalarT>& ypB()
        {
            return ypB_;
        }

        const std::vector<ScalarT>& ypB() const
        {
            return ypB_;
        }

        std::vector<ScalarT>& param()
        {
            return param_;
        }

        const std::vector<ScalarT>& param() const
        {
            return param_;
        }

        std::vector<ScalarT>& param_up()
        {
            return param_up_;
        }

        const std::vector<ScalarT>& param_up() const
        {
            return param_up_;
        }

        std::vector<ScalarT>& param_lo()
        {
            return param_lo_;
        }

        const std::vector<ScalarT>& param_lo() const
        {
            return param_lo_;
        }

        std::vector<ScalarT>& getResidual()
        {
            return f_;
        }

        const std::vector<ScalarT>& getResidual() const
        {
            return f_;
        }

        COO_Matrix<ScalarT, IdxT>& getJacobian()
        {
            return J_;
        }

        const COO_Matrix<ScalarT, IdxT>& getJacobian() const
        {
            return J_;
        }

        std::vector<ScalarT>& getIntegrand()
        {
            return g_;
        }

        const std::vector<ScalarT>& getIntegrand() const
        {
            return g_;
        }

        std::vector<ScalarT>& getAdjointResidual()
        {
            return fB_;
        }

        const std::vector<ScalarT>& getAdjointResidual() const
        {
            return fB_;
        }

        std::vector<ScalarT>& getAdjointIntegrand()
        {
            return gB_;
        }

        const std::vector<ScalarT>& getAdjointIntegrand() const
        {
            return gB_;
        }

        //@todo Fix ID naming
        IdxT getComponentID() const
        {
            return component_id_;
        }



    protected:
        IdxT size_;
        IdxT nnz_;
        IdxT size_quad_;
        IdxT size_param_;

        std::vector<ScalarT> y_;
        std::vector<ScalarT> yp_;
        std::vector<bool> tag_;
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

        IdxT component_id_;
    };

} // namespace PhasorDynamics
} // namespace GridKit
