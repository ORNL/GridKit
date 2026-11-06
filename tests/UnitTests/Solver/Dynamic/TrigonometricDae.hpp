#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

#include <GridKit/LinearAlgebra/SparseMatrix/CsrMatrix.hpp>
#include <GridKit/Model/Evaluator.hpp>

namespace GridKit
{
  namespace Model
  {
    /**
     * @brief A test DAE (referred to as the "trigonometric DAE") that is useful for evaluating order
     * of DAE integrators. The system is 2-dimensional, with one differential variable and one algebraic variable.
     * The DAE is index 1 and the derivatives of the model are non-vanishing.
     *
     * The valid simulation time interval is \f([0.5,2]\f) and the model is initialized at \f(t = 0.5\f).
     *
     */
    template <class ScalarT, typename IdxT>
    class TrigonometricDaeEvaluator : public Model::Evaluator<ScalarT, IdxT>
    {
    public:
      using RealT   = typename Model::Evaluator<ScalarT, IdxT>::RealT;
      using VectorT = typename Model::Evaluator<ScalarT, IdxT>::VectorT;

      constexpr static size_t SIZE = 2;

      /**
       * @brief Evaluate the analytic solution at a specified time.
       *
       * @param[in] time Simulation time.
       * @return The differential and algebraic solution components, respectively.
       */
      static std::array<ScalarT, SIZE> analyticSolution(RealT time)
      {
        return {std::sinh(time), std::tanh(time)};
      }

      /**
       * @brief Evaluate the derivative of the analytic solution at a specified time.
       *
       * @param[in] time Simulation time.
       * @return The derivatives of the differential and algebraic solution components.
       */
      static std::array<ScalarT, SIZE> analyticDerivative(RealT time)
      {
        return {std::cosh(time), 1.0 / std::pow(std::cosh(time), 2)};
      }

      TrigonometricDaeEvaluator()
      {
      }

      int allocate() override
      {
        constexpr size_t NNZ = SIZE * SIZE;

        y_.resize(SIZE);
        yp_.resize(SIZE);
        f_.resize(SIZE);

        IdxT*  row_ptrs = new IdxT[SIZE + 1];
        IdxT*  cols     = new IdxT[NNZ];
        RealT* vals     = new RealT[NNZ];

        for (size_t i = 0; i < SIZE + 1; i++)
        {
          row_ptrs[i] = static_cast<IdxT>(i * SIZE);
        }

        for (size_t i = 0; i < SIZE; i++)
        {
          for (size_t j = 0; j < SIZE; j++)
          {
            cols[i * SIZE + j] = static_cast<IdxT>(j);
          }
        }

        csr_jac_ = std::make_unique<GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>>(SIZE, SIZE, NNZ, &row_ptrs, &cols, &vals);

        return 0;
      }

      int initialize() override
      {
        auto* y  = y_.getData();
        auto* yp = yp_.getData();
        auto* f  = f_.getData();

        const auto solution   = analyticSolution(0.5);
        const auto derivative = analyticDerivative(0.5);
        y[0]                  = solution[0];
        y[1]                  = solution[1];
        yp[0]                 = derivative[0];
        yp[1]                 = derivative[1];

        tag_ = {true, false};

        f[0] = 0.0;
        f[1] = 0.0;

        y_.setDataUpdated();
        yp_.setDataUpdated();
        f_.setDataUpdated();

        return 0;
      }

      IdxT size() override
      {
        return SIZE;
      }

      IdxT nnz() override
      {
        return SIZE * SIZE;
      }

      bool hasJacobian() override
      {
        return true;
      }

      IdxT sizeQuadrature() override
      {
        return 0;
      }

      IdxT sizeParams() override
      {
        return 0;
      }

      int setAbsoluteTolerance([[maybe_unused]] RealT rel_tol) override
      {
        return 0;
      }

      VectorT& absoluteTolerance() override
      {
        return abs_tol_;
      }

      const VectorT& absoluteTolerance() const override
      {
        return abs_tol_;
      }

      int tagDifferentiable() override
      {
        return 0;
      }

      int evaluateResidual() override
      {
        const auto* y  = y_.getData();
        const auto* yp = yp_.getData();
        auto*       f  = f_.getData();

        ScalarT y02 = y[0] * y[0];
        ScalarT y12 = y[1] * y[1];

        f[0] = -yp[0] + y02 / (y[1] * std::sqrt(std::pow(y[0] / y[1], 2) - 1));
        f[1] = y12 + 1 / (1 + y02) - (y02 / y12 - y02);

        f_.setDataUpdated();

        return 0;
      }

      int evaluateJacobian() override
      {
        RealT* vals = csr_jac_->getValues();

        const auto* y = y_.getData();

        ScalarT y1 = y[0];
        ScalarT y2 = y[1];

        ScalarT y12 = y1 * y1;
        ScalarT y13 = y12 * y1;
        ScalarT y22 = y2 * y2;
        ScalarT y23 = y22 * y2;
        ScalarT y24 = y22 * y22;

        ScalarT tmp  = std::pow(y12 / y22 - 1, 1.5);
        ScalarT tmp2 = std::pow(y12 + 1, 2);

        vals[0] = static_cast<RealT>(-alpha_ + -(-y13 + 2 * y1 * y22) / (y23 * tmp));
        vals[1] = static_cast<RealT>(y12 / (y22 * tmp));
        vals[2] = static_cast<RealT>(-2 * y1 * (1 / y22 - 1) - (2 * y1) / tmp2);
        vals[3] = static_cast<RealT>((2 * (y12 + y24)) / y23);

        return 0;
      }

      int evaluateIntegrand() override
      {
        return 0;
      }

      int initializeAdjoint() override
      {
        return 0;
      }

      int evaluateAdjointResidual() override
      {
        return 0;
      }

      int evaluateAdjointIntegrand() override
      {
        return 0;
      }

      void updateTime([[maybe_unused]] RealT t, RealT a) override
      {
        alpha_ = a;
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

      std::vector<bool>& tag() override
      {
        return tag_;
      }

      const std::vector<bool>& tag() const override
      {
        return tag_;
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

      IdxT getIDcomponent()
      {
        return 0;
      }

      GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>* getCsrJacobian() const override
      {
        return csr_jac_.get();
      }

    protected:
      VectorT           y_;
      VectorT           yp_;
      std::vector<bool> tag_;
      VectorT           f_;
      VectorT           g_;

      VectorT yB_;
      VectorT ypB_;
      VectorT fB_;
      VectorT gB_;

      VectorT abs_tol_;

      std::unique_ptr<GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>> csr_jac_;

      RealT alpha_;

      VectorT param_;
      VectorT param_up_;
      VectorT param_lo_;
    };
  } // namespace Model
} // namespace GridKit
