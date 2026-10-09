#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Solver/Dynamic/Partitioned.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Model
  {
    template <class ScalarT, typename IdxT>
    class PartitionedDAEEvaluator : public Evaluator<ScalarT, IdxT>
    {
      using Base       = Evaluator<ScalarT, IdxT>;
      using RealT      = typename Base::RealT;
      using VectorT    = typename Base::VectorT;
      using CsrMatrixT = typename Base::CsrMatrixT;

    public:
      ~PartitionedDAEEvaluator() override
      {
        delete csr_jac_;
      }

      IdxT size() override
      {
        return 5;
      }

      IdxT nnz() override
      {
        return size() * size();
      }

      bool hasJacobian() override
      {
        return true;
      }

      int allocate() override
      {
        if (allocated_)
        {
          return 0;
        }
        y_.resize(size());
        yp_.resize(size());
        f_.resize(size());
        abs_tol_.resize(size());
        csr_jac_ = new CsrMatrixT(size(), size(), nnz());
        csr_jac_->allocateMatrixData(memory::HOST);
        auto* rows = csr_jac_->getRowData();
        auto* cols = csr_jac_->getColData();
        for (IdxT row = 0; row <= size(); ++row)
        {
          rows[row] = size() * row;
        }
        for (IdxT entry = 0; entry < nnz(); ++entry)
        {
          cols[entry] = entry % size();
        }
        allocated_ = true;
        return 0;
      }

      int initialize() override
      {
        allocate();
        auto* y        = y_.getData();
        auto* yp       = yp_.getData();
        y[0]           = 2.0;
        y[1]           = 1.0;
        y[2]           = 1.0;
        y[3]           = 2.0;
        y[4]           = 0.0;
        const RealT s1 = std::sin(RealT(1.0));
        const RealT s2 = std::sin(RealT(2.0));
        yp[0]          = -2.0 - s1;
        yp[1]          = -1.0 - s2;
        yp[2]          = 2.0 * (2.0 + s1);
        yp[3]          = (1.0 + s2) / 2.0;
        yp[4]          = 1.5 * (5.0 + 2.0 * s1 + s2);
        tag_           = {true, false, true, false, false};
        y_.setDataUpdated();
        yp_.setDataUpdated();
        return 0;
      }

      int setAbsoluteTolerance(RealT rel_tol) override
      {
        abs_tol_.setToConst(rel_tol);
        return 0;
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
        f[0]           = std::sin(y[4] - y[1]) - y[0] - yp[0];
        f[1]           = y[0] * y[0] + y[1] * y[1] + y[4] * y[4] - 5.0;
        f[2]           = std::sin(y[4] - y[3]) - y[2] - yp[2];
        f[3]           = y[2] * y[2] + y[3] * y[3] + y[4] * y[4] - 5.0;
        f[4]           = y[0] - y[1] + y[2] - y[3] + y[4];
        f_.setDataUpdated();
        return 0;
      }

      int evaluateJacobian() override
      {
        const auto* y      = y_.getData();
        auto*       values = csr_jac_->getValues();
        std::fill(values, values + nnz(), 0.0);
        values[0]  = -1.0 - alpha_;
        values[1]  = -std::cos(y[4] - y[1]);
        values[4]  = std::cos(y[4] - y[1]);
        values[5]  = 2.0 * y[0];
        values[6]  = 2.0 * y[1];
        values[9]  = 2.0 * y[4];
        values[12] = -1.0 - alpha_;
        values[13] = -std::cos(y[4] - y[3]);
        values[14] = std::cos(y[4] - y[3]);
        values[17] = 2.0 * y[2];
        values[18] = 2.0 * y[3];
        values[19] = 2.0 * y[4];
        values[20] = 1.0;
        values[21] = -1.0;
        values[22] = 1.0;
        values[23] = -1.0;
        values[24] = 1.0;
        csr_jac_->setUpdated(memory::HOST);
        return 0;
      }

      void updateTime(RealT, RealT alpha) override
      {
        alpha_ = alpha;
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

      IdxT sizeQuadrature() override
      {
        return 0;
      }

      IdxT sizeParams() override
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
        return integrand_;
      }

      const VectorT& getIntegrand() const override
      {
        return integrand_;
      }

      VectorT& getAdjointResidual() override
      {
        return adjoint_residual_;
      }

      const VectorT& getAdjointResidual() const override
      {
        return adjoint_residual_;
      }

      VectorT& getAdjointIntegrand() override
      {
        return adjoint_integrand_;
      }

      const VectorT& getAdjointIntegrand() const override
      {
        return adjoint_integrand_;
      }

      CsrMatrixT* getCsrJacobian() const override
      {
        return csr_jac_;
      }

    private:
      VectorT           y_;
      VectorT           yp_;
      VectorT           f_;
      VectorT           abs_tol_;
      std::vector<bool> tag_;
      VectorT           yB_;
      VectorT           ypB_;
      VectorT           param_;
      VectorT           param_up_;
      VectorT           param_lo_;
      VectorT           integrand_;
      VectorT           adjoint_residual_;
      VectorT           adjoint_integrand_;
      CsrMatrixT*       csr_jac_{nullptr};
      RealT             alpha_{};
      bool              allocated_{false};
    };
  } // namespace Model

  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class PartitionedTests
    {
      using Solver = AnalysisManager::Sundials::Partitioned<ScalarT, IdxT>;
      using Mask   = AnalysisManager::Sundials::Mask;

      static std::vector<Mask> componentMasks()
      {
        return {{true, true, false, false, false},
                {false, false, true, true, false}};
      }

      static Mask couplingMask()
      {
        return {false, false, false, false, true};
      }

    public:
      TestOutcome integration()
      {
        TestStatus                                    success = true;
        Model::PartitionedDAEEvaluator<ScalarT, IdxT> model;
        Solver                                        solver(&model, componentMasks(), couplingMask());
        solver.setFixedStep(0.001);
        solver.setTolerance(1.0e-8);
        solver.configureSimulation();
        solver.initializeSimulation(0.0);

        int callbacks = 0;
        solver.runSimulation(1.0, 0.25, [&](ScalarT)
                             { ++callbacks; });

        // Compare against reference solution generated with Mathematica
        const auto*           y    = model.y().getData();
        static constexpr auto tol  = 5.0e-3;
        success                   *= std::abs(y[0] - 0.8181325137273604) < tol;
        success                   *= std::abs(y[1] - 1.2467607076859357) < tol;
        success                   *= std::abs(y[2] - 0.2350113335864333) < tol;
        success                   *= std::abs(y[3] - 1.4725904879949860) < tol;
        success                   *= std::abs(y[4] - 1.6662073483671279) < tol;
        success                   *= callbacks == 4;
        success                   *= solver.getStats().num_steps_ == 1000;
        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
