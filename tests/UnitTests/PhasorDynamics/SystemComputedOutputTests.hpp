#pragma once

#include <cstddef>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/SignalIn.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SignalOut.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

#ifdef GRIDKIT_ENABLE_ENZYME
#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>
#endif

namespace GridKit
{
  namespace Testing
  {
    /// Publishes the squared voltage magnitude of a bus as a computed output.
    template <class ScalarT, typename IdxT>
    class SquaredVoltage final : public PhasorDynamics::Component<ScalarT, IdxT>
    {
    public:
      using ComponentT = PhasorDynamics::Component<ScalarT, IdxT>;
      using RealT      = typename ComponentT::RealT;
      using BusT       = PhasorDynamics::BusBase<ScalarT, IdxT>;

      SquaredVoltage(BusT* bus)
        : bus_(bus)
      {
      }

      PhasorDynamics::SignalOut<ScalarT, IdxT>& output()
      {
        return output_;
      }

      int setGridKitComponentID(IdxT component_id) override final
      {
        this->gridkit_component_id_ = component_id;
        return 0;
      }

      int allocate() override final
      {
        if (!this->allocated_)
        {
          this->allocateVectors(this->size_);
        }

        this->wb_.resize(2);
        this->wb_.setToZero();
        this->output_indices_.assign(1, INVALID_INDEX<IdxT>);
        output_.link(&value_, &this->output_indices_[0]);

        this->allocated_ = true;
        return 0;
      }

      Model::ConfigurationChecks verify() const override final
      {
        return {};
      }

      int initialize() override final
      {
        return 0;
      }

      int tagDifferentiable() override final
      {
        return 0;
      }

      int setAbsoluteTolerance(RealT) override final
      {
        return 0;
      }

      int evaluateResidual() override final
      {
        return 0;
      }

      int evaluateOutputs() override final
      {
        auto* wb = this->wb_.getData();
        wb[0]    = bus_->Vr();
        wb[1]    = bus_->Vi();

        return evaluateOutputs(nullptr, nullptr, wb, nullptr, &value_);
      }

      __attribute__((always_inline)) inline int evaluateOutputs([[maybe_unused]] const ScalarT* y,
                                                                [[maybe_unused]] const ScalarT* yp,
                                                                const ScalarT*                  wb,
                                                                [[maybe_unused]] const ScalarT* ws,
                                                                ScalarT*                        o)
      {
        o[0] = wb[0] * wb[0] + wb[1] * wb[1];
        return 0;
      }

      int evaluateJacobian() override final
      {
#ifdef GRIDKIT_ENABLE_ENZYME
        if (this->J_rows_buffer_ == nullptr)
        {
          this->J_rows_buffer_ = new IdxT[2];
          this->J_cols_buffer_ = new IdxT[2];
          this->J_vals_buffer_ = new RealT[2];
        }

        this->nnz_ = 0;

        Enzyme::Sparse::DfDwb<SquaredVoltage, Enzyme::Sparse::MemberFunctions::Outputs>::eval(this,
                                                                                              1,
                                                                                              2,
                                                                                              this->output_indices_.data(),
                                                                                              bus_->getVariableIndices().data(),
                                                                                              nullptr,
                                                                                              nullptr,
                                                                                              this->wb_.getData(),
                                                                                              nullptr,
                                                                                              this->J_rows_buffer_,
                                                                                              this->J_cols_buffer_,
                                                                                              this->J_vals_buffer_,
                                                                                              this->nnz_);
#endif
        return this->constructCoo();
      }

    private:
      BusT*                                    bus_;
      PhasorDynamics::SignalOut<ScalarT, IdxT> output_;
      ScalarT                                  value_{0.0};
    };

    /// Algebraic consumer of a signal, 0 = -x + s^2.
    template <class ScalarT, typename IdxT>
    class SquaredSignal final : public PhasorDynamics::Component<ScalarT, IdxT>
    {
    public:
      using ComponentT = PhasorDynamics::Component<ScalarT, IdxT>;
      using RealT      = typename ComponentT::RealT;

      SquaredSignal()
      {
        this->size_ = 1;
      }

      PhasorDynamics::SignalIn<ScalarT, IdxT>& input()
      {
        return input_;
      }

      int setGridKitComponentID(IdxT component_id) override final
      {
        this->gridkit_component_id_ = component_id;
        return 0;
      }

      int allocate() override final
      {
        if (!this->allocated_)
        {
          this->allocateVectors(this->size_);
        }

        const auto size = static_cast<std::size_t>(this->size_);
        this->tag_.assign(size, false);
        this->variable_indices_.resize(size);
        this->residual_indices_.resize(size);

        this->ws_.resize(1);
        this->ws_.setToZero();
        this->ws_indices_.assign(1, INVALID_INDEX<IdxT>);

        this->allocated_ = true;
        return 0;
      }

      Model::ConfigurationChecks verify() const override final
      {
        Model::ConfigurationChecks checks("SquaredSignal");
        checks.check(input_.linked(), "input is not linked");
        return checks;
      }

      int initialize() override final
      {
        this->y_.getData()[0]  = ScalarT{0.0};
        this->yp_.getData()[0] = ScalarT{0.0};
        return 0;
      }

      int tagDifferentiable() override final
      {
        return 0;
      }

      int setAbsoluteTolerance(RealT) override final
      {
        return 0;
      }

      int evaluateResidual() override final
      {
        auto* ws             = this->ws_.getData();
        ws[0]                = input_.readSignal();
        this->ws_indices_[0] = input_.signalVariableIndex();

        return evaluateInternalResidual(this->y_.getData(), this->yp_.getData(), nullptr, ws, this->f_.getData());
      }

      __attribute__((always_inline)) inline int evaluateInternalResidual(const ScalarT*                  y,
                                                                         [[maybe_unused]] const ScalarT* yp,
                                                                         [[maybe_unused]] const ScalarT* wb,
                                                                         const ScalarT*                  ws,
                                                                         ScalarT*                        f)
      {
        f[0] = -y[0] + ws[0] * ws[0];
        return 0;
      }

      int evaluateJacobian() override final
      {
#ifdef GRIDKIT_ENABLE_ENZYME
        using Fn = Enzyme::Sparse::MemberFunctions;

        if (this->J_rows_buffer_ == nullptr)
        {
          this->J_rows_buffer_ = new IdxT[2];
          this->J_cols_buffer_ = new IdxT[2];
          this->J_vals_buffer_ = new RealT[2];
        }

        this->nnz_ = 0;

        Enzyme::Sparse::DfDy<SquaredSignal, Fn::InternalResidualWithSignal>::eval(this,
                                                                                  1,
                                                                                  1,
                                                                                  this->residual_indices_.data(),
                                                                                  this->variable_indices_.data(),
                                                                                  this->y_.getData(),
                                                                                  this->yp_.getData(),
                                                                                  nullptr,
                                                                                  this->ws_.getData(),
                                                                                  this->J_rows_buffer_,
                                                                                  this->J_cols_buffer_,
                                                                                  this->J_vals_buffer_,
                                                                                  this->nnz_);

        Enzyme::Sparse::DfDws<SquaredSignal, Fn::InternalResidualWithSignal>::eval(this,
                                                                                   1,
                                                                                   1,
                                                                                   this->residual_indices_.data(),
                                                                                   this->ws_indices_.data(),
                                                                                   this->y_.getData(),
                                                                                   this->yp_.getData(),
                                                                                   nullptr,
                                                                                   this->ws_.getData(),
                                                                                   this->J_rows_buffer_,
                                                                                   this->J_cols_buffer_,
                                                                                   this->J_vals_buffer_,
                                                                                   this->nnz_);
#endif
        return this->constructCoo();
      }

    private:
      PhasorDynamics::SignalIn<ScalarT, IdxT> input_;
    };

    /// Computed outputs published through signal ports without solver variables
    template <class ScalarT, typename IdxT>
    class SystemComputedOutputTests
    {
      using RealT = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

    public:
      SystemComputedOutputTests()  = default;
      ~SystemComputedOutputTests() = default;

      /// Outputs are indexed after the variables, and a consumer evaluated
      /// before its producer reads the output of the current state.
      TestOutcome outputsPrecedeResiduals()
      {
        TestStatus success = true;

        PhasorDynamics::SignalNode<ScalarT, IdxT> node;
        PhasorDynamics::Bus<ScalarT, IdxT>        bus(0.6, 0.8);
        SquaredVoltage<ScalarT, IdxT>             producer(&bus);
        SquaredSignal<ScalarT, IdxT>              consumer;
        producer.output().connect(&node);
        consumer.input().connect(&node);

        PhasorDynamics::SystemModel<ScalarT, IdxT> system;
        system.addBus(&bus);
        system.addComponent(&consumer);
        system.addComponent(&producer);

        success *= system.allocate() == 0;
        success *= system.size() == bus.size() + consumer.size();
        success *= consumer.input().signalVariableIndex() == system.size();

        success *= system.initialize() == 0;
        success *= system.evaluateResidual() == 0;
        success *= isEqual(consumer.getResidual().getData()[0], ScalarT{1.0});

        bus.Vr()  = 1.2;
        bus.Vi()  = 1.6;
        success  *= system.evaluateResidual() == 0;
        success  *= isEqual(consumer.getResidual().getData()[0], ScalarT{16.0});

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// The Enzyme Jacobian applies the chain rule through a computed output.
      TestOutcome chainRuleJacobian()
      {
        TestStatus success = true;

        PhasorDynamics::SignalNode<ScalarT, IdxT> node;
        PhasorDynamics::Bus<ScalarT, IdxT>        bus(0.6, 0.8);
        SquaredVoltage<ScalarT, IdxT>             producer(&bus);
        SquaredSignal<ScalarT, IdxT>              consumer;
        producer.output().connect(&node);
        consumer.input().connect(&node);

        PhasorDynamics::SystemModel<ScalarT, IdxT> system;
        system.addBus(&bus);
        system.addComponent(&producer);
        system.addComponent(&consumer);

        success *= system.allocate() == 0;
        success *= system.initialize() == 0;
        success *= system.evaluateResidual() == 0;
        success *= system.evaluateJacobian() == 0;

        const auto enzyme_jacobian = MapFromCsr(system.getCsrJacobian());

        const auto x  = static_cast<std::size_t>(consumer.getVariableIndex(0));
        const auto vr = static_cast<std::size_t>(bus.getVariableIndices()[0]);
        const auto vi = static_cast<std::size_t>(bus.getVariableIndices()[1]);

        // d(s^2)/dV = 2 s (2 V) with s = |V|^2 = 1
        success *= isEqual(entry(enzyme_jacobian, x, x), RealT{-1.0});
        success *= isEqual(entry(enzyme_jacobian, x, vr), RealT{2.4});
        success *= isEqual(entry(enzyme_jacobian, x, vi), RealT{3.2});

        // Repeated evaluations update the chain-rule products in place; s = 4
        bus.Vr()  = 1.2;
        bus.Vi()  = 1.6;
        success  *= system.evaluateResidual() == 0;
        success  *= system.evaluateJacobian() == 0;

        const auto updated_jacobian = MapFromCsr(system.getCsrJacobian());

        success *= isEqual(entry(updated_jacobian, x, x), RealT{-1.0});
        success *= isEqual(entry(updated_jacobian, x, vr), RealT{19.2});
        success *= isEqual(entry(updated_jacobian, x, vi), RealT{25.6});

        return success.report(__func__);
      }

    private:
      using JacobianT = std::vector<DependencyTracking::Variable::DependencyMap>;

      /// Jacobian entry, zero when the entry is not stored.
      static RealT entry(const JacobianT& jacobian, std::size_t row, std::size_t column)
      {
        const auto& dependencies = jacobian[row];
        const auto  found        = dependencies.find(column);
        if (found == dependencies.end())
        {
          return RealT{0.0};
        }
        return static_cast<RealT>(found->second);
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
