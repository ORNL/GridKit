/**
 * @file Component.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Base class of EMT components.
 */

#pragma once

#include <cassert>
#include <concepts>
#include <cstddef>
#include <memory>
#include <numeric>
#include <span>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <GridKit/Constants.hpp>
#include <GridKit/LinearAlgebra/Vector/Vector.hpp>
#include <GridKit/Model/ConfigurationChecks.hpp>
#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Utilities/Enum.hpp>
#include <GridKit/Utilities/Errors.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace EMT
  {
    using Log = ::GridKit::Utilities::Logger;

    template <typename scalar_type, typename index_type>
    class Component;

    template <typename scalar_type, typename index_type>
    class SystemModel;

    /// Constant sizes and a final type, so Enzyme sees the layout and calls the kernels directly.
    /// A model satisfies it, and so does a set of kernels that a model selects among.
    template <typename ModelT, typename ScalarT, typename IdxT>
    concept ComponentKernels =
        std::is_final_v<ModelT>
        && ModelT::NUM_INTERNAL_VARIABLES.size() == Utilities::enum_size<typename ModelT::InternalVariablesT>()
        && ModelT::NUM_TERMINAL_VARIABLES.size() == Utilities::enum_size<typename ModelT::TerminalsT>()
        && ModelT::NUM_INPUT_VARIABLES.size() == Utilities::enum_size<typename ModelT::InputsT>()
        && ModelT::NUM_OUTPUT_VARIABLES.size() == Utilities::enum_size<typename ModelT::OutputsT>();

    /**
     * @brief EMT component, which views its blocks of the system vectors
     *
     * A terminal reads the variables it is connected to and adds to their
     * residuals. An input reads the variables it is connected to. An output
     * is a value that the model computes, added to the residuals of the
     * variables it is connected to. Each output is evaluated with the
     * residual only when it is connected, and otherwise only for a monitor.
     */
    template <typename scalar_type, typename index_type>
    class Component : public Model::Evaluator<scalar_type, index_type>
    {
      friend class SystemModel<scalar_type, index_type>;

      using EvaluatorT = Model::Evaluator<scalar_type, index_type>;

    protected:
      using ScalarT      = scalar_type;
      using IdxT         = index_type;
      using RealT        = typename EvaluatorT::RealT;
      using VectorT      = typename EvaluatorT::VectorT;
      using CsrMatrixT   = typename EvaluatorT::CsrMatrixT;
      using RealVectorT  = LinearAlgebra::Vector<RealT, IdxT>;
      using IndexVectorT = LinearAlgebra::Vector<IdxT, IdxT>;

    private:
      struct Time
      {
        RealT t{ZERO<RealT>};
        RealT alpha{ONE<RealT>};
      };

      Time  local_time_;
      Time* clock_{&local_time_}; ///< Bound components share the system time and Jacobian coefficient.

    public:
      /// Local variables of `component`, at `indices`, that a terminal, an input or an output is connected to
      struct Connection
      {
        Component*        component{nullptr};
        std::vector<IdxT> indices;
      };

      /// Derivative entries in global indices
      struct CooBuffer
      {
        int resize(IdxT capacity)
        {
          if (rows.resize(capacity) != 0 || cols.resize(capacity) != 0 || values.resize(capacity) != 0)
          {
            return 1;
          }
          return 0;
        }

        IndexVectorT rows;
        IndexVectorT cols;
        RealVectorT  values;
        IdxT         nnz{0};
      };

      virtual ~Component() = default;

      virtual Model::ConfigurationChecks verify() const = 0;

      template <Utilities::SizedEnum InputT>
      int connectInput(InputT input, Component& component, std::vector<IdxT> indices);

      template <Utilities::SizedEnum OutputT>
      int connectOutput(OutputT output, Component& component, std::vector<IdxT> indices);

      /// Size of internal variable `variable`
      template <Utilities::SizedEnum VariableT>
      IdxT size(VariableT variable) const
      {
        const auto k = static_cast<IdxT>(variable);
        assert(k + 1 < offsets_.size());
        return offsets_[k + 1] - offsets_[k];
      }

      /// Values of internal variable `variable`, once the system is allocated
      template <Utilities::SizedEnum VariableT>
      std::span<ScalarT> y(VariableT variable)
      {
        return internal(y_, variable);
      }

      /// Residuals of internal variable `variable`, once the system is allocated
      template <Utilities::SizedEnum VariableT>
      std::span<ScalarT> getResidual(VariableT variable)
      {
        return internal(f_, variable);
      }

      /// Absolute tolerances of internal variable `variable`, once the system is allocated
      template <Utilities::SizedEnum VariableT>
      std::span<ScalarT> absoluteTolerance(VariableT variable)
      {
        return internal(abs_tol_, variable);
      }

      /// Number of entries of output `output`
      template <Utilities::SizedEnum OutputT>
      IdxT outputSize(OutputT output) const
      {
        return outputs_[static_cast<IdxT>(output)].indices.size();
      }

      /// True if `output` is added to a residual
      template <Utilities::SizedEnum OutputT>
      bool outputConnected(OutputT output) const
      {
        return residual_indices_[size_ + terminals_size_ + outputOffset(output)] != INVALID_INDEX<IdxT>;
      }

      /// Entries of output `output` that `evaluateOutputs` evaluated last, zero before the first evaluation
      template <Utilities::SizedEnum OutputT>
      std::span<const ScalarT> getOutput(OutputT output) const
      {
        return {output_values_.data() + outputOffset(output), outputSize(output)};
      }

      /// Evaluate every output at the current values, whether or not it is connected
      virtual int evaluateOutputs()
      {
        return 0;
      }

      /// Case id, or `id.symbol` for a submodel
      const std::string& id() const
      {
        return id_;
      }

      /// The system binds each component to its blocks
      int allocate() override
      {
        return 0;
      }

      int initialize() override
      {
        return 0;
      }

      int tagDifferentiable() override
      {
        return 0;
      }

      /// The relative tolerance, unless a model scales its own
      int setAbsoluteTolerance(RealT rel_tol) override
      {
        return abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      }

      int evaluateResidual() override;

      int evaluateJacobian() override;

      void updateTime(RealT t, RealT alpha) override
      {
        *clock_ = {t, alpha};
      }

      IdxT size() override
      {
        return size_;
      }

      IdxT nnz() override
      {
        return jacobian_entries_.nnz;
      }

      bool hasJacobian() override
      {
        return true;
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

      VectorT& getResidual() override
      {
        return f_;
      }

      const VectorT& getResidual() const override
      {
        return f_;
      }

      VectorT& absoluteTolerance() override
      {
        return abs_tol_;
      }

      const VectorT& absoluteTolerance() const override
      {
        return abs_tol_;
      }

      std::vector<bool>& tag() override
      {
        return tag_;
      }

      const std::vector<bool>& tag() const override
      {
        return tag_;
      }

      [[noreturn]] IdxT sizeQuadrature() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] IdxT sizeParams() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& yB() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& yB() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& ypB() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& ypB() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& param() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& param() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& param_up() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& param_up() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& param_lo() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& param_lo() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] int evaluateIntegrand() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] int initializeAdjoint() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] int evaluateAdjointResidual() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] int evaluateAdjointIntegrand() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& getIntegrand() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& getIntegrand() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& getAdjointResidual() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& getAdjointResidual() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] VectorT& getAdjointIntegrand() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const VectorT& getAdjointIntegrand() const override
      {
        throw NotImplementedError(__func__);
      }

    protected:
      using NotImplementedError = GridKit::Utilities::NotImplementedError;

      Component();

      template <ComponentKernels<scalar_type, index_type> ModelT>
      explicit Component(std::type_identity<ModelT>);

      void setLayout(const std::vector<IdxT>& variables,
                     const std::vector<IdxT>& terminals,
                     const std::vector<IdxT>& inputs,
                     const std::vector<IdxT>& outputs);

      /// Local index of the first entry of internal variable `variable`
      template <Utilities::SizedEnum VariableT>
      IdxT local(VariableT variable) const
      {
        return offsets_[static_cast<IdxT>(variable)];
      }

      template <Utilities::SizedEnum TerminalT>
      IdxT local(TerminalT terminal, IdxT k) const;

      /// Index of the first external variable of `terminal`, among the external values and the residuals at the terminals
      template <Utilities::SizedEnum TerminalT>
      IdxT terminalOffset(TerminalT terminal) const
      {
        return terminal_offsets_[static_cast<IdxT>(terminal)];
      }

      /// Index of the first external value of `input`, after those of the terminals
      template <Utilities::SizedEnum InputT>
      IdxT inputOffset(InputT input) const
      {
        return input_offsets_[static_cast<IdxT>(input)];
      }

      /// Index of the first entry of `output`, among the outputs
      template <Utilities::SizedEnum OutputT>
      IdxT outputOffset(OutputT output) const
      {
        return output_offsets_[static_cast<IdxT>(output)];
      }

      template <typename ModelT, Utilities::SizedEnum SymbolT, typename... Args>
      ModelT& addSubmodel(SymbolT symbol, Args&&... args);

      template <typename ModelT, Utilities::SizedEnum SymbolT>
      ModelT& submodel(SymbolT symbol);

      template <typename ModelT, Utilities::SizedEnum SymbolT>
      const ModelT& submodel(SymbolT symbol) const;

      template <typename ModelT, typename T>
      void evaluateKernels(const ModelT& model, const T* y, const T* yp, const T* w, const T* wp, T* f) const;

      template <typename ModelT>
      int evaluateResidual(const ModelT& model);

      template <typename ModelT>
      int evaluateOutputs(const ModelT& model);

      std::pair<const ScalarT*, const ScalarT*> gather();

      /// Partials of a model whose parameters size its layout
      template <typename ModelT>
      int tagDifferentiable(const ModelT& model);

      /// Partials at the current values, then dF/dy + alpha dF/dy', of such a model that is not linear
      template <typename ModelT>
      int evaluateJacobian(const ModelT& model);

      ScalarT externalTolerance(IdxT k) const;

      /// Absolute tolerances of output `output`, which a signal that follows it takes
      template <Utilities::SizedEnum OutputT>
      std::span<ScalarT> outputTolerance(OutputT output)
      {
        return {output_tol_.data() + outputOffset(output), outputSize(output)};
      }

      /// Enzyme dF/dy' of a fixed-layout model
      template <ComponentKernels<scalar_type, index_type> KernelsT>
      int tagDifferentiable(const KernelsT& kernels);

      /// Enzyme dF/dy + alpha dF/dy' of a fixed-layout model
      template <ComponentKernels<scalar_type, index_type> KernelsT>
      int evaluateJacobian(const KernelsT& kernels);

      /// Simulation time [s], shared with the monitor and time-dependent models.
      const RealT& time() const
      {
        return clock_->t;
      }

      VectorT           y_;
      VectorT           yp_;
      VectorT           f_;
      VectorT           abs_tol_;
      std::vector<bool> tag_;
      RealT             freq_system_base_{ZERO<RealT>}; ///< System frequency base [Hz]
      RealT             va_system_base_{ZERO<RealT>};   ///< System power base [VA]

    private:
      /// Entries of internal variable `variable` in `vector`, a view of the internal block
      template <Utilities::SizedEnum VariableT>
      std::span<ScalarT> internal(VectorT& vector, VariableT variable)
      {
        return {vector.getData() + local(variable), size(variable)};
      }

      /// Local values of a fixed layout: internal variables, terminals and inputs
      template <typename ModelT>
      static constexpr IdxT valueSize()
      {
        return std::accumulate(ModelT::NUM_INTERNAL_VARIABLES.begin(), ModelT::NUM_INTERNAL_VARIABLES.end(), IdxT{0})
               + std::accumulate(ModelT::NUM_TERMINAL_VARIABLES.begin(), ModelT::NUM_TERMINAL_VARIABLES.end(), IdxT{0})
               + std::accumulate(ModelT::NUM_INPUT_VARIABLES.begin(), ModelT::NUM_INPUT_VARIABLES.end(), IdxT{0});
      }

      /// Local residuals of a fixed layout: internal variables, terminals and outputs
      template <typename ModelT>
      static constexpr IdxT residualSize()
      {
        return std::accumulate(ModelT::NUM_INTERNAL_VARIABLES.begin(), ModelT::NUM_INTERNAL_VARIABLES.end(), IdxT{0})
               + std::accumulate(ModelT::NUM_TERMINAL_VARIABLES.begin(), ModelT::NUM_TERMINAL_VARIABLES.end(), IdxT{0})
               + std::accumulate(ModelT::NUM_OUTPUT_VARIABLES.begin(), ModelT::NUM_OUTPUT_VARIABLES.end(), IdxT{0});
      }

      /// Enzyme entries of beta dF/dy + alpha dF/dy'
      template <ComponentKernels<scalar_type, index_type> ModelT>
      int evaluateDfDy(const ModelT& model, RealT beta, RealT alpha);

      int  connect(Connection& connection, Component& component, std::vector<IdxT> indices, IdxT limit);
      void scatter(const ScalarT* f);
      void setIndices();

      std::string                             id_;
      IdxT                                    size_{0};
      std::vector<IdxT>                       offsets_;          ///< First local index of each internal variable, then the internal size
      std::vector<IdxT>                       terminal_offsets_; ///< First external value of each terminal
      std::vector<IdxT>                       input_offsets_;    ///< First external value of each input, after those of the terminals
      std::vector<IdxT>                       output_offsets_;   ///< First entry of each output
      std::vector<Connection>                 terminals_;
      std::vector<Connection>                 inputs_;
      std::vector<Connection>                 outputs_;
      IdxT                                    terminals_size_{0}; ///< Number of external variables of the terminals
      std::vector<std::unique_ptr<Component>> submodels_;         ///< In the order of the parent's submodel enumeration

      /// External residual and the system residual it is added to
      struct Added
      {
        IdxT external;
        IdxT residual;
      };

      IndexVectorT         variable_indices_; ///< System variable of each local value: internal, terminals, inputs
      std::vector<IdxT>    residual_indices_; ///< System residual of each local residual: internal, terminals, outputs
      std::vector<IdxT>    external_indices_; ///< System variable of each external value
      std::vector<Added>   added_;            ///< Residuals of the terminals and of the outputs that are connected
      std::vector<ScalarT> output_values_;    ///< Entries of each output, as evaluated last
      std::vector<ScalarT> output_tol_;       ///< Absolute tolerance of each output entry, the relative tolerance unless the model sets it
      VectorT              local_y_;
      VectorT              local_yp_;
      VectorT              local_f_;
      CooBuffer            jacobian_entries_;
      RealVectorT          dfdy_;
      RealVectorT          dfdyp_;

      ScalarT* system_y_{nullptr};
      ScalarT* system_yp_{nullptr};
      ScalarT* system_f_{nullptr};
      ScalarT* system_abs_tol_{nullptr};
    };
  } // namespace EMT
} // namespace GridKit
