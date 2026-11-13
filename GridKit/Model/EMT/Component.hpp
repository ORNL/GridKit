/**
 * @file Component.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Base class of EMT components.
 */

#pragma once

#include <cassert>
#include <cstddef>
#include <memory>
#include <span>
#include <string>
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
    class SystemModel;

    /**
     * @brief EMT component, which views its blocks of the system vectors
     *
     * The local variables are the internal variables, then the external
     * variables of each port. Each has a value, a derivative and a residual.
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

    public:
      /// External variables read from and added to the local variables of `component` at `indices`
      struct Port
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

      /// Connect `port` to the local variables of `component` from `offset`
      template <typename PortT>
      int connect(PortT port, Component& component, IdxT offset = 0);

      /// Connect `port` to the local variables of `component` at `indices`
      template <typename PortT>
      int connect(PortT port, Component& component, std::vector<IdxT> indices);

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

      int evaluateResidual() override
      {
        return 0;
      }

      /// dF/dy + alpha dF/dy' from the partials of the last tagging
      int evaluateJacobian() override;

      void updateTime(RealT t, RealT alpha) override
      {
        time_  = t;
        alpha_ = alpha;
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

      /// Size of each internal variable, then of each port
      void setLayout(const std::vector<IdxT>& variables, const std::vector<IdxT>& ports);

      /// Local index of the first entry of internal variable `variable`
      template <Utilities::SizedEnum VariableT>
      IdxT local(VariableT variable) const
      {
        return offsets_[static_cast<IdxT>(variable)];
      }

      /// Local index of external variable `k` of `port`
      template <Utilities::SizedEnum PortT>
      IdxT local(PortT port, IdxT k) const;

      /// Index of the first external variable of `port` among the external variables
      template <Utilities::SizedEnum PortT>
      IdxT external(PortT port) const;

      /// Submodel `symbol`, owned by this component and evaluated by the system
      template <typename ModelT, Utilities::SizedEnum SymbolT, typename... Args>
      ModelT& addSubmodel(SymbolT symbol, Args&&... args);

      template <typename ModelT, Utilities::SizedEnum SymbolT>
      ModelT& submodel(SymbolT symbol);

      template <typename ModelT>
      int evaluateResidual(const ModelT& model);

      /// Partials of a model that is linear between events
      template <typename ModelT>
      int tagDifferentiable(const ModelT& model);

      /// Absolute tolerance of external variable `k` of `port`
      template <typename PortT>
      ScalarT absoluteTolerance(PortT port, IdxT k) const;

      VectorT           y_;
      VectorT           yp_;
      VectorT           f_;
      VectorT           abs_tol_;
      std::vector<bool> tag_;
      RealT             time_{ZERO<RealT>};
      RealT             alpha_{ONE<RealT>};
      RealT             va_system_base_{ZERO<RealT>}; ///< System power base [VA]

    private:
      /// Entries of internal variable `variable` in `vector`, a view of the internal block
      template <Utilities::SizedEnum VariableT>
      std::span<ScalarT> internal(VectorT& vector, VariableT variable)
      {
        return {vector.getData() + local(variable), size(variable)};
      }

      void gather(ScalarT* y, ScalarT* yp) const;
      void scatter(const ScalarT* f);
      void setIndices();

      std::string                             id_;
      IdxT                                    size_{0};
      std::vector<IdxT>                       offsets_; ///< First local index of each internal variable, then the internal size
      std::vector<Port>                       ports_;
      std::vector<std::unique_ptr<Component>> submodels_; ///< In the order of the parent's submodel enumeration

      IndexVectorT variable_indices_; ///< Global variable and residual of each local variable
      VectorT      local_y_;
      VectorT      local_yp_;
      VectorT      local_f_;
      CooBuffer    jacobian_entries_;
      RealVectorT  dfdy_;
      RealVectorT  dfdyp_;

      ScalarT* system_y_{nullptr};
      ScalarT* system_yp_{nullptr};
      ScalarT* system_f_{nullptr};
      ScalarT* system_abs_tol_{nullptr};
    };
  } // namespace EMT
} // namespace GridKit
