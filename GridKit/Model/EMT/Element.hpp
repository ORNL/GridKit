/**
 * @file Element.hpp
 *
 * @brief Shared EMT model contract. An Element owns a contiguous block of
 * variables and residual rows, binds to the aggregate global arrays once, reads
 * borrowed inputs through Signals, and supplies its own sparse-jet Jacobian. Only
 * aggregates implement Model::Evaluator.
 *
 */

#pragma once

#include <algorithm>
#include <numeric>
#include <span>
#include <tuple>
#include <type_traits>
#include <vector>

#include <GridKit/Constants.hpp>
#include <GridKit/LinearAlgebra/SparseJet.hpp>
#include <GridKit/Model/EMT/Signal.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    class Element
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename GridKit::ScalarTraits<ScalarT>::RealT;

      /**
       * @brief Accumulating COO triplet store for element Jacobians.
       *
       * Entries are appended by setValues() and combined by deduplicate(),
       * which sorts by (row, col) and sums duplicates.
       */
      struct JacobianTriplets
      {
        std::vector<IdxT>  rows;
        std::vector<IdxT>  cols;
        std::vector<RealT> vals;

        void zeroMatrix()
        {
          rows.clear();
          cols.clear();
          vals.clear();
        }

        void setValues(RealT alpha, const IdxT* r, const IdxT* c, const RealT* v, IdxT n)
        {
          for (IdxT i = 0; i < n; ++i)
          {
            rows.push_back(r[i]);
            cols.push_back(c[i]);
            vals.push_back(alpha * v[i]);
          }
        }

        void deduplicate()
        {
          if (vals.empty())
          {
            return;
          }
          std::vector<size_t> order(vals.size());
          std::iota(order.begin(), order.end(), size_t{0});
          std::sort(order.begin(), order.end(), [this](size_t a, size_t b)
                    { return rows[a] != rows[b] ? rows[a] < rows[b] : cols[a] < cols[b]; });
          std::vector<IdxT>  sorted_rows;
          std::vector<IdxT>  sorted_cols;
          std::vector<RealT> sorted_vals;
          sorted_rows.reserve(order.size());
          sorted_cols.reserve(order.size());
          sorted_vals.reserve(order.size());
          for (const size_t i : order)
          {
            if (!sorted_vals.empty() && sorted_rows.back() == rows[i] && sorted_cols.back() == cols[i])
            {
              sorted_vals.back() += vals[i];
            }
            else
            {
              sorted_rows.push_back(rows[i]);
              sorted_cols.push_back(cols[i]);
              sorted_vals.push_back(vals[i]);
            }
          }
          rows = std::move(sorted_rows);
          cols = std::move(sorted_cols);
          vals = std::move(sorted_vals);
        }

        std::tuple<const std::vector<IdxT>&, const std::vector<IdxT>&, const std::vector<RealT>&>
        getEntries(bool) const
        {
          return {rows, cols, vals};
        }

        std::tuple<std::vector<IdxT>, std::vector<IdxT>, std::vector<RealT>> getEntryCopies() const
        {
          return {rows, cols, vals};
        }

        IdxT nnz() const
        {
          return static_cast<IdxT>(vals.size());
        }
      };

      using MatrixT        = JacobianTriplets;
      using SignalT        = Signal<IdxT>;
      using SignalStorageT = typename SignalT::Storage;
      using SliceT         = Slice<ScalarT, IdxT>;
      using JetT           = GridKit::LinearAlgebra::SparseJet<RealT, IdxT>;

      virtual ~Element() = default;

      virtual IdxT                     size() const                                    = 0;
      virtual std::span<const SignalT> inputs() const                                  = 0;
      virtual int                      initialize()                                    = 0;
      virtual void                     tagDifferentiable(std::vector<bool>& tag) const = 0;
      virtual int                      evaluateResidual()                              = 0;
      virtual int                      evaluateJacobian()                              = 0;

      IdxT end() const
      {
        return base_ + size();
      }

      void bind(ScalarT* y, ScalarT* yp, ScalarT* f, IdxT* index, IdxT base)
      {
        base_   = base;
        gy_     = y;
        gyp_    = yp;
        gindex_ = index;
        y_      = y + base;
        yp_     = yp + base;
        f_      = f + base;
        index_  = index + base;
      }

      void setCoordinate(RealT coordinate, RealT alpha)
      {
        coordinate_ = coordinate;
        alpha_      = alpha;
      }

      MatrixT& jacobian()
      {
        return J_;
      }

    protected:
      RealT coordinate() const
      {
        return coordinate_;
      }

      ScalarT* inputValues(SignalT signal) const
      {
        return inputBase(signal) + signal.offset;
      }

      IdxT* inputIndices(SignalT signal) const
      {
        return gindex_ + signal.offset;
      }

      SliceT state(IdxT offset, IdxT rows, IdxT cols = 1) const
      {
        return {y_, index_, rows, cols, offset};
      }

      SliceT derivative(IdxT offset, IdxT rows, IdxT cols = 1) const
      {
        return {yp_, index_, rows, cols, offset};
      }

      SliceT residual(IdxT offset, IdxT rows, IdxT cols = 1) const
      {
        return {f_, index_, rows, cols, offset};
      }

      SliceT input(SignalT signal) const
      {
        return {inputBase(signal), gindex_, signal.rows, signal.cols, signal.offset};
      }

      void setActiveInput(SignalT signal) const
      {
        active_input_     = signal;
        has_active_input_ = true;
      }

      void clearActiveInput() const
      {
        has_active_input_ = false;
      }

      template <typename ValueT>
      Slice<ValueT, IdxT> input(ValueT* u, SignalT signal) const
      {
        if constexpr (std::is_same_v<ValueT, ScalarT>)
        {
          if (u != nullptr && isActiveInput(signal))
          {
            return {u, inputIndices(signal), signal.rows, signal.cols};
          }
          return {inputBase(signal), gindex_, signal.rows, signal.cols, signal.offset};
        }
        else
        {
          (void) u;
          return {jetInputBase(signal), gindex_, signal.rows, signal.cols, signal.offset};
        }
      }

      RealT inputJacobianScale(SignalT signal) const
      {
        return signal.storage == SignalStorageT::Derivative ? alpha_ : RealT{1.0};
      }

      template <typename ValueT>
      Slice<ValueT, IdxT> slice(ValueT* values, IdxT offset, IdxT rows, IdxT cols = 1) const
      {
        return {values, index_, rows, cols, offset};
      }

      template <typename DerivedT>
      int evaluateLocalAlgebraicJacobian()
      {
        return evaluateSparseJetJacobian<DerivedT>();
      }

      template <typename DerivedT>
      int evaluateInputAlgebraicJacobian()
      {
        return evaluateSparseJetJacobian<DerivedT>();
      }

      template <typename DerivedT>
      int evaluateLocalDifferentialJacobian()
      {
        return evaluateSparseJetJacobian<DerivedT>();
      }

      template <typename DerivedT>
      int evaluateInputDifferentialJacobian()
      {
        return evaluateSparseJetJacobian<DerivedT>();
      }

      template <typename DerivedT>
      int evaluateSparseJetJacobian()
      {
        auto* model = static_cast<DerivedT*>(this);

        J_.zeroMatrix();
        prepareJetInputs();

        std::vector<JetT> yj(static_cast<size_t>(size()));
        std::vector<JetT> ypj(static_cast<size_t>(size()));
        std::vector<JetT> fj(static_cast<size_t>(size()));

        for (IdxT i = 0; i < size(); ++i)
        {
          const auto n = static_cast<size_t>(i);
          yj[n]        = makeJetVariable(y_[n], index_[n], RealT{1.0});
          ypj[n]       = makeJetVariable(yp_[n], index_[n], alpha_);
        }

        model->template evaluateResidualKernel<JetT>(yj.data(), ypj.data(), nullptr, fj.data());

        std::vector<IdxT>  rows;
        std::vector<IdxT>  cols;
        std::vector<RealT> vals;
        for (IdxT row = 0; row < size(); ++row)
        {
          const IdxT global_row = index_[row];
          for (const auto& [col, value] : fj[static_cast<size_t>(row)].tangent)
          {
            if (col == INVALID_INDEX<IdxT>)
            {
              continue;
            }
            rows.push_back(global_row);
            cols.push_back(col);
            vals.push_back(value);
          }
        }

        if (!vals.empty())
        {
          J_.setValues(RealT{1.0},
                       rows.data(),
                       cols.data(),
                       vals.data(),
                       static_cast<IdxT>(vals.size()));
          J_.deduplicate();
        }

        clearJetInputs();
        return 0;
      }

      void prepareJetInputs()
      {
        IdxT state_size      = 0;
        IdxT derivative_size = 0;
        for (const auto signal : inputs())
        {
          const IdxT end = signal.offset + signal.size();
          if (signal.storage == SignalStorageT::Derivative)
          {
            derivative_size = std::max(derivative_size, end);
          }
          else
          {
            state_size = std::max(state_size, end);
          }
        }

        jet_gy_.assign(static_cast<size_t>(state_size), JetT{});
        jet_gyp_.assign(static_cast<size_t>(derivative_size), JetT{});

        for (const auto signal : inputs())
        {
          auto*       jet_values    = mutableJetInputBase(signal);
          const auto* scalar_values = inputValues(signal);
          const auto* indices       = inputIndices(signal);
          const RealT scale         = inputJacobianScale(signal);

          for (IdxT i = 0; i < signal.size(); ++i)
          {
            jet_values[signal.offset + i] =
                makeJetVariable(scalar_values[i], indices[i], scale);
          }
        }
      }

      void clearJetInputs()
      {
        jet_gy_.clear();
        jet_gyp_.clear();
      }

      JetT makeJetVariable(ScalarT value, IdxT col, RealT scale) const
      {
        if (col == INVALID_INDEX<IdxT>)
        {
          return JetT(static_cast<RealT>(value));
        }
        return JetT::variable(static_cast<RealT>(value), col, scale);
      }

      IdxT  base_{0};
      RealT coordinate_{0.0};
      RealT alpha_{0.0};

      ScalarT* gy_{nullptr};
      ScalarT* gyp_{nullptr};
      IdxT*    gindex_{nullptr};
      ScalarT* y_{nullptr};
      ScalarT* yp_{nullptr};
      ScalarT* f_{nullptr};
      IdxT*    index_{nullptr};

      MatrixT J_;

    private:
      ScalarT* inputBase(SignalT signal) const
      {
        return signal.storage == SignalStorageT::Derivative ? gyp_ : gy_;
      }

      JetT* mutableJetInputBase(SignalT signal)
      {
        return signal.storage == SignalStorageT::Derivative ? jet_gyp_.data() : jet_gy_.data();
      }

      JetT* jetInputBase(SignalT signal) const
      {
        return signal.storage == SignalStorageT::Derivative ? jet_gyp_.data() : jet_gy_.data();
      }

      bool isActiveInput(SignalT signal) const
      {
        return has_active_input_
               && active_input_.offset == signal.offset
               && active_input_.rows == signal.rows
               && active_input_.cols == signal.cols
               && active_input_.storage == signal.storage;
      }

      mutable SignalT           active_input_{};
      mutable bool              has_active_input_{false};
      mutable std::vector<JetT> jet_gy_;
      mutable std::vector<JetT> jet_gyp_;
    };
  } // namespace EMT
} // namespace GridKit
