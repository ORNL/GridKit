#pragma once

#include <array>
#include <cmath>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <typename RealT, typename InternalT, typename ExternalT>
    struct ComponentTestPoint
    {
      using StateValues = std::vector<std::pair<InternalT, RealT>>;
      using InputValues = std::vector<std::pair<ExternalT, RealT>>;

      InputValues inputs;
      StateValues state;
      StateValues derivative;
    };

    /// Owns a signal-connected component and its signal storage.
    template <template <typename, typename> class model_type,
              typename scalar_type,
              typename index_type>
    class ComponentTestFixture
    {
    public:
      using ModelT    = model_type<scalar_type, index_type>;
      using ScalarT   = scalar_type;
      using IdxT      = index_type;
      using RealT     = typename ModelT::RealT;
      using InternalT = typename ModelT::InternalVariablesT;
      using ExternalT = typename ModelT::ExternalVariablesT;
      using Data      = typename ModelT::ModelDataT;
      using Point     = ComponentTestPoint<RealT, InternalT, ExternalT>;
      using Values    = typename Point::StateValues;
      using Inputs    = typename Point::InputValues;

      struct Snapshot
      {
        Values state;
        Values derivative;
        Inputs inputs;
        Values outputs;
      };

      ComponentTestFixture(const Data& data, const char* context, RealT tolerance)
        : model_(data), context_(context), tolerance_(tolerance)
      {
      }

      ComponentTestFixture(const ComponentTestFixture&)            = delete;
      ComponentTestFixture& operator=(const ComponentTestFixture&) = delete;
      ComponentTestFixture(ComponentTestFixture&&)                 = delete;
      ComponentTestFixture& operator=(ComponentTestFixture&&)      = delete;

      /// Connect signals before prepare().
      template <ExternalT variable>
      void attachInput(RealT value)
      {
        static_assert(variable < ExternalT::MAXIMUM);
        auto& input = inputs_[static_cast<size_t>(variable)];
        input.value = value;
        input.index = rowCount() + static_cast<IdxT>(variable);
        input.node.set(&input.value, &input.index);
        model_.getSignals().template attachSignalNode<variable>(&input.node);
      }

      template <InternalT variable>
      void assignOutput()
      {
        static_assert(variable < InternalT::MAXIMUM);
        model_.getSignals().template assignSignalNode<variable>(
            &outputs_[static_cast<size_t>(variable)]);
      }

      bool prepare()
      {
        if (!allocated_)
        {
          if (!checkStatus(model_.allocate(), "allocate"))
            return false;
          if (model_.size() != rowCount())
            return fail("internal enum does not match the component size");
          for (IdxT row = 0; row < rowCount(); ++row)
          {
            if (!checkStatus(model_.setVariableIndex(row, row), "setVariableIndex")
                || !checkStatus(model_.setResidualIndex(row, row), "setResidualIndex"))
              return false;
          }
          allocated_ = true;
        }
        return checkStatus(model_.verify(), "verify");
      }

      /// Preparation and initialization only; residual checks evaluate explicitly.
      bool initialize(const Values& seed = {})
      {
        initialized_ = prepare() && setState(seed)
                       && checkStatus(model_.initialize(), "initialize");
        return initialized_;
      }

      bool setInput(ExternalT variable, RealT value)
      {
        const auto port = static_cast<size_t>(variable);
        if (port >= inputs_.size() || !inputs_[port].node.linked())
          return fail("input is not attached");
        inputs_[port].value = value;
        return true;
      }

      bool setState(const Values& values)
      {
        return setValues(model_.y(), values);
      }

      bool setDerivative(const Values& values)
      {
        return setValues(model_.yp(), values);
      }

      /// Omitted entries retain their values; inputs must already be attached.
      bool setPoint(const Point& point)
      {
        for (const auto& [port, value] : point.inputs)
        {
          if (!setInput(port, value))
            return false;
        }
        return setState(point.state) && setDerivative(point.derivative);
      }

      RealT state(InternalT variable) const
      {
        return static_cast<RealT>(model_.y().getData()[static_cast<size_t>(variable)]);
      }

      RealT residual(InternalT variable) const
      {
        return static_cast<RealT>(model_.getResidual().getData()[static_cast<size_t>(variable)]);
      }

      RealT input(ExternalT variable) const
      {
        return static_cast<RealT>(inputs_.at(static_cast<size_t>(variable)).node.read());
      }

      RealT output(InternalT variable) const
      {
        return static_cast<RealT>(outputs_.at(static_cast<size_t>(variable)).read());
      }

      bool evaluateResidual()
      {
        if (!initialized_)
          return fail("component is not initialized");
        return checkStatus(model_.evaluateResidual(), "evaluateResidual");
      }

      /// Require every residual row exactly once, including zero rows.
      bool checkResiduals(const Values& expected, const char* label = "")
      {
        return evaluateResidual() && checkRows(model_.getResidual(), expected, label, true);
      }

      bool checkResiduals(const Point& point, const Values& expected, const char* label = "")
      {
        return setPoint(point) && checkResiduals(expected, label);
      }

      bool checkResidualRows(const Values& expected, const char* label = "")
      {
        return evaluateResidual() && checkRows(model_.getResidual(), expected, label);
      }

      bool checkResidualRows(const Point& point, const Values& expected, const char* label = "")
      {
        return setPoint(point) && checkResidualRows(expected, label);
      }

      bool checkStateRows(const Values& expected, const char* label = "") const
      {
        return checkRows(model_.y(), expected, label);
      }

      bool checkDerivativeRows(const Values& expected, const char* label = "") const
      {
        return checkRows(model_.yp(), expected, label);
      }

      bool checkSteadyState()
      {
        if (!evaluateResidual())
          return false;
        bool success = true;
        for (IdxT row = 0; row < rowCount(); ++row)
        {
          success &= checkValue(model_.getResidual().getData()[row], 0.0, "residual at rest", row);
          success &= checkValue(model_.yp().getData()[row], 0.0, "derivative at rest", row);
        }
        return success;
      }

      Snapshot snapshot() const
      {
        Snapshot result{values(model_.y()), values(model_.yp()), {}, {}};
        for (size_t port = 0; port < inputs_.size(); ++port)
        {
          if (inputs_[port].node.linked())
            result.inputs.emplace_back(static_cast<ExternalT>(port), input(static_cast<ExternalT>(port)));
        }
        for (size_t row = 0; row < outputs_.size(); ++row)
        {
          if (outputs_[row].linked())
            result.outputs.emplace_back(static_cast<InternalT>(row), output(static_cast<InternalT>(row)));
        }
        return result;
      }

      /// Exact value preservation, treating a retained NaN as unchanged.
      bool checkUnchanged(const Snapshot& before) const
      {
        bool success  = checkRows(model_.y(), before.state, "preserved state", true, true);
        success      &= checkRows(model_.yp(), before.derivative, "preserved derivative", true, true);
        for (const auto& [port, expected] : before.inputs)
          success &= checkValue(input(port), expected, "preserved input", static_cast<size_t>(port), true);
        for (const auto& [row, expected] : before.outputs)
          success &= checkValue(output(row), expected, "preserved output", static_cast<size_t>(row), true);
        return success;
      }

      Values residuals() const
      {
        return values(model_.getResidual());
      }

      /// Evaluate F_y + alpha F_yp once; component CSR may be cached.
      bool evaluateJacobian(RealT time, RealT alpha)
      {
        if (!initialized_ || jacobian_used_)
          return fail("Jacobian requires a fresh initialized fixture");
        jacobian_used_ = true;
        model_.updateTime(time, alpha);
        if constexpr (std::is_same_v<ScalarT, DependencyTracking::Variable>)
        {
          if (!checkStatus(model_.initializeDependencyTrackingVariableNumbers(), "number dependencies"))
            return false;
          for (auto& input : inputs_)
          {
            if (input.node.linked())
              input.value.setVariableNumber(2 * static_cast<size_t>(input.index));
          }
        }
#ifndef GRIDKIT_ENABLE_ENZYME
        else
          return fail("Enzyme is disabled");
#endif
        if (!evaluateResidual() || !checkStatus(model_.evaluateJacobian(), "evaluateJacobian"))
          return false;
        if (model_.getCsrJacobian() == nullptr && !checkStatus(model_.constructCsr(), "constructCsr"))
          return false;
        const auto* matrix = model_.getCsrJacobian();
        if (matrix == nullptr || matrix->getNumRows() != rowCount()
            || matrix->getNumColumns() != columnCount())
          return fail("unexpected Jacobian dimensions");
        return true;
      }

      bool checkJacobianRow(InternalT row, const Values& expected, RealT alpha = 0.0)
      {
        if (static_cast<size_t>(row) >= INTERNAL_COUNT)
          return fail("invalid Jacobian row");
        if (!evaluateJacobian(0.0, alpha) || !validRows(expected))
          return false;
        DependencyMap entries;
        for (const auto& [column, value] : expected)
          entries.emplace(static_cast<size_t>(column), value);
        return compareJacobianRow(MapFromCsr(model_.getCsrJacobian())[static_cast<size_t>(row)],
                                  entries,
                                  static_cast<size_t>(row),
                                  alpha);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      template <typename Setup>
      static bool checkJacobian(
          const Data&                                            data,
          Setup                                                  setup,
          std::initializer_list<RealT>                           alphas,
          const char*                                            context,
          RealT                                                  tolerance,
          std::initializer_list<std::pair<InternalT, InternalT>> required_entries = {})
      {
        if (alphas.size() == 0)
        {
          std::cout << context << ": no Jacobian coefficients specified\n";
          return false;
        }
        bool success = true;
        for (const RealT alpha : alphas)
        {
          const std::string    numeric_context  = std::string(context) + " [Enzyme]";
          const std::string    tracking_context = std::string(context) + " [DependencyTracking]";
          ComponentTestFixture numeric(data, numeric_context.c_str(), tolerance);
          ComponentTestFixture<model_type, DependencyTracking::Variable, IdxT>
              tracking(data, tracking_context.c_str(), tolerance);
          if (!setup(numeric) || !setup(tracking)
              || !numeric.evaluateJacobian(0.0, alpha)
              || !tracking.evaluateJacobian(0.0, alpha))
          {
            std::cout << context << " alpha=" << alpha << ": Jacobian setup or evaluation failed\n";
            success = false;
            continue;
          }
          const auto actual   = MapFromCsr(numeric.model().getCsrJacobian());
          const auto expected = MapFromCsr(tracking.model().getCsrJacobian());
          if (numeric.columnCount() != tracking.columnCount())
          {
            success = numeric.fail("Jacobian column counts differ");
            continue;
          }
          for (size_t row = 0; row < actual.size(); ++row)
            success &= numeric.compareJacobianRow(actual[row], expected[row], row, alpha);
          for (const auto& [row, column] : required_entries)
          {
            const auto row_index    = static_cast<size_t>(row);
            const auto column_index = static_cast<size_t>(column);
            if (row_index >= expected.size() || expected[row_index].count(column_index) == 0)
            {
              std::cout << context << " alpha=" << alpha << " row=" << row_index
                        << " column=" << column_index << ": required entry is missing\n";
              success = false;
            }
          }
        }
        return success;
      }
#endif

      ModelT& model()
      {
        return model_;
      }

      const ModelT& model() const
      {
        return model_;
      }

    private:
      template <template <typename, typename> class, typename, typename>
      friend class ComponentTestFixture;

      using VectorT                          = typename ModelT::VectorT;
      using SignalNodeT                      = PhasorDynamics::SignalNode<ScalarT, IdxT>;
      using DependencyMap                    = DependencyTracking::Variable::DependencyMap;
      static constexpr size_t INTERNAL_COUNT = static_cast<size_t>(InternalT::MAXIMUM);
      static constexpr size_t EXTERNAL_COUNT = static_cast<size_t>(ExternalT::MAXIMUM);

      static constexpr IdxT rowCount()
      {
        return static_cast<IdxT>(INTERNAL_COUNT);
      }

      IdxT columnCount() const
      {
        IdxT columns = rowCount();
        for (const auto& input : inputs_)
          if (input.node.linked())
            columns = input.index + 1;
        return columns;
      }

      bool fail(const char* message) const
      {
        std::cout << context_ << ": " << message << '\n';
        return false;
      }

      bool checkStatus(int status, const char* operation) const
      {
        if (status == 0)
          return true;
        std::cout << context_ << ": " << operation << " returned " << status << '\n';
        return false;
      }

      bool validRows(const Values& rows, bool complete = false) const
      {
        if (complete && rows.size() != INTERNAL_COUNT)
          return fail("incomplete vector rows");
        std::array<bool, INTERNAL_COUNT> seen{};
        for (const auto& [variable, value] : rows)
        {
          const auto row = static_cast<size_t>(variable);
          if (row >= INTERNAL_COUNT || seen[row])
            return fail("invalid or duplicate vector row");
          seen[row] = true;
        }
        return true;
      }

      bool setValues(VectorT& vector, const Values& rows)
      {
        if (!allocated_ || !validRows(rows))
          return fail("cannot set vector rows");
        for (const auto& [variable, value] : rows)
          vector.getData()[static_cast<size_t>(variable)] = value;
        return checkStatus(vector.setDataUpdated(), "setDataUpdated");
      }

      static Values values(const VectorT& vector)
      {
        Values result;
        for (IdxT row = 0; row < vector.getSize(); ++row)
          result.emplace_back(static_cast<InternalT>(row), static_cast<RealT>(vector.getData()[row]));
        return result;
      }

      bool checkValue(RealT actual, RealT expected, const char* label, size_t row, bool exact = false) const
      {
        const bool equal = exact ? (actual == expected || (std::isnan(actual) && std::isnan(expected)))
                                 : isEqual(actual, expected, tolerance_);
        if (equal)
          return true;
        std::cout << context_ << ' ' << label << " row=" << row << ": "
                  << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                  << actual << " != " << expected << '\n';
        return false;
      }

      bool checkRows(const VectorT& vector, const Values& expected, const char* label, bool complete = false, bool exact = false) const
      {
        if (!allocated_ || vector.getSize() != rowCount() || !validRows(expected, complete))
          return fail("cannot check vector rows");
        bool success = true;
        for (const auto& [variable, value] : expected)
          success &= checkValue(static_cast<RealT>(vector.getData()[static_cast<size_t>(variable)]),
                                value,
                                label,
                                static_cast<size_t>(variable),
                                exact);
        return success;
      }

      bool compareJacobianRow(const DependencyMap& actual,
                              const DependencyMap& expected,
                              size_t               row,
                              RealT                alpha) const
      {
        bool success        = true;
        auto actual_entry   = actual.begin();
        auto expected_entry = expected.begin();
        while (actual_entry != actual.end() || expected_entry != expected.end())
        {
          if (expected_entry == expected.end()
              || (actual_entry != actual.end() && actual_entry->first < expected_entry->first))
          {
            std::cout << context_ << " alpha=" << alpha << " row=" << row
                      << " column=" << actual_entry->first << ": unexpected entry " << actual_entry->second << '\n';
            ++actual_entry;
            success = false;
          }
          else if (actual_entry == actual.end() || expected_entry->first < actual_entry->first)
          {
            std::cout << context_ << " alpha=" << alpha << " row=" << row
                      << " column=" << expected_entry->first << ": missing entry, expected " << expected_entry->second << '\n';
            ++expected_entry;
            success = false;
          }
          else
          {
            if (!isEqual(static_cast<RealT>(actual_entry->second),
                         static_cast<RealT>(expected_entry->second),
                         tolerance_))
            {
              std::cout << context_ << " alpha=" << alpha << " row=" << row
                        << " column=" << actual_entry->first << ": "
                        << std::setprecision(std::numeric_limits<RealT>::max_digits10)
                        << actual_entry->second << " != " << expected_entry->second << '\n';
              success = false;
            }
            ++actual_entry;
            ++expected_entry;
          }
        }
        return success;
      }

      struct Input
      {
        ScalarT     value{};
        IdxT        index{};
        SignalNodeT node;
      };

      // Signal storage must outlive the model.
      std::array<Input, EXTERNAL_COUNT>       inputs_{};
      std::array<SignalNodeT, INTERNAL_COUNT> outputs_{};
      ModelT                                  model_;
      std::string                             context_;
      RealT                                   tolerance_;
      bool                                    allocated_{false};
      bool                                    initialized_{false};
      bool                                    jacobian_used_{false};
    };
  } // namespace Testing
} // namespace GridKit
