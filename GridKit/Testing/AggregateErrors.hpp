#pragma once

#include <cmath>
#include <format>
#include <limits>
#include <memory>
#include <vector>

#include <GridKit/Testing/OutputAtTime.hpp>

namespace GridKit
{
  namespace Testing
  {

    inline constexpr double DEFAULT_ABS_ERROR_THRESHOLD =
        std::numeric_limits<double>::epsilon();

    /**
     * @brief Aggregate norm data for a single variable over time
     */
    struct TemporalNormAggregate
    {
      /// Name of variable
      std::string label;
      /// L-inf norm
      double      max_value{0.0};
      /// Time at which max value occurred
      double      max_value_time{0.0};
      /// L2 norm
      double      L2{0.0};

      /**
       * @brief Add a new value for the variable to be aggregated
       */
      void push(double val, double time)
      {
        if (val > max_value)
        {
          max_value      = val;
          max_value_time = time;
        }
        L2 += val * val;
      }

      /**
       * @brief Finalize the calculation(s)
       */
      void wrap()
      {
        L2 = std::sqrt(L2);
      }

      /**
       * @brief Scale by a reference value (only if the current value is above
       * the given threshold (used internally for relative errors)
       */
      void scale(const TemporalNormAggregate& ref, double threshold)
      {
        if (ref.max_value > threshold)
        {
          max_value /= ref.max_value;
        }
        if (ref.L2 > threshold)
        {
          L2 /= ref.L2;
        }
      }

      /**
       * @brief Pretty-print label and values to output stream
       */
      std::ostream& display(
          std::ostream& os = std::cout, const std::string& indent = "  ") const
      {
        os << indent << label << ":\n"
           << indent << indent << "max     : "
           << std::format("{:.6e} (at time {:.3e})", max_value, max_value_time)
           << '\n'
           << indent << indent << "L2-norm : "
           << std::format("{:.6e}", L2) << '\n';
        return os;
      }
    };

    /**
     * @brief A set of aggregate norms (represented with TemporalNormAggregate)
     * for the error in each variable plus one for the total (combined) error
     *
     * @note The "total_error" aggregate is based on the L-infinity norm of the
     * local error of variables at a given time step.
     */
    struct ErrorSet
    {
      /// Aggregate of the combined error value at each time step
      TemporalNormAggregate              total_error{"Total"};
      /// Aggregate error for each variable
      std::vector<TemporalNormAggregate> var_errors{};

      /**
       * @brief Construct with variable labels
       */
      template <typename C>
      explicit ErrorSet(const C& labels)
        : var_errors(std::size(labels))
      {
        for (std::size_t i = 0; i < std::size(labels); ++i)
        {
          var_errors[i].label = labels[i];
        }
      }

      virtual ~ErrorSet()
      {
      }

      /**
       * @brief Finalize the calculations
       */
      virtual void wrap()
      {
        for (auto& agg : var_errors)
        {
          agg.wrap();
        }
        total_error.wrap();
      }

      /**
       * @brief Take the error between the two output parameters for each
       * variable and add to the aggregate
       */
      virtual void push(const OutputAtTime&, const OutputAtTime&) = 0;

      /**
       * @brief Pretty-print the set of errors for each variable and total
       */
      std::ostream& display(std::ostream& os = std::cout) const
      {
        std::string indent{"  "};
        os << "Error Set:\n";
        for (const auto& var : var_errors)
        {
          var.display(os, indent);
        }
        total_error.display(os, indent);
        return os;
      }
    };

    enum class ErrorType
    {
      RELATIVE,
      ABSOLUTE
    };

    struct RelativeError;
    struct AbsoluteError;

    template <typename error_type = RelativeError>
    struct ErrorSetImpl;

    template <>
    struct ErrorSetImpl<RelativeError> : ErrorSet
    {
      template <typename C>
      explicit ErrorSetImpl(
          const C& labels,
          double   abs_threshold = DEFAULT_ABS_ERROR_THRESHOLD)
        : ErrorSet(labels),
          abs_threshold_(abs_threshold),
          ref_norms_(std::size(labels))
      {
      }

      void push(const OutputAtTime& test, const OutputAtTime& ref) override
      {
        auto err = test - ref;
        for (std::size_t i = 0; i < var_errors.size(); ++i)
        {
          var_errors[i].push(std::abs(err.data[i]), err.t);
          ref_norms_[i].push(std::abs(ref.data[i]), ref.t);
        }
        total_error.push(lInfNorm(err), err.t);
        ref_total_norm_.push(lInfNorm(ref), ref.t);
      }

      void wrap() override
      {
        ErrorSet::wrap();
        for (std::size_t i = 0; i < var_errors.size(); ++i)
        {
          ref_norms_[i].wrap();
          var_errors[i].scale(ref_norms_[i], abs_threshold_);
        }
        ref_total_norm_.wrap();
        total_error.scale(ref_total_norm_, abs_threshold_);
      }

    private:
      double                             abs_threshold_{};
      TemporalNormAggregate              ref_total_norm_{};
      std::vector<TemporalNormAggregate> ref_norms_{};
    };

    template <>
    struct ErrorSetImpl<AbsoluteError> : ErrorSet
    {
      using ErrorSet::ErrorSet;

      void push(const OutputAtTime& test, const OutputAtTime& ref) override
      {
        auto err = test - ref;
        for (std::size_t i = 0; i < var_errors.size(); ++i)
        {
          var_errors[i].push(std::abs(err.data[i]), err.t);
        }
        total_error.push(lInfNorm(err), err.t);
      }
    };

    template <typename C>
    std::unique_ptr<ErrorSet> makeErrorSet(
        ErrorType type,
        const C&  labels,
        double    abs_threshold = DEFAULT_ABS_ERROR_THRESHOLD)
    {
      switch (type)
      {
      case ErrorType::RELATIVE:
        return std::make_unique<ErrorSetImpl<RelativeError>>(labels, abs_threshold);
      case ErrorType::ABSOLUTE:
        return std::make_unique<ErrorSetImpl<AbsoluteError>>(labels);
      default:
        throw std::runtime_error("Invalid error type");
      }
    }

  } // namespace Testing
} // namespace GridKit
