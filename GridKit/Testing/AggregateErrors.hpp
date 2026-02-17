#pragma once

#include <cmath>
#include <format>
#include <vector>

#include <GridKit/Testing/OutputAtTime.hpp>

namespace GridKit
{
  namespace Testing
  {

    /**
     * @brief Aggregate error data for a single variable
     */
    struct ErrorAggregate
    {
      std::string label;
      double      max_error{0.0};
      double      max_error_time{0.0};
      double      error_norm_L2{0.0};

      void push(double err, double time)
      {
        if (err > max_error)
        {
          max_error      = err;
          max_error_time = time;
        }
        error_norm_L2 += err * err;
      }

      void wrap()
      {
        error_norm_L2 = std::sqrt(error_norm_L2);
      }

      std::ostream& display(
          std::ostream& os = std::cout, const std::string& indent = "  ") const
      {
        os << indent << label << ":\n"
           << indent << indent << "(max, time): "
           << std::format("{:.6e}, {:.3e}", max_error, max_error_time) << '\n'
           << indent << indent << "L2-norm    : "
           << std::format("{:.6e}", error_norm_L2) << '\n';
        return os;
      }
    };

    /**
     * @brief A set of ErrorAggregate for each variable plus a total (combined)
     * ErrorAggregate
     *
     * @note The "total" aggregate is based on the L-infinity norm of the local
     * error of variables at a given time step.
     */
    struct ErrorSet
    {
      ErrorAggregate              total{"Total"};
      std::vector<ErrorAggregate> vars{};

      template <typename C>
      explicit ErrorSet(const C& labels)
        : vars(std::size(labels))
      {
        for (std::size_t i = 0; i < std::size(labels); ++i)
        {
          vars[i].label = labels[i];
        }
      }

      void push(const OutputAtTime& err)
      {
        for (std::size_t i = 0; i < vars.size(); ++i)
        {
          vars[i].push(std::abs(err.data[i]), err.t);
        }
        total.push(lInfNorm(err), err.t);
      }

      void wrap()
      {
        for (auto& agg : vars)
        {
          agg.wrap();
        }
        total.wrap();
      }

      std::ostream& display(std::ostream& os = std::cout) const
      {
        std::string indent{"  "};
        os << "Error Set:\n";
        for (const auto& var : vars)
        {
          var.display(os, indent);
        }
        total.display(os, indent);
        return os;
      }
    };

  } // namespace Testing
} // namespace GridKit
