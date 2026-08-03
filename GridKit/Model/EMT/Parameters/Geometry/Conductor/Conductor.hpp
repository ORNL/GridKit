/**
 * @file Conductor.hpp
 *
 * @brief Static conductor geometry, material, weight, and phase parameter model.
 *
 */

#pragma once

#include <cmath>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <GridKit/Model/EMT/Element.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Conductor/ConductorData.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename scalar_type, typename index_type>
      class Conductor : public Element<scalar_type, index_type>
      {
      public:
        using ScalarT = scalar_type;
        using IdxT    = index_type;
        using SignalT = typename Element<ScalarT, IdxT>::SignalT;

        explicit Conductor(ConductorData<ScalarT, IdxT> data)
          : data_(std::move(data))
        {
        }

        IdxT size() const override
        {
          return 0;
        }

        std::span<const SignalT> inputs() const override
        {
          return {};
        }

        IdxT conductorCount() const
        {
          return data_.K;
        }

        const std::vector<ScalarT>& radius() const
        {
          return data_.radius;
        }

        ScalarT radius(IdxT i) const
        {
          return data_.radius[static_cast<size_t>(i)];
        }

        ScalarT innerRadius(IdxT i) const
        {
          return data_.inner_radius[static_cast<size_t>(i)];
        }

        ScalarT sigma(IdxT i) const
        {
          return data_.sigma[static_cast<size_t>(i)];
        }

        ScalarT mu(IdxT i) const
        {
          return data_.mu[static_cast<size_t>(i)];
        }

        const std::vector<ScalarT>& weight() const
        {
          return data_.weight;
        }

        ScalarT weight(IdxT i) const
        {
          return data_.weight[static_cast<size_t>(i)];
        }

        const std::vector<std::string>& phase() const
        {
          return data_.phase;
        }

        const std::string& phase(IdxT i) const
        {
          return data_.phase[static_cast<size_t>(i)];
        }

        int initialize() override
        {
          validateInputs();
          return 0;
        }

        void tagDifferentiable(std::vector<bool>&) const override
        {
        }

        int evaluateResidual() override
        {
          return 0;
        }

        int evaluateJacobian() override
        {
          J_.zeroMatrix();
          return 0;
        }

      private:
        using Element<ScalarT, IdxT>::J_;

        static bool finite(ScalarT value)
        {
          return std::isfinite(static_cast<double>(value));
        }

        static bool positiveFinite(ScalarT value)
        {
          return finite(value) && value > static_cast<ScalarT>(0.0);
        }

        static bool validPhase(const std::string& phase)
        {
          return phase == "a" || phase == "b" || phase == "c" || phase == "n" || phase == "g";
        }

        void validateInputs() const
        {
          const auto K = static_cast<size_t>(data_.K);
          if (data_.K == 0)
          {
            throw std::runtime_error("Conductor requires at least one conductor");
          }
          if (data_.radius.size() != K || data_.inner_radius.size() != K
              || data_.sigma.size() != K || data_.mu.size() != K || data_.weight.size() != K)
          {
            throw std::runtime_error("Conductor vector sizes must match conductor count");
          }

          for (IdxT i = 0; i < data_.K; ++i)
          {
            const auto n = static_cast<size_t>(i);
            if (!positiveFinite(data_.radius[n]))
            {
              throw std::runtime_error("Conductor radii must be positive and finite");
            }
            if (!finite(data_.inner_radius[n]) || data_.inner_radius[n] < static_cast<ScalarT>(0.0)
                || data_.inner_radius[n] >= data_.radius[n])
            {
              throw std::runtime_error("Conductor inner radii must be finite and less than outer radii");
            }
            if (!positiveFinite(data_.sigma[n]))
            {
              throw std::runtime_error("Conductor conductivities must be positive and finite");
            }
            if (!positiveFinite(data_.mu[n]))
            {
              throw std::runtime_error("Conductor permeabilities must be positive and finite");
            }
            if (!positiveFinite(data_.weight[n]))
            {
              throw std::runtime_error("Conductor weights must be positive and finite");
            }
          }

          if (!data_.phase.empty())
          {
            if (data_.phase.size() != K)
            {
              throw std::runtime_error("Conductor phase size must match conductor count");
            }
            for (const auto& phase : data_.phase)
            {
              if (!validPhase(phase))
              {
                throw std::runtime_error("Conductor phase labels must be one of a, b, c, n, or g");
              }
            }
          }
        }

        ConductorData<ScalarT, IdxT> data_;
      };
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
