/**
 * @file Path.hpp
 *
 * @brief Zero-state EMT element that owns static path-length geometry.
 *
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <GridKit/Model/EMT/Constants.hpp>
#include <GridKit/Model/EMT/Element.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Path/PathData.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Tower/Tower.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename scalar_type, typename index_type>
      class Path : public Element<scalar_type, index_type>
      {
      public:
        using ScalarT = scalar_type;
        using IdxT    = index_type;
        using SignalT = typename Element<ScalarT, IdxT>::SignalT;

        Path(PathData<ScalarT, IdxT> data, const Tower<ScalarT, IdxT>& tower)
          : data_(std::move(data)),
            tower_(tower)
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

        ScalarT length() const
        {
          return length_;
        }

        static constexpr ScalarT earthRadius()
        {
          return static_cast<ScalarT>(6371008.8);
        }

        int initialize() override
        {
          length_ = computeLength(data_, tower_);
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

        static void validatePoint(const PathPoint<ScalarT>& point, size_t point_index)
        {
          const auto context = "Path point " + std::to_string(point_index);
          if (!std::isfinite(static_cast<double>(point.latitude)))
          {
            throw std::runtime_error(context + ": latitude must be finite");
          }
          if (!std::isfinite(static_cast<double>(point.longitude)))
          {
            throw std::runtime_error(context + ": longitude must be finite");
          }
          if (point.latitude < static_cast<ScalarT>(-90.0)
              || point.latitude > static_cast<ScalarT>(90.0))
          {
            throw std::runtime_error(context + ": latitude must be in [-90, 90]");
          }
          if (point.longitude < static_cast<ScalarT>(-180.0)
              || point.longitude > static_cast<ScalarT>(180.0))
          {
            throw std::runtime_error(context + ": longitude must be in [-180, 180]");
          }
        }

        static ScalarT greatCircleDistance(const PathPoint<ScalarT>& point0,
                                           const PathPoint<ScalarT>& point1)
        {
          const ScalarT pi       = Constants::pi<ScalarT>();
          const ScalarT to_rad   = pi / static_cast<ScalarT>(180.0);
          const ScalarT phi0     = point0.latitude * to_rad;
          const ScalarT phi1     = point1.latitude * to_rad;
          const ScalarT lambda0  = point0.longitude * to_rad;
          const ScalarT lambda1  = point1.longitude * to_rad;
          const ScalarT dphi     = phi1 - phi0;
          const ScalarT dlambda  = lambda1 - lambda0;
          const ScalarT sdphi    = std::sin(dphi / static_cast<ScalarT>(2.0));
          const ScalarT sdlambda = std::sin(dlambda / static_cast<ScalarT>(2.0));
          const ScalarT a =
              sdphi * sdphi + std::cos(phi0) * std::cos(phi1) * sdlambda * sdlambda;
          const ScalarT ac = std::clamp(a, static_cast<ScalarT>(0.0), static_cast<ScalarT>(1.0));
          const ScalarT c  = static_cast<ScalarT>(2.0)
                            * std::atan2(std::sqrt(ac),
                                         std::sqrt(static_cast<ScalarT>(1.0) - ac));
          return earthRadius() * c;
        }

        static ScalarT computePathLength(const std::vector<PathPoint<ScalarT>>& path)
        {
          if (path.size() < 2)
          {
            throw std::runtime_error("Path requires at least two GIS points when length is not supplied");
          }

          ScalarT length = 0.0;
          for (size_t i = 0; i < path.size(); ++i)
          {
            validatePoint(path[i], i);
            if (i > 0)
            {
              length += greatCircleDistance(path[i - 1], path[i]);
            }
          }
          return length;
        }

        static ScalarT computeLength(const PathData<ScalarT, IdxT>& data,
                                     const Tower<ScalarT, IdxT>&    tower)
        {
          if (!data.path.empty() && data.path.size() < 2)
          {
            throw std::runtime_error("Path GIS data must contain at least two points");
          }
          for (size_t i = 0; i < data.path.size(); ++i)
          {
            validatePoint(data.path[i], i);
          }

          if (data.length.has_value())
          {
            const ScalarT length = data.length.value();
            if (!std::isfinite(static_cast<double>(length)) || length <= static_cast<ScalarT>(0.0))
            {
              throw std::runtime_error("Path length must be positive and finite");
            }
            return length;
          }

          const ScalarT length = computePathLength(data.path) * tower.saggedToSpanRatio();
          if (!std::isfinite(static_cast<double>(length)) || length <= static_cast<ScalarT>(0.0))
          {
            throw std::runtime_error("Path derived length must be positive and finite");
          }
          return length;
        }

        PathData<ScalarT, IdxT>     data_;
        const Tower<ScalarT, IdxT>& tower_;
        ScalarT                     length_{0.0};
      };
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
