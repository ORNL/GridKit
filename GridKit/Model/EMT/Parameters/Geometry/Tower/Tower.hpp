/**
 * @file Tower.hpp
 *
 * @brief Zero-state EMT element that owns static conductor tower geometry.
 *
 */

#pragma once

#include <cmath>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include <GridKit/Model/EMT/Element.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Conductor/Conductor.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Tower/TowerData.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename scalar_type, typename index_type>
      class Tower : public Element<scalar_type, index_type>
      {
      public:
        using ScalarT = scalar_type;
        using IdxT    = index_type;
        using SignalT = typename Element<ScalarT, IdxT>::SignalT;

        Tower(TowerData<ScalarT, IdxT> data, const Conductor<ScalarT, IdxT>& conductor)
          : data_(std::move(data)),
            conductor_(conductor)
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

        ScalarT span() const
        {
          return data_.span;
        }

        ScalarT height(IdxT i) const
        {
          return effective_height_[static_cast<size_t>(i)];
        }

        ScalarT position(IdxT i) const
        {
          return data_.position[static_cast<size_t>(i)];
        }

        ScalarT horizontalDistance(IdxT i, IdxT j) const
        {
          return std::abs(position(i) - position(j));
        }

        ScalarT distance(IdxT i, IdxT j) const
        {
          return D_[static_cast<size_t>(idx(i, j))];
        }

        ScalarT imageDistance(IdxT i, IdxT j) const
        {
          return Dp_[static_cast<size_t>(idx(i, j))];
        }

        ScalarT sag(IdxT i) const
        {
          return sag_[static_cast<size_t>(i)];
        }

        ScalarT spanLength(IdxT i) const
        {
          return span_length_[static_cast<size_t>(i)];
        }

        ScalarT minimumHeight(IdxT i) const
        {
          return minimum_height_[static_cast<size_t>(i)];
        }

        ScalarT saggedToSpanRatio() const
        {
          return sagged_to_span_ratio_;
        }

        int initialize() override
        {
          validateInputs();
          computeSpanGeometry();
          computeDistances();

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

        IdxT idx(IdxT i, IdxT j) const
        {
          return i * data_.K + j;
        }

        static bool finite(ScalarT value)
        {
          return std::isfinite(static_cast<double>(value));
        }

        static bool positiveFinite(ScalarT value)
        {
          return finite(value) && value > static_cast<ScalarT>(0.0);
        }

        ScalarT attachmentHeight(IdxT i) const
        {
          return data_.height[static_cast<size_t>(i)];
        }

        void computeDistances()
        {
          const size_t data_size = static_cast<size_t>(data_.K) * static_cast<size_t>(data_.K);
          D_.assign(data_size, 0.0);
          Dp_.assign(data_size, 0.0);

          for (IdxT i = 0; i < data_.K; ++i)
          {
            for (IdxT j = 0; j < data_.K; ++j)
            {
              const size_t  si = static_cast<size_t>(i);
              const size_t  sj = static_cast<size_t>(j);
              const size_t  s  = static_cast<size_t>(idx(i, j));
              const ScalarT dx = data_.position[si] - data_.position[sj];
              const ScalarT dh = effective_height_[si] - effective_height_[sj];
              const ScalarT sh = effective_height_[si] + effective_height_[sj];

              D_[s]  = std::sqrt(dx * dx + dh * dh);
              Dp_[s] = std::sqrt(dx * dx + sh * sh);
              if (i != j && (!positiveFinite(D_[s]) || !(Dp_[s] > D_[s])))
              {
                throw std::runtime_error("Tower conductor attachment coordinates must be distinct");
              }
            }
          }
        }

        void validateInputs() const
        {
          const IdxT K = conductor_.conductorCount();
          if (data_.K == 0 || K == 0)
          {
            throw std::runtime_error("Tower requires at least one conductor");
          }
          if (data_.K != K)
          {
            throw std::runtime_error("Tower conductor count must match conductor data");
          }
          if (data_.height.size() != static_cast<size_t>(data_.K)
              || data_.position.size() != static_cast<size_t>(data_.K))
          {
            throw std::runtime_error("Tower coordinate sizes must match conductor count");
          }
          if (!positiveFinite(data_.span))
          {
            throw std::runtime_error("Tower span must be positive and finite");
          }

          for (IdxT i = 0; i < data_.K; ++i)
          {
            if (!positiveFinite(attachmentHeight(i)))
            {
              throw std::runtime_error("Tower heights must be positive and finite");
            }
            if (!finite(position(i)))
            {
              throw std::runtime_error("Tower positions must be finite");
            }
          }

          if (!data_.tension.empty())
          {
            if (data_.tension.size() != static_cast<size_t>(data_.K))
            {
              throw std::runtime_error("Tower tension size must match conductor count");
            }
            for (const auto& tension : data_.tension)
            {
              if (tension.has_value() && !positiveFinite(tension.value()))
              {
                throw std::runtime_error("Tower tension entries must be positive and finite");
              }
            }
          }
        }

        void computeSpanGeometry()
        {
          sag_.assign(static_cast<size_t>(data_.K), static_cast<ScalarT>(0.0));
          span_length_.assign(static_cast<size_t>(data_.K), data_.span);
          minimum_height_       = data_.height;
          effective_height_     = data_.height;
          sagged_to_span_ratio_ = static_cast<ScalarT>(1.0);

          if (data_.tension.empty())
          {
            return;
          }

          ScalarT ratio_sum{0.0};
          for (IdxT i = 0; i < data_.K; ++i)
          {
            const auto& tension = data_.tension[static_cast<size_t>(i)];
            if (!tension.has_value())
            {
              ratio_sum += static_cast<ScalarT>(1.0);
              continue;
            }

            const ScalarT weight = conductor_.weight(i);
            if (!positiveFinite(weight))
            {
              throw std::runtime_error("Tower conductor weights must be positive and finite");
            }

            const ScalarT a        = tension.value() / weight;
            const ScalarT eta      = data_.span / (static_cast<ScalarT>(2.0) * a);
            const ScalarT sinh_eta = std::sinh(eta);
            const ScalarT sag      = a * (std::cosh(eta) - static_cast<ScalarT>(1.0));
            const ScalarT length   = static_cast<ScalarT>(2.0) * a * sinh_eta;
            const ScalarT effective_height =
                attachmentHeight(i) - sag
                + (static_cast<ScalarT>(2.0) * a * a / data_.span) * sinh_eta - a;

            if (!finite(sag) || sag >= attachmentHeight(i))
            {
              throw std::runtime_error("Tower sag must be finite and below conductor height");
            }
            if (!positiveFinite(length))
            {
              throw std::runtime_error("Tower span length must be positive and finite");
            }
            if (!finite(effective_height))
            {
              throw std::runtime_error("Tower height must be finite");
            }

            sag_[static_cast<size_t>(i)]               = sag;
            span_length_[static_cast<size_t>(i)]       = length;
            minimum_height_[static_cast<size_t>(i)]    = attachmentHeight(i) - sag;
            effective_height_[static_cast<size_t>(i)]  = effective_height;
            ratio_sum                                 += length / data_.span;
          }

          sagged_to_span_ratio_ = ratio_sum / static_cast<ScalarT>(data_.K);
        }

        TowerData<ScalarT, IdxT>        data_;
        const Conductor<ScalarT, IdxT>& conductor_;
        std::vector<ScalarT>            D_;
        std::vector<ScalarT>            Dp_;
        std::vector<ScalarT>            sag_;
        std::vector<ScalarT>            span_length_;
        std::vector<ScalarT>            minimum_height_;
        std::vector<ScalarT>            effective_height_;
        ScalarT                         sagged_to_span_ratio_{1.0};
      };
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
