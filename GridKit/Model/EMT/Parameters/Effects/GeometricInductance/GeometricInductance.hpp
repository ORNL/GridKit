/**
 * @file GeometricInductance.hpp
 *
 * @brief External geometric series inductance for overhead conductors.
 *
 */

#pragma once

#include <cmath>
#include <span>
#include <vector>

#include <GridKit/LinearAlgebra/DenseExpression.hpp>
#include <GridKit/Model/EMT/Constants.hpp>
#include <GridKit/Model/EMT/Element.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Conductor/Conductor.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Tower/Tower.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename scalar_type, typename index_type>
      class GeometricInductance : public Element<scalar_type, index_type>
      {
      public:
        using ScalarT = scalar_type;
        using IdxT    = index_type;
        using SignalT = typename Element<ScalarT, IdxT>::SignalT;

        GeometricInductance(const Tower<ScalarT, IdxT>&     tower,
                            const Conductor<ScalarT, IdxT>& conductor)
          : tower_(tower),
            conductor_(conductor)
        {
        }

        IdxT size() const override
        {
          const IdxT K = conductor_.conductorCount();
          return K * K;
        }

        std::span<const SignalT> inputs() const override
        {
          return {};
        }

        SignalT Lgeo() const
        {
          const IdxT K = conductor_.conductorCount();
          return {base_, K, K};
        }

        int initialize() override
        {
          const IdxT    K     = conductor_.conductorCount();
          const ScalarT scale = Constants::mu0<ScalarT>() / (2.0 * Constants::pi<ScalarT>());
          auto          Lgeo  = state(0, K, K);
          for (IdxT i = 0; i < K; ++i)
          {
            for (IdxT j = 0; j < K; ++j)
            {
              Lgeo(i, j) = scale * lambda(i, j);
            }
          }
          return 0;
        }

        template <typename V>
        __attribute__((always_inline)) inline void
        evaluateResidualKernel(V* y, V*, V*, V* f) const
        {
          using namespace GridKit::LinearAlgebra;

          const IdxT K      = conductor_.conductorCount();
          const auto two_pi = 2.0 * Constants::pi<ScalarT>();
          const auto mu0    = Constants::mu0<ScalarT>();
          auto       Lgeo   = slice(y, 0, K, K);
          auto       F      = slice(f, 0, K, K);
          auto       Lambda = matrix(K,
                               K,
                               [&](IdxT i, IdxT j)
                               {
                                 return lambda(i, j);
                               });

          equation(F) = -two_pi * Lgeo + mu0 * Lambda;
        }

        void evaluateInternalResidual(ScalarT* y, ScalarT* yp, ScalarT* u, ScalarT* f) const
        {
          evaluateResidualKernel<ScalarT>(y, yp, u, f);
        }

        void tagDifferentiable(std::vector<bool>&) const override
        {
        }

        int evaluateResidual() override
        {
          evaluateInternalResidual(y_, yp_, nullptr, f_);
          return 0;
        }

        int evaluateJacobian() override;

      private:
        using Element<ScalarT, IdxT>::base_;
        using Element<ScalarT, IdxT>::y_;
        using Element<ScalarT, IdxT>::yp_;
        using Element<ScalarT, IdxT>::f_;
        using Element<ScalarT, IdxT>::state;
        using Element<ScalarT, IdxT>::slice;

        ScalarT lambda(IdxT i, IdxT j) const
        {
          if (i == j)
          {
            return std::log(2.0 * tower_.height(i) / conductor_.radius(i));
          }
          return std::log(tower_.imageDistance(i, j) / tower_.distance(i, j));
        }

        const Tower<ScalarT, IdxT>&     tower_;
        const Conductor<ScalarT, IdxT>& conductor_;
      };
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
