/**
 * @file ShuntPotential.hpp
 *
 * @brief Potential-derived zero conductance and capacitance for overhead
 * conductors.
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
      class ShuntPotential : public Element<scalar_type, index_type>
      {
      public:
        using ScalarT = scalar_type;
        using IdxT    = index_type;
        using SignalT = typename Element<ScalarT, IdxT>::SignalT;

        ShuntPotential(const Tower<ScalarT, IdxT>&     tower,
                       const Conductor<ScalarT, IdxT>& conductor);

        IdxT size() const override
        {
          return layout_.n;
        }

        std::span<const SignalT> inputs() const override
        {
          return {};
        }

        SignalT Gpot() const
        {
          return {base_ + layout_.g, layout_.K, layout_.K};
        }

        SignalT Cpot() const
        {
          return {base_ + layout_.c, layout_.K, layout_.K};
        }

        int initialize() override;

        template <typename V>
        void evaluateResidualKernel(V* y, V*, V*, V* f) const;

        void evaluateInternalResidual(ScalarT* y, ScalarT*, ScalarT*, ScalarT* f) const;

        void tagDifferentiable(std::vector<bool>&) const override
        {
        }

        int evaluateResidual() override;
        int evaluateJacobian() override;

      private:
        struct Layout
        {
          explicit Layout(IdxT conductor_count)
            : K(conductor_count),
              KK(K * K),
              p(0),
              g(KK),
              c(2 * KK),
              n(3 * KK)
          {
          }

          IdxT ij(IdxT i, IdxT j) const
          {
            return i * K + j;
          }

          const IdxT K;
          const IdxT KK;
          const IdxT p;
          const IdxT g;
          const IdxT c;
          const IdxT n;
        };

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

        void solveCapacitance(ScalarT* y) const;

        const Tower<ScalarT, IdxT>&     tower_;
        const Conductor<ScalarT, IdxT>& conductor_;
        const Layout                    layout_;
      };
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
