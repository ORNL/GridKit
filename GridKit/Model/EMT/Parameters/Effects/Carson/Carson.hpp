/**
 * @file Carson.hpp
 *
 * @brief Carson earth-return series resistance and inductance for overhead
 * conductors.
 *
 */

#pragma once

#include <span>
#include <vector>

#include <GridKit/Model/EMT/Element.hpp>
#include <GridKit/Model/EMT/Parameters/Effects/Carson/CarsonData.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Tower/Tower.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename scalar_type, typename index_type>
      class Carson : public Element<scalar_type, index_type>
      {
      public:
        using ScalarT = scalar_type;
        using IdxT    = index_type;
        using SignalT = typename Element<ScalarT, IdxT>::SignalT;

        Carson(const Tower<ScalarT, IdxT>& tower,
               CarsonData<ScalarT, IdxT>   data = {});

        IdxT size() const override
        {
          return layout_.n;
        }

        std::span<const SignalT> inputs() const override
        {
          return {};
        }

        ScalarT sigma() const
        {
          return data_.earth_sigma;
        }

        ScalarT epsilon() const
        {
          return data_.earth_eps;
        }

        SignalT eta() const
        {
          return {base_ + layout_.eta, 1, 1};
        }

        SignalT gammaR() const
        {
          return {base_ + layout_.gamma_r, 1, 1};
        }

        SignalT gammaI() const
        {
          return {base_ + layout_.gamma_i, 1, 1};
        }

        SignalT Rcarson() const
        {
          return {base_ + layout_.r, layout_.K, layout_.K};
        }

        SignalT Lcarson() const
        {
          return {base_ + layout_.l, layout_.K, layout_.K};
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
              eta(0),
              gamma_r(1),
              gamma_i(2),
              r(3),
              l(r + KK),
              n(l + KK)
          {
          }

          IdxT ij(IdxT i, IdxT j) const
          {
            return i * K + j;
          }

          const IdxT K;
          const IdxT KK;
          const IdxT eta;
          const IdxT gamma_r;
          const IdxT gamma_i;
          const IdxT r;
          const IdxT l;
          const IdxT n;
        };

        template <typename V>
        struct PairT
        {
          V r{0.0};
          V l{0.0};
        };

        using Pair = PairT<ScalarT>;

        using Element<ScalarT, IdxT>::base_;
        using Element<ScalarT, IdxT>::y_;
        using Element<ScalarT, IdxT>::yp_;
        using Element<ScalarT, IdxT>::f_;
        using Element<ScalarT, IdxT>::state;
        using Element<ScalarT, IdxT>::slice;
        using Element<ScalarT, IdxT>::coordinate;

        template <typename V>
        PairT<V> evaluatePair(IdxT i, IdxT j, V gr, V gi) const;

        const Tower<ScalarT, IdxT>& tower_;
        CarsonData<ScalarT, IdxT>   data_;
        const Layout                layout_;
      };
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
