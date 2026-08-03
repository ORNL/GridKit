/**
 * @file Zc.hpp
 *
 * @brief Characteristic impedance assembly for overhead conductors.
 *
 */

#pragma once

#include <array>
#include <span>
#include <vector>

#include <GridKit/Model/EMT/Element.hpp>
#include <GridKit/Model/EMT/Parameters/Views.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename scalar_type, typename index_type>
      class Zc : public Element<scalar_type, index_type>
      {
      public:
        using ScalarT = scalar_type;
        using IdxT    = index_type;
        using SignalT = typename Element<ScalarT, IdxT>::SignalT;

        Zc(const SeriesView<ScalarT, IdxT>& series,
           const ShuntView<ScalarT, IdxT>&  shunt);

        IdxT size() const override
        {
          return layout_.n;
        }

        std::span<const SignalT> inputs() const override
        {
          inputs_[0] = series_.R();
          inputs_[1] = series_.L();
          inputs_[2] = shunt_.G();
          inputs_[3] = shunt_.C();
          return inputs_;
        }

        SignalT Rc() const
        {
          return {base_ + layout_.rc, layout_.K, layout_.K};
        }

        SignalT Xc() const
        {
          return {base_ + layout_.xc, layout_.K, layout_.K};
        }

        int initialize() override;

        template <typename V>
        void evaluateResidualKernel(V* y, V*, V* u, V* f) const;

        void evaluateInternalResidual(ScalarT* y, ScalarT*, ScalarT* u, ScalarT* f) const;

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
              rc(0),
              xc(KK),
              n(2 * KK)
          {
          }

          const IdxT K;
          const IdxT KK;
          const IdxT rc;
          const IdxT xc;
          const IdxT n;
        };

        using Element<ScalarT, IdxT>::base_;
        using Element<ScalarT, IdxT>::y_;
        using Element<ScalarT, IdxT>::yp_;
        using Element<ScalarT, IdxT>::f_;
        using Element<ScalarT, IdxT>::state;
        using Element<ScalarT, IdxT>::slice;
        using Element<ScalarT, IdxT>::input;
        using Element<ScalarT, IdxT>::coordinate;

        const SeriesView<ScalarT, IdxT>& series_;
        const ShuntView<ScalarT, IdxT>&  shunt_;
        const Layout                     layout_;
        mutable std::array<SignalT, 4>   inputs_{};
      };
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
