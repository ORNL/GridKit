/**
 * @file Yc.hpp
 *
 * @brief Characteristic admittance assembly for overhead conductors.
 *
 */

#pragma once

#include <array>
#include <span>
#include <vector>

#include <GridKit/Model/EMT/Element.hpp>
#include <GridKit/Model/EMT/Parameters/Effects/SeriesImpedance/SeriesImpedance.hpp>
#include <GridKit/Model/EMT/Parameters/Effects/ShuntAdmittance/ShuntAdmittance.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename scalar_type, typename index_type>
      class Yc : public Element<scalar_type, index_type>
      {
      public:
        using ScalarT = scalar_type;
        using IdxT    = index_type;
        using SignalT = typename Element<ScalarT, IdxT>::SignalT;

        Yc(const SeriesImpedance<ScalarT, IdxT>& series,
           const ShuntAdmittance<ScalarT, IdxT>& shunt);

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

        SignalT Gc() const
        {
          return {base_ + layout_.gc, layout_.K, layout_.K};
        }

        SignalT Bc() const
        {
          return {base_ + layout_.bc, layout_.K, layout_.K};
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
              gc(0),
              bc(KK),
              n(2 * KK)
          {
          }

          const IdxT K;
          const IdxT KK;
          const IdxT gc;
          const IdxT bc;
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

        const SeriesImpedance<ScalarT, IdxT>& series_;
        const ShuntAdmittance<ScalarT, IdxT>& shunt_;
        const Layout                          layout_;
        mutable std::array<SignalT, 4>        inputs_{};
      };
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
