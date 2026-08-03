/**
 * @file Tau.hpp
 *
 * @brief Modal phase-delay assembly from Gamma modal phase constants.
 *
 */

#pragma once

#include <array>
#include <span>
#include <vector>

#include <GridKit/Model/EMT/Element.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Path/Path.hpp>
#include <GridKit/Model/EMT/Parameters/Response/Gamma/Gamma.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename scalar_type, typename index_type>
      class Tau : public Element<scalar_type, index_type>
      {
      public:
        using ScalarT = scalar_type;
        using IdxT    = index_type;
        using SignalT = typename Element<ScalarT, IdxT>::SignalT;

        Tau(const Gamma<ScalarT, IdxT>& gamma, const Path<ScalarT, IdxT>& path);

        IdxT size() const override
        {
          return layout_.n;
        }

        std::span<const SignalT> inputs() const override
        {
          inputs_[0] = gamma_.betaM();
          return inputs_;
        }

        SignalT tau() const
        {
          return {base_ + layout_.tau, layout_.K, 1};
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
          explicit Layout(IdxT mode_count)
            : K(mode_count),
              tau(0),
              n(K)
          {
          }

          const IdxT K;
          const IdxT tau;
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

        const Gamma<ScalarT, IdxT>&    gamma_;
        const Path<ScalarT, IdxT>&     path_;
        const Layout                   layout_;
        mutable std::array<SignalT, 1> inputs_{};
      };
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
