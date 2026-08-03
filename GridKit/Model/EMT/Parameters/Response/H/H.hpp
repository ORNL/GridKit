/**
 * @file H.hpp
 *
 * @brief Modal propagation function from Gamma modal propagation constants.
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
      class H : public Element<scalar_type, index_type>
      {
      public:
        using ScalarT = scalar_type;
        using IdxT    = index_type;
        using SignalT = typename Element<ScalarT, IdxT>::SignalT;

        H(const Gamma<ScalarT, IdxT>& gamma, const Path<ScalarT, IdxT>& path);

        IdxT size() const override
        {
          return layout_.n;
        }

        std::span<const SignalT> inputs() const override
        {
          inputs_[0] = gamma_.alphaMDot();
          inputs_[1] = gamma_.betaMDot();
          return inputs_;
        }

        SignalT Hr() const
        {
          return {base_ + layout_.hr, layout_.K, 1};
        }

        SignalT Hi() const
        {
          return {base_ + layout_.hi, layout_.K, 1};
        }

        int initialize() override;

        template <typename V>
        void evaluateResidualKernel(V* y, V* yp, V* u, V* f) const;

        void evaluateInternalResidual(ScalarT* y, ScalarT* yp, ScalarT* u, ScalarT* f) const;

        void tagDifferentiable(std::vector<bool>& tag) const override
        {
          for (IdxT i = 0; i < size(); ++i)
          {
            tag[static_cast<size_t>(base_ + i)] = true;
          }
        }

        int evaluateResidual() override;
        int evaluateJacobian() override;

      private:
        struct Layout
        {
          explicit Layout(IdxT mode_count)
            : K(mode_count),
              hr(0),
              hi(K),
              n(2 * K)
          {
          }

          const IdxT K;
          const IdxT hr;
          const IdxT hi;
          const IdxT n;
        };

        using Element<ScalarT, IdxT>::base_;
        using Element<ScalarT, IdxT>::y_;
        using Element<ScalarT, IdxT>::yp_;
        using Element<ScalarT, IdxT>::f_;
        using Element<ScalarT, IdxT>::state;
        using Element<ScalarT, IdxT>::derivative;
        using Element<ScalarT, IdxT>::slice;
        using Element<ScalarT, IdxT>::input;

        const Gamma<ScalarT, IdxT>&    gamma_;
        const Path<ScalarT, IdxT>&     path_;
        const Layout                   layout_;
        mutable std::array<SignalT, 2> inputs_{};
      };
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
