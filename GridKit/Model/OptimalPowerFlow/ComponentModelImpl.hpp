/**
 * @file ComponentModelImpl.hpp
 * @brief Definition of `ComponentModel`, included by model translation units.
 */

#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <type_traits>
#include <utility>

#include <GridKit/Model/OptimalPowerFlow/ComponentModel.hpp>

namespace GridKit
{
  namespace Enzyme
  {
    namespace Sparse
    {
      // Defined in the Enzyme drivers, which only Enzyme translation units include
      template <typename ModelT>
      struct DfDx;

      template <typename ModelT>
      struct DgDx;

      template <typename ModelT>
      struct D2LDx2;
    } // namespace Sparse
  } // namespace Enzyme

  namespace OptimalPowerFlow
  {
    /**
     * @brief Enzyme derivatives exist for double scalars of models that
     * depend on the variables
     */
    template <typename model_type, typename scalar_type, typename index_type>
    constexpr bool ComponentModel<model_type, scalar_type, index_type>::hasEnzymeDerivatives()
    {
      return std::is_same_v<ScalarT, double> && !model_type::CONSTANT;
    }

    template <typename model_type, typename scalar_type, typename index_type>
    const model_type& ComponentModel<model_type, scalar_type, index_type>::model() const
    {
      return static_cast<const model_type&>(*this);
    }

    /**
     * @brief Local variables gathered from the system variables
     */
    template <typename model_type, typename scalar_type, typename index_type>
    auto ComponentModel<model_type, scalar_type, index_type>::gather(const ScalarT* x) const
    {
      std::array<ScalarT, model_type::VARIABLE_SIZE> local{};
      for (IdxT j = 0; j < model_type::VARIABLE_SIZE; ++j)
      {
        local[j] = x[this->variable_indices_[j]];
      }
      return local;
    }

    template <typename model_type, typename scalar_type, typename index_type>
    ComponentModel<model_type, scalar_type, index_type>::ComponentModel(std::string       id,
                                                                        std::vector<IdxT> terminals)
      : Component<scalar_type, index_type>(std::move(id),
                                           std::move(terminals),
                                           model_type::INTERNAL_SIZE,
                                           model_type::INTERNAL_CONSTRAINT_SIZE)
    {
      assert(this->size() == model_type::VARIABLE_SIZE);
      assert(this->sizeConstraints() == model_type::CONSTRAINT_SIZE);
    }

    template <typename model_type, typename scalar_type, typename index_type>
    int ComponentModel<model_type, scalar_type, index_type>::evaluateObjective(const ScalarT* x, ScalarT& f)
    {
      const auto local  = gather(x);
      f                += model().objective(local.data());
      return 0;
    }

    template <typename model_type, typename scalar_type, typename index_type>
    int ComponentModel<model_type, scalar_type, index_type>::evaluateGradient(const ScalarT* x, RealT* gradient)
    {
      if constexpr (hasEnzymeDerivatives())
      {
        const auto                                   local = gather(x);
        std::array<RealT, model_type::VARIABLE_SIZE> local_gradient{};
        Enzyme::Sparse::DfDx<model_type>::eval(&model(), local.data(), local_gradient.data());

        for (IdxT j = 0; j < model_type::VARIABLE_SIZE; ++j)
        {
          gradient[this->variable_indices_[j]] += local_gradient[j];
        }
      }
      return 0;
    }

    template <typename model_type, typename scalar_type, typename index_type>
    int ComponentModel<model_type, scalar_type, index_type>::evaluateConstraints(const ScalarT* x, ScalarT* g)
    {
      const auto                                       local = gather(x);
      std::array<ScalarT, model_type::CONSTRAINT_SIZE> local_g{};
      model().constraints(local.data(), local_g.data());

      for (IdxT i = 0; i < model_type::CONSTRAINT_SIZE; ++i)
      {
        const IdxT row = this->constraint_indices_[i];
        if (row != INVALID_INDEX<IdxT>)
        {
          g[row] += local_g[i];
        }
      }
      return 0;
    }

    template <typename model_type, typename scalar_type, typename index_type>
    int ComponentModel<model_type, scalar_type, index_type>::evaluateJacobian(const ScalarT* x)
    {
      this->jacobian_.clear();
      if constexpr (hasEnzymeDerivatives())
      {
        const auto local = gather(x);
        Enzyme::Sparse::DgDx<model_type>::eval(&model(),
                                               this->variable_indices_.size(),
                                               this->constraint_indices_.data(),
                                               this->variable_indices_.data(),
                                               local.data(),
                                               &this->jacobian_);
      }
      return 0;
    }

    template <typename model_type, typename scalar_type, typename index_type>
    int ComponentModel<model_type, scalar_type, index_type>::evaluateHessian(const ScalarT* x,
                                                                             RealT          sigma,
                                                                             const RealT*   lambda)
    {
      this->hessian_.clear();
      if constexpr (hasEnzymeDerivatives())
      {
        const auto                                     local = gather(x);
        std::array<RealT, model_type::CONSTRAINT_SIZE> multipliers{};
        for (IdxT i = 0; i < model_type::CONSTRAINT_SIZE; ++i)
        {
          const IdxT row = this->constraint_indices_[i];
          if (row != INVALID_INDEX<IdxT>)
          {
            multipliers[i] = lambda[row];
          }
        }

        Enzyme::Sparse::D2LDx2<model_type>::eval(&model(),
                                                 this->variable_indices_.size(),
                                                 this->variable_indices_.data(),
                                                 local.data(),
                                                 sigma,
                                                 multipliers.data(),
                                                 &this->hessian_);
      }
      return 0;
    }

    template <typename model_type, typename scalar_type, typename index_type>
    int ComponentModel<model_type, scalar_type, index_type>::evaluateTerminalPower(const ScalarT* x, ScalarT* power)
    {
      const auto                                       local = gather(x);
      std::array<ScalarT, model_type::CONSTRAINT_SIZE> local_g{};
      model().constraints(local.data(), local_g.data());

      std::copy(local_g.begin() + model_type::INTERNAL_CONSTRAINT_SIZE, local_g.end(), power);
      return 0;
    }
  } // namespace OptimalPowerFlow
} // namespace GridKit
