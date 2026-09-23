/**
 * @file Generator.hpp
 * @brief Declaration of an optimal power flow generator.
 */

#pragma once

#include <GridKit/Model/OptimalPowerFlow/ComponentModel.hpp>
#include <GridKit/Model/OptimalPowerFlow/Generator/GeneratorData.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /**
     * @brief Dispatchable active and reactive power injection with a
     * quadratic cost
     */
    template <typename scalar_type, typename index_type>
    class Generator : public ComponentModel<Generator<scalar_type, index_type>, scalar_type, index_type>
    {
      using ComponentModelT = ComponentModel<Generator<scalar_type, index_type>, scalar_type, index_type>;

      using ComponentModelT::id_;
      using ComponentModelT::variable_indices_;

    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename ComponentModelT::RealT;
      using ModelDataT = GeneratorData<RealT, IdxT>;

      /// Local variables: active and reactive power, then terminal real and imaginary voltage
      static constexpr IdxT VARIABLE_SIZE            = 4;
      /// Local constraints: power into the bus
      static constexpr IdxT CONSTRAINT_SIZE          = 2;
      static constexpr IdxT INTERNAL_SIZE            = 2;
      static constexpr IdxT INTERNAL_CONSTRAINT_SIZE = 0;
      static constexpr bool CONSTANT                 = false;

      explicit Generator(const ModelDataT& data);

      int verify() const override;
      int initialize(const Model::StateData& state, ScalarT* x, RealT* x_lower, RealT* x_upper) override;

      __attribute__((always_inline)) inline ScalarT objective(const ScalarT* x) const;
      __attribute__((always_inline)) inline void    constraints(const ScalarT* x, ScalarT* g) const;

    private:
      RealT Pmin_{0.0};
      RealT Pmax_{0.0};
      RealT Qmin_{0.0};
      RealT Qmax_{0.0};
      RealT c0_{0.0};
      RealT c1_{0.0};
      RealT c2_{0.0};
      RealT online_{1.0};
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
