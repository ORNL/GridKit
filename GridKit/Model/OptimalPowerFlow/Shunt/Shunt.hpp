/**
 * @file Shunt.hpp
 * @brief Declaration of an optimal power flow shunt.
 */

#pragma once

#include <GridKit/Model/OptimalPowerFlow/ComponentModel.hpp>
#include <GridKit/Model/OptimalPowerFlow/Shunt/ShuntData.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /**
     * @brief Constant admittance \f$Y = G + jB\f$ to ground
     *
     * Power into the bus is \f$-Y^* (V_\mathrm{r}^2 + V_\mathrm{i}^2)\f$.
     */
    template <typename scalar_type, typename index_type>
    class Shunt : public ComponentModel<Shunt<scalar_type, index_type>, scalar_type, index_type>
    {
      using ComponentModelT = ComponentModel<Shunt<scalar_type, index_type>, scalar_type, index_type>;

      using ComponentModelT::in_service_;

    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename ComponentModelT::RealT;
      using ModelDataT = ShuntData<RealT, IdxT>;

      /// Local variables: terminal real and imaginary voltage
      static constexpr IdxT VARIABLE_SIZE            = 2;
      /// Local constraint functions: power into the bus
      static constexpr IdxT CONSTRAINT_SIZE          = 2;
      /// Internal variables: none
      static constexpr IdxT INTERNAL_VARIABLE_SIZE   = 0;
      /// Internal constraint functions: none
      static constexpr IdxT INTERNAL_CONSTRAINT_SIZE = 0;
      static constexpr bool CONSTANT                 = false;

      explicit Shunt(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        initialize(ScalarT* x) override;

      __attribute__((always_inline)) inline ScalarT objective(const ScalarT* x) const;
      __attribute__((always_inline)) inline void    constraints(const ScalarT* x, ScalarT* g) const;

    private:
      RealT G_{0.0};
      RealT B_{0.0};
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
