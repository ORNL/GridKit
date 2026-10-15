/**
 * @file Load.hpp
 * @brief Declaration of an optimal power flow load.
 */

#pragma once

#include <GridKit/Model/OptimalPowerFlow/ComponentModel.hpp>
#include <GridKit/Model/OptimalPowerFlow/Load/LoadData.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /**
     * @brief Constant power demand
     */
    template <typename scalar_type, typename index_type>
    class Load : public ComponentModel<Load<scalar_type, index_type>, scalar_type, index_type>
    {
      using ComponentModelT = ComponentModel<Load<scalar_type, index_type>, scalar_type, index_type>;

      using ComponentModelT::in_service_;

    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename ComponentModelT::RealT;
      using ModelDataT = LoadData<RealT, IdxT>;

      /// Local variables: terminal real and imaginary voltage
      static constexpr IdxT VARIABLE_SIZE            = 2;
      /// Local constraint functions: power into the bus
      static constexpr IdxT CONSTRAINT_SIZE          = 2;
      /// Internal variables: none
      static constexpr IdxT INTERNAL_VARIABLE_SIZE   = 0;
      /// Internal constraint functions: none
      static constexpr IdxT INTERNAL_CONSTRAINT_SIZE = 0;
      static constexpr bool CONSTANT                 = true;

      explicit Load(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        initialize(ScalarT* x) override;

      __attribute__((always_inline)) inline ScalarT objective(const ScalarT* x) const;
      __attribute__((always_inline)) inline void    constraints(const ScalarT* x, ScalarT* g) const;

    private:
      RealT P_{0.0};
      RealT Q_{0.0};
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
