#pragma once

#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/PhasorDynamics/GridElement.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using Log = ::GridKit::Utilities::Logger;

    /**
     * @brief Component model implementation base class.
     */
    template <class ScalarT, typename IdxT>
    class Component : public GridElement<ScalarT, IdxT>
    {
    public:
      using RealT   = typename GridElement<ScalarT, IdxT>::RealT;
      using MatrixT = typename GridElement<ScalarT, IdxT>::MatrixT;

      Component() = default;

      virtual ~Component()
      {
      }

      /// @todo Remove this method. It should be part of DynamicSolver class.
      bool hasJacobian() override
      {
        return true;
      }

      void updateTime(RealT t, RealT a) override
      {
        time_  = t;
        alpha_ = a;
      }

      /**
       * @brief Set system frequency and power bases.
       *
       * @param[in] freq_system_base - System frequency base in Hz.
       * @param[in] va_system_base - System power base in VA.
       */
      void setSystemBase(RealT freq_system_base, RealT va_system_base)
      {
        freq_system_base_ = freq_system_base;
        va_system_base_   = va_system_base;
      }

      virtual int setGridKitComponentID(IdxT) = 0;

      IdxT getGridKitComponentID() const
      {
        return gridkit_component_id_;
      }

    protected:
      IdxT gridkit_component_id_{0};

      std::vector<ScalarT> wb_;
      std::vector<ScalarT> h_;

      RealT time_;
      RealT alpha_;

      RealT freq_system_base_{60.0};
      RealT va_system_base_{100.0e6};
    };

  } // namespace PhasorDynamics
} // namespace GridKit
