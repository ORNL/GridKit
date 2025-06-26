#pragma once

#include <vector>


#include <Model/PhasorDynamics/Bus/BusControl/BusControlData.hpp>
#include <Model/PhasorDynamics/BusBase.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /*!
     * @brief BusControl model implementation base class.
     *
     */
    template <class ScalarT, typename IdxT>
    class BusControl : public BusBase<ScalarT, IdxT>
    {

    public:
      using real_type = typename BusBase<ScalarT, IdxT>::real_type;
      using DataT     = BusControlData<real_type, IdxT>;

      // Use all inherited constructors
      using BusBase<ScalarT, IdxT>::BusBase;

      BusControl()
      {
      }

      BusControl(const DataT& data)
      {
      }

      virtual ~BusControl()
      {
      }

      // Pure virtual function, returns bus type (DEFAULT or SLACK).
      virtual int BusType() const = 0;

      // Assign source of signal via function pointer
      virtual void set_source(ScalarT (*callback)()) = 0;

      // Accessor function of signal
      virtual ScalarT poll() = 0;

    private:
     
      // Private callback to source
      ScalarT (*callback_)();

    };

  } // namespace PhasorDynamics
} // namespace GridKit
