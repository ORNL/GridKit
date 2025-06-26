#pragma once

#include <vector>

#include <Model/PhasorDynamics/Bus/BusBase.hpp>

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

      using BusBase<ScalarT, IdxT>::size_;
      using BusBase<ScalarT, IdxT>::y_;
      using BusBase<ScalarT, IdxT>::yp_;
      using BusBase<ScalarT, IdxT>::yB_;
      using BusBase<ScalarT, IdxT>::ypB_;
      using BusBase<ScalarT, IdxT>::f_;
      using BusBase<ScalarT, IdxT>::fB_;
      using BusBase<ScalarT, IdxT>::tag_;

    public:
      using real_type = typename BusBase<ScalarT, IdxT>::real_type;

      BusControl();
      BusControl(const DataT& data);
      virtual ~BusControl();

      /// Pure virtual function, returns bus type (DEFAULT or SLACK).
      virtual int BusType() const = 0;

    private:
      ScalarT channel;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
