#pragma once

#include <vector>

#include <Model/PhasorDynamics/Bus/BusBase.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /*!
     * @brief BusElectric model implementation base class.
     *
     */
    template <class ScalarT, typename IdxT>
    class BusElectric: public BusBase<ScalarT, IdxT>
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

      BusElectric();
      BusElectric(ScalarT Vr, ScalarT Vi);
      BusElectric(const DataT& data);
      virtual ~BusElectric();

      /// Pure virtual function, returns bus type (DEFAULT or SLACK).
      virtual int BusType() const = 0;

      // Voltage and Current Accessors
      virtual ScalarT&       Vr()       = 0;
      virtual const ScalarT& Vr() const = 0;
      virtual ScalarT&       Vi()       = 0;
      virtual const ScalarT& Vi() const = 0;
      virtual ScalarT&       Ir()       = 0;
      virtual const ScalarT& Ir() const = 0;
      virtual ScalarT&       Ii()       = 0;
      virtual const ScalarT& Ii() const = 0;


      virtual int BusType() const override
      {
        return BusData<real_type, IdxT>::BusType::DEFAULT;
      }

      virtual ScalarT& Vr() override
      {
        return y_[0];
      }

      virtual const ScalarT& Vr() const override
      {
        return y_[0];
      }

      virtual ScalarT& Vi() override
      {
        return y_[1];
      }

      virtual const ScalarT& Vi() const override
      {
        return y_[1];
      }

      virtual ScalarT& Ir() override
      {
        return f_[0];
      }

      virtual const ScalarT& Ir() const override
      {
        return f_[0];
      }

      virtual ScalarT& Ii() override
      {
        return f_[1];
      }

      virtual const ScalarT& Ii() const override
      {
        return f_[1];
      }

    private:
      ScalarT Vr0_{0.0};
      ScalarT Vi0_{0.0};
    };

  } // namespace PhasorDynamics
} // namespace GridKit
