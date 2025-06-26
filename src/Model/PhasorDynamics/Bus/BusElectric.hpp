#pragma once

#include <vector>

#include <Model/PhasorDynamics/BusBase.hpp>

// Forward declaration of BusData structure
namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename RealT, typename IdxT>
    struct BusElectricData;
  }
} // namespace GridKit

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

    public:
    
      using real_type = typename BusBase<ScalarT, IdxT>::real_type;
      using DataT     = BusElectricData<real_type, IdxT>;

      BusElectric()
      {
      }

      BusElectric(ScalarT Vr, ScalarT Vi)
      {
      }

      BusElectric(const DataT& data)
      {
      }

      virtual ~BusElectric()
      {
      }

      /// Pure virtual function, returns bus type (DEFAULT or SLACK).
      virtual int BusType() const = 0;

      virtual ScalarT& Vr() = 0;
      virtual ScalarT& Vi() = 0;
      virtual ScalarT& Ir() = 0; 
      virtual ScalarT& Ii() = 0;
      
      virtual const ScalarT& Vr() const = 0;
      virtual const ScalarT& Vi() const = 0;
      virtual const ScalarT& Ir() const = 0;
      virtual const ScalarT& Ii() const = 0;

    private:
      ScalarT Vr0_{0.0};
      ScalarT Vi0_{0.0};

    };

  } // namespace PhasorDynamics
} // namespace GridKit
