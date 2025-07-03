
#pragma once

#include <Model/PhasorDynamics/BusBase.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    
    // TODO change base class so there is no Vr and Vi

    template <class ScalarT, typename IdxT>
    class BusSignal : public BusBase<ScalarT, IdxT>
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
      using DataT     = BusData<real_type, IdxT>;

      BusSignal(const DataT& data);
      
      ~BusSignal()
      {
      }

      virtual int allocate() override;
      virtual int tagDifferentiable() override;
      virtual int initialize() override;
      virtual int evaluateResidual() override;
      virtual int evaluateIntegrand() override;
      virtual int evaluateJacobian() override;
      virtual int initializeAdjoint() override;
      virtual int evaluateAdjointIntegrand() override;
      virtual int evaluateAdjointResidual() override;

      
       /**
       * @brief A one-time initialization function
       * that can be called by any.
       * 
       * @return state value of signal
       */
      void initial_value(ScalarT value) override
      {
        if(is_initialized_)
        {
          std::cout << "ERROR!";
        }
        else{
          y_[0] = value;
          is_initialized_ = true;
        }
      }

       /**
       * @brief A read only accessor to the signal state
       * 
       * @return state value of signal
       */
      ScalarT& read() override
      {
        return y_[0];
      }

       /**
       * @brief A write operation only available
       * to the source signal
       * 
       * @param value The value of the signal
       */
      void send(ScalarT& value) override
      {
        f_[0] += value;
      }


      // TODO remove jsut junk so it compiles
      ScalarT& Vr() override
      {
        return y_[0];
      }
      const ScalarT& Vr() const override
      {
        return y_[0];
      }
      ScalarT& Vi() override
      {
        return y_[0];
      }

      const ScalarT& Vi() const override
      {
        return y_[0];
      }

      ScalarT& Ir() override
      {
        return y_[0];
      }

      const ScalarT& Ir() const override
      {
        return y_[0];
      }

      ScalarT& Ii() override
      {
        return y_[0];
      }

      const ScalarT& Ii() const override
      {
        return y_[0];
      }
    
    private:

      bool is_initialized_{false};

    };

  } // namespace PhasorDynamics
} // namespace GridKit