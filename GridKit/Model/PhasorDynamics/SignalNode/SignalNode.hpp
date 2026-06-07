#pragma once

#include <functional>
#include <type_traits>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalHistory.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename real_type, typename index_type>
    struct SignalNodeData;
  } // namespace PhasorDynamics
} // namespace GridKit

namespace GridKit
{
  namespace PhasorDynamics
  {
    /*!
     * @brief SignalNode model implementation base class.
     *
     */
    template <typename scalar_type, typename index_type>
    class SignalNode
    {
    public:
      using ScalarT      = scalar_type;
      using IdxT         = index_type;
      using RealT        = typename GridKit::ScalarTraits<ScalarT>::RealT;
      using PrehistoryFn = std::function<RealT(RealT)>;

      SignalNode();
      SignalNode(const SignalNodeData<RealT, IdxT>& data);

      virtual ~SignalNode() = default;

      void    set(ScalarT* signal_in, IdxT* global_index, PrehistoryFn prehistory = {});
      bool    linked() const;
      ScalarT read() const;
      ScalarT readWithDelay(RealT delay) const;
      void    init(ScalarT signal_in);
      void    requireHistoryWindow(RealT delay);
      void    setReadTime(RealT t);
      void    resetHistory(RealT t0);
      int     stepAccepted(RealT t);

      const IdxT signalId() const
      {
        return signal_id_;
      }

      IdxT getVariableIndex() const
      {
        return variableIndexForDelay(0.0);
      }

      IdxT variableIndexForDelay(RealT delay) const
      {
        if (delay > 0.0 || variable_index_ == nullptr)
        {
          return INVALID_INDEX<IdxT>;
        }
        return *variable_index_;
      }

      RealT historyWindow() const
      {
        return history_window_;
      }

      virtual const IdxT busID() const
      {
        return bus_id_;
      }

    private:
      RealT realValue() const;

      ScalarT*             signal_{nullptr};
      IdxT                 signal_id_{0};
      SignalHistory<RealT> history_;
      PrehistoryFn         prehistory_;
      RealT                history_window_{0.0};
      RealT                initial_value_{0.0};
      RealT                read_time_{0.0};

    protected:
      const IdxT bus_id_{INVALID_INDEX<IdxT>};

      IdxT* variable_index_{nullptr};
    };

  } // namespace PhasorDynamics
} // namespace GridKit
