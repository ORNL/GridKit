#pragma once

#include <GridKit/Constants.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename RealT, typename IdxT>
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
    template <class ScalarT, typename IdxT>
    class SignalNode
    {
    public:
      using RealT = typename GridKit::ScalarTraits<ScalarT>::RealT;

      SignalNode();
      SignalNode(const SignalNodeData<RealT, IdxT>& data);

      virtual ~SignalNode() = default;

      void    set(ScalarT* signal_in, IdxT* global_index);
      bool    linked() const;
      ScalarT read() const;
      void    init(ScalarT signal_in);

      const IdxT signalId() const
      {
        return signal_id_;
      }

      IdxT getVariableIndex() const
      {
        return *variable_index_;
      }

      virtual const IdxT busID() const
      {
        return bus_id_;
      }

    private:
      ScalarT* signal_{nullptr};
      IdxT     signal_id_{0};

    protected:
      const IdxT bus_id_{INVALID_INDEX<IdxT>};

      IdxT* variable_index_{nullptr};
    };

  } // namespace PhasorDynamics
} // namespace GridKit
