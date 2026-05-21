/**
 * @file BusToSignalAdapterImpl.cpp
 * @author Philip Fackler (facklerpw@ornl.gov)
 *
 * @brief Definition of a BusToSignalAdapter component.
 *
 */

#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/BusToSignalAdapter/BusToSignalAdapter.hpp>
#include <GridKit/Model/PhasorDynamics/BusToSignalAdapter/BusToSignalAdapterData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using Log = ::GridKit::Utilities::Logger;

    /**
     * @brief  Constructor for BusToSignalAdapter component
     *
     * @param bus   Signal used for voltage
     * @param data  Data object
     */
    template <typename ScalarT, typename IdxT>
    BusToSignalAdapter<ScalarT, IdxT>::BusToSignalAdapter(bus_type* bus)
      : bus_(bus)
    {
      size_ = 0;
    }

    template <class ScalarT, typename IdxT>
    BusToSignalAdapter<ScalarT, IdxT>::~BusToSignalAdapter()
    {
    }

    /**
     * @brief Set the component ID
     */
    template <class ScalarT, typename IdxT>
    int BusToSignalAdapter<ScalarT, IdxT>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /**
     * @brief Allocate memory for model
     */
    template <class ScalarT, typename IdxT>
    int BusToSignalAdapter<ScalarT, IdxT>::allocate()
    {
      static constexpr auto VREAL = BusToSignalAdapterInternalVariables::VREAL;
      static constexpr auto VIMAG = BusToSignalAdapterInternalVariables::VIMAG;
      // FIXME: using a placeholder for variable index here is wrong. But it
      //        also seems like nothing ever uses the index from the signal.
      //        So we should either...
      //        - add methods to BusBase to get variable indices (satisfy the
      //          interface correctly)
      //        - remove variable index from SignalNode (because maybe it is
      //          not actually needed)
      if (signals_.template isAssigned<VREAL>())
      {
        signals_.template getSignalNode<VREAL>()->set(&bus_->Vr(), &vr_index_);
      }
      if (signals_.template isAssigned<VIMAG>())
      {
        signals_.template getSignalNode<VIMAG>()->set(&bus_->Vi(), &vi_index_);
      }

      return 0;
    }

    /**
     * @brief verify method checks that attached signals are also linked
     */
    template <class ScalarT, typename IdxT>
    int BusToSignalAdapter<ScalarT, IdxT>::verify() const
    {
      static constexpr auto IREAL = BusToSignalAdapterExternalVariables::IREAL;
      static constexpr auto IIMAG = BusToSignalAdapterExternalVariables::IIMAG;

      int ret = 0;

      if (signals_.template isAttached<IREAL>())
      {
        if (!signals_.template isLinked<IREAL>())
        {
          Log::error() << "BusToSignalAdapter: Ir signal attached with no linked source\n";
          ret += 1;
        }
      }

      if (signals_.template isAttached<IIMAG>())
      {
        if (!signals_.template isLinked<IIMAG>())
        {
          Log::error() << "BusToSignalAdapter: Ii signal attached with no linked source\n";
          ret += 1;
        }
      }

      return ret;
    }

    /**
     * @brief Initialize the adapter
     */
    template <class ScalarT, typename IdxT>
    int BusToSignalAdapter<ScalarT, IdxT>::initialize()
    {
      static constexpr auto IREAL = BusToSignalAdapterExternalVariables::IREAL;
      static constexpr auto IIMAG = BusToSignalAdapterExternalVariables::IIMAG;

      if (signals_.template isAttached<IREAL>())
      {
        bus_->Ir() += signals_.template readExternalVariable<IREAL>();
      }
      if (signals_.template isAttached<IIMAG>())
      {
        bus_->Ii() += signals_.template readExternalVariable<IIMAG>();
      }

      return 0;
    }

    /**
     * @brief No variables to differentiate
     */
    template <class ScalarT, typename IdxT>
    int BusToSignalAdapter<ScalarT, IdxT>::tagDifferentiable()
    {
      return 0;
    }

    /**
     * @brief Residual evaluation
     */
    template <class ScalarT, typename IdxT>
    int BusToSignalAdapter<ScalarT, IdxT>::evaluateResidual()
    {
      static constexpr auto IREAL = BusToSignalAdapterExternalVariables::IREAL;
      static constexpr auto IIMAG = BusToSignalAdapterExternalVariables::IIMAG;

      if (signals_.template isAttached<IREAL>())
      {
        bus_->Ir() += signals_.template readExternalVariable<IREAL>();
      }
      if (signals_.template isAttached<IIMAG>())
      {
        bus_->Ii() += signals_.template readExternalVariable<IIMAG>();
      }

      return 0;
    }

    /**
     * @brief Jacobian evaluation
     */
    template <class ScalarT, typename IdxT>
    int BusToSignalAdapter<ScalarT, IdxT>::evaluateJacobian()
    {
      return 0;
    }

  } // namespace PhasorDynamics
} // namespace GridKit
