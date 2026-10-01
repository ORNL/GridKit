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
#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeSet.hpp>
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
    template <typename scalar_type, typename index_type>
    BusToSignalAdapter<scalar_type, index_type>::BusToSignalAdapter(BusT* bus)
      : bus_(bus)
    {
      size_ = 0;
    }

    template <typename scalar_type, typename index_type>
    BusToSignalAdapter<scalar_type, index_type>::BusToSignalAdapter(
        BusT* bus, const ModelDataT&)
      : bus_(bus)
    {
      size_ = 0;
    }

    template <typename scalar_type, typename index_type>
    BusToSignalAdapter<scalar_type, index_type>::~BusToSignalAdapter()
    {
    }

    /**
     * @brief Set the component ID
     */
    template <typename scalar_type, typename index_type>
    int BusToSignalAdapter<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /**
     * @brief Allocate memory for model
     */
    template <typename scalar_type, typename index_type>
    int BusToSignalAdapter<scalar_type, index_type>::allocate()
    {
      using SignalOut = BusToSignalAdapterSignalOutputs;

      if (!allocated_)
      {
        this->allocateVectors(size_);
      }
      auto size = static_cast<std::size_t>(size_);

      variable_indices_.resize(size);
      residual_indices_.resize(size);

      // vr_index_ and vi_index_ are both set to INVALID_INDEX. This component
      // simply passes voltage from bus to output signal, so indices are
      // ignored.
      if (auto vr_port = ports_.out.template port<SignalOut::vr>())
      {
        vr_port.link(&bus_->Vr(), &vr_index_);
      }
      if (auto vi_port = ports_.out.template port<SignalOut::vi>())
      {
        vi_port.link(&bus_->Vi(), &vi_index_);
      }

      allocated_ = true;
      return 0;
    }

    /**
     * @brief verify method checks that attached signals are also linked
     */
    template <typename scalar_type, typename index_type>
    int BusToSignalAdapter<scalar_type, index_type>::verify() const
    {
      using SignalIn = BusToSignalAdapterSignalInputs;

      int ret = 0;

      auto ir_port = ports_.in.template port<SignalIn::ir>();
      if (ir_port.connected() && !ir_port.linked())
      {
        Log::error() << "BusToSignalAdapter: Ir signal attached with no linked source\n";
        ret += 1;
      }

      auto ii_port = ports_.in.template port<SignalIn::ii>();
      if (ii_port.connected() && !ii_port.linked())
      {
        Log::error() << "BusToSignalAdapter: Ii signal attached with no linked source\n";
        ret += 1;
      }

      return ret;
    }

    /**
     * @brief Initialize the adapter
     */
    template <typename scalar_type, typename index_type>
    int BusToSignalAdapter<scalar_type, index_type>::initialize()
    {
      return 0;
    }

    /**
     * @brief No variables to differentiate
     */
    template <typename scalar_type, typename index_type>
    int BusToSignalAdapter<scalar_type, index_type>::tagDifferentiable()
    {
      return 0;
    }

    /**
     * @brief Compute the absolute tolerance for each variable in the model
     *
     * @param rel_tol The relative tolerance which can be used to pick the
     *        absolute tolerance.
     * @tparam scalar_type Scalar data type
     * @tparam index_type Index data type
     * @return int 0 if successful, non-zero otherwise.
     *
     * This represents a "noise" level close to zero for which pure relative
     * error cannot be used.
     */
    template <typename scalar_type, typename index_type>
    int BusToSignalAdapter<scalar_type, index_type>::setAbsoluteTolerance(RealT)
    {
      return 0;
    }

    /**
     * @brief Residual evaluation
     */
    template <typename scalar_type, typename index_type>
    int BusToSignalAdapter<scalar_type, index_type>::evaluateResidual()
    {
      using SignalIn = BusToSignalAdapterSignalInputs;

      if (auto ir_port = ports_.in.template port<SignalIn::ir>())
      {
        bus_->Ir() += ir_port.readSignal();
        if (bus_->size() > 0)
        {
          bus_->getResidual().setDataUpdated();
        }
      }
      if (auto ii_port = ports_.in.template port<SignalIn::ii>())
      {
        bus_->Ii() += ii_port.readSignal();
        if (bus_->size() > 0)
        {
          bus_->getResidual().setDataUpdated();
        }
      }

      return 0;
    }

    /**
     * @brief Jacobian evaluation
     */
    template <typename scalar_type, typename index_type>
    int BusToSignalAdapter<scalar_type, index_type>::evaluateJacobian()
    {
      return 0;
    }

  } // namespace PhasorDynamics
} // namespace GridKit
