#include <cassert>
#include <iostream>

#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusFactory.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusSignalVoltageIn/BusSignalVoltageIn.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusSignalVoltageOut/BusSignalVoltageOut.hpp>
#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Model/VariableMonitorController.hpp>

// Include all components
#include <GridKit/Model/PhasorDynamics/ComponentLibrary.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Constructor for the system model
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::SystemModel()
      : monitor_(std::make_unique<MonitorT>())
    {
    }

    /**
     * @brief Construct a new System Model object
     *
     * @param[in] data - Data structure with complete system data
     *
     * @pre SystemModelData contains consistent connectivity information
     * and physically meaningful model parameters.
     *
     * @post All component models in SystemModelData are created, and
     * correctly connected into the system model.
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::SystemModel(const SystemModelData<RealT, IdxT>& data)
      : monitor_(std::make_unique<MonitorT>(time_))
    {
      using namespace Governor;
      using namespace Exciter;
      using namespace Stabilizer;
      using namespace Controller;
      using namespace Converter;

      owns_components_ = true;

      // Store parsed system bases before constructing data-driven components.
      this->setSystemBase(data.freq_base, data.va_base);

      // Add signal nodes
      for (const auto& signaldata : data.signal)
      {
        signal_nodes_.add(signaldata);
      }

      // Add electrical buses
      for (const auto& busdata : data.bus)
      {
        BusBase<ScalarT, IdxT>* bus = BusFactory<ScalarT, IdxT>::create(busdata, signal_nodes_);
        addBus(bus);
      }

      // Components initialize in insertion order, so each producer precedes
      // the consumers that read it during initialization.
      addDevices<Regca<ScalarT, IdxT>>(data.regca);
      addDevices<Branch<ScalarT, IdxT>>(data.branch);
      addDevices<LoadZ<ScalarT, IdxT>>(data.loadz);
      addDevices<LoadZIP<ScalarT, IdxT>>(data.loadzip);
      addDevices<Genrou<ScalarT, IdxT>>(data.genrou);
      addDevices<Gensal<ScalarT, IdxT>>(data.gensal);
      addDevices<GenClassical<ScalarT, IdxT>>(data.genclassical);
      // REECB follows its current-command and feedback producers.
      addDevices<Reecb<ScalarT, IdxT>>(data.reecb);
      addDevices<Tgov1<ScalarT, IdxT>>(data.gov);
      addDevices<GastPti<ScalarT, IdxT>>(data.gastpti);
      addDevices<Hygov<ScalarT, IdxT>>(data.hygov);
      // IEEEST precedes the exciters that read its output.
      addDevices<Ieeest<ScalarT, IdxT>>(data.stabilizer);
      addDevices<Ieeet1<ScalarT, IdxT>>(data.exciter);
      addDevices<Esdc1a<ScalarT, IdxT>>(data.esdc1a);
      addDevices<SexsPti<ScalarT, IdxT>>(data.sexspti);
      // REPCA follows the signal producers it reads.
      addDevices<Repca<ScalarT, IdxT>>(data.repca);
      addDevices<ConstantSignalSource<ScalarT, IdxT>>(data.constant_source);
      addDevices<FunctionSignalSource<ScalarT, IdxT>>(data.function_source);
      addDevices<BusFault<ScalarT, IdxT>>(data.bus_fault);

      for (const auto& sink : data.monitor_sink)
      {
        monitor_->addSink(sink);
      }
    }

    /**
     * @brief Destructor for the system model
     *
     * If the SystemModel owns the components, it needs to delete them upon
     * destructor call.
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::~SystemModel()
    {
      if (owns_components_)
      {
        for (auto component : components_)
        {
          delete component;
        }

        for (auto bus : buses_)
        {
          delete bus;
        }
      }
    }

    /**
     * @brief Set component ID
     *
     * @note Should default to 0. Nested system models are not currently
     * supported.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /**
     * @brief Allocate system storage and bind buses and components to it.
     *
     * First sum the bus and component sizes, then allocate the system vectors
     * and bind each bus and component to its portion of those vectors.
     *
     * @pre Buses and components with nonzero size are unallocated or already
     * bound to external system storage.
     *
     * @note System model composition is flat; nested systems are not supported.
     *
     * @throws std::runtime_error if storage allocation, child binding, model
     * verification, or initialization for sparse Jacobian discovery fails.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::allocate()
    {
      size_ = 0;

      for (const auto& bus : buses_)
      {
        size_ += bus->size();
      }

      for (const auto& component : components_)
      {
        size_ += component->size();
      }

      // Allocate global vectors
      if (!allocated_)
      {
        // Topology changes invalidate the Jacobian sparsity pattern and COO-to-CSR map.
        delete csr_jac_;
        csr_jac_ = nullptr;

        delete[] map_to_csr_;
        map_to_csr_ = nullptr;

        nnz_ = 0;
        this->allocateVectors(size_);
      }

      if (y_.getSize() != size_
          || yp_.getSize() != size_
          || f_.getSize() != size_
          || abs_tol_.getSize() != size_)
      {
        Log::error() << "SystemModel vector sizes do not match the system size\n";
        throw std::runtime_error("SystemModel allocation failed");
      }

      tag_.resize(size_);
      variable_indices_.resize(size_);
      residual_indices_.resize(size_);

      // Default variable and residual index mapping to local index
      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, j);
        this->setResidualIndex(j, j);
      }

      IdxT offset = 0;

      for (const auto& bus : buses_)
      {
        const int bind_status = bus->bind(y_, yp_, f_, abs_tol_, offset);
        if (bind_status != 0)
        {
          Log::error() << "Failed to bind bus vectors to system storage\n";
          throw std::runtime_error("SystemModel allocation failed");
        }

        if (bus->allocate() != 0)
        {
          Log::error() << "Failed to allocate bus data\n";
          throw std::runtime_error("SystemModel allocation failed");
        }

        for (IdxT j = 0; j < bus->size(); ++j)
        {
          bus->setVariableIndex(j, offset + j);
          bus->setResidualIndex(j, offset + j);
        }

        offset += bus->size();
      }

      for (const auto& component : components_)
      {
        const int bind_status = component->bind(y_, yp_, f_, abs_tol_, offset);
        if (bind_status != 0)
        {
          Log::error() << "Failed to bind component vectors to system storage\n";
          throw std::runtime_error("SystemModel allocation failed");
        }

        if (component->allocate() != 0)
        {
          Log::error() << "Failed to allocate component data\n";
          throw std::runtime_error("SystemModel allocation failed");
        }

        for (IdxT j = 0; j < component->size(); ++j)
        {
          component->setVariableIndex(j, offset + j);
          component->setResidualIndex(j, offset + j);
        }

        offset += component->size();
      }

      if (offset != size_)
      {
        Log::error() << "Bound vector sizes do not match the system size\n";
        throw std::runtime_error("SystemModel allocation failed");
      }

      // Verify component configuration
      int errorCount = this->verify();
      if (errorCount > 0)
      {
        Log::error() << "Component errors: " << errorCount << std::endl;
        throw std::runtime_error("SystemModel allocation failed");
      }

      // Sparse-pattern discovery requires an initialized operating point. A failed
      // initialization aborts allocation before residual/Jacobian evaluation or
      // monitor startup.
      // @todo Replace with a sparsity analysis that sets the NNZ and allocates
      // the Jacobian without needing the Jacobian values.
      if (hasJacobian())
      {
        const int status = initialize();
        if (status != 0)
        {
          Log::error() << "System model initialization failed with status "
                       << status << '\n';
          throw std::runtime_error("SystemModel allocation failed");
        }
        evaluateResidual();
        evaluateJacobian();
      }

      initializeMonitor();
      startMonitor();

      allocated_ = true;
      return 0;
    }

    /**
     * @brief Verify all components are configured correctly
     *
     * This method accumulates and returns the number of errors given by
     * components. It should return 0 when all is well.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::verify() const
    {
      int ret = 0;

      // Verify components
      for (const auto& component : components_)
      {
        ret += component->verify();
      }

      return ret;
    }

    /**
     * @brief Initialize buses first, then all the other components.
     *
     * @pre All buses and components must be allocated at this point.
     * @pre Bus variables are written before component variables in the
     * system variable vector.
     *
     * Buses must be initialized before other components, because other
     * components may write to buses during the initialization.
     *
     * Also, generators may write to control devices (e.g. governors,
     * exciters, etc.) during the initialization.
     *
     * @todo Implement writing to system vectors in a thread-safe way.
     *
     * @note Currently assuming each component stores variables contiguously in memory and
     * that these are simply concatenated in the global system.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::initialize()
    {
      int status = 0;

      for (const auto& bus : buses_)
      {
        status += bus->initialize();
      }

      for (const auto& component : components_)
      {
        status += component->initialize();
      }

      y_.setDataUpdated();
      yp_.setDataUpdated();

      // For DependencyTracking::Variable, set variable numbers
      if constexpr (std::is_same_v<scalar_type, DependencyTracking::Variable>)
      {
        this->initializeDependencyTrackingVariableNumbers();
      }

      return status;
    }

    /**
     * @brief Add monitors from buses and components and start monitor
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::initializeMonitor()
    {
      for (const auto* bus : buses_)
      {
        auto* mon = bus->getMonitor();
        if (mon && !mon->empty())
        {
          monitor_->addMonitor(mon);
        }
      }

      for (const auto* component : components_)
      {
        auto* mon = component->getMonitor();
        if (mon && !mon->empty())
        {
          monitor_->addMonitor(mon);
        }
      }
    }

    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::startMonitor()
    {
      monitor_->start();
    }

    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::stopMonitor()
    {
      monitor_->stop();
    }

    template <typename scalar_type, typename index_type>
    bool SystemModel<scalar_type, index_type>::monitoring() const
    {
      return !monitor_->empty();
    }

    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::printMonitoredVariables() const
    {
      monitor_->print();
    }

    /**
     * @todo Tagging differential variables
     *
     * Identify what variables in the system of differential-algebraic
     * equations are differential variables, i.e. their derivatives
     * appear in the equations.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::tagDifferentiable()
    {
      // Set initial values for global solution vectors
      for (const auto& bus : buses_)
      {
        bus->tagDifferentiable();
        for (IdxT j = 0; j < bus->size(); ++j)
        {
          tag_[bus->getVariableIndex(j)] = bus->tag()[j];
        }
      }

      for (const auto& component : components_)
      {
        component->tagDifferentiable();
        for (IdxT j = 0; j < component->size(); ++j)
        {
          tag_[component->getVariableIndex(j)] = component->tag()[j];
        }
      }

      return 0;
    }

    /**
     * @brief Compute the absolute tolerance for each variable in the model
     *
     * @param rel_tol The relative tolerance which can be used to pick the
     *        absolute tolerance.
     * @return int 0 if successful, non-zero otherwise.
     *
     * This represents a "noise" level close to zero for which pure relative
     * error cannot be used.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      for (const auto& bus : buses_)
      {
        bus->setAbsoluteTolerance(rel_tol);
      }

      for (const auto& component : components_)
      {
        component->setAbsoluteTolerance(rel_tol);
      }

      abs_tol_.setDataUpdated();

      return 0;
    }

    /**
     * @brief Compute system residual vector
     *
     * Buses and components read and write their bound system-vector slices
     * directly.
     *
     * @warning Residuals must be computed for buses, before component
     * residuals are computed. Buses own residuals for currents
     * Ir and Ii, but the contributions to these residuals come
     * from components. Buses assign their residual values, while components
     * add to those values by in-place addition. This is why (for now) bus
     * residuals need to be computed first.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateResidual()
    {
      for (const auto& bus : buses_)
      {
        bus->evaluateResidual();
      }

      for (const auto& component : components_)
      {
        component->evaluateResidual();
      }

      f_.setDataUpdated();

      return 0;
    }

    /**
     * @brief Update time
     *
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::updateTime(RealT t, RealT a)
    {
      time_  = t;
      alpha_ = a;
      for (const auto& component : components_)
      {
        component->updateTime(t, a);
      }
    }

    /**
     * @brief Add bus
     *
     * Add bus at the end of the bus array and map bus ID with GridKit's ID for the bus
     *
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::addBus(BusT* bus)
    {
      IdxT gridkit_bus_id                = static_cast<IdxT>(buses_.size());
      gridkit_bus_indices_[bus->busID()] = gridkit_bus_id;
      buses_.push_back(bus);
      allocated_ = false;
    }

    /**
     * @brief Add component
     *
     * Add component at the end of the components array and set GridKit's component ID
     *
     * @todo: No integer user-facing component_id for now, but we could map GridKit's
     * component ID to the disambiguation_string
     *
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::addComponent(ComponentT* component)
    {
      IdxT gridkit_component_id = static_cast<IdxT>(components_.size());
      component->setGridKitComponentID(gridkit_component_id);
      component->setSystemBase(this->freq_system_base_,
                               this->va_system_base_);
      components_.push_back(component);
      allocated_ = false;
    }

    /**
     * @brief Add fault
     *
     * The fault is added to the components array, and we keep a map to its
     * location, so it can easily be accessed.
     *
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::addFault(ComponentT* component)
    {
      IdxT gridkit_component_id                = static_cast<IdxT>(components_.size());
      IdxT gridkit_fault_id                    = static_cast<IdxT>(gridkit_fault_indices_.size());
      gridkit_fault_indices_[gridkit_fault_id] = gridkit_component_id;
      addComponent(component);
    }

    /**
     * @brief Construct, connect, and add one device per data entry
     *
     * The device constructor takes one bus per enumerator of its `Buses`
     * enum, followed by its model data.
     *
     * @throws std::out_of_range if the data omits a bus or names an unknown one.
     */
    template <typename scalar_type, typename index_type>
    template <typename DeviceT>
    void SystemModel<scalar_type, index_type>::addDevices(
        const std::vector<typename DeviceT::ModelDataT>& device_data)
    {
      using BusesT = typename DeviceT::ModelDataT::Buses;

      constexpr std::size_t bus_count = Utilities::enum_size<BusesT>();
      static_assert(bus_count <= 2, "Devices attach to at most two buses");

      for (const auto& data : device_data)
      {
        DeviceT* device = nullptr;
        if constexpr (bus_count == 0)
        {
          device = new DeviceT(data);
        }
        else if constexpr (bus_count == 1)
        {
          device = new DeviceT(getBus(data.buses.at(BusesT::bus)), data);
        }
        else
        {
          device = new DeviceT(getBus(data.buses.at(BusesT::bus1)),
                               getBus(data.buses.at(BusesT::bus2)),
                               data);
        }

        if constexpr (requires { device->getPorts(); })
        {
          device->getPorts().connect(data, signal_nodes_);
        }

        if constexpr (std::is_same_v<DeviceT, BusFault<ScalarT, IdxT>>)
        {
          addFault(device);
        }
        else
        {
          addComponent(device);
        }
      }
    }

    /**
     * @brief Set system bases and propagate them to existing components.
     *
     * @param[in] freq_system_base - System frequency base in Hz.
     * @param[in] va_system_base - System power base in VA.
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::setSystemBase(
        RealT freq_system_base, RealT va_system_base)
    {
      ComponentT::setSystemBase(freq_system_base, va_system_base);

      for (auto* component : components_)
      {
        component->setSystemBase(freq_system_base, va_system_base);
      }
    }

    /**
     * @brief Return pointer to a bus
     *
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::BusT*
    SystemModel<scalar_type, index_type>::getBus(IdxT bus_id)
    {
      // Should fail if user-provided bus_id is incorrect
      IdxT gridkit_bus_id = gridkit_bus_indices_.at(bus_id);
      assert((buses_[gridkit_bus_id])->busID() == bus_id);
      return buses_[gridkit_bus_id];
    }

    /**
     * @brief Return pointer to a signal
     *
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::SignalNodeT*
    SystemModel<scalar_type, index_type>::getSignalNode(IdxT signal_id)
    {
      return signal_nodes_[signal_id];
    }

    /**
     * @brief Return pointer to a component
     *
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::ComponentT*
    SystemModel<scalar_type, index_type>::getComponent(IdxT gridkit_component_id)
    {
      // gridkit_component_id_ is set by System model and guaranteed to be unique
      return components_[gridkit_component_id];
    }

    /**
     * @brief Return pointer to a bus fault model
     *
     * This function is used to provide easier access to setting and
     * clearing faults from the SystemModel interface.
     *
     */
    template <typename scalar_type, typename index_type>
    BusFault<scalar_type, index_type>*
    SystemModel<scalar_type, index_type>::getBusFault(IdxT fault_id)
    {
      IdxT component_id = gridkit_fault_indices_.at(fault_id);
      return dynamic_cast<BusFault<ScalarT, IdxT>*>(components_[component_id]);
    }

  } // namespace PhasorDynamics
} // namespace GridKit
