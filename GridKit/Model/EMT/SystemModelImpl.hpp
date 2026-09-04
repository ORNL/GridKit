#include <cassert>
#include <iostream>

#include <GridKit/Definitions.hpp>
#include <GridKit/Model/EMT/SystemModel.hpp>
#include <GridKit/Model/EMT/SystemModelData.hpp>
#include <GridKit/Model/VariableMonitorController.hpp>

// Include all components
#include <GridKit/Model/EMT/ComponentLibrary.hpp>

namespace GridKit
{
  namespace EMT
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
     *
     * Buses are constructed first so subsequent devices can attach their
     * electrical ports; buses initialize first for the same reason.
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::SystemModel(const SystemModelData<RealT, IdxT>& data)
      : monitor_(std::make_unique<MonitorT>(time_))
    {
      owns_components_ = true;

      // Add electrical buses
      for (const auto& busdata : data.bus)
      {
        auto* bus = new Bus<ScalarT, IdxT>(busdata);
        addBus(bus);
      }

      // Add signals
      for (const auto& signaldata : data.signal)
      {
        Signal<ScalarT, IdxT>* signal = new Signal<ScalarT, IdxT>(signaldata);
        addSignal(signal);
      }

      // Add voltage sources
      for (const auto& sourcedata : data.voltage_source)
      {
        IdxT bus_index = 0;
        if (sourcedata.buses.contains(VoltageSourceBuses::bus))
        {
          bus_index = sourcedata.buses.at(VoltageSourceBuses::bus);
        }

        auto* source = new VoltageSource<ScalarT, IdxT>(sourcedata);

        constexpr auto VA = VoltageSourceExternalVariables::VA;
        source->getSignals().template attachPort<VA>(&(getBus(bus_index)->voltagePort()));

        addComponent(source);
      }

      // Add dependent voltage sources
      for (const auto& sourcedata : data.dependent_voltage_source)
      {
        IdxT bus_index = 0;
        if (sourcedata.buses.contains(DependentVoltageSourceBuses::bus))
        {
          bus_index = sourcedata.buses.at(DependentVoltageSourceBuses::bus);
        }

        auto* source = new DependentVoltageSource<ScalarT, IdxT>(sourcedata);

        constexpr auto VA = DependentVoltageSourceExternalVariables::VA;
        source->getSignals().template attachPort<VA>(&(getBus(bus_index)->voltagePort()));

        if (sourcedata.signal_inputs.contains(DependentVoltageSourceSignalInputs::ea))
        {
          IdxT           ea = sourcedata.signal_inputs.at(DependentVoltageSourceSignalInputs::ea);
          constexpr auto EA = DependentVoltageSourceExternalVariables::EA;
          source->getSignals().template attachSignal<EA>(getSignal(ea));
        }

        if (sourcedata.signal_inputs.contains(DependentVoltageSourceSignalInputs::eb))
        {
          IdxT           eb = sourcedata.signal_inputs.at(DependentVoltageSourceSignalInputs::eb);
          constexpr auto EB = DependentVoltageSourceExternalVariables::EB;
          source->getSignals().template attachSignal<EB>(getSignal(eb));
        }

        if (sourcedata.signal_inputs.contains(DependentVoltageSourceSignalInputs::ec))
        {
          IdxT           ec = sourcedata.signal_inputs.at(DependentVoltageSourceSignalInputs::ec);
          constexpr auto EC = DependentVoltageSourceExternalVariables::EC;
          source->getSignals().template attachSignal<EC>(getSignal(ec));
        }

        addComponent(source);
      }

      // Add synchronous machines
      for (const auto& machinedata : data.machine)
      {
        IdxT bus_index = 0;
        if (machinedata.buses.contains(MachineBuses::bus))
        {
          bus_index = machinedata.buses.at(MachineBuses::bus);
        }

        auto* machine = new Machine<ScalarT, IdxT>(machinedata);

        constexpr auto VA = MachineExternalVariables::VA;
        machine->getSignals().template attachPort<VA>(&(getBus(bus_index)->voltagePort()));

        if (machinedata.signal_outputs.contains(MachineSignalOutputs::speed))
        {
          IdxT           speed = machinedata.signal_outputs.at(MachineSignalOutputs::speed);
          constexpr auto OMEGA = MachineInternalVariables::OMEGA;
          machine->getSignals().template assignSignal<OMEGA>(getSignal(speed));
        }

        if (machinedata.signal_inputs.contains(MachineSignalInputs::pm))
        {
          IdxT           pm = machinedata.signal_inputs.at(MachineSignalInputs::pm);
          constexpr auto PM = MachineExternalVariables::PM;
          machine->getSignals().template attachSignal<PM>(getSignal(pm));
        }

        if (machinedata.signal_inputs.contains(MachineSignalInputs::efd))
        {
          IdxT           efd = machinedata.signal_inputs.at(MachineSignalInputs::efd);
          constexpr auto EFD = MachineExternalVariables::EFD;
          machine->getSignals().template attachSignal<EFD>(getSignal(efd));
        }

        addComponent(machine);
      }

      // Add lumped lines
      for (const auto& linedata : data.line_lumped)
      {
        IdxT bus1_index = 0;
        if (linedata.buses.contains(LineLumpedBuses::bus1))
        {
          bus1_index = linedata.buses.at(LineLumpedBuses::bus1);
        }

        IdxT bus2_index = 0;
        if (linedata.buses.contains(LineLumpedBuses::bus2))
        {
          bus2_index = linedata.buses.at(LineLumpedBuses::bus2);
        }

        auto* line = new LineLumped<ScalarT, IdxT>(linedata);

        constexpr auto V1A = LineLumpedExternalVariables::V1A;
        constexpr auto V2A = LineLumpedExternalVariables::V2A;
        line->getSignals().template attachPort<V1A>(&(getBus(bus1_index)->voltagePort()));
        line->getSignals().template attachPort<V2A>(&(getBus(bus2_index)->voltagePort()));

        addComponent(line);
      }

      // Add loads
      for (const auto& loaddata : data.loadz)
      {
        IdxT bus_index = 0;
        if (loaddata.buses.contains(LoadZBuses::bus))
        {
          bus_index = loaddata.buses.at(LoadZBuses::bus);
        }

        auto* load = new LoadZ<ScalarT, IdxT>(loaddata);

        constexpr auto VA = LoadZExternalVariables::VA;
        load->getSignals().template attachPort<VA>(&(getBus(bus_index)->voltagePort()));

        addComponent(load);
      }

      // Add Tgov1 governors after the machines they read at initialization
      for (const auto& govdata : data.gov)
      {
        auto* gov = new Controller::Tgov1<ScalarT, IdxT>(govdata);

        if (govdata.signal_inputs.contains(Controller::Tgov1SignalInputs::speed))
        {
          IdxT           speed = govdata.signal_inputs.at(Controller::Tgov1SignalInputs::speed);
          constexpr auto OMEGA = Controller::Tgov1ExternalVariables::OMEGA;
          gov->getSignals().template attachSignal<OMEGA>(getSignal(speed));
        }

        if (govdata.signal_inputs.contains(Controller::Tgov1SignalInputs::pref))
        {
          IdxT           pref = govdata.signal_inputs.at(Controller::Tgov1SignalInputs::pref);
          constexpr auto PREF = Controller::Tgov1ExternalVariables::PREF;
          gov->getSignals().template attachSignal<PREF>(getSignal(pref));
        }

        if (govdata.signal_outputs.contains(Controller::Tgov1SignalOutputs::pmech))
        {
          IdxT           pmech = govdata.signal_outputs.at(Controller::Tgov1SignalOutputs::pmech);
          constexpr auto PM    = Controller::Tgov1InternalVariables::PM;
          gov->getSignals().template assignSignal<PM>(getSignal(pmech));
        }

        addComponent(gov);
      }

      // Add switches
      for (const auto& switchdata : data.sw)
      {
        IdxT bus1_index = 0;
        if (switchdata.buses.contains(SwitchBuses::bus1))
        {
          bus1_index = switchdata.buses.at(SwitchBuses::bus1);
        }

        IdxT bus2_index = 0;
        if (switchdata.buses.contains(SwitchBuses::bus2))
        {
          bus2_index = switchdata.buses.at(SwitchBuses::bus2);
        }

        auto* sw = new Switch<ScalarT, IdxT>(switchdata);

        constexpr auto V1A = SwitchExternalVariables::V1A;
        constexpr auto V2A = SwitchExternalVariables::V2A;
        sw->getSignals().template attachPort<V1A>(&(getBus(bus1_index)->voltagePort()));
        sw->getSignals().template attachPort<V2A>(&(getBus(bus2_index)->voltagePort()));

        addSwitch(sw);
      }

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

        for (auto signal : signals_)
        {
          delete signal;
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
     * @brief Allocate system storage and bind components to it.
     *
     * First sum the component sizes, then allocate the system vectors and
     * bind each component to its portion of those vectors.
     *
     * @pre Components with nonzero size are unallocated or already bound to
     * external system storage.
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

      for (const auto& component : components_)
      {
        size_ += component->size();
      }

      // Allocate global vectors
      if (!allocated_)
      {
        // Topology changes invalidate the Jacobian sparsity pattern and COO-to-CSR map.
        resetJacobianStructure();
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
     * @brief Initialize all components in insertion order.
     *
     * @pre All components must be allocated at this point.
     *
     * Buses are inserted first by the data constructor, so they initialize
     * before the components that read or write bus voltages during their
     * initialization.
     *
     * @todo Implement writing to system vectors in a thread-safe way.
     *
     * @note Currently assuming each component stores variables contiguously in memory and
     * that these are simply concateneted in the global system.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::initialize()
    {
      int status = 0;

      for (const auto& component : components_)
      {
        status += component->initialize();
      }

      y_.setDataUpdated();
      yp_.setDataUpdated();

      return status;
    }

    /**
     * @brief Add monitors from components and start monitor
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::initializeMonitor()
    {
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
     * @brief Tag differential variables.
     *
     * Identify what variables in the system of differential-algebraic
     * equations are differential variables, i.e. their derivatives
     * appear in the equations.
     *
     * @pre All components must be allocated, so connection-dependent
     * derivative-coupling flags are final.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::tagDifferentiable()
    {
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
      for (const auto& component : components_)
      {
        component->setAbsoluteTolerance(rel_tol);
      }

      abs_tol_.setDataUpdated();

      return 0;
    }

    /**
     * @brief Compute the residuals each component owns.
     *
     * Components read and write their bound system-vector slices directly.
     * The internal phase assigns every owned residual row over the whole
     * component list before any external accumulation begins.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateInternalResidual()
    {
      for (const auto& component : components_)
      {
        component->evaluateInternalResidual();
      }

      return 0;
    }

    /**
     * @brief Accumulate component contributions to residuals owned elsewhere,
     * e.g. bus current balances.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateExternalResidual()
    {
      for (const auto& component : components_)
      {
        component->evaluateExternalResidual();
      }

      return 0;
    }

    /**
     * @brief Compute system residual vector
     *
     * Internal residuals assign every owned entry of the residual vector,
     * then external residuals accumulate the remaining contributions.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateResidual()
    {
      evaluateInternalResidual();
      evaluateExternalResidual();

      f_.setDataUpdated();

      return 0;
    }

    /**
     * @brief Evaluate system Jacobian.
     *
     * Evaluate component Jacobians, then concatenate them into the system COO
     * and deduplicate into the CSR Jacobian. On subsequent calls only the CSR
     * values are refreshed through the COO-to-CSR map.
     *
     * A component whose rows are constant, such as a bus, owns no Jacobian
     * entries and returns a null COO Jacobian; it is skipped in the same
     * order in both passes.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateJacobian()
    {
      // Evaluate component Jacobians, including contributions to rows owned
      // by other components
      for (const auto& component : components_)
      {
        component->evaluateJacobian();
      }

      // Build or update system CSR Jacobian
      if (csr_jac_ == nullptr)
      {
        // Count the number of non-zeros
        IdxT nnz_dup = 0;
        for (const auto& component : components_)
        {
          auto component_jacobian = component->getCooJacobian();

          if (component_jacobian != nullptr)
          {
            nnz_dup += component_jacobian->getNnz();
          }
        }

        // Allocate COO triplet arrays (we own these until we hand off to CsrMatrix)
        IdxT*  rows_dup = new IdxT[static_cast<size_t>(nnz_dup)];
        IdxT*  cols_dup = new IdxT[static_cast<size_t>(nnz_dup)];
        RealT* vals_dup = new RealT[static_cast<size_t>(nnz_dup)];

        IdxT counter = 0;
        for (const auto& component : components_)
        {
          auto component_jacobian = component->getCooJacobian();

          if (component_jacobian != nullptr)
          {
            const IdxT*  rows    = component_jacobian->getRowData();
            const IdxT*  columns = component_jacobian->getColData();
            const RealT* values  = component_jacobian->getValues();
            for (IdxT i = 0; i < component_jacobian->getNnz(); ++i)
            {
              rows_dup[counter] = rows[i];
              cols_dup[counter] = columns[i];
              vals_dup[counter] = values[i];
              counter++;
            }
          }
        }

        // Build the system COO Jacobian
        CooMatrixT jac(size_, size_, nnz_dup, &rows_dup, &cols_dup, &vals_dup);

        // Populate CSR data with sort and deduplicate
        IdxT* row_ptrs = jac.getCsrRowData();

        // Deduplicated nnz
        nnz_ = jac.getNnz();

        // Allocate cols/vals with deduplicated nnz
        IdxT*  cols = new IdxT[static_cast<size_t>(nnz_)];
        RealT* vals = new RealT[static_cast<size_t>(nnz_)];

        std::copy(jac.getColData(), jac.getColData() + nnz_, cols);
        std::copy(jac.getValues(), jac.getValues() + nnz_, vals);

        // Create the CSR Jacobian
        csr_jac_ = new CsrMatrixT(size_, size_, nnz_, &row_ptrs, &cols, &vals);

        const IdxT* map_to_sorted = jac.getMapToSorted();
        const IdxT* map_to_dedup  = jac.getMapToDeduplicated();

        // Build a mappping from original COO index to CSR index
        map_to_csr_ = new IdxT[static_cast<size_t>(nnz_dup)];
        for (IdxT i = 0; i < nnz_dup; ++i)
        {
          map_to_csr_[map_to_sorted[i]] = map_to_dedup[i];
        }
      }
      else
      {
        // Zero out values
        RealT* vals = csr_jac_->getValues();
        for (IdxT i = 0; i < csr_jac_->getNnz(); ++i)
        {
          vals[i] = 0.0;
        }

        // Update CSR values from component Jacobians
        IdxT counter = 0;
        for (const auto& component : components_)
        {
          auto component_jacobian = component->getCooJacobian();

          if (component_jacobian != nullptr)
          {
            const RealT* values = component_jacobian->getValues();
            for (IdxT i = 0; i < component_jacobian->getNnz(); ++i)
            {
              vals[map_to_csr_[counter]] += values[i];
              counter++;
            }
          }
        }
      }

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
     * @brief Invalidate the assembled Jacobian structure.
     *
     * Every component's cached structure is invalidated along with the
     * system CSR Jacobian and the COO-to-CSR map.
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::resetJacobianStructure()
    {
      delete csr_jac_;
      csr_jac_ = nullptr;

      delete[] map_to_csr_;
      map_to_csr_ = nullptr;

      nnz_ = 0;

      for (const auto& component : components_)
      {
        component->resetJacobianStructure();
      }
    }

    /**
     * @brief Add bus
     *
     * The bus is added to the components array like any other component, and
     * we keep a map from its user-facing bus ID to its component index so
     * device ports can be wired to it while building the system.
     *
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::addBus(BusT* bus)
    {
      gridkit_bus_indices_[bus->busID()] = static_cast<IdxT>(components_.size());
      addComponent(bus);
    }

    /**
     * @brief Add signal
     *
     * Add signal at the end of the signals array and map signal ID with
     * GridKit's ID for the signal
     *
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::addSignal(SignalT* signal)
    {
      IdxT gridkit_signal_id                      = static_cast<IdxT>(signals_.size());
      gridkit_signal_indices_[signal->signalId()] = gridkit_signal_id;
      signals_.push_back(signal);
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
      components_.push_back(component);
      allocated_ = false;
    }

    /**
     * @brief Add switch
     *
     * The switch is added to the components array, and we keep a map to its
     * location, so it can easily be accessed for discrete events.
     *
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::addSwitch(ComponentT* component)
    {
      IdxT gridkit_component_id                  = static_cast<IdxT>(components_.size());
      IdxT gridkit_switch_id                     = static_cast<IdxT>(gridkit_switch_indices_.size());
      gridkit_switch_indices_[gridkit_switch_id] = gridkit_component_id;
      addComponent(component);
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
      IdxT  gridkit_component_id = gridkit_bus_indices_.at(bus_id);
      auto* bus                  = dynamic_cast<BusT*>(components_[gridkit_component_id]);
      assert(bus != nullptr && bus->busID() == bus_id);
      return bus;
    }

    /**
     * @brief Return pointer to a signal
     *
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::SignalT*
    SystemModel<scalar_type, index_type>::getSignal(IdxT signal_id)
    {
      // Should fail if user-provided signal_id is incorrect
      IdxT gridkit_signal_id = gridkit_signal_indices_.at(signal_id);
      assert((signals_[gridkit_signal_id])->signalId() == signal_id);
      return signals_[gridkit_signal_id];
    }

    /**
     * @brief Return pointer to a component
     *
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::ComponentT*
    SystemModel<scalar_type, index_type>::getComponent(IdxT gridkit_component_id)
    {
      // gridkit_component_id_ is set by System model and guarantied to be unique
      return components_[gridkit_component_id];
    }

    /**
     * @brief Return pointer to a switch model
     *
     * This function is used to provide easier access to opening and closing
     * switches from the SystemModel interface.
     *
     */
    template <typename scalar_type, typename index_type>
    Switch<scalar_type, index_type>*
    SystemModel<scalar_type, index_type>::getSwitch(IdxT switch_id)
    {
      IdxT component_id = gridkit_switch_indices_.at(switch_id);
      return dynamic_cast<Switch<ScalarT, IdxT>*>(components_[component_id]);
    }

  } // namespace EMT
} // namespace GridKit
