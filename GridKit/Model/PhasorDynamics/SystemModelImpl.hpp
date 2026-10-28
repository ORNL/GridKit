#include <cassert>
#include <format>
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
     * @post All component models in SystemModelData, including its
     * subsystems, are created, and correctly connected into the system model.
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::SystemModel(const SystemModelData<RealT, IdxT>& data)
      : SystemModel(data, {})
    {
    }

    /**
     * @brief Construct a subsystem connected to its parent's buses (external buses)
     *
     * @param[in] data           - Subsystem data with the buses it owns
     * @param[in] external_buses - Parent's buses that subsystem devices may
     * use
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::SystemModel(const SystemModelData<RealT, IdxT>& data,
                                                      const std::vector<BusT*>&           external_buses)
      : SystemModel()
    {
      using namespace Governor;
      using namespace Exciter;
      using namespace Stabilizer;
      using namespace Controller;
      using namespace Converter;

      owns_components_ = true;
      monitor_         = std::make_unique<MonitorT>(time_);

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

      for (auto* bus : external_buses)
      {
        if (gridkit_bus_indices_.contains(bus->busID())
            || !external_buses_.emplace(bus->busID(), bus).second)
        {
          Log::error() << "External bus " << bus->busID()
                       << " is already in the subsystem\n";
          throw std::runtime_error("SystemModel construction failed");
        }
        ordered_external_buses_.push_back(bus);
      }

      // Add REGCA converters
      for (const auto& regcadata : data.regca)
      {
        IdxT bus_index = 0;
        if (regcadata.buses.contains(RegcaBuses::bus))
        {
          bus_index = regcadata.buses.at(RegcaBuses::bus);
        }

        auto* regca = new Regca<ScalarT, IdxT>(getBus(bus_index), regcadata);
        addComponent(regca, regcadata.label());

        regca->getPorts().connect(regcadata, signal_nodes_);
      }

      // Add branches
      for (const auto& branchdata : data.branch)
      {
        IdxT bus1_index = 0;
        if (branchdata.buses.contains(BranchBuses::bus1))
        {
          bus1_index = branchdata.buses.at(BranchBuses::bus1);
        }

        IdxT bus2_index = 0;
        if (branchdata.buses.contains(BranchBuses::bus2))
        {
          bus2_index = branchdata.buses.at(BranchBuses::bus2);
        }

        auto* branch = new BranchT(getBus(bus1_index),
                                   getBus(bus2_index),
                                   branchdata);

        addComponent(branch, branchdata.label());
      }

      // Add loads
      /// @todo Add loads to JSON parser
      for (const auto& loaddata : data.loadz)
      {
        IdxT bus_index = 0;
        if (loaddata.buses.contains(LoadZBuses::bus))
        {
          bus_index = loaddata.buses.at(LoadZBuses::bus);
        }
        auto* load = new LoadZ<ScalarT, IdxT>(getBus(bus_index), loaddata);
        addComponent(load, loaddata.label());
      }

      // Add zip loads
      /// @todo Add zip loads to JSON parser
      for (const auto& loadzipdata : data.loadzip)
      {
        IdxT bus_index = 0;
        if (loadzipdata.buses.contains(LoadZIPBuses::bus))
        {
          bus_index = loadzipdata.buses.at(LoadZIPBuses::bus);
        }
        auto* loadzip = new LoadZIP<ScalarT, IdxT>(getBus(bus_index),
                                                   loadzipdata);
        addComponent(loadzip, loadzipdata.label());
      }

      // Add GENROU generators
      for (const auto& gendata : data.genrou)
      {
        IdxT bus_index = 0;
        if (gendata.buses.contains(GenrouBuses::bus))
        {
          bus_index = gendata.buses.at(GenrouBuses::bus);
        }
        auto* gen = new Genrou<ScalarT, IdxT>(getBus(bus_index), gendata);
        addComponent(gen, gendata.label());
        gen->getPorts().connect(gendata, signal_nodes_);
      }

      // Add GENSAL generators
      for (const auto& gendata : data.gensal)
      {
        IdxT bus_index = 0;
        if (gendata.buses.contains(GensalBuses::bus))
        {
          bus_index = gendata.buses.at(GensalBuses::bus);
        }

        auto* gen = new Gensal<ScalarT, IdxT>(getBus(bus_index), gendata);
        addComponent(gen, gendata.label());
        gen->getPorts().connect(gendata, signal_nodes_);
      }

      // Add classical generators
      for (const auto& gendata : data.genclassical)
      {
        IdxT bus_index = 0;
        if (gendata.buses.contains(GenClassicalBuses::bus))
        {
          bus_index = gendata.buses.at(GenClassicalBuses::bus);
        }
        auto* gen = new GenClassical<ScalarT, IdxT>(getBus(bus_index), gendata);
        addComponent(gen, gendata.label());
        gen->getPorts().connect(gendata, signal_nodes_);
      }

      // Add REECB after its current-command and feedback producers because
      // components initialize in insertion order.
      for (const auto& reecbdata : data.reecb)
      {
        BusT* bus = nullptr;
        if (reecbdata.buses.contains(ReecbBuses::bus))
        {
          bus = getBus(reecbdata.buses.at(ReecbBuses::bus));
        }

        auto* reecb = new Reecb<ScalarT, IdxT>(bus, reecbdata);
        addComponent(reecb, reecbdata.label());

        reecb->getPorts().connect(reecbdata, signal_nodes_);
      }

      // Add Tgov1 governors
      for (const auto& govdata : data.gov)
      {
        auto* gov = new Tgov1<ScalarT, IdxT>(govdata);
        addComponent(gov, govdata.label());

        gov->getPorts().connect(govdata, signal_nodes_);
      }

      // Add GASTPTI governors
      for (const auto& gastptidata : data.gastpti)
      {
        auto* gastpti = new GastPti<ScalarT, IdxT>(gastptidata);
        addComponent(gastpti, gastptidata.label());

        gastpti->getPorts().connect(gastptidata, signal_nodes_);
      }

      // Add HYGOV governors
      for (const auto& hygovdata : data.hygov)
      {
        auto* hygov = new Hygov<ScalarT, IdxT>(hygovdata);
        addComponent(hygov, hygovdata.label());

        hygov->getPorts().connect(hygovdata, signal_nodes_);
      }

      // Add IEEEST stabilizers before exciters that consume their output during
      // initialization.
      for (const auto& stabdata : data.stabilizer)
      {
        auto* stabilizer = new Ieeest<ScalarT, IdxT>(stabdata);
        addComponent(stabilizer, stabdata.label());
        stabilizer->getPorts().connect(stabdata, signal_nodes_);
      }

      for (const auto& excitedata : data.exciter)
      {
        IdxT bus_index = 0;
        if (excitedata.buses.contains(Ieeet1Buses::bus))
        {
          bus_index = excitedata.buses.at(Ieeet1Buses::bus);
        }

        auto* exciter = new Ieeet1<ScalarT, IdxT>(getBus(bus_index), excitedata);
        addComponent(exciter, excitedata.label());

        exciter->getPorts().connect(excitedata, signal_nodes_);
      }

      for (const auto& excitedata : data.esdc1a)
      {
        BusT* bus = nullptr;
        if (excitedata.buses.contains(Esdc1aBuses::bus))
        {
          bus = getBus(excitedata.buses.at(Esdc1aBuses::bus));
        }

        auto* exciter = new Esdc1a<ScalarT, IdxT>(bus, excitedata);
        addComponent(exciter, excitedata.label());

        exciter->getPorts().connect(excitedata, signal_nodes_);
      }

      for (const auto& excitedata : data.sexspti)
      {
        IdxT bus_index = 0;
        if (excitedata.buses.contains(SexsPtiBuses::bus))
        {
          bus_index = excitedata.buses.at(SexsPtiBuses::bus);
        }

        auto* exciter = new SexsPti<ScalarT, IdxT>(getBus(bus_index), excitedata);
        addComponent(exciter, excitedata.label());

        exciter->getPorts().connect(excitedata, signal_nodes_);
      }

      // Add REPCA plant controllers after the signal producers they read at
      // initialization
      for (const auto& repcadata : data.repca)
      {
        BusT* bus = nullptr;
        if (repcadata.buses.contains(RepcaBuses::bus))
        {
          bus = getBus(repcadata.buses.at(RepcaBuses::bus));
        }

        auto* repca = new Repca<ScalarT, IdxT>(bus, repcadata);
        addComponent(repca, repcadata.label());

        repca->getPorts().connect(repcadata, signal_nodes_);
      }

      // Add constant signal sources
      for (const auto& srcdata : data.constant_source)
      {
        auto* source = new ConstantSignalSource<ScalarT, IdxT>(srcdata);
        addComponent(source, srcdata.label());
        source->getPorts().connect(srcdata, signal_nodes_);
      }

      // Add constant signal sources
      for (const auto& srcdata : data.function_source)
      {
        auto* source = new FunctionSignalSource<ScalarT, IdxT>(srcdata);
        addComponent(source, srcdata.label());
        source->getPorts().connect(srcdata, signal_nodes_);
      }

      // Add subsystems; they may use this system's buses and its external
      // buses. Added after devices, so devices above cannot use subsystem buses.
      if (!data.subsystems.empty())
      {
        std::vector<BusT*> visible_buses(external_buses);
        visible_buses.insert(visible_buses.end(), buses_.begin(), buses_.end());
        for (const auto& subsystem_data : data.subsystems)
        {
          auto subsystem = std::make_unique<SystemModel>(subsystem_data, visible_buses);
          addSubsystem(subsystem.get());
          subsystem.release();
        }
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
        // Subsystem devices may hold this system's buses; delete them first.
        for (auto subsystem : subsystems_)
        {
          delete subsystem;
        }

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
     * @note Should default to 0.
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
     * @pre This system is not a subsystem.
     *
     * @throws std::runtime_error if storage allocation, child binding, model
     * verification, or initialization for sparse Jacobian discovery fails.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::allocate()
    {
      if (is_subsystem_)
      {
        Log::error() << "A subsystem is allocated by its top-level system\n";
        throw std::runtime_error("SystemModel allocation failed");
      }

      computeSize();

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

      bindChildren(0);

      // Verify component configuration
      if (!this->verifyAndLog())
      {
        throw std::runtime_error("SystemModel allocation failed");
      }

      // Signal input indices need all variable indices and a verified configuration
      if (this->setSignalIndices() != 0)
      {
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
     * @brief Set and return the size of this system, including subsystems
     */
    template <typename scalar_type, typename index_type>
    typename SystemModel<scalar_type, index_type>::IdxT
    SystemModel<scalar_type, index_type>::computeSize()
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

      for (const auto& subsystem : subsystems_)
      {
        size_ += subsystem->computeSize();
      }

      return size_;
    }

    /**
     * @brief Bind buses, components and subsystems to this system's vectors
     *
     * @param[in] global_offset - Global index of this system's first variable
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::bindChildren(IdxT global_offset)
    {
      offset_ = global_offset;

      // New indices invalidate a subsystem's Jacobian.
      if (is_subsystem_)
      {
        delete csr_jac_;
        csr_jac_ = nullptr;

        delete[] map_to_csr_;
        map_to_csr_ = nullptr;

        nnz_ = 0;
      }

      tag_.resize(size_);
      variable_indices_.resize(size_);
      residual_indices_.resize(size_);

      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, global_offset + j);
        this->setResidualIndex(j, global_offset + j);
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
          bus->setVariableIndex(j, global_offset + offset + j);
          bus->setResidualIndex(j, global_offset + offset + j);
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
          component->setVariableIndex(j, global_offset + offset + j);
          component->setResidualIndex(j, global_offset + offset + j);
        }

        offset += component->size();
      }

      for (const auto& subsystem : subsystems_)
      {
        const int bind_status = subsystem->bind(y_, yp_, f_, abs_tol_, offset);
        if (bind_status != 0)
        {
          Log::error() << "Failed to bind subsystem vectors to system storage\n";
          throw std::runtime_error("SystemModel allocation failed");
        }

        subsystem->bindChildren(global_offset + offset);

        offset += subsystem->size();
      }

      if (offset != size_)
      {
        Log::error() << "Bound vector sizes do not match the system size\n";
        throw std::runtime_error("SystemModel allocation failed");
      }
    }

    /**
     * @brief Verify all components are configured correctly
     *
     * Collects every problem reported by the components, each prefixed with
     * the component ID. passed() is true when all is well.
     */
    template <typename scalar_type, typename index_type>
    Model::ConfigurationChecks SystemModel<scalar_type, index_type>::verify() const
    {
      Model::ConfigurationChecks checks;

      // Verify components
      for (const auto& component : components_)
      {
        const auto component_checks = component->verify();
        for (const auto& error : component_checks.errors())
        {
          checks.fail(std::format("component {}: {}", component->getGridKitComponentID(), error));
        }
      }

      for (const auto& subsystem : subsystems_)
      {
        const auto subsystem_checks = subsystem->verify();
        for (const auto& error : subsystem_checks.errors())
        {
          checks.fail(std::format("subsystem {}: {}", subsystem->getGridKitComponentID(), error));
        }
      }

      return checks;
    }

    /**
     * @brief Set signal input indices of all components, including those
     * in subsystems.
     *
     * @pre All variable indices are assigned and the configuration is verified.
     *
     * Signal nodes belong to one system, so each signal source and its
     * inputs are in the same subsystem.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::setSignalIndices()
    {
      int status = 0;

      for (const auto& component : components_)
      {
        if (component->setSignalIndices() != 0)
        {
          Log::error() << "Failed to set signal indices of component "
                       << component->getGridKitComponentID() << '\n';
          status = 1;
        }
      }

      for (const auto& subsystem : subsystems_)
      {
        if (subsystem->setSignalIndices() != 0)
        {
          Log::error() << "Failed to set signal indices in subsystem "
                       << subsystem->getGridKitComponentID() << '\n';
          status = 1;
        }
      }

      return status;
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

      // Subsystem devices need the external buses initialized above.
      for (const auto& subsystem : subsystems_)
      {
        status += subsystem->initialize();
      }

      y_.setDataUpdated();
      yp_.setDataUpdated();

      // For DependencyTracking::Variable, set variable numbers at the top level
      if constexpr (std::is_same_v<scalar_type, DependencyTracking::Variable>)
      {
        if (!is_subsystem_)
        {
          this->initializeDependencyTrackingVariableNumbers();
        }
      }

      return status;
    }

    /**
     * @brief Add monitors from buses and components and start monitor
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::initializeMonitor()
    {
      addMonitorsTo(*monitor_);
    }

    /**
     * @brief Add monitors from buses, components and subsystems to `monitor`
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::addMonitorsTo(MonitorT& monitor) const
    {
      for (const auto* bus : buses_)
      {
        auto* mon = bus->getMonitor();
        if (mon && !mon->empty())
        {
          monitor.addMonitor(mon);
        }
      }

      for (const auto* component : components_)
      {
        auto* mon = component->getMonitor();
        if (mon && !mon->empty())
        {
          monitor.addMonitor(mon);
        }
      }

      for (const auto* subsystem : subsystems_)
      {
        subsystem->addMonitorsTo(monitor);
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
      // Set initial values for global solution vectors; tag_ starts at offset_
      for (const auto& bus : buses_)
      {
        bus->tagDifferentiable();
        for (IdxT j = 0; j < bus->size(); ++j)
        {
          tag_[bus->getVariableIndex(j) - offset_] = bus->tag()[j];
        }
      }

      for (const auto& component : components_)
      {
        component->tagDifferentiable();
        for (IdxT j = 0; j < component->size(); ++j)
        {
          tag_[component->getVariableIndex(j) - offset_] = component->tag()[j];
        }
      }

      for (const auto& subsystem : subsystems_)
      {
        subsystem->tagDifferentiable();
        for (IdxT j = 0; j < subsystem->size(); ++j)
        {
          tag_[subsystem->getVariableIndex(j) - offset_] = subsystem->tag()[j];
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

      for (const auto& subsystem : subsystems_)
      {
        subsystem->setAbsoluteTolerance(rel_tol);
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
     * For the same reason, external buses are evaluated before subsystems.
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

      for (const auto& subsystem : subsystems_)
      {
        subsystem->evaluateResidual();
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

      for (const auto& subsystem : subsystems_)
      {
        subsystem->updateTime(t, a);
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
     * @param[in] component - Component to add
     * @param[in] label     - ComponentData::label() of the device (e.g. "Genrou_1"),
     *                        or empty for components without case-file data
     *
     * @pre This system has no subsystems.
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::addComponent(ComponentT* component, const std::string& label)
    {
      if (!subsystems_.empty())
      {
        Log::error() << "A system with subsystems cannot contain devices\n";
        throw std::runtime_error("SystemModel::addComponent failed");
      }

      IdxT gridkit_component_id = static_cast<IdxT>(components_.size());
      component->setGridKitComponentID(gridkit_component_id);
      component->setSystemBase(this->freq_system_base_,
                               this->va_system_base_);
      components_.push_back(component);
      component_labels_.push_back(label);
      if (!label.empty())
      {
        gridkit_component_indices_[label] = gridkit_component_id;
      }
      allocated_ = false;
    }

    /**
     * @brief Add subsystem; deleted with this system if it owns its components
     *
     * @pre The subsystem is not allocated and not part of another system, and
     * has the same system bases. This system has no devices.
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::addSubsystem(SystemModel* subsystem)
    {
      if (subsystem == nullptr || subsystem == this || subsystem->is_subsystem_)
      {
        Log::error() << "Invalid subsystem or subsystem already in a system\n";
        throw std::runtime_error("SystemModel::addSubsystem failed");
      }

      if (!components_.empty())
      {
        Log::error() << "A system with devices cannot contain subsystems\n";
        throw std::runtime_error("SystemModel::addSubsystem failed");
      }

      if (subsystem->freq_system_base_ != this->freq_system_base_
          || subsystem->va_system_base_ != this->va_system_base_)
      {
        Log::error() << "Subsystem bases differ from the parent's\n";
        throw std::runtime_error("SystemModel::addSubsystem failed");
      }

      subsystem->is_subsystem_ = true;
      subsystem->setGridKitComponentID(static_cast<IdxT>(subsystems_.size()));
      subsystems_.push_back(subsystem);
      allocated_ = false;
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

      for (auto* subsystem : subsystems_)
      {
        subsystem->setSystemBase(freq_system_base, va_system_base);
      }
    }

    /**
     * @brief Return pointer to a bus of this system, its subsystems, or an
     * external bus
     *
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::BusT*
    SystemModel<scalar_type, index_type>::getBus(IdxT bus_id)
    {
      BusT* bus = findBus(bus_id);
      if (bus == nullptr && external_buses_.contains(bus_id))
      {
        bus = external_buses_.at(bus_id);
      }
      if (bus == nullptr)
      {
        Log::error() << "Unknown bus " << bus_id << "\n";
        throw std::runtime_error("SystemModel bus lookup failed");
      }
      return bus;
    }

    /**
     * @brief Return pointer to a bus of this system or its subsystems, or
     * nullptr if there is none
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::BusT*
    SystemModel<scalar_type, index_type>::findBus(IdxT bus_id)
    {
      BusT*      bus = nullptr;
      const auto it  = gridkit_bus_indices_.find(bus_id);
      if (it != gridkit_bus_indices_.end())
      {
        bus = buses_[it->second];
        assert(bus->busID() == bus_id);
      }

      for (IdxT subsystem = 0; bus == nullptr && subsystem < numSubsystems(); ++subsystem)
      {
        bus = subsystems_[subsystem]->findBus(bus_id);
      }
      return bus;
    }

    template <typename scalar_type, typename index_type>
    typename SystemModel<scalar_type, index_type>::IdxT
    SystemModel<scalar_type, index_type>::numComponents() const
    {
      return static_cast<IdxT>(components_.size());
    }

    /**
     * @brief Return the label (device class and case-file ID) of a component of
     * this system, empty if it has none
     */
    template <typename scalar_type, typename index_type>
    const std::string&
    SystemModel<scalar_type, index_type>::getComponentLabel(IdxT gridkit_component_id) const
    {
      return component_labels_.at(gridkit_component_id);
    }

    template <typename scalar_type, typename index_type>
    typename SystemModel<scalar_type, index_type>::IdxT
    SystemModel<scalar_type, index_type>::numSubsystems() const
    {
      return static_cast<IdxT>(subsystems_.size());
    }

    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>*
    SystemModel<scalar_type, index_type>::getSubsystem(IdxT index)
    {
      return subsystems_.at(index);
    }

    /**
     * @brief Global index of this system's first variable and residual
     */
    template <typename scalar_type, typename index_type>
    typename SystemModel<scalar_type, index_type>::IdxT
    SystemModel<scalar_type, index_type>::offset() const
    {
      return offset_;
    }

    /**
     * @brief External buses of this subsystem, in construction order
     */
    template <typename scalar_type, typename index_type>
    const std::vector<typename SystemModel<scalar_type, index_type>::BusT*>&
    SystemModel<scalar_type, index_type>::getExternalBuses() const
    {
      return ordered_external_buses_;
    }

    /**
     * @brief Return pointer to a signal of this system or its subsystems
     *
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::SignalNodeT*
    SystemModel<scalar_type, index_type>::getSignalNode(IdxT signal_id)
    {
      SignalNodeT* node = findSignalNode(signal_id);
      if (node == nullptr)
      {
        Log::error() << "Unknown signal " << signal_id << "\n";
        throw std::runtime_error("SystemModel signal lookup failed");
      }
      return node;
    }

    /**
     * @brief Return pointer to a signal of this system or its subsystems, or
     * nullptr if there is none
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::SignalNodeT*
    SystemModel<scalar_type, index_type>::findSignalNode(IdxT signal_id)
    {
      SignalNodeT* node = signal_nodes_.find(signal_id);
      for (IdxT subsystem = 0; node == nullptr && subsystem < numSubsystems(); ++subsystem)
      {
        node = subsystems_[subsystem]->findSignalNode(signal_id);
      }
      return node;
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
     * @brief Return pointer to a branch of this system or its subsystems by
     * case-file `id`
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::BranchT*
    SystemModel<scalar_type, index_type>::getBranch(const std::string& id)
    {
      BranchT* branch = dynamic_cast<BranchT*>(findComponent("Branch_" + id));
      if (branch == nullptr)
      {
        Log::error() << "Unknown branch " << id << "\n";
        throw std::runtime_error("SystemModel branch lookup failed");
      }
      return branch;
    }

    /**
     * @brief Return pointer to a component of this system or its subsystems by
     * label (device class and case-file ID, e.g. "Genrou_1")
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::ComponentT*
    SystemModel<scalar_type, index_type>::getComponent(const std::string& label)
    {
      ComponentT* component = findComponent(label);
      if (component == nullptr)
      {
        Log::error() << "Unknown component " << label << "\n";
        throw std::runtime_error("SystemModel component lookup failed");
      }
      return component;
    }

    /**
     * @brief Return pointer to a component of this system or its subsystems by
     * label, or nullptr if there is none
     */
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::ComponentT*
    SystemModel<scalar_type, index_type>::findComponent(const std::string& label)
    {
      ComponentT* component = nullptr;
      const auto  it        = gridkit_component_indices_.find(label);
      if (it != gridkit_component_indices_.end())
      {
        component = components_[it->second];
      }

      for (IdxT subsystem = 0; component == nullptr && subsystem < numSubsystems(); ++subsystem)
      {
        component = subsystems_[subsystem]->findComponent(label);
      }
      return component;
    }

  } // namespace PhasorDynamics
} // namespace GridKit
