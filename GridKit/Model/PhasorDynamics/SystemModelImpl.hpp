#include <algorithm>
#include <cassert>
#include <iostream>

#include <GridKit/Constants.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusFactory.hpp>
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

      // Add electrical buses. Each bus publishes its voltage on two signal
      // nodes the system creates and assigns here.
      for (const auto& busdata : data.bus)
      {
        BusBase<ScalarT, IdxT>* bus = BusFactory<ScalarT, IdxT>::create(busdata);

        auto vr = std::make_unique<SignalNode<ScalarT, IdxT>>();
        auto vi = std::make_unique<SignalNode<ScalarT, IdxT>>();
        bus->getSignals().template assignSignalNode<BusInternalVariables::VR>(vr.get());
        bus->getSignals().template assignSignalNode<BusInternalVariables::VI>(vi.get());
        bus_signals_.push_back(std::move(vr));
        bus_signals_.push_back(std::move(vi));

        addBus(bus);
      }

      // Add signal nodes
      for (const auto& signaldata : data.signal)
      {
        signal_nodes_.add(signaldata);
      }

      // Add bus-to-signal adapters
      for (const auto& adapterdata : data.adapter)
      {
        IdxT bus_index = 0;
        if (adapterdata.buses.contains(BusToSignalAdapterBuses::bus))
        {
          bus_index = adapterdata.buses.at(BusToSignalAdapterBuses::bus);
        }

        auto* adapter = new BusToSignalAdapter<ScalarT, IdxT>(getBus(bus_index), adapterdata);
        adapter->getPorts().connect(adapterdata, signal_nodes_);
        addComponent(adapter);
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
        regca->getPorts().connect(regcadata, signal_nodes_);

        addComponent(regca);
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

        auto* branch = new Branch<ScalarT, IdxT>(branchdata);

        auto& bus1_signals = getBus(bus1_index)->getSignals();
        auto& bus2_signals = getBus(bus2_index)->getSignals();
        branch->getSignals().template attachSignalNode<BranchExternalVariables::VR1>(
            bus1_signals.template getSignalNode<BusInternalVariables::VR>());
        branch->getSignals().template attachSignalNode<BranchExternalVariables::VI1>(
            bus1_signals.template getSignalNode<BusInternalVariables::VI>());
        branch->getSignals().template attachSignalNode<BranchExternalVariables::VR2>(
            bus2_signals.template getSignalNode<BusInternalVariables::VR>());
        branch->getSignals().template attachSignalNode<BranchExternalVariables::VI2>(
            bus2_signals.template getSignalNode<BusInternalVariables::VI>());
        addComponent(branch);
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
        auto* load = new LoadZ<ScalarT, IdxT>(loaddata);

        auto& bus_signals = getBus(bus_index)->getSignals();
        load->getSignals().template attachSignalNode<LoadZExternalVariables::VR>(
            bus_signals.template getSignalNode<BusInternalVariables::VR>());
        load->getSignals().template attachSignalNode<LoadZExternalVariables::VI>(
            bus_signals.template getSignalNode<BusInternalVariables::VI>());
        addComponent(load);
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
        auto* loadzip = new LoadZIP<ScalarT, IdxT>(loadzipdata);

        auto& bus_signals = getBus(bus_index)->getSignals();
        loadzip->getSignals().template attachSignalNode<LoadZIPExternalVariables::VR>(
            bus_signals.template getSignalNode<BusInternalVariables::VR>());
        loadzip->getSignals().template attachSignalNode<LoadZIPExternalVariables::VI>(
            bus_signals.template getSignalNode<BusInternalVariables::VI>());
        addComponent(loadzip);
      }

      // Add GENROU generators
      for (const auto& gendata : data.genrou)
      {
        IdxT bus_index = 0;
        if (gendata.buses.contains(GenrouBuses::bus))
        {
          bus_index = gendata.buses.at(GenrouBuses::bus);
        }

        auto* gen = new Genrou<ScalarT, IdxT>(gendata);

        auto& bus_signals = getBus(bus_index)->getSignals();
        gen->getSignals().template attachSignalNode<GenrouExternalVariables::VR>(
            bus_signals.template getSignalNode<BusInternalVariables::VR>());
        gen->getSignals().template attachSignalNode<GenrouExternalVariables::VI>(
            bus_signals.template getSignalNode<BusInternalVariables::VI>());

        gen->getPorts().connect(gendata, signal_nodes_);
        addComponent(gen);
      }

      // Add GENSAL generators
      for (const auto& gendata : data.gensal)
      {
        IdxT bus_index = 0;
        if (gendata.buses.contains(GensalBuses::bus))
        {
          bus_index = gendata.buses.at(GensalBuses::bus);
        }

        auto* gen = new Gensal<ScalarT, IdxT>(gendata);

        auto& bus_signals = getBus(bus_index)->getSignals();
        gen->getSignals().template attachSignalNode<GensalExternalVariables::VR>(
            bus_signals.template getSignalNode<BusInternalVariables::VR>());
        gen->getSignals().template attachSignalNode<GensalExternalVariables::VI>(
            bus_signals.template getSignalNode<BusInternalVariables::VI>());

        gen->getPorts().connect(gendata, signal_nodes_);
        addComponent(gen);
      }

      // Add classical generators
      for (const auto& gendata : data.genclassical)
      {
        IdxT bus_index = 0;
        if (gendata.buses.contains(GenClassicalBuses::bus))
        {
          bus_index = gendata.buses.at(GenClassicalBuses::bus);
        }
        auto* gen = new GenClassical<ScalarT, IdxT>(gendata);

        auto& bus_signals = getBus(bus_index)->getSignals();
        gen->getSignals().template attachSignalNode<GenClassicalExternalVariables::VR>(
            bus_signals.template getSignalNode<BusInternalVariables::VR>());
        gen->getSignals().template attachSignalNode<GenClassicalExternalVariables::VI>(
            bus_signals.template getSignalNode<BusInternalVariables::VI>());
        gen->getPorts().connect(gendata, signal_nodes_);
        addComponent(gen);
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
        reecb->getPorts().connect(reecbdata, signal_nodes_);

        addComponent(reecb);
      }

      // Add Tgov1 governors
      for (const auto& govdata : data.gov)
      {
        auto* gov = new Tgov1<ScalarT, IdxT>(govdata);
        gov->getPorts().connect(govdata, signal_nodes_);

        addComponent(gov);
      }

      // Add GASTPTI governors
      for (const auto& gastptidata : data.gastpti)
      {
        auto* gastpti = new GastPti<ScalarT, IdxT>(gastptidata);
        gastpti->getPorts().connect(gastptidata, signal_nodes_);

        addComponent(gastpti);
      }

      // Add HYGOV governors
      for (const auto& hygovdata : data.hygov)
      {
        auto* hygov = new Hygov<ScalarT, IdxT>(hygovdata);
        hygov->getPorts().connect(hygovdata, signal_nodes_);

        addComponent(hygov);
      }

      // Add IEEEST stabilizers before exciters that consume their output during
      // initialization.
      for (const auto& stabdata : data.stabilizer)
      {
        auto* stabilizer = new Ieeest<ScalarT, IdxT>(stabdata);
        stabilizer->getPorts().connect(stabdata, signal_nodes_);
        addComponent(stabilizer);
      }

      for (const auto& excitedata : data.exciter)
      {
        IdxT bus_index = 0;
        if (excitedata.buses.contains(Ieeet1Buses::bus))
        {
          bus_index = excitedata.buses.at(Ieeet1Buses::bus);
        }

        auto* exciter = new Ieeet1<ScalarT, IdxT>(excitedata);

        auto& bus_signals = getBus(bus_index)->getSignals();
        exciter->getSignals().template attachSignalNode<Ieeet1ExternalVariables::VREAL>(
            bus_signals.template getSignalNode<BusInternalVariables::VR>());
        exciter->getSignals().template attachSignalNode<Ieeet1ExternalVariables::VIMAG>(
            bus_signals.template getSignalNode<BusInternalVariables::VI>());
        exciter->getPorts().connect(excitedata, signal_nodes_);

        addComponent(exciter);
      }

      for (const auto& excitedata : data.esdc1a)
      {
        BusT* bus = nullptr;
        if (excitedata.buses.contains(Esdc1aBuses::bus))
        {
          bus = getBus(excitedata.buses.at(Esdc1aBuses::bus));
        }

        auto* exciter = new Esdc1a<ScalarT, IdxT>(bus, excitedata);
        exciter->getPorts().connect(excitedata, signal_nodes_);

        addComponent(exciter);
      }

      for (const auto& excitedata : data.sexspti)
      {
        IdxT bus_index = 0;
        if (excitedata.buses.contains(SexsPtiBuses::bus))
        {
          bus_index = excitedata.buses.at(SexsPtiBuses::bus);
        }

        auto* exciter = new SexsPti<ScalarT, IdxT>(excitedata);

        auto& bus_signals = getBus(bus_index)->getSignals();
        exciter->getSignals().template attachSignalNode<SexsPtiExternalVariables::VREAL>(
            bus_signals.template getSignalNode<BusInternalVariables::VR>());
        exciter->getSignals().template attachSignalNode<SexsPtiExternalVariables::VIMAG>(
            bus_signals.template getSignalNode<BusInternalVariables::VI>());

        if (excitedata.signal_outputs.contains(SexsPtiSignalOutputs::efd))
        {
          IdxT           efd = excitedata.signal_outputs.at(SexsPtiSignalOutputs::efd);
          constexpr auto EFD = SexsPtiInternalVariables::EFD;
          exciter->getSignals().template assignSignalNode<EFD>(getSignal(efd));
        }

        if (excitedata.signal_inputs.contains(SexsPtiSignalInputs::vs))
        {
          IdxT           vs = excitedata.signal_inputs.at(SexsPtiSignalInputs::vs);
          constexpr auto VS = SexsPtiExternalVariables::VS;
          exciter->getSignals().template attachSignalNode<VS>(getSignal(vs));
        }
        exciter->getPorts().connect(excitedata, signal_nodes_);

        addComponent(exciter);
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
        repca->getPorts().connect(repcadata, signal_nodes_);

        addComponent(repca);
      }

      // Add constant signal sources
      for (const auto& srcdata : data.constant_source)
      {
        auto* source = new ConstantSignalSource<ScalarT, IdxT>(srcdata);
        source->getPorts().connect(srcdata, signal_nodes_);
        addComponent(source);
      }

      // Add faults
      for (const auto& faultdata : data.bus_fault)
      {
        IdxT bus_index = 0;
        if (faultdata.buses.contains(BusFaultBuses::bus))
        {
          bus_index = faultdata.buses.at(BusFaultBuses::bus);
        }
        auto* fault = new BusFault<ScalarT, IdxT>(faultdata);

        auto& bus_signals = getBus(bus_index)->getSignals();
        fault->getSignals().template attachSignalNode<BusFaultExternalVariables::VR>(
            bus_signals.template getSignalNode<BusInternalVariables::VR>());
        fault->getSignals().template attachSignalNode<BusFaultExternalVariables::VI>(
            bus_signals.template getSignalNode<BusInternalVariables::VI>());
        addFault(fault);
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
      }
    }

    /**
     * @brief Set component ID
     *
     * @note Defaults to 0 for a root model; a parent model assigns the ID
     * of a nested system like any other component.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /**
     * @brief Find the component that owns a system-local variable index.
     *
     * Component offsets are monotonically increasing, so the owner is the
     * last component whose offset is not past the index. Zero-size
     * components share their successor's offset and are never selected.
     *
     * @pre allocate() has populated component_offsets_ and
     * 0 <= local_index < size().
     */
    template <typename scalar_type, typename index_type>
    size_t SystemModel<scalar_type, index_type>::findOwningComponent(IdxT local_index) const
    {
      auto it = std::upper_bound(component_offsets_.begin(), component_offsets_.end(), local_index);
      return static_cast<size_t>(std::distance(component_offsets_.begin(), it)) - 1;
    }

    /**
     * @brief Assign a global variable index to a system-local variable.
     *
     * The write is routed to the component that owns the variable, so the
     * global index lands in the component's index-map slot, where published
     * signal nodes point. The system's own index map is kept in sync.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::setVariableIndex(IdxT local_index, IdxT global_index)
    {
      ComponentT::setVariableIndex(local_index, global_index);
      const size_t c = findOwningComponent(local_index);
      return components_[c]->setVariableIndex(local_index - component_offsets_[c], global_index);
    }

    /**
     * @brief Assign a global residual index to a system-local residual.
     *
     * The write is routed to the component that owns the residual, so the
     * global index lands in the component's index-map slot, where published
     * signal nodes point. The system's own index map is kept in sync.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::setResidualIndex(IdxT local_index, IdxT global_index)
    {
      ComponentT::setResidualIndex(local_index, global_index);
      const size_t c = findOwningComponent(local_index);
      return components_[c]->setResidualIndex(local_index - component_offsets_[c], global_index);
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
     * @note System models compose hierarchically: a nested system binds to
     * its parent's storage and allocates like any other component. Residual
     * assembly and Jacobian priming run only on the root model.
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

      component_offsets_.clear();

      IdxT offset = 0;

      for (const auto& component : components_)
      {
        component_offsets_.push_back(offset);

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

        offset += component->size();
      }

      if (offset != size_)
      {
        Log::error() << "Bound vector sizes do not match the system size\n";
        throw std::runtime_error("SystemModel allocation failed");
      }

      // Default variable and residual index mapping to local index. The
      // routed writes land in the owning component's index-map slots, where
      // published signal nodes point. A parent model can reassign global
      // indices later through the same setters.
      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, j);
        this->setResidualIndex(j, j);
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
      if (!bound_ && hasJacobian())
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
     * @brief Compute the residuals each bus and component owns.
     *
     * Components read and write their bound system-vector slices directly.
     * Every owned residual row, including the bus current balance rows, is
     * assigned here before external contributions accumulate in phase two.
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
     * @brief Fill component contributions to residuals owned elsewhere,
     * e.g. bus current balances.
     *
     * Every component fills its external residual buffer; composite
     * children recurse. Contributions are not scattered here; the root
     * model assembles them in evaluateResidual().
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateExternalResidual()
    {
      for (const auto& component : components_)
      {
        component->evaluateExternalResidual();
      }

      for (const auto& component : components_)
      {
        component->getResidual().setDataUpdated();
      }

      return 0;
    }

    /**
     * @brief Forward the root residual vector to every child.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::scatterExternalResidual(ScalarT* f_root)
    {
      for (const auto& component : components_)
      {
        component->scatterExternalResidual(f_root);
      }

      return 0;
    }

    /**
     * @brief Compute system residual vector
     *
     * Internal residuals assign every owned entry of the residual vector,
     * external residuals fill the component contribution buffers, and the
     * scatter assembles every contribution at its global row. Only a root
     * model assembles; a bound system is evaluated by its root.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateResidual()
    {
      if (bound_)
      {
        Log::error() << "A bound system model is evaluated by its root model\n";
        return 1;
      }

      evaluateInternalResidual();
      evaluateExternalResidual();
      this->scatterExternalResidual(f_.getData());

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
     * Register the bus for ID lookup and add it to the system as an ordinary
     * component.
     *
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::addBus(BusT* bus)
    {
      IdxT gridkit_bus_id                = static_cast<IdxT>(bus_lookup_.size());
      gridkit_bus_indices_[bus->busID()] = gridkit_bus_id;
      bus_lookup_.push_back(bus);
      addComponent(bus);
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

      // Keep the size current so a parent model can read it before this
      // model allocates. allocate() recomputes and checks the sum.
      size_      += component->size();
      allocated_  = false;
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
      assert((bus_lookup_[gridkit_bus_id])->busID() == bus_id);
      return bus_lookup_[gridkit_bus_id];
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
      // gridkit_component_id_ is set by System model and guarantied to be unique
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
