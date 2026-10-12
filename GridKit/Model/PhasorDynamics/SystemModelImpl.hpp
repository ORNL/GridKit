#include <algorithm>
#include <cassert>
#include <iostream>

#include <GridKit/Constants.hpp>
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

      // Add electrical buses. Each bus publishes its voltage on outlets `vr`
      // and `vi`, connected here to unnamed nodes unless the bus data names them.
      for (const auto& busdata : data.bus)
      {
        BusBase<ScalarT, IdxT>* bus = BusFactory<ScalarT, IdxT>::create(busdata, signal_nodes_);
        for (auto output : {BusSignalOutputs::vr, BusSignalOutputs::vi})
        {
          auto& port = bus->getPorts().out[output];
          if (!port.connected())
          {
            port.connect(signal_nodes_.add());
          }
        }
        addBus(bus);
      }

      // Add REGCA converters
      for (const auto& regcadata : data.regca)
      {
        auto* regca = new Regca<ScalarT, IdxT>(regcadata);
        connectTerminals(*regca, regcadata);
        regca->getPorts().connect(regcadata, signal_nodes_);

        addComponent(regca);
      }

      // Add branches
      for (const auto& branchdata : data.branch)
      {
        auto* branch = new Branch<ScalarT, IdxT>(branchdata);
        connectTerminals(*branch, branchdata);
        addComponent(branch);
      }

      // Add loads
      /// @todo Add loads to JSON parser
      for (const auto& loaddata : data.loadz)
      {
        auto* load = new LoadZ<ScalarT, IdxT>(loaddata);
        connectTerminals(*load, loaddata);
        addComponent(load);
      }

      // Add zip loads
      /// @todo Add zip loads to JSON parser
      for (const auto& loadzipdata : data.loadzip)
      {
        auto* loadzip = new LoadZIP<ScalarT, IdxT>(loadzipdata);
        connectTerminals(*loadzip, loadzipdata);
        addComponent(loadzip);
      }

      // Add GENROU generators
      for (const auto& gendata : data.genrou)
      {
        auto* gen = new Genrou<ScalarT, IdxT>(gendata);
        connectTerminals(*gen, gendata);
        gen->getPorts().connect(gendata, signal_nodes_);
        addComponent(gen);
      }

      // Add GENSAL generators
      for (const auto& gendata : data.gensal)
      {
        auto* gen = new Gensal<ScalarT, IdxT>(gendata);
        connectTerminals(*gen, gendata);
        gen->getPorts().connect(gendata, signal_nodes_);
        addComponent(gen);
      }

      // Add classical generators
      for (const auto& gendata : data.genclassical)
      {
        auto* gen = new GenClassical<ScalarT, IdxT>(gendata);
        connectTerminals(*gen, gendata);
        gen->getPorts().connect(gendata, signal_nodes_);
        addComponent(gen);
      }

      // Add REECB after its current-command and feedback producers because
      // components initialize in insertion order.
      for (const auto& reecbdata : data.reecb)
      {
        auto* reecb = new Reecb<ScalarT, IdxT>(reecbdata);
        connectTerminals(*reecb, reecbdata);
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
        auto* exciter = new Ieeet1<ScalarT, IdxT>(excitedata);
        connectTerminals(*exciter, excitedata);
        exciter->getPorts().connect(excitedata, signal_nodes_);

        addComponent(exciter);
      }

      for (const auto& excitedata : data.esdc1a)
      {
        auto* exciter = new Esdc1a<ScalarT, IdxT>(excitedata);
        connectTerminals(*exciter, excitedata);
        exciter->getPorts().connect(excitedata, signal_nodes_);

        addComponent(exciter);
      }

      for (const auto& excitedata : data.sexspti)
      {
        auto* exciter = new SexsPti<ScalarT, IdxT>(excitedata);
        connectTerminals(*exciter, excitedata);
        exciter->getPorts().connect(excitedata, signal_nodes_);

        addComponent(exciter);
      }

      // Add REPCA plant controllers after the signal producers they read at
      // initialization
      for (const auto& repcadata : data.repca)
      {
        auto* repca = new Repca<ScalarT, IdxT>(repcadata);
        connectTerminals(*repca, repcadata);
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

      // Add constant signal sources
      for (const auto& srcdata : data.function_source)
      {
        auto* source = new FunctionSignalSource<ScalarT, IdxT>(srcdata);
        source->getPorts().connect(srcdata, signal_nodes_);
        addComponent(source);
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
     * its parent's storage and allocates like any other component. Index
     * assignment, verification, Jacobian priming, and monitors run only on
     * the root model, once every model in the tree is allocated.
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

      // External slots mirror the slots of the components in insertion order.
      IdxT size_ext = 0;
      for (const auto& component : components_)
      {
        size_ext += static_cast<IdxT>(component->getExternalVariableIndices().size());
      }
      this->allocateExternalVectors(size_ext);

      if (bound_)
      {
        allocated_ = true;
        return 0;
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
     * @pre All buses and components must be allocated at this point.
     *
     * Components read their terminal bus voltages during initialization,
     * so a bus must be added before the components attached to it.
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
     * @brief Forward the monitor controller to every child.
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::collectMonitors(MonitorT& controller)
    {
      for (const auto& component : components_)
      {
        component->collectMonitors(controller);
      }
    }

    /**
     * @brief Register every monitor in the component tree with this model's
     * controller.
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::initializeMonitor()
    {
      this->collectMonitors(*monitor_);
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
      // Copy component tags into this model's tag vector. Local storage is
      // indexed locally; global indices live only in the routed index maps.
      for (size_t c = 0; c < components_.size(); ++c)
      {
        components_[c]->tagDifferentiable();
        for (IdxT j = 0; j < components_[c]->size(); ++j)
        {
          tag_[static_cast<size_t>(component_offsets_[c] + j)] = components_[c]->tag()[j];
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
     * @brief Evaluate every component, then this system's own residual.
     *
     * Each component assigns the rows it owns and fills its external
     * residual. The internal residual completes the rows this system owns
     * with those contributions; the external residual holds the
     * contributions to rows outside it, which its parent assembles in turn.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateResidual()
    {
      // A parent gathers its components; a root model gathers itself.
      if (!bound_)
      {
        gatherExternalVariables();
      }

      for (const auto& component : components_)
      {
        component->evaluateResidual();
      }

      evaluateInternalResidual();
      evaluateExternalResidual();
      f_.setDataUpdated();

      return 0;
    }

    /**
     * @brief Gather the external variables of every component, then collect
     * those that lie outside this system, i.e. its coupling to the rest of
     * the model.
     *
     * Slot k of each component, in insertion order, is one slot of this
     * system. Slots this system owns stay at INVALID_INDEX.
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::gatherExternalVariables()
    {
      for (const auto& component : components_)
      {
        component->gatherExternalVariables();
      }

      auto*  y_ext = y_ext_.getData();
      size_t slot  = 0;

      for (const auto& component : components_)
      {
        const auto  values           = component->getExternalVariables();
        const auto& variable_indices = component->getExternalVariableIndices();
        const auto& residual_indices = component->getExternalResidualIndices();

        for (size_t k = 0; k < variable_indices.size(); ++k, ++slot)
        {
          variable_indices_ext_[slot] = INVALID_INDEX<IdxT>;
          residual_indices_ext_[slot] = INVALID_INDEX<IdxT>;
          if (!owns(variable_indices_, variable_indices[k]))
          {
            y_ext[slot]                 = values[k];
            variable_indices_ext_[slot] = variable_indices[k];
          }
          if (!owns(residual_indices_, residual_indices[k]))
          {
            residual_indices_ext_[slot] = residual_indices[k];
          }
        }
      }
    }

    /**
     * @brief Complete the rows this system owns with the external residuals
     * of its components.
     *
     * Every component has been evaluated, so each row is assigned by its
     * owner before any contribution is added.
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::evaluateInternalResidual()
    {
      auto* f = f_.getData();

      for (const auto& component : components_)
      {
        const auto  values           = component->getExternalResidual();
        const auto& residual_indices = component->getExternalResidualIndices();

        for (size_t k = 0; k < residual_indices.size(); ++k)
        {
          if (owns(residual_indices_, residual_indices[k]))
          {
            f[static_cast<size_t>(residual_indices[k] - residual_indices_.front())] += values[k];
          }
        }
      }
    }

    /**
     * @brief Copy the contributions to rows outside this system into its
     * external residual.
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::evaluateExternalResidual()
    {
      auto*  f_ext = f_ext_.getData();
      size_t slot  = 0;

      for (const auto& component : components_)
      {
        const auto values = component->getExternalResidual();

        for (size_t k = 0; k < values.size(); ++k, ++slot)
        {
          if (residual_indices_ext_[slot] != INVALID_INDEX<IdxT>)
          {
            f_ext[slot] = values[k];
          }
        }
      }
    }

    /**
     * @brief Whether a global index falls in the range of an index map.
     *
     * @pre The map holds contiguous global indices, as a bound system does.
     */
    template <typename scalar_type, typename index_type>
    bool SystemModel<scalar_type, index_type>::owns(const std::vector<IdxT>& indices, IdxT index)
    {
      if (indices.empty() || index == INVALID_INDEX<IdxT> || index < indices.front())
      {
        return false;
      }
      return index - indices.front() < static_cast<IdxT>(indices.size());
    }

    /**
     * @brief Raw COO entry count of the component tree.
     */
    template <typename scalar_type, typename index_type>
    index_type SystemModel<scalar_type, index_type>::jacobianNnz()
    {
      IdxT nnz_dup = 0;
      for (const auto& component : components_)
      {
        nnz_dup += component->jacobianNnz();
      }

      return nnz_dup;
    }

    /**
     * @brief Forward the root COO arrays to every child.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::scatterJacobian(IdxT* rows, IdxT* cols, RealT* vals, IdxT& counter)
    {
      for (const auto& component : components_)
      {
        component->scatterJacobian(rows, cols, vals, counter);
      }

      return 0;
    }

    /**
     * @brief Forward the root CSR values and map to every child.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::scatterJacobianValues(RealT* vals_csr, const IdxT* map_to_csr, IdxT& counter)
    {
      for (const auto& component : components_)
      {
        component->scatterJacobianValues(vals_csr, map_to_csr, counter);
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
     * @brief Connect the bus terminals of a component to the buses named in
     * its data
     *
     */
    template <typename scalar_type, typename index_type>
    template <typename ModelT>
    void SystemModel<scalar_type, index_type>::connectTerminals(
        ModelT& model, const typename ModelT::ModelDataT& data)
    {
      for (const auto& [terminal, bus_id] : data.buses)
      {
        model.getPorts().connectTerminal(terminal, getBus(bus_id)->getPorts());
      }
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

  } // namespace PhasorDynamics
} // namespace GridKit
