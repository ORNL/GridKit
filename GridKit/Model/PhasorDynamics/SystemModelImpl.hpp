#include <cassert>
#include <iostream>
#include <set>

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

      // Resolve declarations once. Device feedback is independent of this alias graph.
      std::map<IdxT, const SignalNodeData<RealT, IdxT>*> signals;
      std::map<IdxT, std::pair<IdxT, RealT>>             resolved_delays;
      for (const auto& signal : data.signal)
      {
        signals.emplace(signal.signal_id, &signal);
      }
      for (const auto& signal : data.signal)
      {
        if (!signal.source)
        {
          if (signal.delay != 0)
          {
            throw std::invalid_argument("Signal delay requires a source");
          }
          continue;
        }
        IdxT           root  = signal.signal_id;
        RealT          delay = 0;
        std::set<IdxT> visited;
        while (signals.at(root)->source)
        {
          if (!visited.insert(root).second)
          {
            throw std::invalid_argument("Cyclic delayed signal declaration");
          }
          const auto& declaration = *signals.at(root);
          if (!std::isfinite(declaration.delay) || declaration.delay < 0)
          {
            throw std::invalid_argument("Signal delay must be finite and nonnegative");
          }
          delay += declaration.delay;
          root   = *declaration.source;
          if (!signals.contains(root))
          {
            throw std::invalid_argument("Delayed signal source does not exist");
          }
        }
        if (!std::isfinite(delay))
        {
          throw std::invalid_argument("Accumulated signal delay is not finite");
        }
        if (delay == 0 && signal.monitor)
        {
          throw std::invalid_argument("Signal monitor z requires a positive total delay");
        }
        resolved_delays.emplace(signal.signal_id, std::pair{root, delay});
      }

      // Add electrical buses
      for (const auto& busdata : data.bus)
      {
        BusBase<ScalarT, IdxT>* bus = BusFactory<ScalarT, IdxT>::create(busdata, signal_nodes_);
        addBus(bus);
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

        auto* branch = new Branch<ScalarT, IdxT>(getBus(bus1_index),
                                                 getBus(bus2_index),
                                                 branchdata);
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
        auto* load = new LoadZ<ScalarT, IdxT>(getBus(bus_index), loaddata);
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
        auto* loadzip = new LoadZIP<ScalarT, IdxT>(getBus(bus_index),
                                                   loadzipdata);
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
        auto* gen = new Genrou<ScalarT, IdxT>(getBus(bus_index), gendata);
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

        auto* gen = new Gensal<ScalarT, IdxT>(getBus(bus_index), gendata);
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
        auto* gen = new GenClassical<ScalarT, IdxT>(getBus(bus_index), gendata);
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

        auto* exciter = new Ieeet1<ScalarT, IdxT>(getBus(bus_index), excitedata);
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

        auto* exciter = new SexsPti<ScalarT, IdxT>(getBus(bus_index), excitedata);
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

      // Add constant signal sources
      for (const auto& srcdata : data.function_source)
      {
        auto* source = new FunctionSignalSource<ScalarT, IdxT>(srcdata);
        source->getPorts().connect(srcdata, signal_nodes_);
        addComponent(source);
      }

      for (const auto& signal : data.signal)
      {
        if (!signal.source)
        {
          continue;
        }
        const auto [root, delay] = resolved_delays.at(signal.signal_id);
        auto* output             = signal_nodes_[signal.signal_id];
        auto* source             = signal_nodes_[root];
        if (output->assigned())
        {
          throw std::invalid_argument("Delayed signal also has a device producer");
        }
        if (delay == 0)
        {
          output->setAssigned();
          delayed_connections_.push_back({output, source, nullptr, root, delay});
          continue;
        }
        auto [entry, inserted] = signal_histories_.try_emplace(root, SourceHistory{source, {}, {}});
        entry->second.delays.push_back(delay);
        SignalDelayData<RealT, IdxT> delay_data;
        delay_data.device_class                           = "SignalDelay";
        delay_data.disambiguation_string                  = signal.name;
        delay_data.parameters[SignalDelayParameters::tau] = delay;
        delay_data.signal_inputs[SignalDelayInputs::u]    = root;
        delay_data.signal_outputs[SignalDelayOutputs::z]  = signal.signal_id;
        if (signal.monitor)
        {
          delay_data.monitored_variables.insert(SignalDelayMonitors::z);
        }
        auto component = std::make_unique<SignalDelay<ScalarT, IdxT>>(delay_data, entry->second.history);
        component->getPorts().connect(delay_data, signal_nodes_);
        delayed_connections_.push_back({output, source, component.get(), root, delay});
        addComponent(component.release());
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

      for (const auto& connection : delayed_connections_)
      {
        if (!connection.source->linked())
        {
          throw std::invalid_argument("Delayed signal source has no linked producer");
        }
        if (!connection.component)
        {
          connection.output->link(*connection.source);
        }
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
      int                   status = 0;
      std::set<ComponentT*> deferred_delays;

      for (auto& [id, source] : signal_histories_)
      {
        source.history.reset();
        source.history.begin(time_);
      }
      for (const auto& connection : delayed_connections_)
      {
        if (connection.component)
        {
          deferred_delays.insert(connection.component);
        }
        if (!connection.component || !signal_histories_.at(connection.source_id).history.hasPrehistory())
        {
          connection.output->link(*connection.source);
        }
        else
        {
          connection.component->initialize();
        }
      }

      for (const auto& bus : buses_)
      {
        status += bus->initialize();
      }

      for (const auto& component : components_)
      {
        if (!deferred_delays.contains(component))
        {
          status += component->initialize();
        }
      }
      // Keep equilibrium aliases active for initialization writes even when
      // ordinary devices were appended after the generated delay components.
      for (auto* component : deferred_delays)
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

    template <typename ScalarT, typename IdxT>
    void SystemModel<ScalarT, IdxT>::setSignalHistoryOptions(SignalHistoryOptions<RealT> options)
    {
      for (auto& [id, source] : signal_histories_)
      {
        source.history.setOptions(options);
      }
    }

    template <typename ScalarT, typename IdxT>
    void SystemModel<ScalarT, IdxT>::initializeSignalHistory(RealT t, bool consistent)
    {
      for (auto& [id, source] : signal_histories_)
      {
        if (consistent)
        {
          source.history.initialize(t, static_cast<RealT>(source.source->read()));
        }
        else
        {
          source.history.begin(t);
        }
      }
    }

    template <typename ScalarT, typename IdxT>
    typename SystemModel<ScalarT, IdxT>::RealT
    SystemModel<ScalarT, IdxT>::assessSignalHistory(RealT begin, RealT end, const std::function<void(RealT)>& sample_state)
    {
      // The solver has restored the candidate endpoint before entering this hook.
      for (auto& [id, source] : signal_histories_)
      {
        source.history.prepare(end, static_cast<RealT>(source.source->read()));
      }
      std::vector<RealT>                    points{begin + (end - begin) / 4, begin + (end - begin) / 2, begin + 3 * (end - begin) / 4};
      std::vector<const DelayedConnection*> delayed_checks;
      for (const auto& connection : delayed_connections_)
      {
        if (!connection.component)
        {
          continue;
        }
        const RealT delay       = connection.delay;
        const bool  within_step = end - delay > begin;
        if (within_step)
        {
          points.push_back(end - delay);
        }
        // IDA controls ordinary delayed-state interpolation. Also check it
        // against history when crossing finer stored intervals or using the
        // current source iterate, where a step can hide a resolved feature.
        const auto& segments      = signal_histories_.at(connection.source_id).history.segments();
        auto        first         = std::upper_bound(segments.begin(), segments.end(), begin - delay, [](RealT q, const auto& segment)
                                      { return q < segment.end; });
        auto        last          = std::lower_bound(first, segments.end(), end - delay, [](const auto& segment, RealT q)
                                     { return segment.begin < q; });
        const bool  finer_history = std::distance(first, last) > 2;
        if (finer_history)
        {
          for (auto it = first; it != last; ++it)
          {
            const RealT left  = std::max(begin, it->begin + delay);
            const RealT right = std::min(end, it->end + delay);
            points.push_back(left + (right - left) / 2);
          }
        }
        if (within_step || finer_history)
        {
          delayed_checks.push_back(&connection);
        }
      }
      std::sort(points.begin(), points.end());
      points.erase(std::unique(points.begin(), points.end()), points.end());
      RealT error = 0;
      for (RealT t : points)
      {
        sample_state(t);
        for (const auto& [id, source] : signal_histories_)
        {
          error = std::max(error, source.history.estimate(t, static_cast<RealT>(source.source->read())));
        }
        for (const auto* connection : delayed_checks)
        {
          const auto& history = signal_histories_.at(connection->source_id).history;
          error               = std::max(error, history.estimateDelayed(t, connection->delay, static_cast<RealT>(connection->output->read())));
        }
        if (error > 1)
        {
          return error;
        }
      }
      return error;
    }

    template <typename ScalarT, typename IdxT>
    void SystemModel<ScalarT, IdxT>::acceptSignalHistory()
    {
      for (auto& [id, source] : signal_histories_)
      {
        source.history.accept();
      }
    }

    template <typename ScalarT, typename IdxT>
    void SystemModel<ScalarT, IdxT>::finishSignalHistory()
    {
      for (auto& [id, source] : signal_histories_)
      {
        source.history.maintain(source.delays);
      }
    }

    template <typename ScalarT, typename IdxT>
    void SystemModel<ScalarT, IdxT>::discardSignalHistory()
    {
      for (auto& [id, source] : signal_histories_)
      {
        source.history.discard();
      }
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

  } // namespace PhasorDynamics
} // namespace GridKit
