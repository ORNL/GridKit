#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

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

      auto parameterValue = [](const auto& component_data,
                               auto        parameter,
                               RealT       fallback)
      {
        auto value = component_data.parameters.find(parameter);
        if (value == component_data.parameters.end())
        {
          return fallback;
        }

        return std::visit(
            [](const auto& raw_value)
            { return static_cast<RealT>(raw_value); },
            value->second);
      };

      auto mappedSignalId = [](const auto& component_data,
                               auto        signal_input) -> std::optional<IdxT>
      {
        auto signal = component_data.signal_inputs.find(signal_input);
        if (signal == component_data.signal_inputs.end())
        {
          return std::nullopt;
        }
        return signal->second;
      };

      auto addMappedInput = [this, &mappedSignalId](const auto& component_data,
                                                    auto        signal_input,
                                                    const char* input_name,
                                                    RealT       initial_value)
      {
        return addInputSignal(component_data.disambiguation_string,
                              input_name,
                              initial_value,
                              mappedSignalId(component_data, signal_input));
      };

      struct DispatchInputs
      {
        RealT p;
        RealT q;
        RealT online;
      };

      auto dispatchInputs = [](const auto& component_data,
                               RealT       p,
                               RealT       q)
      {
        DispatchInputs inputs{p, q, RealT{1.0}};
        if (component_data.initial_state)
        {
          const auto& state = *component_data.initial_state;
          inputs.p          = state.p.value_or(inputs.p);
          inputs.q          = state.q.value_or(inputs.q);
          inputs.online     = static_cast<RealT>(state.online.value_or(true));
        }
        return inputs;
      };

      owns_components_ = true;

      // Store parsed system bases before constructing data-driven components.
      this->setSystemBase(data.freq_base, data.va_base);

      std::map<IdxT, std::pair<RealT, RealT>> initial_bus_voltages;

      // Add electrical buses
      for (const auto& busdata : data.bus)
      {
        BusBase<ScalarT, IdxT>* bus = BusFactory<ScalarT, IdxT>::create(busdata);

        RealT vr = busdata.Vr0;
        RealT vi = busdata.Vi0;

        if (busdata.initial_state)
        {
          const auto& state = *busdata.initial_state;
          if (state.vr)
          {
            vr = *state.vr;
            bus->setVr(vr);
          }
          if (state.vi)
          {
            vi = *state.vi;
            bus->setVi(vi);
          }
        }

        initial_bus_voltages[busdata.bus_id] = {vr, vi};
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

        auto* branch = new Branch<ScalarT, IdxT>(
            getBus(bus1_index), getBus(bus2_index), branchdata);

        RealT tap   = parameterValue(branchdata, BranchParameters::tap, RealT{1.0});
        RealT phase = parameterValue(branchdata, BranchParameters::phase, RealT{0.0});
        RealT open  = RealT{0.0};
        if (branchdata.initial_state)
        {
          const auto& state = *branchdata.initial_state;
          tap               = state.tap.value_or(tap);
          phase             = state.phase.value_or(phase);
          open              = static_cast<RealT>(state.open.value_or(false));
        }

        auto* tap_signal   = addMappedInput(branchdata, BranchSignalInputs::tap, "tap", tap);
        auto* phase_signal = addMappedInput(branchdata, BranchSignalInputs::phase, "phase", phase);
        auto* open_signal  = addMappedInput(branchdata, BranchSignalInputs::open, "open", open);

        branch->getSignals().template attachSignalNode<BranchExternalVariables::TAP>(tap_signal);
        branch->getSignals().template attachSignalNode<BranchExternalVariables::PHASE>(phase_signal);
        branch->getSignals().template attachSignalNode<BranchExternalVariables::OPEN>(open_signal);
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

        RealT online = RealT{1.0};
        if (loaddata.initial_state)
        {
          online = static_cast<RealT>(loaddata.initial_state->online.value_or(true));
        }

        auto* online_signal = addMappedInput(
            loaddata, LoadZSignalInputs::online, "online", online);
        load->getSignals().template attachSignalNode<LoadZExternalVariables::ONLINE>(online_signal);
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
        auto* loadzip = new LoadZIP<ScalarT, IdxT>(getBus(bus_index), loadzipdata);

        const auto [vr, vi] = initial_bus_voltages.at(bus_index);
        const RealT V       = std::sqrt(vr * vr + vi * vi);
        const RealT Vnom    = parameterValue(
            loadzipdata, LoadZIPParameters::Vnom, RealT{1.0});
        const RealT alphaI = parameterValue(
            loadzipdata, LoadZIPParameters::alphaI, RealT{0.0});
        const RealT alphaP = parameterValue(
            loadzipdata, LoadZIPParameters::alphaP, RealT{0.0});
        const RealT alphaZ = RealT{1.0} - alphaI - alphaP;
        const RealT ratio  = V / Vnom;
        const RealT legacy_dispatch_factor =
            alphaZ * ratio * ratio + alphaI * ratio + alphaP;

        auto inputs = dispatchInputs(
            loadzipdata,
            -parameterValue(loadzipdata, LoadZIPParameters::Pnom, RealT{0.0})
                * legacy_dispatch_factor,
            -parameterValue(loadzipdata, LoadZIPParameters::Qnom, RealT{0.0})
                * legacy_dispatch_factor);

        auto* p_signal = addMappedInput(
            loadzipdata, LoadZIPSignalInputs::p, "p", inputs.p);
        auto* q_signal = addMappedInput(
            loadzipdata, LoadZIPSignalInputs::q, "q", inputs.q);
        auto* online_signal = addMappedInput(
            loadzipdata, LoadZIPSignalInputs::online, "online", inputs.online);

        loadzip->getSignals().template attachSignalNode<LoadZIPExternalVariables::P>(p_signal);
        loadzip->getSignals().template attachSignalNode<LoadZIPExternalVariables::Q>(q_signal);
        loadzip->getSignals().template attachSignalNode<LoadZIPExternalVariables::ONLINE>(online_signal);
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

        /// @todo Genrou (and likely other components) would need to name multiple
        /// signal inlets and outlets. For now we have only speed out and mechanical
        /// power in.
        if (gendata.signal_outputs.contains(GenrouSignalOutputs::speed))
        {
          IdxT           speed = gendata.signal_outputs.at(GenrouSignalOutputs::speed);
          constexpr auto OMEGA = GenrouInternalVariables::OMEGA;
          gen->getSignals().template assignSignalNode<OMEGA>(getSignal(speed));
        }

        if (gendata.signal_inputs.contains(GenrouSignalInputs::pmech))
        {
          IdxT           pmech = gendata.signal_inputs.at(GenrouSignalInputs::pmech);
          constexpr auto PM    = GenrouExternalVariables::PM;
          gen->getSignals().template attachSignalNode<PM>(getSignal(pmech));
        }

        if (gendata.signal_inputs.contains(GenrouSignalInputs::efd))
        {
          IdxT           efd = gendata.signal_inputs.at(GenrouSignalInputs::efd);
          constexpr auto EFD = GenrouExternalVariables::EFD;
          gen->getSignals().template attachSignalNode<EFD>(getSignal(efd));
        }

        auto inputs = dispatchInputs(
            gendata,
            parameterValue(gendata, GenrouParameters::p0, RealT{0.0}),
            parameterValue(gendata, GenrouParameters::q0, RealT{0.0}));

        auto* p_signal = addMappedInput(
            gendata, GenrouSignalInputs::p, "p", inputs.p);
        auto* q_signal = addMappedInput(
            gendata, GenrouSignalInputs::q, "q", inputs.q);
        auto* online_signal = addMappedInput(
            gendata, GenrouSignalInputs::online, "online", inputs.online);

        gen->getSignals().template attachSignalNode<GenrouExternalVariables::P>(p_signal);
        gen->getSignals().template attachSignalNode<GenrouExternalVariables::Q>(q_signal);
        gen->getSignals().template attachSignalNode<GenrouExternalVariables::ONLINE>(online_signal);

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

        if (gendata.signal_outputs.contains(GensalSignalOutputs::speed))
        {
          IdxT           speed = gendata.signal_outputs.at(GensalSignalOutputs::speed);
          constexpr auto OMEGA = GensalInternalVariables::OMEGA;
          gen->getSignals().template assignSignalNode<OMEGA>(getSignal(speed));
        }

        if (gendata.signal_inputs.contains(GensalSignalInputs::pmech))
        {
          IdxT           pmech = gendata.signal_inputs.at(GensalSignalInputs::pmech);
          constexpr auto PM    = GensalExternalVariables::PM;
          gen->getSignals().template attachSignalNode<PM>(getSignal(pmech));
        }

        if (gendata.signal_inputs.contains(GensalSignalInputs::efd))
        {
          IdxT           efd = gendata.signal_inputs.at(GensalSignalInputs::efd);
          constexpr auto EFD = GensalExternalVariables::EFD;
          gen->getSignals().template attachSignalNode<EFD>(getSignal(efd));
        }

        auto inputs = dispatchInputs(
            gendata,
            parameterValue(gendata, GensalParameters::p0, RealT{0.0}),
            parameterValue(gendata, GensalParameters::q0, RealT{0.0}));

        auto* p_signal = addMappedInput(
            gendata, GensalSignalInputs::p, "p", inputs.p);
        auto* q_signal = addMappedInput(
            gendata, GensalSignalInputs::q, "q", inputs.q);
        auto* online_signal = addMappedInput(
            gendata, GensalSignalInputs::online, "online", inputs.online);

        gen->getSignals().template attachSignalNode<GensalExternalVariables::P>(p_signal);
        gen->getSignals().template attachSignalNode<GensalExternalVariables::Q>(q_signal);
        gen->getSignals().template attachSignalNode<GensalExternalVariables::ONLINE>(online_signal);

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

        auto inputs = dispatchInputs(
            gendata,
            parameterValue(gendata, GenClassicalParameters::p0, RealT{0.0}),
            parameterValue(gendata, GenClassicalParameters::q0, RealT{0.0}));

        auto* p_signal = addMappedInput(
            gendata, GenClassicalSignalInputs::p, "p", inputs.p);
        auto* q_signal = addMappedInput(
            gendata, GenClassicalSignalInputs::q, "q", inputs.q);
        auto* online_signal = addMappedInput(
            gendata, GenClassicalSignalInputs::online, "online", inputs.online);

        gen->getSignals().template attachSignalNode<GenClassicalExternalVariables::P>(p_signal);
        gen->getSignals().template attachSignalNode<GenClassicalExternalVariables::Q>(q_signal);
        gen->getSignals().template attachSignalNode<GenClassicalExternalVariables::ONLINE>(online_signal);
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

      // Add faults
      for (const auto& faultdata : data.bus_fault)
      {
        IdxT bus_index = 0;
        if (faultdata.buses.contains(BusFaultBuses::bus))
        {
          bus_index = faultdata.buses.at(BusFaultBuses::bus);
        }
        auto* fault = new BusFault<ScalarT, IdxT>(getBus(bus_index), faultdata);
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

      for (const auto& [key, signals] : input_signals_)
      {
        for (const auto* signal : signals)
        {
          if (!signal->linked())
          {
            throw std::runtime_error(
                "SystemModel input signal has no backing storage: "
                + key.first + "." + key.second);
          }
        }
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
     * that these are simply concateneted in the global system.
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
     * add to those values by in-place adition. This is why (for now) bus
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
     * @brief Add a model input signal and register its semantic endpoint.
     */
    template <typename scalar_type, typename index_type>
    typename SystemModel<scalar_type, index_type>::SignalT*
    SystemModel<scalar_type, index_type>::addInputSignal(
        const std::string&  device_id,
        const std::string&  input_name,
        RealT               initial_value,
        std::optional<IdxT> signal_id)
    {
      InputKey key{device_id, input_name};
      SignalT* signal = nullptr;
      if (signal_id)
      {
        // An explicit graph connection supplies the value and takes
        // precedence over state and legacy fallbacks.
        signal = getSignal(*signal_id);
      }
      else
      {
        auto input = std::make_unique<InputSignal>(
            static_cast<ScalarT>(initial_value));
        signal = &input->node;
        owned_input_signals_.push_back(std::move(input));
      }

      input_signals_[std::move(key)].push_back(signal);
      return signal;
    }

    /**
     * @brief Set a named model input signal.
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::setInput(
        const std::string& device_id,
        const std::string& input_name,
        RealT              value)
    {
      auto& signals = input_signals_.at({device_id, input_name});
      for (auto* signal : signals)
      {
        if (!signal->linked())
        {
          throw std::logic_error("SystemModel input signal has no backing storage");
        }
      }
      for (auto* signal : signals)
      {
        signal->init(static_cast<ScalarT>(value));
      }
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
