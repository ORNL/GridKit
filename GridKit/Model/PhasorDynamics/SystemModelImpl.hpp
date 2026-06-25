#include <cassert>
#include <iostream>

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

      owns_components_ = true;

      // Store parsed system bases before constructing data-driven components.
      this->setSystemBase(data.freq_base, data.va_base);

      // Add electrical buses
      for (const auto& busdata : data.bus)
      {
        BusBase<ScalarT, IdxT>* bus = BusFactory<ScalarT, IdxT>::create(busdata);
        addBus(bus);
      }

      /// @todo Signal data needs to be populated by the parser.
      /// See TODOs in SystemModelDataJSONParser
      for (const auto& signaldata : data.signal)
      {
        SignalNode<ScalarT, IdxT>* signal = new SignalNode<ScalarT, IdxT>(signaldata);
        addSignal(signal);
      }

      // Add bus-to-signal adapters
      for (const auto& adapterdata : data.adapter)
      {
        using AdapterPorts = BusToSignalAdapterPorts;
        IdxT bus_index     = 0;
        if (adapterdata.ports.contains(AdapterPorts::bus))
        {
          bus_index = adapterdata.ports.at(AdapterPorts::bus);
        }

        auto* adapter = new BusToSignalAdapter<ScalarT, IdxT>(getBus(bus_index));

        if (adapterdata.ports.contains(AdapterPorts::vr))
        {
          IdxT           vr    = adapterdata.ports.at(AdapterPorts::vr);
          constexpr auto VREAL = BusToSignalAdapterInternalVariables::VREAL;
          adapter->getSignals().template assignSignalNode<VREAL>(getSignal(vr));
        }

        if (adapterdata.ports.contains(AdapterPorts::vi))
        {
          IdxT           vi    = adapterdata.ports.at(AdapterPorts::vi);
          constexpr auto VIMAG = BusToSignalAdapterInternalVariables::VIMAG;
          adapter->getSignals().template assignSignalNode<VIMAG>(getSignal(vi));
        }

        if (adapterdata.ports.contains(AdapterPorts::ir))
        {
          IdxT           ir    = adapterdata.ports.at(AdapterPorts::ir);
          constexpr auto IREAL = BusToSignalAdapterExternalVariables::IREAL;
          adapter->getSignals().template attachSignalNode<IREAL>(getSignal(ir));
        }

        if (adapterdata.ports.contains(AdapterPorts::ii))
        {
          IdxT           ii    = adapterdata.ports.at(AdapterPorts::ii);
          constexpr auto IIMAG = BusToSignalAdapterExternalVariables::IIMAG;
          adapter->getSignals().template attachSignalNode<IIMAG>(getSignal(ii));
        }

        addComponent(adapter);
      }

      // Add branches
      for (const auto& branchdata : data.branch)
      {
        IdxT bus1_index = 0;
        if (branchdata.ports.contains(BranchPorts::bus1))
        {
          bus1_index = branchdata.ports.at(BranchPorts::bus1);
        }

        IdxT bus2_index = 0;
        if (branchdata.ports.contains(BranchPorts::bus2))
        {
          bus2_index = branchdata.ports.at(BranchPorts::bus2);
        }

        auto* branch = new Branch<ScalarT, IdxT>(
            getBus(bus1_index), getBus(bus2_index), branchdata);
        addComponent(branch);
      }

      // Add loads
      /// @todo Add loads to JSON parser
      for (const auto& loaddata : data.loadz)
      {
        IdxT bus_index = 0;
        if (loaddata.ports.contains(LoadZPorts::bus))
        {
          bus_index = loaddata.ports.at(LoadZPorts::bus);
        }
        auto* load = new LoadZ<ScalarT, IdxT>(getBus(bus_index), loaddata);
        addComponent(load);
      }

      // Add zip loads
      /// @todo Add zip loads to JSON parser
      for (const auto& loadzipdata : data.loadzip)
      {
        IdxT bus_index = 0;
        if (loadzipdata.ports.contains(LoadZIPPorts::bus))
        {
          bus_index = loadzipdata.ports.at(LoadZIPPorts::bus);
        }
        auto* loadzip = new LoadZIP<ScalarT, IdxT>(getBus(bus_index), loadzipdata);
        addComponent(loadzip);
      }

      // Add GENROU generators
      for (const auto& gendata : data.genrou)
      {
        IdxT bus_index = 0;
        if (gendata.ports.contains(GenrouPorts::bus))
        {
          bus_index = gendata.ports.at(GenrouPorts::bus);
        }

        auto* gen = new Genrou<ScalarT, IdxT>(getBus(bus_index), gendata);

        /// @todo Genrou (and likely other components) would need to name multiple
        /// signal inlets and outlets. For now we have only speed out and mechanical
        /// power in.
        if (gendata.ports.contains(GenrouPorts::speed))
        {
          IdxT           speed = gendata.ports.at(GenrouPorts::speed);
          constexpr auto OMEGA = GenrouInternalVariables::OMEGA;
          gen->getSignals().template assignSignalNode<OMEGA>(getSignal(speed));
        }

        if (gendata.ports.contains(GenrouPorts::pmech))
        {
          IdxT           pmech = gendata.ports.at(GenrouPorts::pmech);
          constexpr auto PM    = GenrouExternalVariables::PM;
          gen->getSignals().template attachSignalNode<PM>(getSignal(pmech));
        }

        if (gendata.ports.contains(GenrouPorts::efd))
        {
          IdxT           efd = gendata.ports.at(GenrouPorts::efd);
          constexpr auto EFD = GenrouExternalVariables::EFD;
          gen->getSignals().template attachSignalNode<EFD>(getSignal(efd));
        }

        addComponent(gen);
      }

      // Add GENSAL generators
      for (const auto& gendata : data.gensal)
      {
        IdxT bus_index = 0;
        if (gendata.ports.contains(GensalPorts::bus))
        {
          bus_index = gendata.ports.at(GensalPorts::bus);
        }

        auto* gen = new Gensal<ScalarT, IdxT>(getBus(bus_index), gendata);

        if (gendata.ports.contains(GensalPorts::speed))
        {
          IdxT           speed = gendata.ports.at(GensalPorts::speed);
          constexpr auto OMEGA = GensalInternalVariables::OMEGA;
          gen->getSignals().template assignSignalNode<OMEGA>(getSignal(speed));
        }

        if (gendata.ports.contains(GensalPorts::pmech))
        {
          IdxT           pmech = gendata.ports.at(GensalPorts::pmech);
          constexpr auto PM    = GensalExternalVariables::PM;
          gen->getSignals().template attachSignalNode<PM>(getSignal(pmech));
        }

        if (gendata.ports.contains(GensalPorts::efd))
        {
          IdxT           efd = gendata.ports.at(GensalPorts::efd);
          constexpr auto EFD = GensalExternalVariables::EFD;
          gen->getSignals().template attachSignalNode<EFD>(getSignal(efd));
        }

        addComponent(gen);
      }

      // Add classical generators
      for (const auto& gendata : data.genclassical)
      {
        IdxT bus_index = 0;
        if (gendata.ports.contains(GenClassicalPorts::bus))
        {
          bus_index = gendata.ports.at(GenClassicalPorts::bus);
        }
        auto* gen = new GenClassical<ScalarT, IdxT>(getBus(bus_index), gendata);
        addComponent(gen);
      }

      // Add Tgov1 governors
      for (const auto& govdata : data.gov)
      {
        auto* gov = new Tgov1<ScalarT, IdxT>(govdata);

        if (govdata.ports.contains(Tgov1Ports::speed))
        {
          IdxT           speed      = govdata.ports.at(Tgov1Ports::speed);
          constexpr auto DELTAOMEGA = Tgov1ExternalVariables::DELTAOMEGA;
          gov->getSignals().template attachSignalNode<DELTAOMEGA>(getSignal(speed));
        }

        if (govdata.ports.contains(Tgov1Ports::pmech))
        {
          IdxT           pmech = govdata.ports.at(Tgov1Ports::pmech);
          constexpr auto PM    = Tgov1InternalVariables::PM;
          gov->getSignals().template assignSignalNode<PM>(getSignal(pmech));
        }

        addComponent(gov);
      }

      for (const auto& excitedata : data.exciter)
      {
        IdxT bus_index = 0;
        if (excitedata.ports.contains(Ieeet1Ports::bus))
        {
          bus_index = excitedata.ports.at(Ieeet1Ports::bus);
        }

        auto* exciter = new Ieeet1<ScalarT, IdxT>(getBus(bus_index), excitedata);

        if (excitedata.ports.contains(Ieeet1Ports::speed))
        {
          IdxT           speed = excitedata.ports.at(Ieeet1Ports::speed);
          constexpr auto OMEGA = Ieeet1ExternalVariables::OMEGA;
          exciter->getSignals().template attachSignalNode<OMEGA>(getSignal(speed));
        }

        if (excitedata.ports.contains(Ieeet1Ports::efd))
        {
          IdxT           efd = excitedata.ports.at(Ieeet1Ports::efd);
          constexpr auto EFD = Ieeet1InternalVariables::EFD;
          exciter->getSignals().template assignSignalNode<EFD>(getSignal(efd));
        }

        if (excitedata.ports.contains(Ieeet1Ports::vs))
        {
          IdxT           vs = excitedata.ports.at(Ieeet1Ports::vs);
          constexpr auto VS = Ieeet1ExternalVariables::VS;
          exciter->getSignals().template attachSignalNode<VS>(getSignal(vs));
        }

        addComponent(exciter);
      }

      for (const auto& excitedata : data.sexspti)
      {
        IdxT bus_index = 0;
        if (excitedata.ports.contains(SexsPtiPorts::bus))
        {
          bus_index = excitedata.ports.at(SexsPtiPorts::bus);
        }

        auto* exciter = new SexsPti<ScalarT, IdxT>(getBus(bus_index), excitedata);

        if (excitedata.ports.contains(SexsPtiPorts::efd))
        {
          IdxT           efd = excitedata.ports.at(SexsPtiPorts::efd);
          constexpr auto EFD = SexsPtiInternalVariables::EFD;
          exciter->getSignals().template assignSignalNode<EFD>(getSignal(efd));
        }

        if (excitedata.ports.contains(SexsPtiPorts::vs))
        {
          IdxT           vs = excitedata.ports.at(SexsPtiPorts::vs);
          constexpr auto VS = SexsPtiExternalVariables::VS;
          exciter->getSignals().template attachSignalNode<VS>(getSignal(vs));
        }

        addComponent(exciter);
      }

      // Add IEEEST stabilizers
      for (const auto& stabdata : data.stabilizer)
      {
        auto* stabilizer = new Ieeest<ScalarT, IdxT>(stabdata);

        if (stabdata.ports.contains(IeeestPorts::input))
        {
          IdxT           input = stabdata.ports.at(IeeestPorts::input);
          constexpr auto U     = IeeestExternalVariables::U;
          stabilizer->getSignals().template attachSignalNode<U>(getSignal(input));
        }

        if (stabdata.ports.contains(IeeestPorts::output))
        {
          IdxT           output = stabdata.ports.at(IeeestPorts::output);
          constexpr auto VSS    = IeeestInternalVariables::VSS;
          stabilizer->getSignals().template assignSignalNode<VSS>(getSignal(output));
        }

        addComponent(stabilizer);
      }

      // Add faults
      for (const auto& faultdata : data.bus_fault)
      {
        IdxT bus_index = 0;
        if (faultdata.ports.contains(BusFaultPorts::bus))
        {
          bus_index = faultdata.ports.at(BusFaultPorts::bus);
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

        for (auto signal : signals_)
        {
          delete signal;
        }
      }
      if (csr_jac_ != nullptr)
      {
        delete csr_jac_;
        csr_jac_ = nullptr;
      }
      if (map_to_csr_ != nullptr)
      {
        delete[] map_to_csr_;
        map_to_csr_ = nullptr;
      }
    }

    /**
     * @brief Set component ID
     *
     * @note Should default to 0. The system model could be used as a
     * component in a larger system that would need to set this value.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /**
     * @brief Allocate buses, components, and system objects.
     *
     * This method first allocates bus objects, then component objects,
     * and computes system size (number of unknowns). Once the size is
     * computed, system global objects are allocated.
     *
     * @post size_ >= 1
     *
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::allocate()
    {
      size_ = 0;

      // Allocate all buses
      for (const auto& bus : buses_)
      {
        bus->allocate();
        for (IdxT j = 0; j < bus->size(); ++j)
        {
          bus->setVariableIndex(j, size_ + j);
          bus->setResidualIndex(j, size_ + j);
        }
        size_ += bus->size();
      }

      // Allocate all components
      for (const auto& component : components_)
      {
        component->allocate();
        for (IdxT j = 0; j < component->size(); ++j)
        {
          component->setVariableIndex(j, size_ + j);
          component->setResidualIndex(j, size_ + j);
        }
        size_ += component->size();
      }

      // Allocate global vectors
      y_.resize(size_);
      yp_.resize(size_);
      f_.resize(size_);
      tag_.resize(size_);
      abs_tol_.resize(size_);
      variable_indices_.resize(size_);
      residual_indices_.resize(size_);

      // Default variable and residual index mapping to local index
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

      // Start variable monitors
      initializeMonitor();
      startMonitor();

      // Perform an initial Jacobian evaluation for sparse Jacobians, such that
      // the dynamic solver can querry the NNZ value when it is configured.
      // @todo Replace with a sparsity analysis that sets the NNZ and allocates the Jacobian
      // without needing the Jacobian values.
      if (hasJacobian())
      {
        initialize();
        evaluateResidual();
        evaluateJacobian();
      }

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
     * @brief Check components for Jacobian availability
     *
     * @return true
     * @return false
     */
    template <typename scalar_type, typename index_type>
    bool SystemModel<scalar_type, index_type>::hasJacobian()
    {
      bool has_jacobian = false;
#ifdef GRIDKIT_ENABLE_ENZYME
      has_jacobian = true;
      for (const auto& component : components_)
      {
        has_jacobian = has_jacobian && component->hasJacobian();
      }

      if (!has_jacobian)
      {
        Log::warning() << "GridKit was built with Enzyme, but some models "
                          "don't have Jacobians available. "
                          "Falling back to dense Jacobians in PhasorDynamics.\n";
      }
#else
      Log::warning() << "GridKit was not built with Enzyme. "
                     << "Falling back to dense Jacobians in PhasorDynamics.\n";
#endif
      return has_jacobian;
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
     * @todo Implement writting to system vectors in a thread-safe way.
     *
     * @note Currently assuming each component stores variables contiguously in memory and
     * that these are simply concateneted in the global system.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::initialize()
    {
      for (const auto& bus : buses_)
      {
        bus->initialize();
      }

      for (const auto& bus : buses_)
      {
        for (IdxT j = 0; j < bus->size(); ++j)
        {
          y_[bus->getVariableIndex(j)]  = bus->y()[j];
          yp_[bus->getVariableIndex(j)] = bus->yp()[j];
        }
      }

      // Initialize components
      for (const auto& component : components_)
      {
        component->initialize();
      }

      for (const auto& component : components_)
      {
        for (IdxT j = 0; j < component->size(); ++j)
        {
          y_[component->getVariableIndex(j)]  = component->y()[j];
          yp_[component->getVariableIndex(j)] = component->yp()[j];
        }
      }

      return 0;
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
      IdxT offset = 0;
      for (const auto& bus : buses_)
      {
        bus->setAbsoluteTolerance(rel_tol);
        for (IdxT j = 0; j < bus->size(); ++j)
        {
          abs_tol_[offset + j] = bus->absoluteTolerance()[j];
        }
        offset += bus->size();
      }

      for (const auto& component : components_)
      {
        component->setAbsoluteTolerance(rel_tol);
        for (IdxT j = 0; j < component->size(); ++j)
        {
          abs_tol_[offset + j] = component->absoluteTolerance()[j];
        }
        offset += component->size();
      }

      return 0;
    }

    /**
     * @brief Compute system residual vector
     *
     * First, update bus and component variables from the system solution
     * vector. Next, evaluate residuals in buses and components, and
     * then copy values to the global residual vector.
     *
     * @warning Residuals must be computed for buses, before component
     * residuals are computed. Buses own residuals for currents
     * Ir and Ii, but the contributions to these residuals come
     * from components. Buses assign their residual values, while components
     * add to those values by in-place adition. This is why (for now) bus
     * residuals need to be computed first.
     *
     * @todo Here, components write to local values, which are then copied
     * to global system vectors. Make components write to the system
     * vectors directly.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateResidual()
    {
      updateVariables();

      for (const auto& bus : buses_)
      {
        bus->evaluateResidual();
      }

      for (const auto& component : components_)
      {
        component->evaluateResidual();
      }

      // Update residual vector
      for (const auto& bus : buses_)
      {
        for (IdxT j = 0; j < bus->size(); ++j)
        {
          f_[bus->getResidualIndex(j)] = bus->getResidual()[j];
        }
      }

      for (const auto& component : components_)
      {
        for (IdxT j = 0; j < component->size(); ++j)
        {
          f_[component->getResidualIndex(j)] = component->getResidual()[j];
        }
      }

      return 0;
    }

    /**
     * @brief Evaluate system Jacobian.
     *
     * First, initialize bus Jacobians to 0.
     * Then, evaluate component Jacobians (internal block and bus Jacobian contributions).
     * Once component Jacobians are evaluated, store the result in the system Jacobian.
     * Finally, store bus Jacobians into the system Jacobian after all component have added their
     * contributions.
     *
     * @todo split the initial assembly from updating values. This will the
     * slow otherwise.
     *
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateJacobian()
    {
      // Initialize bus Jacobians
      for (const auto& bus : buses_)
      {
        bus->evaluateJacobian();
      }

      // Evaluate component Jacobians and update bus Jacobians
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
          auto component_jacobian = component->getJacobian();

          std::tuple<std::vector<IdxT>&, std::vector<IdxT>&, std::vector<RealT>&> component_jacobian_entries = component_jacobian.getEntries(false);
          const auto [rows, columns, values]                                                                 = component_jacobian_entries;
          for (size_t i = 0; i < rows.size(); ++i)
          {
            ++nnz_dup;
          }
        }
        for (const auto& bus : buses_)
        {
          auto bus_jacobian = bus->getJacobian();

          std::tuple<std::vector<IdxT>&, std::vector<IdxT>&, std::vector<RealT>&> bus_jacobian_entries = bus_jacobian.getEntries(false);
          const auto [rows, columns, values]                                                           = bus_jacobian_entries;
          for (size_t i = 0; i < rows.size(); ++i)
          {
            ++nnz_dup;
          }
        }

        // Allocate COO triplet arrays (we own these until we hand off to CsrMatrix)
        IdxT*  rows_dup = new IdxT[nnz_dup];
        IdxT*  cols_dup = new IdxT[nnz_dup];
        RealT* vals_dup = new RealT[nnz_dup];

        IdxT counter = 0;
        for (const auto& component : components_)
        {
          auto component_jacobian = component->getJacobian();

          std::tuple<std::vector<IdxT>&, std::vector<IdxT>&, std::vector<RealT>&> component_jacobian_entries = component_jacobian.getEntries(false);
          const auto [rows, columns, values]                                                                 = component_jacobian_entries;
          for (size_t i = 0; i < rows.size(); ++i)
          {
            rows_dup[counter] = rows[i];
            cols_dup[counter] = columns[i];
            vals_dup[counter] = values[i];
            counter++;
          }
        }
        for (const auto& bus : buses_)
        {
          auto bus_jacobian = bus->getJacobian();

          std::tuple<std::vector<IdxT>&, std::vector<IdxT>&, std::vector<RealT>&> bus_jacobian_entries = bus_jacobian.getEntries(false);
          const auto [rows, columns, values]                                                           = bus_jacobian_entries;
          for (size_t i = 0; i < rows.size(); ++i)
          {
            rows_dup[counter] = rows[i];
            cols_dup[counter] = columns[i];
            vals_dup[counter] = values[i];
            counter++;
          }
        }

        // Build the system COO Jacobian
        LinearAlgebra::CooMatrix<RealT, IdxT> jac(size_, size_, nnz_dup, &rows_dup, &cols_dup, &vals_dup);

        // Populate CSR data with sort and deduplicate
        IdxT* row_ptrs = jac.getCsrRowData();

        // Deduplicated nnz
        nnz_ = jac.getNnz();

        // Allocate cols/vals with deduplicated nnz
        IdxT*  cols = new IdxT[nnz_];
        RealT* vals = new RealT[nnz_];

        std::copy(jac.getColData(), jac.getColData() + nnz_, cols);
        std::copy(jac.getValues(), jac.getValues() + nnz_, vals);

        // Create the CSR Jacobian
        csr_jac_ = new CsrMatrixT(size_, size_, nnz_, &row_ptrs, &cols, &vals);

        const IdxT* map_to_sorted = jac.getMapToSorted();
        const IdxT* map_to_dedup  = jac.getMapToDeduplicated();

        // Build a mappping from original COO index to CSR index
        map_to_csr_ = new IdxT[nnz_dup];
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

        // Update CSR values from component and bus Jacobians
        IdxT counter = 0;
        for (const auto& component : components_)
        {
          auto component_jacobian = component->getJacobian();

          std::tuple<std::vector<IdxT>&, std::vector<IdxT>&, std::vector<RealT>&> component_jacobian_entries = component_jacobian.getEntries(false);
          const auto [rows, columns, values]                                                                 = component_jacobian_entries;
          for (size_t i = 0; i < rows.size(); ++i)
          {
            vals[map_to_csr_[counter]] += values[i];
            ++counter;
          }
        }
        for (const auto& bus : buses_)
        {
          auto bus_jacobian = bus->getJacobian();

          std::tuple<std::vector<IdxT>&, std::vector<IdxT>&, std::vector<RealT>&> bus_jacobian_entries = bus_jacobian.getEntries(false);
          const auto [rows, columns, values]                                                           = bus_jacobian_entries;
          for (size_t i = 0; i < rows.size(); ++i)
          {
            vals[map_to_csr_[counter]] += values[i];
            ++counter;
          }
        }
      }

      // std::cout << "System Jacobian\n";
      // csr_jac_->print(std::cout);

      return 0;
    }

    /**
     * @brief Update variables in buses and components
     */
    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::updateVariables()
    {
      for (const auto& bus : buses_)
      {
        for (IdxT j = 0; j < bus->size(); ++j)
        {
          bus->y()[j]  = y_[bus->getVariableIndex(j)];
          bus->yp()[j] = yp_[bus->getVariableIndex(j)];
        }
      }
      for (const auto& component : components_)
      {
        for (IdxT j = 0; j < component->size(); ++j)
        {
          component->y()[j]  = y_[component->getVariableIndex(j)];
          component->yp()[j] = yp_[component->getVariableIndex(j)];
        }
      }
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

      updateVariables();
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
      component->setSystemBase(this->freq_system_base_,
                               this->va_system_base_);
      components_.push_back(component);
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
