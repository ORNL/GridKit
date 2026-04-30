
#include <GridKit/Model/PhasorDynamics/ConnectedElementImpl.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename ScalarP, typename IdxP>
    SystemModel<ScalarP, IdxP>::SystemModel()
    {
      // Set system model tolerances
      rel_tol_         = 1e-7;
      abs_tol_         = 1e-9;
      this->max_steps_ = 2000;
    }

    template <typename ScalarP, typename IdxP>
    SystemModel<ScalarP, IdxP>::SystemModel(const ModelDataT& data)
      : monitor_(time_)
    {
      using namespace Governor;
      using namespace Exciter;
      using namespace Stabilizer;

      // Set system model tolerances
      rel_tol_         = 1e-7;
      abs_tol_         = 1e-9;
      this->max_steps_ = 2000;

      owns_components_ = true;

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

      // Add branches
      for (const auto& branchdata : data.branch)
      {
        IdxT bus1_index = 0;
        if (branchdata.ports.contains(BranchData<ScalarT, IdxT>::Ports::bus1))
        {
          bus1_index = branchdata.ports.at(BranchData<ScalarT, IdxT>::Ports::bus1);
        }

        IdxT bus2_index = 0;
        if (branchdata.ports.contains(BranchData<ScalarT, IdxT>::Ports::bus2))
        {
          bus2_index = branchdata.ports.at(BranchData<ScalarT, IdxT>::Ports::bus2);
        }

        auto* branch = new Branch<ScalarT, IdxT>(getBus(bus1_index), getBus(bus2_index), branchdata);
        addComponent(branch);
      }

      // Add loads
      /// @todo Add loads to JSON parser
      for (const auto& loaddata : data.load)
      {
        IdxT bus_index = 0;
        if (loaddata.ports.contains(LoadData<ScalarT, IdxT>::Ports::bus))
        {
          bus_index = loaddata.ports.at(LoadData<ScalarT, IdxT>::Ports::bus);
        }
        auto* load = new Load<ScalarT, IdxT>(getBus(bus_index), loaddata);
        addComponent(load);
      }

      // Add zip loads
      /// @todo Add zip loads to JSON parser
      for (const auto& loadzipdata : data.loadzip)
      {
        IdxT bus_index = 0;
        if (loadzipdata.ports.contains(LoadZIPData<ScalarT, IdxT>::Ports::bus))
        {
          bus_index = loadzipdata.ports.at(LoadZIPData<ScalarT, IdxT>::Ports::bus);
        }
        auto* loadzip = new LoadZIP<ScalarT, IdxT>(getBus(bus_index), loadzipdata);
        addComponent(loadzip);
      }

      // Add GENROU generators
      for (const auto& gendata : data.genrou)
      {
        IdxT bus_index = 0;
        if (gendata.ports.contains(GenrouData<ScalarT, IdxT>::Ports::bus))
        {
          bus_index = gendata.ports.at(GenrouData<ScalarT, IdxT>::Ports::bus);
        }

        auto* gen = new Genrou<ScalarT, IdxT>(getBus(bus_index), gendata);

        /// @todo Genrou (and likely other components) would need to name multiple
        /// signal inlets and outlets. For now we have only speed out and mechanical
        /// power in.
        if (gendata.ports.contains(GenrouData<ScalarT, IdxT>::Ports::speed))
        {
          IdxT speed = gendata.ports.at(GenrouData<ScalarT, IdxT>::Ports::speed);
          gen->getSignals().template assignSignalNode<GenrouInternalVariables::OMEGA>(getSignal(speed));
        }

        if (gendata.ports.contains(GenrouData<ScalarT, IdxT>::Ports::pmech))
        {
          IdxT pmech = gendata.ports.at(GenrouData<ScalarT, IdxT>::Ports::pmech);
          gen->getSignals().template attachSignalNode<GenrouExternalVariables::PM>(getSignal(pmech));
        }

        if (gendata.ports.contains(GenrouData<ScalarT, IdxT>::Ports::efd))
        {
          IdxT efd = gendata.ports.at(GenrouData<ScalarT, IdxT>::Ports::efd);
          gen->getSignals().template attachSignalNode<GenrouExternalVariables::EFD>(getSignal(efd));
        }

        addComponent(gen);
      }

      // Add classical generators
      for (const auto& gendata : data.genclassical)
      {
        IdxT bus_index = 0;
        if (gendata.ports.contains(GenClassicalData<ScalarT, IdxT>::Ports::bus))
        {
          bus_index = gendata.ports.at(GenClassicalData<ScalarT, IdxT>::Ports::bus);
        }
        auto* gen = new GenClassical<ScalarT, IdxT>(getBus(bus_index), gendata);
        addComponent(gen);
      }

      // Add Tgov1 governors
      for (const auto& govdata : data.gov)
      {
        auto* gov = new Tgov1<ScalarT, IdxT>(govdata);

        if (govdata.ports.contains(Tgov1Data<ScalarT, IdxT>::Ports::speed))
        {
          IdxT speed = govdata.ports.at(Tgov1Data<ScalarT, IdxT>::Ports::speed);
          gov->getSignals().template attachSignalNode<Tgov1ExternalVariables::DELTAOMEGA>(getSignal(speed));
        }

        if (govdata.ports.contains(Tgov1Data<ScalarT, IdxT>::Ports::pmech))
        {
          IdxT pmech = govdata.ports.at(Tgov1Data<ScalarT, IdxT>::Ports::pmech);
          gov->getSignals().template assignSignalNode<Tgov1InternalVariables::PM>(getSignal(pmech));
        }

        addComponent(gov);
      }

      for (const auto& excitedata : data.exciter)
      {
        IdxT bus_index = 0;
        if (excitedata.ports.contains(Ieeet1Data<ScalarT, IdxT>::Ports::bus))
        {
          bus_index = excitedata.ports.at(Ieeet1Data<ScalarT, IdxT>::Ports::bus);
        }

        auto* exciter = new Ieeet1<ScalarT, IdxT>(getBus(bus_index), excitedata);

        if (excitedata.ports.contains(Ieeet1Data<ScalarT, IdxT>::Ports::speed))
        {
          IdxT speed = excitedata.ports.at(Ieeet1Data<ScalarT, IdxT>::Ports::speed);
          exciter->getSignals().template attachSignalNode<Ieeet1ExternalVariables::OMEGA>(getSignal(speed));
        }

        if (excitedata.ports.contains(Ieeet1Data<ScalarT, IdxT>::Ports::efd))
        {
          IdxT efd = excitedata.ports.at(Ieeet1Data<ScalarT, IdxT>::Ports::efd);
          exciter->getSignals().template assignSignalNode<Ieeet1InternalVariables::EFD>(getSignal(efd));
        }

        if (excitedata.ports.contains(Ieeet1Data<ScalarT, IdxT>::Ports::vs))
        {
          IdxT vs = excitedata.ports.at(Ieeet1Data<ScalarT, IdxT>::Ports::vs);
          exciter->getSignals().template attachSignalNode<Ieeet1ExternalVariables::VS>(getSignal(vs));
        }

        addComponent(exciter);
      }

      for (const auto& excitedata : data.sexspti)
      {
        IdxT bus_index = 0;
        if (excitedata.ports.contains(SexsPtiData<ScalarT, IdxT>::Ports::bus))
        {
          bus_index = excitedata.ports.at(SexsPtiData<ScalarT, IdxT>::Ports::bus);
        }

        auto* exciter = new SexsPti<ScalarT, IdxT>(getBus(bus_index), excitedata);

        if (excitedata.ports.contains(SexsPtiData<ScalarT, IdxT>::Ports::efd))
        {
          IdxT efd = excitedata.ports.at(SexsPtiData<ScalarT, IdxT>::Ports::efd);
          exciter->getSignals().template assignSignalNode<SexsPtiInternalVariables::EFD>(getSignal(efd));
        }

        if (excitedata.ports.contains(SexsPtiData<ScalarT, IdxT>::Ports::vs))
        {
          IdxT vs = excitedata.ports.at(SexsPtiData<ScalarT, IdxT>::Ports::vs);
          exciter->getSignals().template attachSignalNode<SexsPtiExternalVariables::VS>(getSignal(vs));
        }

        addComponent(exciter);
      }

      // Add IEEEST stabilizers
      for (const auto& stabdata : data.stabilizer)
      {
        auto* stabilizer = new Ieeest<ScalarT, IdxT>(stabdata);

        if (stabdata.ports.contains(IeeestPorts::input))
        {
          IdxT input = stabdata.ports.at(IeeestPorts::input);
          stabilizer->getSignals().template attachSignalNode<IeeestExternalVariables::U>(getSignal(input));
        }

        if (stabdata.ports.contains(IeeestPorts::cutout))
        {
          IdxT cutout = stabdata.ports.at(IeeestPorts::cutout);
          stabilizer->getSignals().template attachSignalNode<IeeestExternalVariables::VCT>(getSignal(cutout));
        }

        if (stabdata.ports.contains(IeeestPorts::output))
        {
          IdxT output = stabdata.ports.at(IeeestPorts::output);
          stabilizer->getSignals().template assignSignalNode<IeeestInternalVariables::VS>(getSignal(output));
        }

        addComponent(stabilizer);
      }

      // Add faults
      for (const auto& faultdata : data.bus_fault)
      {
        IdxT bus_index = 0;
        if (faultdata.ports.contains(BusFaultData<ScalarT, IdxT>::Ports::bus))
        {
          bus_index = faultdata.ports.at(BusFaultData<ScalarT, IdxT>::Ports::bus);
        }
        auto* fault = new BusFault<ScalarT, IdxT>(getBus(bus_index), faultdata);
        addFault(fault);
      }

      for (const auto& sink : data.monitor_sink)
      {
        monitor_.addSink(sink);
      }
    }

    template <typename ScalarP, typename IdxP>
    SystemModel<ScalarP, IdxP>::~SystemModel()
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

    template <typename ScalarP, typename IdxP>
    int SystemModel<ScalarP, IdxP>::allocate()
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
      auto size = static_cast<std::size_t>(size_);
      y_.resize(size);
      yp_.resize(size);
      f_.resize(size);
      tag_.resize(size);
      variable_indices_.resize(size);
      residual_indices_.resize(size);

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

    template <typename ScalarP, typename IdxP>
    int SystemModel<ScalarP, IdxP>::verify() const
    {
      int ret = 0;

      // Verify components
      for (const auto& component : components_)
      {
        ret += component->verify();
      }

      return ret;
    }

    template <typename ScalarP, typename IdxP>
    bool SystemModel<ScalarP, IdxP>::hasJacobian()
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
        Log::warning() << "GritKit was built with Enzyme, but some models don't have Jacobians available. "
                       << "Falling back to dense Jacobians in PhasorDynamics.\n";
      }
#else
      Log::warning() << "GritKit was not built with Enzyme. "
                     << "Falling back to dense Jacobians in PhasorDynamics.\n";
#endif
      return has_jacobian;
    }

    template <typename ScalarP, typename IdxP>
    int SystemModel<ScalarP, IdxP>::initialize()
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

    template <typename ScalarP, typename IdxP>
    void SystemModel<ScalarP, IdxP>::initializeMonitor()
    {
      for (const auto* bus : buses_)
      {
        auto* mon = bus->getMonitor();
        if (mon && !mon->empty())
        {
          monitor_.addMonitor(mon);
        }
      }

      for (const auto* component : components_)
      {
        auto* mon = component->getMonitor();
        if (mon && !mon->empty())
        {
          monitor_.addMonitor(mon);
        }
      }
    }

    template <typename ScalarP, typename IdxP>
    int SystemModel<ScalarP, IdxP>::tagDifferentiable()
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

    template <typename ScalarP, typename IdxP>
    int SystemModel<ScalarP, IdxP>::evaluateResidual()
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

    template <typename ScalarP, typename IdxP>
    int SystemModel<ScalarP, IdxP>::evaluateJacobian()
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

      // J_.printMatrix("System Jacobian");

      return 0;
    }

  } // namespace PhasorDynamics
} // namespace GridKit
