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

        auto* adapter = new BusToSignalAdapter<ScalarT, IdxT>(getBus(bus_index),
                                                              adapterdata,
                                                              signal_nodes_);
        addComponent(adapter);
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
        auto* gen = new Genrou<ScalarT, IdxT>(getBus(bus_index),
                                              gendata,
                                              signal_nodes_);
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

        auto* gen = new Gensal<ScalarT, IdxT>(getBus(bus_index),
                                              gendata,
                                              signal_nodes_);
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
        addComponent(gen);
      }

      // Add Tgov1 governors
      for (const auto& govdata : data.gov)
      {
        auto* gov = new Tgov1<ScalarT, IdxT>(govdata, signal_nodes_);
        addComponent(gov);
      }

      for (const auto& excitedata : data.exciter)
      {
        IdxT bus_index = 0;
        if (excitedata.buses.contains(Ieeet1Buses::bus))
        {
          bus_index = excitedata.buses.at(Ieeet1Buses::bus);
        }

        auto* exciter = new Ieeet1<ScalarT, IdxT>(getBus(bus_index),
                                                  excitedata,
                                                  signal_nodes_);
        addComponent(exciter);
      }

      for (const auto& excitedata : data.sexspti)
      {
        IdxT bus_index = 0;
        if (excitedata.buses.contains(SexsPtiBuses::bus))
        {
          bus_index = excitedata.buses.at(SexsPtiBuses::bus);
        }

        auto* exciter = new SexsPti<ScalarT, IdxT>(getBus(bus_index),
                                                   excitedata,
                                                   signal_nodes_);
        addComponent(exciter);
      }

      // Add IEEEST stabilizers
      for (const auto& stabdata : data.stabilizer)
      {
        auto* stabilizer = new Ieeest<ScalarT, IdxT>(stabdata, signal_nodes_);
        addComponent(stabilizer);
      }

      // Add constant signal sources
      for (const auto& srcdata : data.constant_source)
      {
        auto* source = new ConstantSignalSource<ScalarT, IdxT>(srcdata,
                                                               signal_nodes_);
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

      for (const auto& bus : buses_)
      {
        has_jacobian = has_jacobian && bus->hasJacobian();
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
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::evaluateJacobian()
    {
      // Initialize bus Jacobians
      for (const auto& bus : buses_)
      {
        bus->evaluateJacobian();
      }

      // Evaluate component Jacobians, including contribution to the bus Jacobians
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
          else
          {
            Log::warning() << "A component has returned a nullptr Jacobian.\n";
          }
        }

        for (const auto& bus : buses_)
        {
          auto bus_jacobian = bus->getCooJacobian();

          if (bus_jacobian != nullptr)
          {
            nnz_dup += bus_jacobian->getNnz();
          }
          else
          {
            Log::warning() << "A bus has returned a nullptr Jacobian.\n";
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
          else
          {
            Log::warning() << "A component has returned a nullptr Jacobian.\n";
          }
        }

        for (const auto& bus : buses_)
        {
          auto bus_jacobian = bus->getCooJacobian();

          if (bus_jacobian != nullptr)
          {
            const IdxT*  rows    = bus_jacobian->getRowData();
            const IdxT*  columns = bus_jacobian->getColData();
            const RealT* values  = bus_jacobian->getValues();
            for (IdxT i = 0; i < bus_jacobian->getNnz(); ++i)
            {
              rows_dup[counter] = rows[i];
              cols_dup[counter] = columns[i];
              vals_dup[counter] = values[i];
              counter++;
            }
          }
          else
          {
            Log::warning() << "A bus has returned a nullptr Jacobian.\n";
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

        // Update CSR values from component and bus Jacobians
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
          else
          {
            Log::warning() << "A component has returned a nullptr Jacobian.\n";
          }
        }

        for (const auto& bus : buses_)
        {
          auto bus_jacobian = bus->getCooJacobian();

          if (bus_jacobian != nullptr)
          {
            const RealT* values = bus_jacobian->getValues();
            for (IdxT i = 0; i < bus_jacobian->getNnz(); ++i)
            {
              vals[map_to_csr_[counter]] += values[i];
              counter++;
            }
          }
          else
          {
            Log::warning() << "A bus has returned a nullptr Jacobian.\n";
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
