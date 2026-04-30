#pragma once

#include <cassert>
#include <iostream>
#include <vector>

#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusFactory.hpp>
#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ConnectedElement.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Model/VariableMonitorController.hpp>
#include <GridKit/ScalarTraits.hpp>

// Include all components
#include <GridKit/Model/PhasorDynamics/ComponentLibrary.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename, typename>
    class SystemModel;

    enum class SystemModelInternalVariables : size_t
    {
      MAXIMUM
    };

    enum class SystemModelExternalVariables : size_t
    {
      MAXIMUM
    };

    template <typename ScalarP, typename IdxP>
    struct ConnectedElementTraits<SystemModel<ScalarP, IdxP>>
    {
      using SystemModelT = SystemModel<ScalarP, IdxP>;

      using ElementT           = SystemModelT;
      using ScalarT            = ScalarP;
      using IdxT               = IdxP;
      using RealT              = typename ScalarTraits<ScalarT>::RealT;
      using ModelDataT         = SystemModelData<RealT, IdxT>;
      using InternalVariablesT = SystemModelInternalVariables;
      using ExternalVariablesT = SystemModelExternalVariables;
      using InterfaceT         = Component<ScalarT, IdxT>;
    };

    /**
     * @brief Prototype for a system model class
     *
     * This class maps component data to system data and implements
     * Model::Evaluator for the system model. This is still work in
     * progress and code is not optimized.
     *
     * @todo Address thread safety for the system model methods.
     *
     */
    template <typename ScalarP, typename IdxP>
    class SystemModel : public ConnectedElement<SystemModel<ScalarP, IdxP>>
    {
      using ConnectedElement<SystemModel>::gridkit_component_id_;
      using ConnectedElement<SystemModel>::size_;
      using ConnectedElement<SystemModel>::nnz_;
      using ConnectedElement<SystemModel>::time_;
      using ConnectedElement<SystemModel>::alpha_;
      using ConnectedElement<SystemModel>::y_;
      using ConnectedElement<SystemModel>::yp_;
      using ConnectedElement<SystemModel>::tag_;
      using ConnectedElement<SystemModel>::f_;
      using ConnectedElement<SystemModel>::J_;
      using ConnectedElement<SystemModel>::rel_tol_;
      using ConnectedElement<SystemModel>::abs_tol_;
      using ConnectedElement<SystemModel>::variable_indices_;
      using ConnectedElement<SystemModel>::residual_indices_;

    public:
      using ScalarT    = typename ConnectedElement<SystemModel>::ScalarT;
      using IdxT       = typename ConnectedElement<SystemModel>::IdxT;
      using RealT      = typename ConnectedElement<SystemModel>::RealT;
      using BusT       = BusBase<ScalarT, IdxT>;
      using SignalT    = SignalNode<ScalarT, IdxT>;
      using ComponentT = Component<ScalarT, IdxT>;
      using CsrMatrixT = typename Model::Evaluator<ScalarT, IdxT>::CsrMatrixT;
      using ModelDataT = SystemModelData<RealT, IdxT>;

      /**
       * @brief Constructor for the system model
       */
      SystemModel();

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
      SystemModel(const ModelDataT& data);

      /**
       * @brief Destructor for the system model
       *
       * If the SystemModel owns the components, it needs to delete them upon
       * destructor call.
       */
      virtual ~SystemModel();

      /**
       * @brief Set component ID
       *
       * @note Should default to 0. The system model could be used as a
       * component in a larger system that would need to set this value.
       */
      int setGridKitComponentID(IdxT component_id) override
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
      int allocate() override;

      /**
       * @brief Verify all components are configured correctly
       *
       * This method accumulates and returns the number of errors given by
       * components. It should return 0 when all is well.
       */
      int verify() const override;

      /**
       * @brief Check components for Jacobian availability
       *
       * @return true
       * @return false
       */
      bool hasJacobian() override;

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
      int initialize() override;

      /**
       * @brief Add monitors from buses and components and start monitor
       */
      void initializeMonitor();

      void startMonitor() override
      {
        monitor_.start();
      }

      void stopMonitor() override
      {
        monitor_.stop();
      }

      /**
       * @todo Tagging differential variables
       *
       * Identify what variables in the system of differential-algebraic
       * equations are differential variables, i.e. their derivatives
       * appear in the equations.
       */
      int tagDifferentiable() override;

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
      int evaluateResidual() override;

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
      int evaluateJacobian() override;

      CsrMatrixT* getCsrJacobian() const override
      {
        return csr_jac_;
      }

      bool monitoring() const override
      {
        return !monitor_.empty();
      }

      void printMonitoredVariables() const override
      {
        monitor_.print();
      }

      /**
       * @brief Update variables in buses and components
       */
      void updateVariables()
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
      void updateTime(RealT t, RealT a) override
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
      void addBus(BusT* bus)
      {
        IdxT gridkit_bus_id                = static_cast<IdxT>(buses_.size());
        gridkit_bus_indices_[bus->busID()] = gridkit_bus_id;
        buses_.push_back(bus);
      }

      /**
       * @brief Add signal
       *
       * Add signal at the end of the signals array and map signal ID with GridKit's ID for the signal
       *
       */
      void addSignal(SignalT* signal)
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
      void addComponent(ComponentT* component)
      {
        IdxT gridkit_component_id = static_cast<IdxT>(components_.size());
        component->setGridKitComponentID(gridkit_component_id);
        components_.push_back(component);
      }

      /**
       * @brief Add fault
       *
       * The fault is added to the components array, and we keep a map to its
       * location, so it can easily be accessed.
       *
       */
      void addFault(ComponentT* component)
      {
        IdxT gridkit_component_id                = static_cast<IdxT>(components_.size());
        IdxT gridkit_fault_id                    = static_cast<IdxT>(gridkit_fault_indices_.size());
        gridkit_fault_indices_[gridkit_fault_id] = gridkit_component_id;
        addComponent(component);
      }

      /**
       * @brief Return pointer to a bus
       *
       */
      BusT* getBus(IdxT bus_id)
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
      SignalT* getSignal(IdxT signal_id)
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
      ComponentT* getComponent(IdxT gridkit_component_id)
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
      BusFault<ScalarT, IdxT>* getBusFault(IdxT fault_id)
      {
        IdxT component_id = gridkit_fault_indices_.at(fault_id);
        return dynamic_cast<BusFault<ScalarT, IdxT>*>(components_[component_id]);
      }

    private:
      std::vector<BusT*>       buses_;
      std::vector<SignalT*>    signals_;
      std::vector<ComponentT*> components_;

      std::map<IdxT, IdxT> gridkit_bus_indices_;    ///< Map between gridkit_bus_id and bus_id
      std::map<IdxT, IdxT> gridkit_signal_indices_; ///< Map between gridkit_signal_id and signal_id
      std::map<IdxT, IdxT> gridkit_fault_indices_;  ///< Map between fault_id and component_id

      bool owns_components_{false};

      IdxT*       map_to_csr_{nullptr};
      CsrMatrixT* csr_jac_{nullptr};

      /// Variable monitor
      Model::VariableMonitorController<ScalarT> monitor_;
    }; // class SystemModel

  } // namespace PhasorDynamics
} // namespace GridKit
