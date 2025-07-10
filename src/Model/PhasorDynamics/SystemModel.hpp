#pragma once

#include <cassert>
#include <iostream>
#include <vector>

#include <Model/PhasorDynamics/BusBase.hpp>
#include <Model/PhasorDynamics/Component.hpp>
#include <Model/PhasorDynamics/SystemModelData.hpp>
#include <ScalarTraits.hpp>

// Temporary
#include <Model/PhasorDynamics/Branch/Branch.hpp>
#include <Model/PhasorDynamics/Bus/BusFactory.hpp>
#include <Model/PhasorDynamics/BusFault/BusFault.hpp>
#include <Model/PhasorDynamics/Governor/Tgov1/Tgov1.hpp>
#include <Model/PhasorDynamics/Load/Load.hpp>
#include <Model/PhasorDynamics/SynchronousMachine/GENROUwS/Genrou.hpp>
#include <Model/PhasorDynamics/SynchronousMachine/GenClassical/GenClassical.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {

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
    template <class ScalarT, typename IdxT>
    class SystemModel : public PhasorDynamics::Component<ScalarT, IdxT>
    {
      using bus_type       = PhasorDynamics::BusBase<ScalarT, IdxT>;
      using component_type = PhasorDynamics::Component<ScalarT, IdxT>;
      using real_type      = typename Model::Evaluator<ScalarT, IdxT>::real_type;

      using PhasorDynamics::Component<ScalarT, IdxT>::size_;
      using PhasorDynamics::Component<ScalarT, IdxT>::nnz_;
      using PhasorDynamics::Component<ScalarT, IdxT>::time_;
      using PhasorDynamics::Component<ScalarT, IdxT>::alpha_;
      using PhasorDynamics::Component<ScalarT, IdxT>::y_;
      using PhasorDynamics::Component<ScalarT, IdxT>::yp_;
      using PhasorDynamics::Component<ScalarT, IdxT>::tag_;
      using PhasorDynamics::Component<ScalarT, IdxT>::f_;
      using PhasorDynamics::Component<ScalarT, IdxT>::rel_tol_;
      using PhasorDynamics::Component<ScalarT, IdxT>::abs_tol_;

    public:
      /**
       * @brief Constructor for the system model
       */
      SystemModel()
      {
        // Set system model tolerances
        rel_tol_         = 1e-7;
        abs_tol_         = 1e-9;
        this->max_steps_ = 2000;
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
      SystemModel(SystemModelData<real_type, IdxT>& data)
      {
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

        // Add GENROU generators
        for (const auto& gendata : data.genrou)
        {
          IdxT bus_index = 0;
          IdxT pmech_index = 0;
          IdxT speed_index = 0;

          bool has_bus = gendata.ports.contains(GenrouData<ScalarT, IdxT>::Ports::bus);
          bool has_pmech_signal = gendata.ports.contains(GenrouData<ScalarT, IdxT>::Ports::pmech_signal);
          bool has_speed_signal = gendata.ports.contains(GenrouData<ScalarT, IdxT>::Ports::speed_signal);

          if (has_bus && has_pmech_signal && has_speed_signal)
          {

            bus_index = gendata.ports.at(GenrouData<ScalarT, IdxT>::Ports::bus);
            pmech_index = gendata.ports.at(GenrouData<ScalarT, IdxT>::Ports::pmech_signal);
            speed_index = gendata.ports.at(GenrouData<ScalarT, IdxT>::Ports::speed_signal);

            auto* gen = new Genrou<ScalarT, IdxT>(
              getBus(bus_index), 
              getBus(pmech_index),
              getBus(speed_index),
              gendata);
            addComponent(gen);
          }
          else{
            bus_index = gendata.ports.at(GenrouData<ScalarT, IdxT>::Ports::bus);
            
            auto* gen = new Genrou<ScalarT, IdxT>(
              getBus(bus_index), 
              nullptr,
              nullptr,
              gendata);
            addComponent(gen);
          }
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

        // Add governors
        for (const auto& govdata : data.gov)
        {

          /* Need to add signals to gendata class */
          auto* gov = new Governor::Tgov1<ScalarT, IdxT>(
              getBus(govdata.signal_pmech),
              getBus(govdata.signal_speed),
              govdata);

          addComponent(gov);
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
      }

      /**
       * @brief Destructor for the system model
       *
       * If the SystemModel owns the components, it needs to delete them upon
       * destructor call.
       */
      virtual ~SystemModel()
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
       * @brief Allocate buses, components, and system objects.
       *
       * This method first allocates bus objects, then component objects,
       * and computes system size (number of unknowns). Once the size is
       * computed, system global objects are allocated.
       *
       * @post size_ >= 1
       *
       */
      int allocate()
      {
        size_ = 0;

        // Allocate all buses
        for (const auto& bus : buses_)
        {
          bus->allocate();
          size_ += bus->size();
        }

        // Allocate all components
        for (const auto& component : components_)
        {
          component->allocate();
          size_ += component->size();
        }

        // Allocate global vectors
        y_.resize(size_);
        yp_.resize(size_);
        f_.resize(size_);
        tag_.resize(size_);

        return 0;
      }

      /**
       * @brief Assume that jacobian is not avalible
       *
       * @return true
       * @return false
       */
      bool hasJacobian()
      {
        return false;
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
       */
      int initialize()
      {
        // Set initial values for global solution vectors
        IdxT varOffset = 0;

        for (const auto& bus : buses_)
        {
          bus->initialize();
        }

        for (const auto& bus : buses_)
        {
          for (IdxT j = 0; j < bus->size(); ++j)
          {
            y_[varOffset + j]  = bus->y()[j];
            yp_[varOffset + j] = bus->yp()[j];
          }
          varOffset += bus->size();
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
            y_[varOffset + j]  = component->y()[j];
            yp_[varOffset + j] = component->yp()[j];
          }
          varOffset += component->size();
        }

        return 0;
      }

      /**
       * @todo Tagging differential variables
       *
       * Identify what variables in the system of differential-algebraic
       * equations are differential variables, i.e. their derivatives
       * appear in the equations.
       */
      int tagDifferentiable()
      {
        // Set initial values for global solution vectors
        IdxT offset = 0;
        for (const auto& bus : buses_)
        {
          bus->tagDifferentiable();
          for (IdxT j = 0; j < bus->size(); ++j)
          {
            tag_[offset + j] = bus->tag()[j];
          }
          offset += bus->size();
        }

        for (const auto& component : components_)
        {
          component->tagDifferentiable();
          for (IdxT j = 0; j < component->size(); ++j)
          {
            tag_[offset + j] = component->tag()[j];
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
      int evaluateResidual()
      {
        // Update variables
        IdxT varOffset = 0;
        for (const auto& bus : buses_)
        {
          for (IdxT j = 0; j < bus->size(); ++j)
          {
            bus->y()[j]  = y_[varOffset + j];
            bus->yp()[j] = yp_[varOffset + j];
          }
          varOffset += bus->size();

          bus->evaluateResidual();
        }

        for (const auto& component : components_)
        {
          for (IdxT j = 0; j < component->size(); ++j)
          {
            component->y()[j]  = y_[varOffset + j];
            component->yp()[j] = yp_[varOffset + j];
          }
          varOffset += component->size();

          component->evaluateResidual();
        }

        // Update residual vector
        IdxT resOffset = 0;
        for (const auto& bus : buses_)
        {
          for (IdxT j = 0; j < bus->size(); ++j)
          {
            f_[resOffset + j] = bus->getResidual()[j];
          }
          resOffset += bus->size();
        }

        for (const auto& component : components_)
        {
          for (IdxT j = 0; j < component->size(); ++j)
          {
            f_[resOffset + j] = component->getResidual()[j];
          }
          resOffset += component->size();
        }

        return 0;
      }

      /**
       * @brief Evaluate system Jacobian.
       *
       * @todo Need to implement Jacobian. For now, using finite difference
       * approximation provided by IDA. This works for dense Jacobian matrix
       * only.
       *
       */
      int evaluateJacobian()
      {
        return 0;
      }

      void updateTime(real_type t, real_type a)
      {
        for (const auto& component : components_)
        {
          component->updateTime(t, a);
        }
      }

      void addBus(bus_type* bus)
      {
        buses_.push_back(bus);
      }

      void addComponent(component_type* component)
      {
        components_.push_back(component);
      }

      void addFault(component_type* component)
      {
        components_.push_back(component);
        faults_.push_back(component);
      }

      bus_type* getBus(IdxT busid)
      {
        // Need to implement mapping of bus IDs to buses in the system model
        assert((buses_[busid])->busID() == busid);
        return buses_[busid];
      }

      /**
       * @brief Return pointer to a bus fault model
       *
       * This function is used to provide easier access to setting and
       * clearing faults from the SystemModel interface.
       *
       * @warning This is a hack to get access to bus faults in examples.
       * A more comprehensive solution is needed.
       */
      BusFault<ScalarT, IdxT>* getBusFault(IdxT fault_id)
      {
        return dynamic_cast<BusFault<ScalarT, IdxT>*>(faults_[fault_id]);
      }

    private:
      std::vector<bus_type*>       buses_;
      std::vector<component_type*> components_;
      std::vector<component_type*> faults_;

      bool owns_components_{false};

    }; // class SystemModel

  } // namespace PhasorDynamics
} // namespace GridKit
