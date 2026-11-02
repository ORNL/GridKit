#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeSet.hpp>

namespace GridKit
{
  namespace Model
  {
    template <typename scalar_type>
    class VariableMonitorController;
  }

  namespace PhasorDynamics
  {
    template <typename real_type, typename index_type>
    struct SystemModelData;

    template <typename scalar_type, typename index_type>
    class BusBase;

    template <typename scalar_type, typename index_type>
    class Branch;

    template <typename scalar_type, typename index_type>
    class SignalNode;

    /**
     * @brief Prototype for a system model class
     *
     * This class maps component data to system data and implements
     * Model::Evaluator for the system model. This is still work in
     * progress and code is not optimized.
     *
     * Subsystems own their buses and devices and may connect to the parent's
     * buses (external buses). A system with subsystems contains only buses
     * and subsystems. Only the top-level system allocates storage.
     *
     * evaluateJacobian() on the top-level system assembles the system Jacobian.
     * On a subsystem, it assembles the contributions of the subsystem's buses
     * and devices in global indices, of size offset() + size().
     *
     * @todo Address thread safety for the system model methods.
     *
     */
    template <typename scalar_type, typename index_type>
    class SystemModel : public Component<scalar_type, index_type>
    {
      using Component<scalar_type, index_type>::gridkit_component_id_;
      using Component<scalar_type, index_type>::size_;
      using Component<scalar_type, index_type>::nnz_;
      using Component<scalar_type, index_type>::time_;
      using Component<scalar_type, index_type>::alpha_;
      using Component<scalar_type, index_type>::y_;
      using Component<scalar_type, index_type>::yp_;
      using Component<scalar_type, index_type>::tag_;
      using Component<scalar_type, index_type>::abs_tol_;
      using Component<scalar_type, index_type>::f_;
      using Component<scalar_type, index_type>::variable_indices_;
      using Component<scalar_type, index_type>::residual_indices_;
      using Component<scalar_type, index_type>::csr_jac_;
      using Component<scalar_type, index_type>::map_to_csr_;
      using Component<scalar_type, index_type>::allocated_;

    private:
      using ScalarT        = scalar_type;
      using IdxT           = index_type;
      using RealT          = typename Model::Evaluator<ScalarT, IdxT>::RealT;
      using CsrMatrixT     = typename Model::Evaluator<ScalarT, IdxT>::CsrMatrixT;
      using CooMatrixT     = typename Model::Evaluator<ScalarT, IdxT>::CooMatrixT;
      using BusT           = BusBase<ScalarT, IdxT>;
      using SignalNodeSetT = SignalNodeSet<ScalarT, IdxT>;
      using SignalNodeT    = SignalNodeSetT::SignalNodeT;
      using ComponentT     = Component<ScalarT, IdxT>;
      using BranchT        = Branch<ScalarT, IdxT>;
      using MonitorT       = Model::VariableMonitorController<ScalarT>;

    public:
      SystemModel();

      SystemModel(const SystemModelData<RealT, IdxT>& data);

      SystemModel(const SystemModelData<RealT, IdxT>& data,
                  const std::vector<BusT*>&           external_buses);

      virtual ~SystemModel();

      int setGridKitComponentID(IdxT component_id) override;

      int allocate() override;
      int initialize() override;

      Model::ConfigurationChecks verify() const override;
      int                        setSignalIndices() override;

      bool hasJacobian() override;

      void initializeMonitor();
      void startMonitor() override;
      void stopMonitor() override;
      bool monitoring() const override;
      void printMonitoredVariables() const override;

      int tagDifferentiable() override;
      int setAbsoluteTolerance(RealT rel_tol) override;
      int evaluateResidual() override;
      int evaluateJacobian() override;

      void updateTime(RealT t, RealT a) override;

      void addBus(BusT* bus);
      void addComponent(ComponentT* component, const std::string& label = "");
      void addSubsystem(SystemModel* subsystem);

      void setSystemBase(RealT freq_system_base, RealT va_system_base);

      BusT*        getBus(IdxT bus_id);
      SignalNodeT* getSignalNode(IdxT signal_id);
      ComponentT*  getComponent(IdxT gridkit_component_id);
      ComponentT*  getComponent(const std::string& label);
      BranchT*     getBranch(const std::string& id);

      IdxT               numComponents() const;
      const std::string& getComponentLabel(IdxT gridkit_component_id) const;
      IdxT               numSubsystems() const;
      SystemModel*       getSubsystem(IdxT index);
      IdxT               offset() const;

      const std::vector<BusT*>& getExternalBuses() const;

    private:
      BusT*        findBus(IdxT bus_id);
      SignalNodeT* findSignalNode(IdxT signal_id);
      ComponentT*  findComponent(const std::string& label);
      IdxT         computeSize();
      void         bindChildren(IdxT global_offset);
      void         addMonitorsTo(MonitorT& monitor) const;
      void         evaluateCooJacobians();

      std::vector<BusT*>        buses_;
      SignalNodeSetT            signal_nodes_;
      std::vector<ComponentT*>  components_;
      std::vector<SystemModel*> subsystems_;

      /// Component labels (device class and case-file ID), by gridkit_component_id
      std::vector<std::string> component_labels_;

      /// Map between component label and gridkit_component_id
      std::map<std::string, IdxT> gridkit_component_indices_;

      std::map<IdxT, IdxT>  gridkit_bus_indices_;    ///< Map between gridkit_bus_id and bus_id
      std::map<IdxT, BusT*> external_buses_;         ///< Parent's buses by bus_id, not owned
      std::vector<BusT*>    ordered_external_buses_; ///< External buses in construction order, not owned

      IdxT offset_{0};           ///< Global index of this system's first variable
      bool is_subsystem_{false}; ///< Allocated and bound by a parent system
      bool owns_components_{false};

      /// Variable monitor
      std::unique_ptr<MonitorT> monitor_;
    }; // class SystemModel

  } // namespace PhasorDynamics
} // namespace GridKit
