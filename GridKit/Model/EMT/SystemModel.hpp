#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <GridKit/Model/EMT/Component.hpp>

namespace GridKit
{
  namespace Model
  {
    template <typename scalar_type>
    class VariableMonitorController;
  }

  namespace EMT
  {
    template <typename real_type, typename index_type>
    struct SystemModelData;

    template <typename scalar_type, typename index_type>
    class Bus;

    template <typename scalar_type, typename index_type>
    class Switch;

    template <typename scalar_type, typename index_type>
    class SignalNode;

    /**
     * @brief Prototype for a system model class
     *
     * This class maps component data to system data and implements
     * Model::Evaluator for the system model. Every element of the system,
     * including buses, is an ordinary component held in one list; the bus
     * identifier map exists only to wire ports while building the system.
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

    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename Model::Evaluator<ScalarT, IdxT>::RealT;
      using CsrMatrixT = typename Model::Evaluator<ScalarT, IdxT>::CsrMatrixT;
      using CooMatrixT = typename Model::Evaluator<ScalarT, IdxT>::CooMatrixT;
      using BusT       = Bus<ScalarT, IdxT>;
      using SignalT    = SignalNode<ScalarT, IdxT>;
      using ComponentT = Component<ScalarT, IdxT>;
      using MonitorT   = Model::VariableMonitorController<ScalarT>;

      SystemModel();

      SystemModel(const SystemModelData<RealT, IdxT>& data);

      virtual ~SystemModel();

      int setGridKitComponentID(IdxT component_id) override;

      int allocate() override;
      int verify() const override;
      int initialize() override;

      bool hasJacobian() override;

      void initializeMonitor();
      void startMonitor() override;
      void stopMonitor() override;
      bool monitoring() const override;
      void printMonitoredVariables() const override;

      int tagDifferentiable() override;
      int setAbsoluteTolerance(RealT rel_tol) override;
      int evaluateInternalResidual() override;
      int evaluateExternalResidual() override;
      int evaluateResidual() override;
      int evaluateJacobian() override;

      void updateTime(RealT t, RealT a) override;

      /**
       * @brief Invalidate the assembled Jacobian structure.
       *
       * Discrete events can change which Jacobian entries are exactly zero,
       * and exact zeros are dropped from the discovered sparsity pattern.
       * Drivers call this after applying a discrete event so the next
       * evaluateJacobian() rediscovers the pattern.
       */
      void resetJacobianStructure() override;

      void addBus(BusT* bus);
      void addSignal(SignalT* signal);
      void addComponent(ComponentT* component);
      void addSwitch(ComponentT* component);

      BusT*                  getBus(IdxT bus_id);
      SignalT*               getSignal(IdxT signal_id);
      ComponentT*            getComponent(IdxT gridkit_component_id);
      Switch<ScalarT, IdxT>* getSwitch(IdxT switch_id);

    private:
      std::vector<ComponentT*> components_;
      std::vector<SignalT*>    signals_;

      std::map<IdxT, IdxT> gridkit_bus_indices_;    ///< Map between bus_id and component index, wiring only
      std::map<IdxT, IdxT> gridkit_signal_indices_; ///< Map between gridkit_signal_id and signal_id
      std::map<IdxT, IdxT> gridkit_switch_indices_; ///< Map between switch_id and component index

      bool owns_components_{false};

      /// Variable monitor
      std::unique_ptr<MonitorT> monitor_;
    }; // class SystemModel

  } // namespace EMT
} // namespace GridKit
