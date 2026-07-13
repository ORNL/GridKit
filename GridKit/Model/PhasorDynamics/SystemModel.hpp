#pragma once

#include <map>
#include <memory>
#include <vector>

#include <GridKit/Model/PhasorDynamics/Component.hpp>

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
    class BusFault;

    template <typename scalar_type, typename index_type>
    class SignalNode;

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
      using BusT       = BusBase<ScalarT, IdxT>;
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
      int evaluateResidual() override;
      int evaluateJacobian() override;

      void updateTime(RealT t, RealT a) override;

      void addBus(BusT* bus);
      void addSignal(SignalT* signal);
      void addComponent(ComponentT* component);
      void addFault(ComponentT* component);

      void setSystemBase(RealT freq_system_base, RealT va_system_base);

      BusT*                    getBus(IdxT bus_id);
      SignalT*                 getSignal(IdxT signal_id);
      ComponentT*              getComponent(IdxT gridkit_component_id);
      BusFault<ScalarT, IdxT>* getBusFault(IdxT fault_id);

    private:
      std::vector<BusT*>       buses_;
      std::vector<SignalT*>    signals_;
      std::vector<ComponentT*> components_;

      std::map<IdxT, IdxT> gridkit_bus_indices_;    ///< Map between gridkit_bus_id and bus_id
      std::map<IdxT, IdxT> gridkit_signal_indices_; ///< Map between gridkit_signal_id and signal_id
      std::map<IdxT, IdxT> gridkit_fault_indices_;  ///< Map between fault_id and component_id

      bool owns_components_{false};

      /// Variable monitor
      std::unique_ptr<MonitorT> monitor_;
    }; // class SystemModel

  } // namespace PhasorDynamics
} // namespace GridKit
