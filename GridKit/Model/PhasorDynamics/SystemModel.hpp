#pragma once

#include <array>
#include <map>
#include <memory>
#include <vector>

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/NetworkAdmittance.hpp>

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

      int  tagDifferentiable() override;
      int  setAbsoluteTolerance(RealT rel_tol) override;
      int  evaluateResidual() override;
      int  evaluateJacobian() override;
      void printResidualPerformanceStats() const;

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
      /// Assemble the network admittance matrix and the residual sweep list.
      void assembleNetworkAdmittance();

      std::vector<BusT*>       buses_;
      std::vector<SignalT*>    signals_;
      std::vector<ComponentT*> components_;

      /// The constant linear network, evaluated ahead of the residual sweep
      NetworkAdmittance<ScalarT, IdxT> network_;
      /// Components that evaluate their own residual, in components_ order
      std::vector<ComponentT*>         residual_components_;
      /// Cached HOST storage the network product reads and writes
      ScalarT*                         network_y_data_{nullptr};
      ScalarT*                         network_f_data_{nullptr};

      std::map<IdxT, IdxT> gridkit_bus_indices_;    ///< Map between gridkit_bus_id and bus_id
      std::map<IdxT, IdxT> gridkit_signal_indices_; ///< Map between gridkit_signal_id and signal_id
      std::map<IdxT, IdxT> gridkit_fault_indices_;  ///< Map between fault_id and component_id

      bool owns_components_{false};

      /// Variable monitor
      std::unique_ptr<MonitorT> monitor_;

      static constexpr std::size_t             profile_group_count_ = 10;
      /// Group boundaries in components_, fixed at construction
      std::array<IdxT, profile_group_count_>   profile_component_ends_{};
      /// The same boundaries projected onto residual_components_
      std::array<IdxT, profile_group_count_>   profile_sweep_ends_{};
      std::array<double, profile_group_count_> profile_residual_seconds_{};
      double                                   profile_bus_residual_seconds_{};
      double                                   profile_network_residual_seconds_{};
      long int                                 profile_residual_calls_{};
    }; // class SystemModel

  } // namespace PhasorDynamics
} // namespace GridKit
