#pragma once

#include <memory>
#include <string>

#include <GridKit/Model/EMT/Container.hpp>

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

    /**
     * @brief Executable root of an EMT component hierarchy.
     *
     * Container owns composition, scope, and recursive traversal. SystemModel
     * adds only the root CSR matrix, monitor sinks, and event conveniences
     * required by a simulation driver.
     */
    template <typename scalar_type, typename index_type>
    class SystemModel : public Container<scalar_type, index_type>
    {
      using Container<scalar_type, index_type>::abs_tol_;
      using Container<scalar_type, index_type>::allocated_;
      using Container<scalar_type, index_type>::csr_jac_;
      using Container<scalar_type, index_type>::f_;
      using Container<scalar_type, index_type>::map_to_csr_;
      using Container<scalar_type, index_type>::nnz_;
      using Container<scalar_type, index_type>::size_;
      using Container<scalar_type, index_type>::time_;
      using Container<scalar_type, index_type>::y_;

    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename Container<ScalarT, IdxT>::RealT;
      using CsrMatrixT = typename Container<ScalarT, IdxT>::CsrMatrixT;
      using CooMatrixT = typename Container<ScalarT, IdxT>::CooMatrixT;
      using BusT       = Bus<ScalarT, IdxT>;
      using SignalT    = Signal<ScalarT, IdxT>;
      using ComponentT = Component<ScalarT, IdxT>;
      using ContainerT = Container<ScalarT, IdxT>;
      using MonitorT   = Model::VariableMonitorController<ScalarT>;

      SystemModel();
      explicit SystemModel(const SystemModelData<RealT, IdxT>& data);
      ~SystemModel() override;

      int  allocate() override;
      bool hasJacobian() override;
      int  evaluateJacobian() override;

      void initializeMonitor();
      void startMonitor() override;
      void stopMonitor() override;
      bool monitoring() const override;
      void printMonitoredVariables() const override;

      // Named component(path) and signal(path) are the primary hierarchical
      // interface. These methods preserve the current application API.
      BusT*                  getBus(const std::string& path);
      SignalT*               getSignal(const std::string& path);
      ComponentT*            getComponent(IdxT local_index);
      Switch<ScalarT, IdxT>* getSwitch(const std::string& path);

    private:
      std::unique_ptr<MonitorT> monitor_;
    };
  } // namespace EMT
} // namespace GridKit
