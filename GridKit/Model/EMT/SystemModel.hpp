/**
 * @file SystemModel.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief EMT system model.
 */

#pragma once

#include <format>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <GridKit/Model/ConfigurationChecks.hpp>
#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Monitor.hpp>
#include <GridKit/Model/EMT/SystemModelData.hpp>
#include <GridKit/Model/VariableMonitorController.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief EMT system of buses, devices and their submodels
     *
     * Each component owns a contiguous block of the system vectors, with
     * parents before their submodels.
     *
     * @todo Synchronize the residual scatters and the Jacobian refill, which
     * add to entries shared by the devices at a bus.
     */
    template <typename scalar_type, typename index_type>
    class SystemModel : public Component<scalar_type, index_type>
    {
      using ComponentT   = Component<scalar_type, index_type>;
      using ScalarT      = scalar_type;
      using IdxT         = index_type;
      using RealT        = typename ComponentT::RealT;
      using IndexVectorT = typename ComponentT::IndexVectorT;
      using CsrMatrixT   = typename ComponentT::CsrMatrixT;
      using CooBufferT   = typename ComponentT::CooBuffer;
      using Port         = typename ComponentT::Port;
      using DeviceDataT  = DeviceData<IdxT>;
      using MonitorT     = Model::VariableMonitorController<ScalarT>;

      using ComponentT::abs_tol_;
      using ComponentT::f_;
      using ComponentT::tag_;
      using ComponentT::time_;
      using ComponentT::y_;
      using ComponentT::yp_;

    public:
      /// @throw std::invalid_argument for an unknown class, bus or port, or a duplicate id
      explicit SystemModel(const SystemModelData<RealT, IdxT>& data);

      Model::ConfigurationChecks verify() const override;

      int  allocate() override;
      int  initialize() override;
      int  tagDifferentiable() override;
      int  setAbsoluteTolerance(RealT rel_tol) override;
      int  evaluateResidual() override;
      int  evaluateJacobian() override;
      void updateTime(RealT t, RealT alpha) override;

      IdxT        size() override;
      IdxT        nnz() override;
      CsrMatrixT* getCsrJacobian() const override;

      bool monitoring() const override;
      void printMonitoredVariables() const override;
      void startMonitor() override;
      void stopMonitor() override;

      ComponentT& bus(const std::string& number) const;
      ComponentT& device(const std::string& id) const;

      /// @throw std::invalid_argument if device `id` is not a `ModelT`
      template <typename ModelT>
      ModelT& device(const std::string& id) const
      {
        auto* model = dynamic_cast<ModelT*>(&device(id));
        if (model == nullptr)
        {
          throw std::invalid_argument(std::format("EMT::SystemModel: device {} is not a {}", id, ModelT::CLASS_NAME));
        }
        return *model;
      }

    private:
      /// Component entries of the Jacobian and the CSR entry of each
      struct Assembly
      {
        std::vector<const CooBufferT*> entries;
        IndexVectorT                   map_to_csr;
      };

      ComponentT& add(const DeviceDataT& device);
      void        addBlocks(ComponentT& component);

      template <typename ModelT>
      ModelT& addDevice(const DeviceDataT& device);

      int compose();
      int bind(ComponentT& component, IdxT offset);
      int allocateDerivatives();

      std::unique_ptr<CsrMatrixT> assemble(IdxT rows, IdxT cols, Assembly& assembly) const;
      int                         refill(CsrMatrixT& matrix, const Assembly& assembly) const;

      std::vector<std::unique_ptr<ComponentT>> components_; ///< Buses, then devices
      std::vector<ComponentT*>                 blocks_;     ///< Components and submodels, parents first
      std::map<std::string, ComponentT*>       buses_;
      std::map<std::string, ComponentT*>       devices_;
      Monitor<ScalarT, IdxT>                   monitored_;

      std::unique_ptr<CsrMatrixT> jacobian_;
      Assembly                    jacobian_assembly_;
      std::unique_ptr<MonitorT>   monitor_;
    };
  } // namespace EMT
} // namespace GridKit
