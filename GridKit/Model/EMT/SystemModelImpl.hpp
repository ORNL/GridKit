#include <algorithm>
#include <stdexcept>

#include <GridKit/Definitions.hpp>
#include <GridKit/Model/EMT/ContainerImpl.hpp>
#include <GridKit/Model/EMT/SystemModel.hpp>
#include <GridKit/Model/EMT/SystemModelData.hpp>
#include <GridKit/Model/VariableMonitorController.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::SystemModel()
      : monitor_(std::make_unique<MonitorT>())
    {
    }

    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::SystemModel(
        const SystemModelData<RealT, IdxT>& data)
      : ContainerT(static_cast<const ContainerData<RealT, IdxT>&>(data)),
        monitor_(std::make_unique<MonitorT>(time_))
    {
      for (const auto& sink : data.monitor_sink)
      {
        monitor_->addSink(sink);
      }
    }

    template <typename scalar_type, typename index_type>
    SystemModel<scalar_type, index_type>::~SystemModel() = default;

    /**
     * @brief Allocate the complete hierarchy and prepare the executable root.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::allocate()
    {
      if (this->boundToParent())
      {
        return ContainerT::allocate();
      }

      if (ContainerT::allocate() != 0)
      {
        throw std::runtime_error("SystemModel allocation failed");
      }

      const int errors = this->verify();
      if (errors != 0)
      {
        Log::error() << "Component errors: " << errors << '\n';
        throw std::runtime_error("SystemModel allocation failed");
      }

      initializeMonitor();
      startMonitor();
      allocated_ = true;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::initialize(const std::map<std::string, std::map<std::string, RealT>>& state)
    {
      const int status = ContainerT::initialize(state);
      if (status != 0)
        return status;
      if (hasJacobian())
      {
        this->resetJacobianStructure();
        this->evaluateResidual();
        this->evaluateJacobian();
      }
      return 0;
    }

    /**
     * @brief Classify and validate the assembled DAE before IDA initialization.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::tagDifferentiable()
    {
      if (this->boundToParent())
        return ComponentT::tagDifferentiable();
      if (!hasJacobian())
        throw std::runtime_error("EMT DAE validation requires an exact model Jacobian");

      const auto Fyp      = this->jacobianEntries(ZERO<RealT>, ONE<RealT>);
      const auto Fy       = this->jacobianEntries(ONE<RealT>, ZERO<RealT>);
      const auto analysis = analyzeDae(static_cast<size_t>(size_), Fy, Fyp);
      if (analysis.status != DaeAnalysis::Status::regular)
      {
        std::string message = analysis.status == DaeAnalysis::Status::structurally_singular
                                  ? "EMT DAE is structurally deficient in the current differential/algebraic partition."
                                  : "EMT DAE initial-value Jacobian is numerically singular at the supplied state.";
        for (const auto row : analysis.equations)
          message += "\n  Equation: " + this->describeDaeIndex(static_cast<IdxT>(row));
        for (const auto column : analysis.variables)
          message += "\n  Variable: " + this->describeDaeIndex(static_cast<IdxT>(column));
        throw std::runtime_error(message);
      }
      this->setDifferentialTags(analysis.differential);
      return this->evaluateJacobian();
    }

    /**
     * @brief Deduplicate the hierarchy's aggregate COO into the solver CSR.
     */
    template <typename scalar_type, typename index_type>
    int SystemModel<scalar_type, index_type>::assembleJacobian(RealT y_scale, RealT yp_scale)
    {
      if (this->boundToParent())
      {
        return ContainerT::assembleJacobian(y_scale, yp_scale);
      }

      const int status = ContainerT::assembleJacobian(y_scale, yp_scale);
      if (status != 0)
      {
        return status;
      }

      auto* aggregate = this->getCooJacobian();
      if (aggregate == nullptr)
      {
        return 0;
      }

      const IdxT nnz_dup = aggregate->getNnz();
      if (csr_jac_ == nullptr)
      {
        auto* rows_dup = new IdxT[static_cast<size_t>(nnz_dup)];
        auto* cols_dup = new IdxT[static_cast<size_t>(nnz_dup)];
        auto* vals_dup = new RealT[static_cast<size_t>(nnz_dup)];
        std::copy(aggregate->getRowData(), aggregate->getRowData() + nnz_dup, rows_dup);
        std::copy(aggregate->getColData(), aggregate->getColData() + nnz_dup, cols_dup);
        std::copy(aggregate->getValues(), aggregate->getValues() + nnz_dup, vals_dup);

        CooMatrixT jac(size_, size_, nnz_dup, &rows_dup, &cols_dup, &vals_dup);
        IdxT*      row_ptrs = jac.getCsrRowData();
        nnz_                = jac.getNnz();

        auto* cols = new IdxT[static_cast<size_t>(nnz_)];
        auto* vals = new RealT[static_cast<size_t>(nnz_)];
        std::copy(jac.getColData(), jac.getColData() + nnz_, cols);
        std::copy(jac.getValues(), jac.getValues() + nnz_, vals);
        csr_jac_ = new CsrMatrixT(size_, size_, nnz_, &row_ptrs, &cols, &vals);

        const auto* map_to_sorted = jac.getMapToSorted();
        const auto* map_to_dedup  = jac.getMapToDeduplicated();
        map_to_csr_               = new IdxT[static_cast<size_t>(nnz_dup)];
        for (IdxT i = 0; i < nnz_dup; ++i)
        {
          map_to_csr_[map_to_sorted[i]] = map_to_dedup[i];
        }
      }
      else
      {
        auto* values = csr_jac_->getValues();
        std::fill(values, values + csr_jac_->getNnz(), RealT{});
        for (IdxT i = 0; i < nnz_dup; ++i)
        {
          values[map_to_csr_[i]] += aggregate->getValues()[i];
        }
      }
      nnz_ = csr_jac_->getNnz();
      return 0;
    }

    template <typename scalar_type, typename index_type>
    void SystemModel<scalar_type, index_type>::initializeMonitor()
    {
      this->forEachComponent([this](const ComponentT& component)
                             {
        const auto* child_monitor = component.getMonitor();
        if (child_monitor != nullptr && !child_monitor->empty())
        {
          monitor_->addMonitor(child_monitor);
        } });
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
      typename SignalT::ReadScope scope;
      monitor_->print();
    }

    template <typename scalar_type, typename index_type>
    typename SystemModel<scalar_type, index_type>::BusT*
    SystemModel<scalar_type, index_type>::getBus(const std::string& path)
    {
      return &this->template component<BusT>(path);
    }

    template <typename scalar_type, typename index_type>
    typename SystemModel<scalar_type, index_type>::SignalT*
    SystemModel<scalar_type, index_type>::getSignal(const std::string& path)
    {
      return &this->signal(path);
    }

    template <typename scalar_type, typename index_type>
    typename SystemModel<scalar_type, index_type>::ComponentT*
    SystemModel<scalar_type, index_type>::getComponent(IdxT local_index)
    {
      return &this->component(local_index);
    }

    template <typename scalar_type, typename index_type>
    Switch<scalar_type, index_type>*
    SystemModel<scalar_type, index_type>::getSwitch(const std::string& path)
    {
      return &this->template component<Switch<ScalarT, IdxT>>(path);
    }
  } // namespace EMT
} // namespace GridKit
