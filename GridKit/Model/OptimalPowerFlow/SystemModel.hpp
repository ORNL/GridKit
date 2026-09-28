/**
 * @file SystemModel.hpp
 * @brief Optimal power flow system model.
 */

#pragma once

#include <map>
#include <memory>
#include <vector>

#include <GridKit/Model/OptimalPowerFlow/Component.hpp>
#include <GridKit/Model/OptimalPowerFlow/SystemModelData.hpp>
#include <GridKit/Model/OptimizationEvaluator.hpp>
#include <GridKit/Model/StateData.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    class Bus;

    /**
     * @brief AC optimal power flow in Cartesian coordinates
     *
     * Minimizes the generation cost subject to the power balance at every
     * bus that is not infinite and the bounds of the components. The state
     * gives the starting point, the fixed demand, and the device settings.
     * The reference bus is the lowest-numbered bus unless an infinite bus
     * sets the angle reference.
     */
    template <typename scalar_type, typename index_type>
    class SystemModel : public Model::OptimizationEvaluator<scalar_type, index_type>
    {
      using EvaluatorT = Model::OptimizationEvaluator<scalar_type, index_type>;

      using EvaluatorT::f_;
      using EvaluatorT::g_;
      using EvaluatorT::g_lower_;
      using EvaluatorT::g_upper_;
      using EvaluatorT::gradient_;
      using EvaluatorT::hessian_;
      using EvaluatorT::jacobian_;
      using EvaluatorT::x_;
      using EvaluatorT::x_lower_;
      using EvaluatorT::x_upper_;

    public:
      using ScalarT     = scalar_type;
      using IdxT        = index_type;
      using RealT       = typename EvaluatorT::RealT;
      using CsrMatrixT  = typename EvaluatorT::CsrMatrixT;
      using ComponentT  = Component<ScalarT, IdxT>;
      using BusT        = Bus<ScalarT, IdxT>;
      using CooEntriesT = typename ComponentT::CooEntriesT;

      SystemModel(const SystemModelData<RealT, IdxT>& data, const Model::StateData& state);

      int allocate() override;
      int initialize() override;

      int evaluateObjective() override;
      int evaluateGradient() override;
      int evaluateConstraints() override;
      int evaluateJacobian() override;
      int evaluateHessian(RealT sigma, const RealT* lambda) override;

      /// Input state with the bus voltages and device currents at `x()`
      Model::StateData solutionState() const;

    private:
      using EntriesT = const CooEntriesT& (ComponentT::*) () const;

      int verify() const;
      int connect(ComponentT& component);
      int allocateDerivatives();

      std::unique_ptr<CsrMatrixT> assemble(IdxT rows, IdxT cols, EntriesT entries, std::vector<IdxT>& map_to_csr) const;
      int                         refill(CsrMatrixT& matrix, EntriesT entries, const std::vector<IdxT>& map_to_csr) const;

      Model::StateData state_;

      /// Buses first, then devices
      std::vector<std::unique_ptr<ComponentT>> components_;

      /// Buses by number
      std::map<IdxT, BusT*> buses_;

      /// Map from component entries to CSR entries
      std::vector<IdxT> jacobian_map_;
      std::vector<IdxT> hessian_map_;

      int data_errors_{0};
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
