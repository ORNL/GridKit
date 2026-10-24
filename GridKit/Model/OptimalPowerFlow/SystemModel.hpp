/**
 * @file SystemModel.hpp
 * @brief Optimal power flow system model.
 */

#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <GridKit/Model/ConfigurationChecks.hpp>
#include <GridKit/Model/OptimalPowerFlow/Bus/Bus.hpp>
#include <GridKit/Model/OptimalPowerFlow/Component.hpp>
#include <GridKit/Model/OptimalPowerFlow/Generator/Generator.hpp>
#include <GridKit/Model/OptimalPowerFlow/SystemModelData.hpp>
#include <GridKit/Model/OptimizationEvaluator.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    /**
     * @brief AC optimal power flow in Cartesian coordinates
     *
     * @todo Replace the temporary STL containers of system-level data.
     */
    template <typename scalar_type, typename index_type>
    class SystemModel : public Component<scalar_type, index_type>
    {
      using ComponentT   = Component<scalar_type, index_type>;
      using EvaluatorT   = Model::OptimizationEvaluator<scalar_type, index_type>;
      using ScalarT      = scalar_type;
      using IdxT         = index_type;
      using RealT        = typename EvaluatorT::RealT;
      using RealVectorT  = typename EvaluatorT::RealVectorT;
      using IndexVectorT = LinearAlgebra::Vector<IdxT, IdxT>;
      using CsrMatrixT   = typename EvaluatorT::CsrMatrixT;
      using BusT         = Bus<ScalarT, IdxT>;
      using GeneratorT   = Generator<ScalarT, IdxT>;
      using Block        = typename ComponentT::Block;
      using CooBufferT   = typename ComponentT::CooBuffer;

      using EvaluatorT::f_;
      using EvaluatorT::g_;
      using EvaluatorT::g_lower_;
      using EvaluatorT::g_upper_;
      using EvaluatorT::gradient_;
      using EvaluatorT::hessian_;
      using EvaluatorT::jacobian_;
      using EvaluatorT::lambda_;
      using EvaluatorT::tag_;
      using EvaluatorT::x_;
      using EvaluatorT::x_lower_;
      using EvaluatorT::x_upper_;
      using EvaluatorT::z_lower_;
      using EvaluatorT::z_upper_;

    public:
      explicit SystemModel(const SystemModelData<RealT, IdxT>& data);
      explicit SystemModel(std::vector<std::unique_ptr<ComponentT>> components);

      Model::ConfigurationChecks verify() const override;

      int allocate() override;
      int tagInequality() override;
      int initialize() override;

      int evaluateObjective() override;
      int evaluateGradient() override;
      int evaluateConstraints() override;
      int evaluateJacobian() override;
      int evaluateHessian(RealT sigma) override;

      /// Bus by number, which views its solution after `allocate()`
      BusT* getBus(IdxT number);

      /// Generator by case-file `id`, which views its solution after `allocate()`
      GeneratorT* getGenerator(const std::string& id);

    private:
      /// Component entries of one derivative matrix and the CSR entry of each
      struct Assembly
      {
        std::vector<const CooBufferT*> entries;
        IndexVectorT                   map_to_csr;
      };

      void addBus(const BusData<RealT, IdxT>& data);

      template <typename ModelT, typename DataT>
      ModelT* addDevice(const DataT& data);

      int compose();
      int bind(ComponentT& component, const Block& offset);
      int allocateDerivatives();

      std::unique_ptr<CsrMatrixT> assemble(IdxT rows, IdxT cols, Assembly& assembly) const;
      int                         refill(CsrMatrixT& matrix, const Assembly& assembly) const;

      /// Buses first, then devices
      std::vector<std::unique_ptr<ComponentT>> components_;

      /// Buses by number
      std::map<IdxT, BusT*> buses_;

      /// Generators by id
      std::map<std::string, GeneratorT*> generators_;

      Assembly jacobian_assembly_;
      Assembly hessian_assembly_;
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
