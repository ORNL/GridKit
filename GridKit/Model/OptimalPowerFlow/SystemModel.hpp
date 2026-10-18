/**
 * @file SystemModel.hpp
 * @brief Optimal power flow system model.
 */

#pragma once

#include <map>
#include <memory>
#include <string_view>
#include <vector>

#include <GridKit/Model/ConfigurationChecks.hpp>
#include <GridKit/Model/OptimalPowerFlow/Component.hpp>
#include <GridKit/Model/OptimalPowerFlow/SystemModelData.hpp>
#include <GridKit/Model/OptimizationEvaluator.hpp>

namespace GridKit
{
  namespace OptimalPowerFlow
  {
    template <typename scalar_type, typename index_type>
    class Bus;

    template <typename scalar_type, typename index_type>
    class Branch;

    /**
     * @brief AC optimal power flow in Cartesian coordinates
     *
     * Minimizes the generation cost subject to the power balance at every
     * bus, which are the equality constraints, and the limits of the
     * components, which are the inequality constraints and the variable
     * bounds. Each reference bus has zero imaginary voltage, and every island
     * needs one.
     *
     * @todo Replace the temporary STL containers of system-level data.
     */
    template <typename scalar_type, typename index_type>
    class SystemModel : public Model::OptimizationEvaluator<scalar_type, index_type>
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
      using BranchT      = Branch<ScalarT, IdxT>;
      using Block        = typename ComponentT::Block;
      using CooBufferT   = typename ComponentT::CooBuffer;

      using EvaluatorT::c_;
      using EvaluatorT::d_;
      using EvaluatorT::d_lower_;
      using EvaluatorT::d_upper_;
      using EvaluatorT::equality_jacobian_;
      using EvaluatorT::f_;
      using EvaluatorT::gradient_;
      using EvaluatorT::hessian_;
      using EvaluatorT::inequality_jacobian_;
      using EvaluatorT::x_;
      using EvaluatorT::x_lower_;
      using EvaluatorT::x_upper_;

    public:
      explicit SystemModel(const SystemModelData<RealT, IdxT>& data);
      explicit SystemModel(std::vector<std::unique_ptr<ComponentT>> components);

      Model::ConfigurationChecks verify() const;

      int allocate() override;
      int initialize() override;

      int evaluateObjective() override;
      int evaluateGradient() override;
      int evaluateEqualities() override;
      int evaluateInequalities() override;
      int evaluateEqualityJacobian() override;
      int evaluateInequalityJacobian() override;

      int evaluateHessian(RealT        sigma,
                          const RealT* lambda_c,
                          const RealT* lambda_d) override;

    private:
      /// Component entries of one derivative matrix and the CSR entry of each
      struct Assembly
      {
        std::vector<const CooBufferT*> entries;
        IndexVectorT                   map_to_csr;
      };

      template <typename ModelT, typename DataT>
      static std::unique_ptr<ModelT> build(const DataT& data, std::string_view id);

      void addBus(const BusData<RealT, IdxT>& data);

      template <typename ModelT, typename DataT>
      ModelT& addDevice(const DataT& data);

      int               compose();
      int               bind(ComponentT& component, const Block& offset);
      std::vector<IdxT> unreferencedIslands() const;
      int               allocateDerivatives();

      std::unique_ptr<CsrMatrixT> assemble(IdxT rows, IdxT cols, Assembly& assembly) const;
      int                         refill(CsrMatrixT& matrix, const Assembly& assembly) const;

      /// Buses first, then devices
      std::vector<std::unique_ptr<ComponentT>> components_;

      /// Buses by number
      std::map<IdxT, BusT*> buses_;

      /// Electrical connectivity, independent of component dependencies
      std::vector<BranchT*> branches_;

      Assembly equality_jacobian_assembly_;
      Assembly inequality_jacobian_assembly_;
      Assembly hessian_assembly_;
    };
  } // namespace OptimalPowerFlow
} // namespace GridKit
