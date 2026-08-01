#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Model/OPF/Component.hpp>
#include <GridKit/Model/OPF/SystemData.hpp>
#include <GridKit/Model/StateData.hpp>

namespace GridKit
{
  namespace OPF
  {
    /**
     * @brief System-level nonlinear program exposed through Model::Evaluator.
     *
     * This class is the only OPF object consumed by a numerical solver. It
     * owns global variables, constraints, bounds, objective data, and the CSR
     * constraint Jacobian assembled from local components.
     */
    template <class scalar_type = double, typename index_type = std::size_t>
    class System final : public Model::Evaluator<scalar_type, index_type>
    {
    public:
      using ScalarT     = scalar_type;
      using IdxT        = index_type;
      using EvaluatorT  = Model::Evaluator<ScalarT, IdxT>;
      using RealT       = typename EvaluatorT::RealT;
      using VectorT     = typename EvaluatorT::VectorT;
      using RealVectorT = typename EvaluatorT::RealVectorT;
      using CsrMatrixT  = typename EvaluatorT::CsrMatrixT;
      using ComponentT  = Component<ScalarT, IdxT>;
      using SystemDataT = SystemData<RealT, IdxT>;

      /**
       * @brief Compose an OPF problem from immutable system data and an operating state.
       *
       * @pre system_data satisfies INPUT_FORMAT.md, including at least one
       *      bus, unique bus ids and numbers, unique device ids, and
       *      resolvable bus references.
       * @pre Bus state is keyed as `bus_id_<number>` and device state by device id.
       * @pre Every load has finite p and q values in state.
       *
       * allocate() validates these conditions before assigning global indices.
       * The lowest-numbered bus provides the voltage angle gauge, which fixes
       * the rotational degeneracy of the formulation and has no modeling
       * significance.
       */
      System(const SystemDataT& system_data, const Model::StateData& state);
      ~System() override = default;

      int allocate() override;
      int initialize() override;

      int evaluateResidual() override;
      int evaluateJacobian() override;
      int evaluateObjective() override;

      IdxT size() override;
      IdxT sizeResidual() override;
      IdxT nnz() override;

      int getVariableBounds(RealVectorT& lower, RealVectorT& upper) override;
      int getResidualBounds(RealVectorT& lower, RealVectorT& upper) override;

      bool  hasObjective() const override;
      RealT objective() const override;

      bool        hasJacobian() override;
      CsrMatrixT* getCsrJacobian() const override;

      VectorT&       y() override;
      const VectorT& y() const override;

      VectorT&       getResidual() override;
      const VectorT& getResidual() const override;

      VectorT&       getObjectiveGradient() override;
      const VectorT& getObjectiveGradient() const override;

      const SystemDataT&      systemData() const;
      const Model::StateData& inputState() const;

      /**
       * @brief Create portable state data containing the accepted OPF point.
       *
       * The result starts as a copy of the input state. It updates bus vr/vi,
       * generator p/q, and generator/load/shunt bus-current injections while
       * preserving statuses, taps, phases, fixed load p/q, and unrelated
       * entries. JSON serialization remains a responsibility of the
       * state-data layer.
       */
      Model::StateData solutionState() const;

    private:
      void validate() const;

      SystemDataT      system_data_;
      Model::StateData input_state_;

      std::vector<std::unique_ptr<ComponentT>> components_;

      VectorT     variables_;
      VectorT     constraints_;
      VectorT     objective_gradient_;
      RealVectorT variable_lower_;
      RealVectorT variable_upper_;
      RealVectorT constraint_lower_;
      RealVectorT constraint_upper_;

      std::unique_ptr<CsrMatrixT> jacobian_;

      IdxT  size_variables_{};
      IdxT  size_constraints_{};
      IdxT  nnz_{};
      RealT objective_{};

      bool allocated_{false};
      bool initialized_{false};
    };

  } // namespace OPF
} // namespace GridKit
