#pragma once

#include <cstddef>
#include <vector>

#include <GridKit/Model/OPF/Branch/BranchData.hpp>
#include <GridKit/Model/OPF/Component.hpp>
#include <GridKit/Model/OPF/ComponentConstraints.hpp>
#include <GridKit/Model/OPF/ComponentVariables.hpp>
#include <GridKit/Model/StateData.hpp>

namespace GridKit
{
  namespace OPF
  {
    enum class BranchExternalVariables : std::size_t
    {
      VMF, ///< Voltage magnitude at the from bus
      VAF, ///< Voltage angle at the from bus
      VMT, ///< Voltage magnitude at the to bus
      VAT, ///< Voltage angle at the to bus
      MAXIMUM
    };

    enum class BranchInternalConstraints : std::size_t
    {
      SF2, ///< Squared apparent power at the from terminal
      ST2, ///< Squared apparent power at the to terminal
      MAXIMUM
    };

    // These rows are present only when smax is defined. Otherwise
    // sizeInternalConstraints() is zero and their bindings remain invalid.

    enum class BranchExternalConstraints : std::size_t
    {
      DIVPF, ///< Active-power balance at the from bus
      DIVQF, ///< Reactive-power balance at the from bus
      DIVPT, ///< Active-power balance at the to bus
      DIVQT, ///< Reactive-power balance at the to bus
      MAXIMUM
    };

    /// Pi-model branch contribution with optional terminal power limits.
    template <class scalar_type, typename index_type>
    class Branch final : public Component<scalar_type, index_type>
    {
    public:
      using ScalarT      = scalar_type;
      using IdxT         = index_type;
      using ComponentT   = Component<ScalarT, IdxT>;
      using RealT        = typename ComponentT::RealT;
      using DataT        = BranchData<RealT, IdxT>;
      using StateT       = Model::DeviceState;
      using VariablesT   = ComponentVariables<ScalarT,
                                              IdxT,
                                              NoVariables,
                                              BranchExternalVariables>;
      using ConstraintsT = ComponentConstraints<IdxT,
                                                BranchInternalConstraints,
                                                BranchExternalConstraints>;

      Branch(const DataT& data, const StateT& state);
      ~Branch() override = default;

      VariablesT&         variables();
      const VariablesT&   variables() const;
      ConstraintsT&       constraints();
      const ConstraintsT& constraints() const;

      IdxT sizeInternalVariables() const override;
      IdxT sizeInternalConstraints() const override;
      IdxT nnz() const override;
      void setConstraintOffset(IdxT offset) override;
      void addJacobianSparsity(std::vector<typename ComponentT::JacobianEntry>& entries) const override;
      int  setJacobianSlots(const std::vector<IdxT>& slots) override;
      int  setConstraintBounds(RealT* lower, RealT* upper) const override;
      int  addConstraints(const ScalarT* values, ScalarT* constraints) const override;
      int  addJacobian(const ScalarT* values, RealT* jacobian_values) const override;

    private:
      DataT             data_;
      StateT            state_;
      VariablesT        variables_;
      ConstraintsT      constraints_;
      std::vector<IdxT> jacobian_slots_;
    };

  } // namespace OPF
} // namespace GridKit
