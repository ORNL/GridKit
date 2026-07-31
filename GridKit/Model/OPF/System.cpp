#include "System.hpp"

#include <cstddef>

#include <GridKit/Model/OPF/Branch.hpp>
#include <GridKit/Model/OPF/Bus.hpp>
#include <GridKit/Model/OPF/Generator.hpp>
#include <GridKit/Model/OPF/Load.hpp>
#include <GridKit/Model/OPF/Shunt.hpp>
#include <GridKit/Utilities/Errors.hpp>

namespace GridKit
{
  namespace OPF
  {
    template <class ScalarT, typename IdxT>
    System<ScalarT, IdxT>::System(const SystemDataT&      case_data,
                                  const Model::StateData& state)
      : case_data_(case_data),
        input_state_(state)
    {
    }

    template <class ScalarT, typename IdxT>
    int System<ScalarT, IdxT>::allocate()
    {
      throw Utilities::NotImplementedError(__func__);
    }

    template <class ScalarT, typename IdxT>
    int System<ScalarT, IdxT>::initialize()
    {
      throw Utilities::NotImplementedError(__func__);
    }

    template <class ScalarT, typename IdxT>
    int System<ScalarT, IdxT>::evaluateResidual()
    {
      throw Utilities::NotImplementedError(__func__);
    }

    template <class ScalarT, typename IdxT>
    int System<ScalarT, IdxT>::evaluateJacobian()
    {
      throw Utilities::NotImplementedError(__func__);
    }

    template <class ScalarT, typename IdxT>
    int System<ScalarT, IdxT>::evaluateObjective()
    {
      throw Utilities::NotImplementedError(__func__);
    }

    template <class ScalarT, typename IdxT>
    IdxT System<ScalarT, IdxT>::size()
    {
      return size_variables_;
    }

    template <class ScalarT, typename IdxT>
    IdxT System<ScalarT, IdxT>::sizeResidual()
    {
      return size_constraints_;
    }

    template <class ScalarT, typename IdxT>
    IdxT System<ScalarT, IdxT>::nnz()
    {
      return nnz_;
    }

    template <class ScalarT, typename IdxT>
    int System<ScalarT, IdxT>::getVariableBounds([[maybe_unused]] RealVectorT& lower,
                                                 [[maybe_unused]] RealVectorT& upper)
    {
      throw Utilities::NotImplementedError(__func__);
    }

    template <class ScalarT, typename IdxT>
    int System<ScalarT, IdxT>::getResidualBounds([[maybe_unused]] RealVectorT& lower,
                                                 [[maybe_unused]] RealVectorT& upper)
    {
      throw Utilities::NotImplementedError(__func__);
    }

    template <class ScalarT, typename IdxT>
    bool System<ScalarT, IdxT>::hasObjective() const
    {
      return true;
    }

    template <class ScalarT, typename IdxT>
    typename System<ScalarT, IdxT>::RealT System<ScalarT, IdxT>::objective() const
    {
      return objective_;
    }

    template <class ScalarT, typename IdxT>
    bool System<ScalarT, IdxT>::hasJacobian()
    {
      return jacobian_ != nullptr;
    }

    template <class ScalarT, typename IdxT>
    typename System<ScalarT, IdxT>::CsrMatrixT* System<ScalarT, IdxT>::getCsrJacobian() const
    {
      return jacobian_.get();
    }

    template <class ScalarT, typename IdxT>
    typename System<ScalarT, IdxT>::VectorT& System<ScalarT, IdxT>::y()
    {
      return variables_;
    }

    template <class ScalarT, typename IdxT>
    const typename System<ScalarT, IdxT>::VectorT& System<ScalarT, IdxT>::y() const
    {
      return variables_;
    }

    template <class ScalarT, typename IdxT>
    typename System<ScalarT, IdxT>::VectorT& System<ScalarT, IdxT>::getResidual()
    {
      return constraints_;
    }

    template <class ScalarT, typename IdxT>
    const typename System<ScalarT, IdxT>::VectorT& System<ScalarT, IdxT>::getResidual() const
    {
      return constraints_;
    }

    template <class ScalarT, typename IdxT>
    typename System<ScalarT, IdxT>::VectorT& System<ScalarT, IdxT>::getObjectiveGradient()
    {
      return objective_gradient_;
    }

    template <class ScalarT, typename IdxT>
    const typename System<ScalarT, IdxT>::VectorT& System<ScalarT, IdxT>::getObjectiveGradient() const
    {
      return objective_gradient_;
    }

    template <class ScalarT, typename IdxT>
    const typename System<ScalarT, IdxT>::SystemDataT& System<ScalarT, IdxT>::caseData() const
    {
      return case_data_;
    }

    template <class ScalarT, typename IdxT>
    const Model::StateData& System<ScalarT, IdxT>::inputState() const
    {
      return input_state_;
    }

    template <class ScalarT, typename IdxT>
    Model::StateData System<ScalarT, IdxT>::solutionState() const
    {
      throw Utilities::NotImplementedError(__func__);
    }

    template class System<double, int>;
    template class System<double, long int>;
    template class System<double, std::size_t>;

  } // namespace OPF
} // namespace GridKit
