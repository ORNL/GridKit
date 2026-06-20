#include "ArkStep.hpp"

namespace AnalysisManager
{
  namespace Sundials
  {
    template <class ScalarT, typename IdxT>
    ArkStep<ScalarT, IdxT>::ArkStep(GridKit::Model::Evaluator<ScalarT, IdxT>* model, Type type)
      : Arkode<ScalarT, IdxT>(model),
        type_(type)
    {
    }

    template <class ScalarT, typename IdxT>
    void* ArkStep<ScalarT, IdxT>::createStepper(RealT t0, N_Vector y0)
    {
      if (type_ == Type::Implicit)
        return ARKStepCreate(nullptr, Arkode<ScalarT, IdxT>::Rhs, t0, y0, this->context_);
      return ARKStepCreate(Arkode<ScalarT, IdxT>::Rhs, nullptr, t0, y0, this->context_);
    }

    template class ArkStep<sunrealtype, long int>;
    template class ArkStep<sunrealtype, int>;
    template class ArkStep<sunrealtype, size_t>;

  } // namespace Sundials
} // namespace AnalysisManager
