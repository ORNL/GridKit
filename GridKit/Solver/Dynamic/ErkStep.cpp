#include "ErkStep.hpp"

namespace AnalysisManager
{
  namespace Sundials
  {
    template <class ScalarT, typename IdxT>
    ErkStep<ScalarT, IdxT>::ErkStep(GridKit::Model::Evaluator<ScalarT, IdxT>* model)
      : Arkode<ScalarT, IdxT>(model)
    {
    }

    template <class ScalarT, typename IdxT>
    void* ErkStep<ScalarT, IdxT>::createStepper(RealT t0, N_Vector y0)
    {
      return ERKStepCreate(Arkode<ScalarT, IdxT>::Rhs, t0, y0, this->context_);
    }

    template class ErkStep<sunrealtype, long int>;
    template class ErkStep<sunrealtype, int>;
    template class ErkStep<sunrealtype, size_t>;

  } // namespace Sundials
} // namespace AnalysisManager
