#pragma once

#include <GridKit/Solver/Dynamic/Arkode.hpp>

#include <arkode/arkode_erkstep.h>

namespace AnalysisManager
{
  namespace Sundials
  {
    template <class ScalarT, typename IdxT>
    class ErkStep : public Arkode<ScalarT, IdxT>
    {
      using RealT = typename GridKit::ScalarTraits<ScalarT>::RealT;

    public:
      ErkStep(GridKit::Model::Evaluator<ScalarT, IdxT>* model);

    protected:
      void* createStepper(RealT t0, N_Vector y0) override;

      bool usesImplicitSolver() const override
      {
        return false;
      }

      bool usesMassMatrix() const override
      {
        return false;
      }

      bool requiresIdentityMass() const override
      {
        return true;
      }
    };

  } // namespace Sundials
} // namespace AnalysisManager
