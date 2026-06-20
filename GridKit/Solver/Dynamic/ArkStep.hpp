#pragma once

#include <GridKit/Solver/Dynamic/Arkode.hpp>

#include <arkode/arkode_arkstep.h>

namespace AnalysisManager
{
  namespace Sundials
  {
    template <class ScalarT, typename IdxT>
    class ArkStep : public Arkode<ScalarT, IdxT>
    {
      using RealT = typename GridKit::ScalarTraits<ScalarT>::RealT;

    public:
      enum class Type
      {
        Implicit,
        Explicit
      };

      ArkStep(GridKit::Model::Evaluator<ScalarT, IdxT>* model, Type type = Type::Implicit);

      void setMethod(Type type)
      {
        type_ = type;
      }

    protected:
      void* createStepper(RealT t0, N_Vector y0) override;

      bool usesImplicitSolver() const override
      {
        return type_ == Type::Implicit;
      }

    private:
      Type type_;
    };

  } // namespace Sundials
} // namespace AnalysisManager
