
#pragma once

#include <exception>
#include <functional>
#include <iostream>
#include <optional>

#include <nvector/nvector_serial.h>
#include <sundials/sundials_context.h>
#include <sunlinsol/sunlinsol_dense.h>  /* access to dense linear solver        */
#include <sunmatrix/sunmatrix_sparse.h> /* access to sparse SUNMatrix           */

#include <GridKit/Definitions.hpp>

#ifdef GRIDKIT_ENABLE_SUNDIALS_SPARSE
#include <sunlinsol/sunlinsol_klu.h> /* access to KLU linear solver          */
#endif

#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Solver/Dynamic/DynamicSolver.hpp>

namespace AnalysisManager
{
  namespace Sundials
  {
    template <class ScalarT, typename IdxT>
    class Ida : public DynamicSolver<ScalarT, IdxT>
    {
      using DynamicSolver<ScalarT, IdxT>::model_;

      using RealT = typename GridKit::ScalarTraits<ScalarT>::RealT;

    public:
      Ida(GridKit::Model::Evaluator<ScalarT, IdxT>* model);
      ~Ida();

      int configureSimulation();
      int configureLinearSolver(); // TODO: make private
#ifdef GRIDKIT_ENABLE_SUNDIALS_SPARSE
      int configureLinearSolverSparse(); // TODO: make private
#endif
      int configureLinearSolverDense(); // TODO: make private
      int getDefaultInitialCondition();
      int setIntegrationTime(RealT t_init, RealT t_final, int nout);
      int initializeSimulation(RealT t0, bool findConsistent = false);

      int runSimulation(RealT tf, int nout = 1, std::optional<std::function<void(RealT)>> step_callback = {});
      int deleteSimulation();

      int configureQuadrature();
      int initializeQuadrature();
      int runSimulationQuadrature(RealT tf, int nout = 1);
      int deleteQuadrature();

      int configureAdjoint();
      int configureLinearSolverBackward();
      int initializeAdjoint(IdxT steps = 100);
      int initializeBackwardSimulation(RealT tf);
      int runForwardSimulation(RealT tf, int nout = 1);
      int runBackwardSimulation(RealT t0);
      int deleteAdjoint();
      int deleteBackwardSimulation();

      int saveInitialCondition()
      {
        N_VScale(1.0, yy_, yy0_);
        N_VScale(1.0, yp_, yp0_);
        return 0;
      }

      int getSavedInitialCondition()
      {
        N_VScale(1.0, yy0_, yy_);
        N_VScale(1.0, yp0_, yp_);
        return 0;
      }

      RealT getInitialTime()
      {
        return t_init_;
      }

      RealT getFinalTime()
      {
        return t_final_;
      }

      int getNumberOutputTimes()
      {
        return nout_;
      }

      const RealT* getIntegral() const
      {
        return N_VGetArrayPointer(q_);
      }

      RealT* getIntegral()
      {
        return N_VGetArrayPointer(q_);
      }

      const RealT* getAdjointIntegral() const
      {
        return N_VGetArrayPointer(qB_);
      }

      RealT* getAdjointIntegral()
      {
        return N_VGetArrayPointer(qB_);
      }

      void printOutput(RealT t) const;
      void printSpecial(RealT t, N_Vector x) const;
      void printFinalStats() const;

      void setTimeStep(ScalarT timeStep) override;
      void setBackwardTimeStep(ScalarT timeStep);
      void setFixedStep(ScalarT timeStep, ScalarT relTol, ScalarT absTolFac=1);
      void setBackwardFixedStep(ScalarT timeStep, ScalarT relTol, ScalarT absTolFac=1);
      void setTolerance(ScalarT relTol, ScalarT absTolFac=1) override;
      void setBackwardTolerance(ScalarT relTol, ScalarT absTolFac=1);
      void setQuadratureTolerance(ScalarT relTol, ScalarT);
      void setBackwardQuadratureTolerance(ScalarT relTol, ScalarT);
      void setMaxSteps(IdxT maxSteps) override;
      void setBackwardMaxSteps(IdxT maxSteps);

    private:
      static int Residual(RealT    t,
                          N_Vector yy,
                          N_Vector yp,
                          N_Vector rr,
                          void*    user_data);

      static int Jac(RealT     t,
                     RealT     cj,
                     N_Vector  yy,
                     N_Vector  yp,
                     N_Vector  resvec,
                     SUNMatrix J,
                     void*     user_data,
                     N_Vector  tmp1,
                     N_Vector  tmp2,
                     N_Vector  tmp3);

      static int Integrand(RealT    t,
                           N_Vector yy,
                           N_Vector yp,
                           N_Vector rhsQ,
                           void*    user_data);

      static int adjointResidual(RealT    t,
                                 N_Vector yy,
                                 N_Vector yp,
                                 N_Vector yyB,
                                 N_Vector ypB,
                                 N_Vector rrB,
                                 void*    user_data);

      static int adjointIntegrand(RealT    t,
                                  N_Vector yy,
                                  N_Vector yp,
                                  N_Vector yyB,
                                  N_Vector ypB,
                                  N_Vector rhsQB,
                                  void*    user_data);

    private:
      static constexpr ScalarT DEFAULT_REL_TOL = 1e-5;

      void*           solver_{};
      SUNContext      context_{};
      SUNMatrix       JacobianMat_{};
      SUNMatrix       JacobianMatB_{};
      SUNLinearSolver linearSolver_{};
      SUNLinearSolver linearSolverB_{};

      RealT t_init_{};
      RealT t_final_{};
      int   nout_{}; ///< Number of integration outputs

      N_Vector yy_{};  ///< Solution vector
      N_Vector yp_{};  ///< Solution derivatives vector
      N_Vector tag_{}; ///< Tags differential variables
      N_Vector absTol_{}; ///< Tags differential variables
      N_Vector q_{};   ///< Integrand vector

      N_Vector yy0_{}; ///< Storage for initial values
      N_Vector yp0_{}; ///< Storage for initial derivatives

      N_Vector yyB_{}; ///< Adjoint solution vector
      N_Vector ypB_{}; ///< Adjoint solution derivatives vector
      N_Vector qB_{};  ///< Backward integrand vector

      int backwardID_{};

    private:
      // static void copyMat(Model::Evaluator::Mat& J, SlsMat Jida);
      static void copyVec(const N_Vector x, std::vector<ScalarT>& y);
      static void copyVec(const std::vector<ScalarT>& x, N_Vector y);
      static void copyVec(const std::vector<bool>& x, N_Vector y);

      // int check_flag(void *flagvalue, const char *funcname, int opt);
      static inline void checkAllocation(void* v, const char* functionName);
      static inline void checkOutput(int retval, const char* functionName);

      void setTimeStep(void *mem, ScalarT timeStep);
      void setFixedStep(void *mem, ScalarT timeStep, ScalarT relTol, ScalarT absTolFac);
      void setTolerance(void *mem, ScalarT relTol, ScalarT absTolFac);
      void setMaxSteps(void *mem, IdxT maxSteps);
      void setQuadratureTolerance(void *mem, ScalarT relTol, ScalarT absTolFac);
    };

    /// Simple exception to use within Ida class.
    class SundialsException : public std::exception
    {
      virtual const char* what() const throw()
      {
        return "Method in Ida class failed!\n";
      }
    };

  } // namespace Sundials

} // namespace AnalysisManager
