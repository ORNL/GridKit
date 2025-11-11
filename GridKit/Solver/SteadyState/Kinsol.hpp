
/**
 * @file Kinsol.hpp
 * @author Slaven Peles <slaven.peles@pnnl.gov>
 *
 * Contains declaration of interface to KINSOL nonlinear solver from
 * SUNDIALS library.
 *
 */
#pragma once

#include <exception>
#include <iostream>

#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_sparse.h> /* access to sparse SUNMatrix           */
// #include <sunlinsol/sunlinsol_klu.h>       /* access to KLU linear solver          */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense linear solver        */

#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Solver/SteadyState/SteadyStateSolver.hpp>

namespace AnalysisManager
{
  namespace Sundials
  {
    template <class ScalarT, typename IdxT>
    class Kinsol : public SteadyStateSolver<ScalarT, IdxT>
    {
      using SteadyStateSolver<ScalarT, IdxT>::model_;

      typedef typename GridKit::ScalarTraits<ScalarT>::RealT RealT;

      static_assert(std::is_same_v<RealT, sunrealtype>, "RealT must be the same type as sunrealtype");

    public:
      Kinsol(GridKit::Model::Evaluator<ScalarT, IdxT>* model);
      ~Kinsol();

      int configureSimulation();
      int configureLinearSolver();
      int getDefaultInitialCondition();
      // int setIntegrationTime(RealT t_init, RealT t_final, int nout);
      // int initializeSimulation();
      int runSimulation();
      int deleteSimulation();

      // int configureQuadrature();
      // int initializeQuadrature();
      // int runSimulationQuadrature(RealT tf, int nout=1);
      // int deleteQuadrature();

      // int configureAdjoint();
      // int configureLinearSolverBackward();
      // int initializeAdjoint(IdxT steps = 100);
      // int initializeBackwardSimulation(RealT tf);
      // int runForwardSimulation(RealT tf, int nout=1);
      // int runBackwardSimulation(RealT t0);
      // int deleteAdjoint();

      int saveInitialCondition()
      {
        N_VScale(1.0, yy_, yy0_);
        return 0;
      }

      int getSavedInitialCondition()
      {
        N_VScale(1.0, yy0_, yy_);
        return 0;
      }

      // RealT getInitialTime()
      // {
      //     return t_init_;
      // }

      // RealT getFinalTime()
      // {
      //     return t_final_;
      // }

      // int getNumberOutputTimes()
      // {
      //     return nout_;
      // }

      // const RealT* getIntegral() const
      // {
      //     return N_VGetArrayPointer(q_);
      // }

      // RealT* getIntegral()
      // {
      //     return N_VGetArrayPointer(q_);
      // }

      // const RealT* getAdjointIntegral() const
      // {
      //     return N_VGetArrayPointer(qB_);
      // }

      // RealT* getAdjointIntegral()
      // {
      //     return N_VGetArrayPointer(qB_);
      // }

      void printOutput();
      void printSpecial(sunrealtype t, N_Vector x);
      void printFinalStats();

    private:
      static int Residual(N_Vector yy, N_Vector rr, void* user_data);

      // static int Integrand(sunrealtype t,
      //                      N_Vector yy,   N_Vector yp,
      //                      N_Vector rhsQ, void *user_data);

      // static int adjointResidual(sunrealtype t,
      //                            N_Vector yy,  N_Vector yp,
      //                            N_Vector yyB, N_Vector ypB,
      //                            N_Vector rrB, void *user_data);

      // static int adjointIntegrand(sunrealtype t,
      //                             N_Vector yy,  N_Vector yp,
      //                             N_Vector yyB, N_Vector ypB,
      //                             N_Vector rhsQB, void *user_data);

    private:
      void*           solver_{};
      SUNContext      context_{};
      SUNMatrix       JacobianMat_{};
      SUNLinearSolver linearSolver_{};

      N_Vector yy_{};    ///< Solution vector
      N_Vector scale_{}; ///< Scaling vector
      N_Vector q_{};     ///< Integrand vector

      N_Vector yy0_{}; ///< Storage for initial values

    private:
      // static void copyMat(Model::Evaluator::Mat& J, SlsMat Jida);
      static void copyVec(const N_Vector x, std::vector<ScalarT>& y);
      static void copyVec(const std::vector<ScalarT>& x, N_Vector y);
      static void copyVec(const std::vector<bool>& x, N_Vector y);

      // int check_flag(void *flagvalue, const char *funcname, int opt);
      inline void checkAllocation(void* v, const char* functionName);
      inline void checkOutput(int retval, const char* functionName);
    };

    /// Simple exception to use within Kinsol class.
    class SundialsException : public std::exception
    {
      virtual const char* what() const throw()
      {
        return "Method in Kinsol class failed!\n";
      }
    };

  } // namespace Sundials

} // namespace AnalysisManager
