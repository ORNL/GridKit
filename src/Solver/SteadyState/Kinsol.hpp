
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

#include "Model/Evaluator.hpp"
#include "SteadyStateSolver.hpp"

namespace AnalysisManager
{
  namespace Sundials
  {
    template <class ScalarT, typename IdxT>
    class Kinsol : public SteadyStateSolver<ScalarT, IdxT>
    {
      using SteadyStateSolver<ScalarT, IdxT>::model_;

      typedef typename GridKit::ScalarTraits<ScalarT>::real_type real_type;

      static_assert(std::is_same_v<real_type, sunrealtype>, "real_type must be the same type as sunrealtype");

    public:
      Kinsol(GridKit::Model::Evaluator<ScalarT, IdxT>* model);
      ~Kinsol();

      int configureSimulation();
      int configureLinearSolver();
      int getDefaultInitialCondition();
      // int setIntegrationTime(real_type t_init, real_type t_final, int nout);
      // int initializeSimulation();
      int runSimulation();
      int deleteSimulation();

      // int configureQuadrature();
      // int initializeQuadrature();
      // int runSimulationQuadrature(real_type tf, int nout=1);
      // int deleteQuadrature();

      // int configureAdjoint();
      // int configureLinearSolverBackward();
      // int initializeAdjoint(IdxT steps = 100);
      // int initializeBackwardSimulation(real_type tf);
      // int runForwardSimulation(real_type tf, int nout=1);
      // int runBackwardSimulation(real_type t0);
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

      // real_type getInitialTime()
      // {
      //     return t_init_;
      // }

      // real_type getFinalTime()
      // {
      //     return t_final_;
      // }

      // int getNumberOutputTimes()
      // {
      //     return nout_;
      // }

      // const real_type* getIntegral() const
      // {
      //     return N_VGetArrayPointer(q_);
      // }

      // real_type* getIntegral()
      // {
      //     return N_VGetArrayPointer(q_);
      // }

      // const real_type* getAdjointIntegral() const
      // {
      //     return N_VGetArrayPointer(qB_);
      // }

      // real_type* getAdjointIntegral()
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
      N_Vector tag_{};   ///< Tags differential variables
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
