
#ifndef _IDA_HPP_
#define _IDA_HPP_

#include <exception>
#include <functional>
#include <iostream>
#include <optional>

#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>  /* access to dense linear solver        */
#include <sunlinsol/sunlinsol_klu.h>    /* access to KLU linear solver          */
#include <sunmatrix/sunmatrix_sparse.h> /* access to sparse SUNMatrix           */

#include "DynamicSolver.hpp"
#include "Model/Evaluator.hpp"

namespace AnalysisManager
{
  namespace Sundials
  {
    template <class ScalarT, typename IdxT>
    class Ida : public DynamicSolver<ScalarT, IdxT>
    {
      using DynamicSolver<ScalarT, IdxT>::model_;

      typedef typename GridKit::ScalarTraits<ScalarT>::real_type real_type;

      static_assert(std::is_same_v<real_type, sunrealtype>, "real_type must be the same type as sunrealtype");
      // TODO: can't include this assert because cpp file instantiates a few values of IdxT
      // static_assert(std::is_same_v<IdxT, sunindextype>, "IdxT must be the same type as sunindextype");

    public:
      Ida(GridKit::Model::Evaluator<ScalarT, IdxT>* model);
      ~Ida();

      int configureSimulation();
      int configureLinearSolver();
      int getDefaultInitialCondition();
      int setIntegrationTime(real_type t_init, real_type t_final, int nout);
      int initializeSimulation(real_type t0, bool findConsistent = false);

      /// @brief               Run the IDA solver on the given model and produce a solution at the given final time.
      /// @attention           `Ida::initializeSimulation` must be called with the initial time (and to find consistent
      ///                      initial conditions if needed) before calling `runSimulation`.
      /// @param tf            The final time, used to calculate the size of each step IDA should take. The actual time
      ///                      of the final IDA solution should be somewhat close to `tf`, however due to rounding error
      ///                      the precise final time may be before or after `tf`.
      ///
      /// @param nout          The number of times the IDA integrator should be invoked to calculate the solution at `tf`.
      ///
      /// @param step_callback An optional callback which, if provided, will be called after each time the IDA solver has
      ///                      been invoked with the value of `t` that IDA has calculated the last step at. The provided
      ///                      model will be updated with the latest values of `y` and `yp` before the callback is invoked.
      int runSimulation(real_type tf, int nout = 1, std::optional<std::function<void(real_type)>> step_callback = {});
      int deleteSimulation();

      int configureQuadrature();
      int initializeQuadrature();
      int runSimulationQuadrature(real_type tf, int nout = 1);
      int deleteQuadrature();

      int configureAdjoint();
      int configureLinearSolverBackward();
      int initializeAdjoint(IdxT steps = 100);
      int initializeBackwardSimulation(real_type tf);
      int runForwardSimulation(real_type tf, int nout = 1);
      int runBackwardSimulation(real_type t0);
      int deleteAdjoint();

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

      real_type getInitialTime()
      {
        return t_init_;
      }

      real_type getFinalTime()
      {
        return t_final_;
      }

      int getNumberOutputTimes()
      {
        return nout_;
      }

      const real_type* getIntegral() const
      {
        return N_VGetArrayPointer(q_);
      }

      real_type* getIntegral()
      {
        return N_VGetArrayPointer(q_);
      }

      const real_type* getAdjointIntegral() const
      {
        return N_VGetArrayPointer(qB_);
      }

      real_type* getAdjointIntegral()
      {
        return N_VGetArrayPointer(qB_);
      }

      void printOutput(real_type t);
      void printSpecial(real_type t, N_Vector x);
      void printFinalStats();

    private:
      static int Residual(real_type t,
                          N_Vector  yy,
                          N_Vector  yp,
                          N_Vector  rr,
                          void*     user_data);

      static int Jac(real_type t,
                     real_type cj,
                     N_Vector  yy,
                     N_Vector  yp,
                     N_Vector  resvec,
                     SUNMatrix J,
                     void*     user_data,
                     N_Vector  tmp1,
                     N_Vector  tmp2,
                     N_Vector  tmp3);

      static int Integrand(real_type t,
                           N_Vector  yy,
                           N_Vector  yp,
                           N_Vector  rhsQ,
                           void*     user_data);

      static int adjointResidual(real_type t,
                                 N_Vector  yy,
                                 N_Vector  yp,
                                 N_Vector  yyB,
                                 N_Vector  ypB,
                                 N_Vector  rrB,
                                 void*     user_data);

      static int adjointIntegrand(real_type t,
                                  N_Vector  yy,
                                  N_Vector  yp,
                                  N_Vector  yyB,
                                  N_Vector  ypB,
                                  N_Vector  rhsQB,
                                  void*     user_data);

    private:
      void*           solver_{};
      SUNContext      context_{};
      SUNMatrix       JacobianMat_{};
      SUNMatrix       JacobianMatB_{};
      SUNLinearSolver linearSolver_{};
      SUNLinearSolver linearSolverB_{};

      real_type t_init_{};
      real_type t_final_{};
      int       nout_{}; ///< Number of integration outputs

      N_Vector yy_{};  ///< Solution vector
      N_Vector yp_{};  ///< Solution derivatives vector
      N_Vector tag_{};   ///< Tags differential variables
      N_Vector q_{};   ///< Integrand vector

      N_Vector yy0_{}; ///< Storage for initial values
      N_Vector yp0_{}; ///< Storage for initial derivatives

      N_Vector yyB_{}; ///< Adjoint solution vector
      N_Vector ypB_{}; ///< Adjoint solution derivatives vector
      N_Vector qB_{};  ///< Backward integrand vector

      int backwardID_{};

    private:
      // static void copyMat(Model::Evaluator::Mat& J, SlsMat Jida);
      //TODO: should template type be real_type rather than ScalarT?
      static void copyVec(const N_Vector x, std::vector<ScalarT>& y);
      static void copyVec(const std::vector<ScalarT>& x, N_Vector y);
      static void copyVec(const std::vector<bool>& x, N_Vector y);

      // int check_flag(void *flagvalue, const char *funcname, int opt);
      inline void checkAllocation(void* v, const char* functionName);
      inline void checkOutput(int retval, const char* functionName);
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

#endif // _IDA_HPP_
