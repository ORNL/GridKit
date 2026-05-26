#pragma once

#include <iomanip>
#include <iostream>
#include <limits>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/Ieeest.hpp>
#include <GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/IeeestData.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCOO.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class StabilizerIeeestTests
    {
    public:
      using RealT = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      StabilizerIeeestTests()  = default;
      ~StabilizerIeeestTests() = default;

      TestOutcome constructor()
      {
        TestStatus success = true;

        auto  data = makeTestData();
        auto* stab = new PhasorDynamics::Stabilizer::Ieeest<ScalarT, IdxT>(data);

        success *= (stab != nullptr);
        success *= (stab->getMonitor() != nullptr);

        delete stab;

        return success.report(__func__);
      }

      /**
       * @brief All states initialize to zero (stabilizer at rest).
       * With u = 0, all residuals should be zero.
       */
      TestOutcome zeroInitialResidual()
      {
        TestStatus success = true;

        // Create signal nodes for input (u) and output (Vss)
        PhasorDynamics::SignalNode<ScalarT, IdxT> u_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> vss_node;
        ScalarT                                   u_value{0.0};
        IdxT                                      u_index = 12; // beyond internal variables
        ScalarT                                   vss_value{0.0};
        IdxT                                      vss_index = INVALID_INDEX<IdxT>;

        // Link signal nodes to backing storage
        u_node.set(&u_value, &u_index);
        vss_node.set(&vss_value, &vss_index);

        auto                                              data = makeTestData();
        PhasorDynamics::Stabilizer::Ieeest<ScalarT, IdxT> stab(data);

        // Wire: stabilizer reads u_node as input, writes vss_node as output
        stab.getSignals().template attachSignalNode<PhasorDynamics::Stabilizer::IeeestExternalVariables::U>(&u_node);
        stab.getSignals().template assignSignalNode<PhasorDynamics::Stabilizer::IeeestInternalVariables::VSS>(&vss_node);

        stab.allocate();
        success *= (stab.verify() == 0);
        stab.initialize();
        stab.evaluateResidual();

        auto        tol = 10 * std::numeric_limits<RealT>::epsilon();
        const auto& f   = stab.getResidual();
        for (size_t i = 0; i < f.size(); ++i)
        {
          if (!isEqual(f[i], 0.0, tol))
          {
            std::cout << "Non-zero residual at index " << i << ": " << f[i] << "\n";
            success = false;
          }
        }

        // Verify output signal is linked and reads the correct value
        success *= vss_node.linked();
        success *= (vss_node.getVariableIndex() == 11);
        success *= isEqual(vss_node.read(), static_cast<ScalarT>(0.0), tol);

        return success.report(__func__);
      }

      /**
       * @brief Residual evaluation against hand-computed answer key.
       *
       * Sets specific y/yp values and verifies residuals match
       * pre-computed expected values. See plan for derivation.
       */
      TestOutcome residual()
      {
        TestStatus success = true;

        PhasorDynamics::SignalNode<ScalarT, IdxT> u_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> vss_node;
        ScalarT                                   u_value{0.5};
        IdxT                                      u_index = 12;
        ScalarT                                   vss_value{0.0};
        IdxT                                      vss_index = INVALID_INDEX<IdxT>;

        u_node.set(&u_value, &u_index);
        vss_node.set(&vss_value, &vss_index);

        auto                                              data = makeTestData();
        PhasorDynamics::Stabilizer::Ieeest<ScalarT, IdxT> stab(data);

        stab.getSignals().template attachSignalNode<PhasorDynamics::Stabilizer::IeeestExternalVariables::U>(&u_node);
        stab.getSignals().template assignSignalNode<PhasorDynamics::Stabilizer::IeeestInternalVariables::VSS>(&vss_node);

        stab.allocate();
        stab.initialize();
        setStatePoint(stab);
        stab.evaluateResidual();

        // Hand-computed answer key (see plan for full derivation)
        const std::vector<ScalarT> res_answer = {
            0.19,   // f[0]:  -x1_dot + x2
            0.28,   // f[1]:  -x2_dot + x3
            0.37,   // f[2]:  -x3_dot + x4
            1.0975, // f[3]:  -x4_dot + (-a0*x1 - a1*x2 - a2*x3 - a3*x4 + u) / a4
            0.25,   // f[4]:  -T2*x5_dot - x5 + v4
            0.24,   // f[5]:  -T4*x6_dot - x6 + v5
            -0.05,  // f[6]:  -T6*x7_dot - x7 + v6
            -0.42,  // f[7]:  -v4 + x1 + A5*x2 + A6*x3
            -0.25,  // f[8]:  -T2*(v5 - x5) + T1*(v4 - x5)
            -0.31,  // f[9]:  -T4*(v6 - x6) + T3*(v5 - x6)
            5.75,   // f[10]: -T6*v7 + Ks*T5*(v6 - x7)
            0.0,    // f[11]: limiter (v7=0.05 within [-0.1, 0.1])
        };

        // Looser tolerance for f[11] — Math::clamp is a smooth ramp approximation.
        const auto loose_tol = static_cast<RealT>(1.0e-4);
        auto&      residual  = stab.getResidual();

        for (size_t i = 0; i < res_answer.size(); ++i)
        {
          auto test_tol = (i == 11) ? loose_tol : static_cast<RealT>(10 * std::numeric_limits<ScalarT>::epsilon());
          if (!isEqual(residual[i], res_answer[i], test_tol))
          {
            std::cout << "Incorrect result for residual " << i << ": "
                      << std::setprecision(15) << residual[i]
                      << " != " << res_answer[i] << "\n";
            success = false;
          }
        }

        // Verify output signal reads the stabilizer output
        success *= isEqual(vss_node.read(), static_cast<ScalarT>(0.05), loose_tol);

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /**
       * @brief Compare DependencyTracking Jacobian against Enzyme Jacobian.
       */
      TestOutcome jacobian()
      {
        TestStatus success = true;

        auto data = makeTestData();

        std::vector<DependencyTracking::Variable::DependencyMap>
            dependency_tracking_jacobian = DependencyTrackingJacobian(data);

        std::vector<DependencyTracking::Variable::DependencyMap>
            enzyme_jacobian = EnzymeJacobian(data);

        // Compare DependencyTracking dependencies to Enzyme's
        auto tol = 10 * std::numeric_limits<RealT>::epsilon();
        for (size_t i = 0; i < dependency_tracking_jacobian.size(); ++i)
        {
          success *= (GridKit::Testing::isEqual(dependency_tracking_jacobian[i], enzyme_jacobian[i], tol));
        }

        return success.report(__func__);
      }

    private:
      std::vector<DependencyTracking::Variable::DependencyMap> DependencyTrackingJacobian(
          PhasorDynamics::Stabilizer::IeeestData<RealT, IdxT> ieeestdata)
      {
        using DepVar = DependencyTracking::Variable;

        // Set up signal nodes with DependencyTracking scalar type
        PhasorDynamics::SignalNode<DepVar, IdxT> u_node;
        PhasorDynamics::SignalNode<DepVar, IdxT> vss_node;
        DepVar                                   u_value{0.5};
        IdxT                                     u_index = 12;
        DepVar                                   vss_value{0.0};
        IdxT                                     vss_index = INVALID_INDEX<IdxT>;

        u_node.set(&u_value, &u_index);
        vss_node.set(&vss_value, &vss_index);

        PhasorDynamics::Stabilizer::Ieeest<DepVar, IdxT> stab(ieeestdata);
        stab.getSignals().template attachSignalNode<PhasorDynamics::Stabilizer::IeeestExternalVariables::U>(&u_node);
        stab.getSignals().template assignSignalNode<PhasorDynamics::Stabilizer::IeeestInternalVariables::VSS>(&vss_node);

        stab.allocate();
        stab.initialize();

        // --- d/dy: tag internal variables as independent ---
        for (size_t i = 0; i < stab.size(); ++i)
        {
          stab.y()[i].setVariableNumber(i);
        }
        // Tag external signal u as an additional independent variable
        u_value.setVariableNumber(stab.size());
        u_value.setValue(0.5);

        setStatePointDep(stab);

        stab.evaluateResidual();
        std::vector<DepVar> residual_y = stab.getResidual();

        // --- d/dy': tag derivatives as independent ---
        stab.initialize();
        for (size_t i = 0; i < stab.size(); ++i)
        {
          stab.yp()[i].setVariableNumber(i);
        }

        u_value = 0.5;
        setStatePointDep(stab);

        stab.evaluateResidual();
        std::vector<DepVar> residual_yp = stab.getResidual();

        // Print dependencies for debugging
        for (size_t i = 0; i < residual_y.size(); ++i)
        {
          std::cout << i << "th residual, y: ";
          (residual_y[i]).print(std::cout);
          std::cout << "\n";
          std::cout << i << "th residual, yp: ";
          (residual_yp[i]).print(std::cout);
          std::cout << "\n";
        }

        // Merge d/dy and d/dy' into a single dependency map
        std::vector<DependencyTracking::Variable::DependencyMap> dependencies(residual_y.size());
        for (IdxT i = 0; i < residual_y.size(); ++i)
        {
          auto dependency_y  = (residual_y[i]).getDependencies();
          auto dependency_yp = (residual_yp[i]).getDependencies();

          for (const auto& pair_y : dependency_y)
          {
            auto it_yp = dependency_yp.find(pair_y.first);
            if (it_yp != dependency_yp.end())
            {
              dependencies[i].insert(std::make_pair(pair_y.first, pair_y.second + it_yp->second));
            }
            else
            {
              dependencies[i].insert(std::make_pair(pair_y.first, pair_y.second));
            }
          }

          // Insert yp dependencies that did not exist in the y dependencies
          for (const auto& pair_yp : dependency_yp)
          {
            if (!dependency_y.contains(pair_yp.first))
            {
              dependencies[i].insert(std::make_pair(pair_yp.first, pair_yp.second));
            }
          }
        }

        return dependencies;
      }

      std::vector<DependencyTracking::Variable::DependencyMap> EnzymeJacobian(
          PhasorDynamics::Stabilizer::IeeestData<RealT, IdxT> ieeestdata)
      {
        PhasorDynamics::SignalNode<ScalarT, IdxT> u_node;
        PhasorDynamics::SignalNode<ScalarT, IdxT> vss_node;
        ScalarT                                   u_value{0.5};
        IdxT                                      u_index = 12;
        ScalarT                                   vss_value{0.0};
        IdxT                                      vss_index = INVALID_INDEX<IdxT>;

        u_node.set(&u_value, &u_index);
        vss_node.set(&vss_value, &vss_index);

        PhasorDynamics::Stabilizer::Ieeest<ScalarT, IdxT> stab(ieeestdata);
        stab.getSignals().template attachSignalNode<PhasorDynamics::Stabilizer::IeeestExternalVariables::U>(&u_node);
        stab.getSignals().template assignSignalNode<PhasorDynamics::Stabilizer::IeeestInternalVariables::VSS>(&vss_node);

        stab.allocate();
        stab.initialize();
        setStatePoint(stab);

        stab.updateTime(0.0, 1.0); // alpha = 1.0 to verify d/dy' term

        stab.evaluateResidual();
        stab.evaluateJacobian();

        auto model_jacobian = stab.getJacobian();
        model_jacobian.deduplicate();
        model_jacobian.printMatrix("IEEEST Jacobian");

        return GridKit::Testing::MapFromCOO(model_jacobian);
      }
#endif

    private:
      static constexpr ScalarT tol_ = 10 * std::numeric_limits<ScalarT>::epsilon();

      /**
       * @brief Standard IEEEST parameter set for all tests.
       * Derived: a0=1, a1=0.4, a2=0.63, a3=0.1, a4=0.08
       */
      auto makeTestData() -> PhasorDynamics::Stabilizer::IeeestData<RealT, IdxT>
      {
        using Params = PhasorDynamics::Stabilizer::IeeestParameters;

        PhasorDynamics::Stabilizer::IeeestData<RealT, IdxT> data;
        data.device_class          = "stabilizer";
        data.disambiguation_string = "ieeest_test";
        data.monitored_variables.insert(PhasorDynamics::Stabilizer::IeeestMonitorableVariables::vss);

        data.parameters[Params::A1]     = 0.1;
        data.parameters[Params::A2]     = 0.2;
        data.parameters[Params::A3]     = 0.3;
        data.parameters[Params::A4]     = 0.4;
        data.parameters[Params::A5]     = 0.5;
        data.parameters[Params::A6]     = 0.6;
        data.parameters[Params::T1]     = 0.5;
        data.parameters[Params::T2]     = 1.0;
        data.parameters[Params::T3]     = 0.3;
        data.parameters[Params::T4]     = 1.0;
        data.parameters[Params::T5]     = 2.0;
        data.parameters[Params::T6]     = 5.0;
        data.parameters[Params::Ks]     = 10.0;
        data.parameters[Params::Lsmin]  = -0.1;
        data.parameters[Params::Lsmax]  = 0.1;
        data.parameters[Params::Vcl]    = 0.0;
        data.parameters[Params::Vcu]    = 0.0;
        data.parameters[Params::Tdelay] = 0.0;

        return data;
      }

      /**
       * @brief Set a non-trivial operating point for residual/Jacobian tests.
       * Avoids zeros and ones to catch coefficient errors.
       */
      void setStatePoint(PhasorDynamics::Stabilizer::Ieeest<ScalarT, IdxT>& stab)
      {
        stab.y()[0]  = 0.1;  // x1
        stab.y()[1]  = 0.2;  // x2
        stab.y()[2]  = 0.3;  // x3
        stab.y()[3]  = 0.4;  // x4
        stab.y()[4]  = 0.5;  // x5
        stab.y()[5]  = 0.6;  // x6
        stab.y()[6]  = 0.7;  // x7
        stab.y()[7]  = 0.8;  // v4
        stab.y()[8]  = 0.9;  // v5
        stab.y()[9]  = 1.0;  // v6
        stab.y()[10] = 0.05; // v7  (within limiter range)
        stab.y()[11] = 0.05; // Vss (model output)

        stab.yp()[0] = 0.01; // x1_dot
        stab.yp()[1] = 0.02; // x2_dot
        stab.yp()[2] = 0.03; // x3_dot
        stab.yp()[3] = 0.04; // x4_dot
        stab.yp()[4] = 0.05; // x5_dot
        stab.yp()[5] = 0.06; // x6_dot
        stab.yp()[6] = 0.07; // x7_dot
      }

      /**
       * @brief Set the same operating point for DependencyTracking variables.
       * Uses setValue() to set the numeric value while preserving dependency info.
       */
      void setStatePointDep(PhasorDynamics::Stabilizer::Ieeest<DependencyTracking::Variable, IdxT>& stab)
      {
        stab.y()[0].setValue(0.1);
        stab.y()[1].setValue(0.2);
        stab.y()[2].setValue(0.3);
        stab.y()[3].setValue(0.4);
        stab.y()[4].setValue(0.5);
        stab.y()[5].setValue(0.6);
        stab.y()[6].setValue(0.7);
        stab.y()[7].setValue(0.8);
        stab.y()[8].setValue(0.9);
        stab.y()[9].setValue(1.0);
        stab.y()[10].setValue(0.05);
        stab.y()[11].setValue(0.05);

        stab.yp()[0].setValue(0.01);
        stab.yp()[1].setValue(0.02);
        stab.yp()[2].setValue(0.03);
        stab.yp()[3].setValue(0.04);
        stab.yp()[4].setValue(0.05);
        stab.yp()[5].setValue(0.06);
        stab.yp()[6].setValue(0.07);
      }
    }; // class StabilizerIeeestTests

  } // namespace Testing
} // namespace GridKit
