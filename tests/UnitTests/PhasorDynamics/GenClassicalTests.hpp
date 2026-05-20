/**
 * @file GenClassicalTests.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @author Abdourahman Barry (abdourahman@vt.edu)
 * @brief Tests for classical generator model.
 *
 */
#define _USE_MATH_DEFINES /* need this since directly including GenClassical.cpp for MSVC compiler */
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GenClassical/GenClassical.hpp>
#include <GridKit/Model/VariableMonitorController.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Testing/Tokenizer.hpp>
#include <GridKit/Utilities/MapFromCOO.hpp>

namespace GridKit
{
  namespace Testing
  {

    template <class ScalarT, typename IdxT>
    class GenClassicalTests
    {
    private:
      using RealT                   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;
      using GenClassicalDataT       = PhasorDynamics::GenClassicalData<RealT, IdxT>;
      static constexpr ScalarT tol_ = 10 * std::numeric_limits<ScalarT>::epsilon();

      static GenClassicalDataT makeGenClassicalData()
      {
        using Parameter = typename GenClassicalDataT::Parameters;
        using Port      = typename GenClassicalDataT::Ports;

        GenClassicalDataT data;
        data.device_class               = "GenClassical";
        data.disambiguation_string      = "1";
        data.ports[Port::bus]           = 1;
        data.parameters[Parameter::p0]  = RealT{1.0};
        data.parameters[Parameter::q0]  = RealT{0.0};
        data.parameters[Parameter::H]   = RealT{0.5};
        data.parameters[Parameter::D]   = RealT{0.0};
        data.parameters[Parameter::Ra]  = RealT{0.0};
        data.parameters[Parameter::Xdp] = RealT{0.5};

        return data;
      }

    public:
      GenClassicalTests()  = default;
      ~GenClassicalTests() = default;

      TestOutcome constructor()
      {
        TestStatus success = true;

        auto* bus = new PhasorDynamics::Bus<ScalarT, IdxT>(1.0, 0.0);

        PhasorDynamics::Component<ScalarT, IdxT>* machine =
            new PhasorDynamics::GenClassical<ScalarT, IdxT>(bus, 1);

        success *= (machine != nullptr);

        if (machine)
        {
          delete machine;
        }
        delete bus;

        return success.report(__func__);
      }

      /**
       * A test case to verify residual values
       */
      TestOutcome residual()
      {
        TestStatus success = true;

        // Classical generator parameters
        RealT H{0.5};
        RealT D{-1.0};
        RealT Ra{0.5};
        RealT Xdp{0.5};

        // Classical generator inputs
        RealT Pm{1.0};
        RealT Ep{2.0};

        ScalarT Vr1{1.0}; ///< Bus-1 real voltage
        ScalarT Vi1{1.0}; ///< Bus-1 imaginary voltage

        // Test answer keys
        const std::vector<ScalarT> res_answer = {0.0,
                                                 -0.5,
                                                 -6.0,
                                                 2.0,
                                                 -6.0};

        PhasorDynamics::Bus<ScalarT, IdxT>          bus(Vr1, Vi1);
        PhasorDynamics::GenClassical<ScalarT, IdxT> gen(&bus, 1, 1.0, 1.0, H, D, Ra, Xdp);
        bus.allocate();
        bus.initialize();

        // Allocate but not initialize genrator model
        gen.allocate();
        gen.setPmech(Pm);
        gen.setEp(Ep);

        // Set variable values matching the answer key
        gen.y()[0] = M_PI; // delta
        gen.y()[1] = 1.0;  // omega
        gen.y()[2] = 2.0;  // telec
        gen.y()[3] = -2.0; // ir
        gen.y()[4] = -4.0; // ii

        // Set derivative values matching the answer key
        gen.yp()[0] = 2 * M_PI * 60.0; // delta_dot
        gen.yp()[1] = -1.5;            // omega_dot
        gen.yp()[2] = 0;
        gen.yp()[3] = 0;
        gen.yp()[4] = 0;

        gen.evaluateResidual();
        std::vector<ScalarT>& residual = gen.getResidual();

        for (size_t i = 0; i < res_answer.size(); ++i)
        {
          if (!isEqual(residual[i], res_answer[i], tol_))
          {
            std::cout << "Incorrect result for residual " << i << ": "
                      << residual[i] << " != " << res_answer[i] << "\n";
            success = false;
            break;
          }
        }

        return success.report(__func__);
      }

      /**
       * @brief Checks monitored terminal current and power use system base.
       */
      TestOutcome monitor_system_base()
      {
        TestStatus success = true;

        using Parameter = typename GenClassicalDataT::Parameters;
        using Variable  = typename GenClassicalDataT::MonitorableVariables;

        auto data                            = makeGenClassicalData();
        data.parameters[Parameter::mva_base] = RealT{50.0};
        data.monitored_variables.insert(Variable::ir);
        data.monitored_variables.insert(Variable::p);

        PhasorDynamics::Bus<ScalarT, IdxT>          bus(1.0, 0.0);
        PhasorDynamics::GenClassical<ScalarT, IdxT> gen(&bus, data);

        bus.allocate();
        bus.initialize();
        bus.evaluateResidual();

        gen.setSystemBase(60.0, 100.0e6);
        gen.allocate();
        gen.initialize();
        gen.evaluateResidual();

        RealT                                     time = 0.0;
        Model::VariableMonitorController<ScalarT> monitor(time);
        monitor.addMonitor(gen.getMonitor());

        std::stringstream os;
        monitor.print(os, Model::VariableMonitorBase::Csv{});

        auto values = Tokenizer<RealT>(os.str(), ',')();
        if (values.size() == 3)
        {
          success *= isEqual(values[1], 1.0, tol_);
          success *= isEqual(values[2], 1.0, tol_);
        }
        else
        {
          success = false;
        }

        return success.report(__func__);
      }

      /**
       *
       * Verifies correctness of the system initialization
       */
      TestOutcome initial()
      {
        TestStatus success = true;

        // Classical generator parameters
        RealT p0{3.0};
        RealT q0{-1.0};
        RealT H{1.0};
        RealT D{1.0};
        RealT Ra{0.1};
        RealT Xdp{2.3};

        ScalarT Vr1{1.0}; ///< Bus-1 real voltage
        ScalarT Vi1{1.0}; ///< Bus-1 imaginary voltage

        // Test answer keys
        const std::vector<ScalarT> var_answer = {
            3.0 * M_PI / 4.0, // delta
            0.0,              // omega
            3.5,              // Te
            1.0,              // Ir
            2.0,              // Ii
        };

        PhasorDynamics::Bus<ScalarT, IdxT>          bus(Vr1, Vi1);
        PhasorDynamics::GenClassical<ScalarT, IdxT> gen(&bus, 1, p0, q0, H, D, Ra, Xdp);
        bus.allocate();
        bus.initialize();
        gen.allocate();
        gen.initialize();

        for (size_t i = 0; i < var_answer.size(); ++i)
        {
          if (!isEqual(gen.y()[i], var_answer[i], tol_))
          {
            std::cout << "Incorrect result: "
                      << gen.y()[i] << " != " << var_answer[i] << "\n";
            success = false;
            break;
          }

          if (!isEqual(gen.yp()[i], 0.0, tol_))
          {
            std::cout << "Incorrect result: "
                      << gen.yp()[i] << " != 0\n";
            success = false;
            break;
          }
        }

        return success.report(__func__);
      }

      /**
       * Verifies the residual evaluates to zero for the initial conditions
       */
      TestOutcome zeroInitialResidual()
      {
        TestStatus success = true;

        // Classical generator parameters
        RealT p0{3.0};
        RealT q0{-1.0};
        RealT H{1.0};
        RealT D{1.0};
        RealT Ra{0.6};
        RealT Xdp{0.2};

        ScalarT Vr1{1.0}; ///< Bus real voltage
        ScalarT Vi1{1.0}; ///< Bus imaginary voltage

        PhasorDynamics::Bus<ScalarT, IdxT>          bus(Vr1, Vi1);
        PhasorDynamics::GenClassical<ScalarT, IdxT> gen(&bus, 1, p0, q0, H, D, Ra, Xdp);
        bus.allocate();
        bus.initialize();
        gen.allocate();
        gen.initialize();
        gen.evaluateResidual();
        std::vector<ScalarT> res = gen.getResidual();

        for (size_t i = 0; i < res.size(); ++i)
        {
          if (!isEqual(res[i], 0.0, tol_))
          {
            std::cout << "Incorrect result: "
                      << gen.yp()[i] << " != 0\n";
            success = false;
            break;
          }
        }

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /**
       * A test case to verify Jacobian values
       */
      TestOutcome jacobian()
      {
        TestStatus success = true;

        // Classical generator parameters
        RealT H{0.5};
        RealT D{-1.0};
        RealT Ra{0.5};
        RealT Xdp{0.5};

        // Jacobian via DependencyTracking
        std::vector<DependencyTracking::Variable::DependencyMap> dependency_tracking_jacobian = DependencyTrackingJacobian(H, D, Ra, Xdp);

        // Jacobian via Enzyme
        std::vector<DependencyTracking::Variable::DependencyMap> enzyme_jacobian = EnzymeJacobian(H, D, Ra, Xdp);

        /// Compare DependencyTracking dependencies to Enzyme's
        for (size_t i = 0; i < dependency_tracking_jacobian.size(); ++i)
        {
          success *= (GridKit::Testing::isEqual(dependency_tracking_jacobian[i], enzyme_jacobian[i]));
        }

        return success.report(__func__);
      }

    private:
      std::vector<DependencyTracking::Variable::DependencyMap> DependencyTrackingJacobian(
          const RealT H, const RealT D, const RealT Ra, const RealT Xdp)
      {
        DependencyTracking::Variable Vr1{1.0}; ///< Bus-1 real voltage
        DependencyTracking::Variable Vi1{1.0}; ///< Bus-1 imaginary voltage

        PhasorDynamics::Bus<DependencyTracking::Variable, IdxT>          bus(Vr1, Vi1);
        PhasorDynamics::GenClassical<DependencyTracking::Variable, IdxT> gen(&bus, 1, 1.0, 1.0, H, D, Ra, Xdp);

        bus.allocate();
        gen.allocate();

        // Get d/dy
        bus.initialize();
        gen.initialize();

        for (size_t i = 0; i < gen.size(); ++i)
        {
          gen.y()[i].setVariableNumber(i); ///< Generator independent variables
        }
        for (size_t i = 0; i < bus.size(); ++i)
        {
          bus.y()[i].setVariableNumber(i + gen.size()); // Bus independent variables
        }

        bus.evaluateResidual();
        gen.evaluateResidual(); ///< Computes the residual and the Jacobian values by tracking
                                ///< the dependencies
        std::vector<DependencyTracking::Variable> residual_y = gen.getResidual();

        // Get d/dy'
        bus.initialize();
        gen.initialize();

        for (size_t i = 0; i < gen.size(); ++i)
        {
          gen.yp()[i].setVariableNumber(i); ///< Generator independent variables
        }

        bus.evaluateResidual();
        gen.evaluateResidual(); ///< Computes the residual and the Jacobian values by tracking
                                ///< the dependencies
        std::vector<DependencyTracking::Variable> residual_yp = gen.getResidual();

        // Print the dependencies
        for (size_t i = 0; i < residual_y.size(); ++i)
        {
          std::cout << i << "th residual, y: ";
          (residual_y[i]).print(std::cout);
          std::cout << "\n";
          std::cout << i << "th residual, yp: ";
          (residual_yp[i]).print(std::cout);
          std::cout << "\n";
        }

        // Extract the dependencies and add d/dy' to d/dy
        std::vector<DependencyTracking::Variable::DependencyMap> dependencies(residual_y.size());
        for (IdxT i = 0; i < residual_y.size(); ++i)
        {
          DependencyTracking::Variable::DependencyMap dependency_y  = (residual_y[i]).getDependencies();
          DependencyTracking::Variable::DependencyMap dependency_yp = (residual_yp[i]).getDependencies();

          for (const auto& pair_y : dependency_y)
          {
            auto index_y = pair_y.first;
            auto value_y = pair_y.second;
            auto it_yp   = dependency_yp.find(index_y);
            if (it_yp != dependency_yp.end())
            {
              auto value_yp = it_yp->second;
              dependencies[i].insert(std::make_pair(index_y, value_y + value_yp));
            }
            else
            {
              dependencies[i].insert(std::make_pair(index_y, value_y));
            }
          }

          // Insert yp dependencies that did not exist in the y dependencies
          for (const auto& pair_yp : dependency_yp)
          {
            auto index_yp = pair_yp.first;
            auto value_yp = pair_yp.second;
            auto it_y     = dependency_y.find(index_yp);
            if (it_y == dependency_y.end())
            {
              dependencies[i].insert(std::make_pair(index_yp, value_yp));
            }
          }
        }

        return dependencies;
      }

      std::vector<DependencyTracking::Variable::DependencyMap> EnzymeJacobian(
          const RealT H, const RealT D, const RealT Ra, const RealT Xdp)
      {
        ScalarT Vr1{1.0}; ///< Bus-1 real voltage
        ScalarT Vi1{1.0}; ///< Bus-1 imaginary voltage

        PhasorDynamics::Bus<ScalarT, IdxT>          bus(Vr1, Vi1);
        PhasorDynamics::GenClassical<ScalarT, IdxT> gen(&bus, 1, 1.0, 1.0, H, D, Ra, Xdp);

        bus.allocate();
        gen.allocate();

        bus.initialize();
        gen.initialize();

        gen.updateTime(0.0, 1.0);

        for (size_t i = 0; i < bus.size(); ++i)
        {
          bus.setVariableIndex(i, i + gen.size()); // Reset bus variable indices
          bus.setResidualIndex(i, i + gen.size()); // Reset bus residual indices
        }

        bus.evaluateResidual();
        gen.evaluateResidual();

        bus.evaluateJacobian();
        gen.evaluateJacobian();
        GridKit::LinearAlgebra::COO_Matrix<ScalarT, IdxT>& model_jacobian = gen.getJacobian();
        model_jacobian.deduplicate();
        model_jacobian.printMatrix("Model Jacobian");

        return GridKit::Testing::MapFromCOO(model_jacobian);
      }
#endif

    }; // class GenClassicalTests

  } // namespace Testing
} // namespace GridKit
