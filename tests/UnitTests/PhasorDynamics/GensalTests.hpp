#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENSALwS/Gensal.hpp>
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
    class GensalTests
    {
    private:
      using RealT                   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;
      using GensalDataT             = PhasorDynamics::GensalData<RealT, IdxT>;
      static constexpr ScalarT tol_ = 1000 * std::numeric_limits<ScalarT>::epsilon();

      static GensalDataT makeGensalData()
      {
        using Parameter = typename GensalDataT::Parameters;
        using Port      = typename GensalDataT::Ports;

        GensalDataT data;
        data.device_class                 = "Gensal";
        data.disambiguation_string        = "1";
        data.ports[Port::bus]             = 1;
        data.parameters[Parameter::p0]    = RealT{1.0};
        data.parameters[Parameter::q0]    = RealT{0.05013};
        data.parameters[Parameter::H]     = RealT{3.0};
        data.parameters[Parameter::D]     = RealT{0.0};
        data.parameters[Parameter::Ra]    = RealT{0.0};
        data.parameters[Parameter::Tdop]  = RealT{7.0};
        data.parameters[Parameter::Tdopp] = RealT{0.04};
        data.parameters[Parameter::Tqopp] = RealT{0.05};
        data.parameters[Parameter::Xd]    = RealT{2.1};
        data.parameters[Parameter::Xdp]   = RealT{0.2};
        data.parameters[Parameter::Xdpp]  = RealT{0.18};
        data.parameters[Parameter::Xq]    = RealT{0.5};
        data.parameters[Parameter::Xl]    = RealT{0.15};
        data.parameters[Parameter::S10]   = RealT{0.0};
        data.parameters[Parameter::S12]   = RealT{0.0};

        return data;
      }

    public:
      GensalTests()  = default;
      ~GensalTests() = default;

      TestOutcome constructor()
      {
        TestStatus success = true;

        auto* bus  = new PhasorDynamics::Bus<ScalarT, IdxT>(1.0, 0.0);
        auto  data = makeGensalData();

        PhasorDynamics::Component<ScalarT, IdxT>* machine =
            new PhasorDynamics::Gensal<ScalarT, IdxT>(bus, data);

        success *= (machine != nullptr);

        if (machine)
        {
          delete machine;
        }
        delete bus;

        return success.report(__func__);
      }

      /**
       * @brief Checks residual evaluation at initialized steady state.
       */
      TestOutcome residual()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT>    bus(1.0, 0.0);
        auto                                  data = makeGensalData();
        PhasorDynamics::Gensal<ScalarT, IdxT> gen(&bus, data);

        bus.allocate();
        bus.initialize();
        bus.evaluateResidual();

        gen.allocate();
        gen.initialize();
        gen.evaluateResidual();

        const std::vector<ScalarT>& f = gen.getResidual();
        for (const auto& f_val : f)
        {
          if (!isEqual(f_val, 0.0, tol_))
          {
            success = false;
            break;
          }
        }

        return success.report(__func__);
      }

      /**
       * @brief Checks initialized steady state with nonzero armature resistance.
       */
      TestOutcome residual_nonzero_ra()
      {
        TestStatus success = true;

        using Parameter = typename GensalDataT::Parameters;

        auto data                       = makeGensalData();
        data.parameters[Parameter::p0]  = RealT{0.8};
        data.parameters[Parameter::q0]  = RealT{0.2};
        data.parameters[Parameter::Ra]  = RealT{0.05};
        data.parameters[Parameter::S10] = RealT{0.0};
        data.parameters[Parameter::S12] = RealT{0.0};

        PhasorDynamics::Bus<ScalarT, IdxT>    bus(1.0, 0.1);
        PhasorDynamics::Gensal<ScalarT, IdxT> gen(&bus, data);

        bus.allocate();
        bus.initialize();
        bus.evaluateResidual();

        gen.allocate();
        gen.initialize();
        gen.evaluateResidual();

        const std::vector<ScalarT>& f = gen.getResidual();
        for (const auto& f_val : f)
        {
          if (!isEqual(f_val, 0.0, tol_))
          {
            success = false;
            break;
          }
        }

        return success.report(__func__);
      }

      /**
       * @brief Checks GENSAL uses the configured system frequency base.
       */
      TestOutcome frequency_base()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT>    bus(1.0, 0.0);
        auto                                  data = makeGensalData();
        PhasorDynamics::Gensal<ScalarT, IdxT> gen(&bus, data);

        bus.allocate();
        bus.initialize();

        gen.setSystemBase(50.0, 100.0e6);
        gen.allocate();

        gen.y()[1]  = 1.0;
        gen.yp()[0] = TWO<RealT> * M_PI * 50.0;

        gen.evaluateResidual();

        success *= isEqual(gen.getResidual()[0], 0.0, tol_);

        return success.report(__func__);
      }

      /**
       * @brief Checks monitored terminal current and power use system base.
       */
      TestOutcome monitor_system_base()
      {
        TestStatus success = true;

        using Parameter = typename GensalDataT::Parameters;
        using Variable  = typename GensalDataT::MonitorableVariables;

        auto data                            = makeGensalData();
        data.parameters[Parameter::mva_base] = RealT{50.0};
        data.monitored_variables.insert(Variable::ir);
        data.monitored_variables.insert(Variable::p);

        PhasorDynamics::Bus<ScalarT, IdxT>    bus(1.0, 0.0);
        PhasorDynamics::Gensal<ScalarT, IdxT> gen(&bus, data);

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
       * @brief Verifies residual equations against hard-coded values.
       */
      TestOutcome hard_coded_residual()
      {
        TestStatus success = true;

        using Parameter = typename GensalDataT::Parameters;

        auto data                         = makeGensalData();
        data.parameters[Parameter::p0]    = RealT{1.0};
        data.parameters[Parameter::q0]    = RealT{1.0};
        data.parameters[Parameter::H]     = RealT{0.5};
        data.parameters[Parameter::D]     = RealT{-1.0};
        data.parameters[Parameter::Ra]    = RealT{0.5};
        data.parameters[Parameter::Tdop]  = RealT{2.0};
        data.parameters[Parameter::Tdopp] = RealT{4.0};
        data.parameters[Parameter::Tqopp] = RealT{5.0};
        data.parameters[Parameter::Xd]    = RealT{2.0};
        data.parameters[Parameter::Xdp]   = RealT{1.0};
        data.parameters[Parameter::Xdpp]  = RealT{0.5};
        data.parameters[Parameter::Xq]    = RealT{1.5};
        data.parameters[Parameter::Xl]    = RealT{0.25};
        data.parameters[Parameter::S10]   = RealT{0.0};
        data.parameters[Parameter::S12]   = RealT{0.0};

        ScalarT Vr1{1.0};
        ScalarT Vi1{1.0};

        PhasorDynamics::Bus<ScalarT, IdxT>    bus(Vr1, Vi1);
        PhasorDynamics::Gensal<ScalarT, IdxT> gen(&bus, data);

        const std::vector<ScalarT> res_answer = {
            0.0,
            0.0,
            2.1083333333333334,
            -1.028125,
            0.65,
            0.0,
            0.2,
            -1.1,
            -1.4,
            1.8125,
            0.5,
            0.25,
            2.65,
            -0.05,
            0.3,
            -1.2};

        bus.allocate();
        bus.initialize();

        gen.allocate();

        gen.y()[0]  = M_PI;  // delta
        gen.y()[1]  = 1.0;   // omega
        gen.y()[2]  = 2.0;   // Eqp
        gen.y()[3]  = 0.5;   // psidp
        gen.y()[4]  = -0.75; // psiqpp
        gen.y()[5]  = 1.0;   // psidpp
        gen.y()[6]  = 0.2;   // ksat
        gen.y()[7]  = 0.4;   // vd
        gen.y()[8]  = 0.6;   // vq
        gen.y()[9]  = 1.5;   // telec
        gen.y()[10] = 0.25;  // id
        gen.y()[11] = -0.5;  // iq
        gen.y()[12] = 0.75;  // ir
        gen.y()[13] = -0.25; // ii
        gen.y()[14] = 0.1;   // inr
        gen.y()[15] = -0.2;  // ini

        gen.yp()[0] = 2 * M_PI * 60.0; // delta_dot
        gen.yp()[1] = -1.0;            // omega_dot
        gen.yp()[2] = 0.3;             // Eqp_dot
        gen.yp()[3] = -0.7;            // psidp_dot
        gen.yp()[4] = 0.9;             // psiqpp_dot

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

#ifdef GRIDKIT_ENABLE_ENZYME
      /**
       * @brief Checks Jacobian evaluation.
       */
      TestOutcome jacobian()
      {
        TestStatus success = true;

        auto tol = 10 * std::numeric_limits<RealT>::epsilon();

        std::vector<DependencyTracking::Variable::DependencyMap> dependency_tracking_jacobian = DependencyTrackingJacobian();
        std::vector<DependencyTracking::Variable::DependencyMap> enzyme_jacobian              = EnzymeJacobian();

        for (size_t i = 0; i < dependency_tracking_jacobian.size(); ++i)
        {
          success *= (GridKit::Testing::isEqual(dependency_tracking_jacobian[i], enzyme_jacobian[i], tol));
        }

        return success.report(__func__);
      }

    private:
      std::vector<DependencyTracking::Variable::DependencyMap> DependencyTrackingJacobian()
      {
        DependencyTracking::Variable                               Vr1{1.0};
        DependencyTracking::Variable                               Vi1{0.0};
        PhasorDynamics::Bus<DependencyTracking::Variable, IdxT>    bus(Vr1, Vi1);
        auto                                                       data = makeGensalData();
        PhasorDynamics::Gensal<DependencyTracking::Variable, IdxT> gen(&bus, data);

        bus.allocate();
        gen.allocate();

        bus.initialize();
        gen.initialize();

        for (size_t i = 0; i < gen.size(); ++i)
        {
          gen.y()[i].setVariableNumber(i);
        }
        for (size_t i = 0; i < bus.size(); ++i)
        {
          bus.y()[i].setVariableNumber(i + gen.size());
        }

        bus.evaluateResidual();
        gen.evaluateResidual();
        std::vector<DependencyTracking::Variable> residual_y = gen.getResidual();

        bus.initialize();
        gen.initialize();

        for (size_t i = 0; i < gen.size(); ++i)
        {
          gen.yp()[i].setVariableNumber(i);
        }

        bus.evaluateResidual();
        gen.evaluateResidual();
        std::vector<DependencyTracking::Variable> residual_yp = gen.getResidual();

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

      std::vector<DependencyTracking::Variable::DependencyMap> EnzymeJacobian()
      {
        ScalarT                               Vr1{1.0};
        ScalarT                               Vi1{0.0};
        PhasorDynamics::Bus<ScalarT, IdxT>    bus(Vr1, Vi1);
        auto                                  data = makeGensalData();
        PhasorDynamics::Gensal<ScalarT, IdxT> gen(&bus, data);

        bus.allocate();
        gen.allocate();

        bus.initialize();
        gen.initialize();

        gen.updateTime(0.0, 1.0);

        for (size_t i = 0; i < bus.size(); ++i)
        {
          bus.setVariableIndex(i, i + gen.size());
          bus.setResidualIndex(i, i + gen.size());
        }

        bus.evaluateResidual();
        gen.evaluateResidual();

        bus.evaluateJacobian();
        gen.evaluateJacobian();
        GridKit::LinearAlgebra::COO_Matrix<ScalarT, IdxT>& model_jacobian = gen.getJacobian();
        model_jacobian.deduplicate();

        return GridKit::Testing::MapFromCOO(model_jacobian);
      }
#endif
    }; // class GensalTests

  } // namespace Testing
} // namespace GridKit
