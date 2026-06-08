#pragma once

#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

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
      using DataT = PhasorDynamics::Stabilizer::IeeestData<RealT, IdxT>;

      StabilizerIeeestTests()  = default;
      ~StabilizerIeeestTests() = default;

      TestOutcome init()
      {
        TestStatus success = true;

        using Params = PhasorDynamics::Stabilizer::IeeestParameters;

        struct InitCase
        {
          RealT   T6;
          RealT   Ks;
          ScalarT u;
          RealT   Lsmin;
          RealT   Lsmax;
          ScalarT raw_v7;
          ScalarT expected_vss;
        };

        const auto                  loose_tol = static_cast<RealT>(1.0e-4);
        const std::vector<InitCase> cases     = {
            {0.0, 0.0, 0.25, -1.0, 1.0, 0.0, 0.0},
            {0.0, 1.0, 0.25, -1.0, 1.0, 0.25, 0.25},
            {0.0, 2.0, 0.25, -1.0, 1.0, 0.50, 0.50},
            {0.0, 4.0, 0.25, -1.0, 0.6, 1.00, 0.60},
            {5.0, 3.0, 0.25, -1.0, 1.0, 0.00, 0.00},
        };

        for (const auto& test : cases)
        {
          PhasorDynamics::SignalNode<ScalarT, IdxT> u_node;
          PhasorDynamics::SignalNode<ScalarT, IdxT> vss_node;
          ScalarT                                   u_value{test.u};
          IdxT                                      u_index{12};
          ScalarT                                   vss_value{0.0};
          IdxT                                      vss_index{INVALID_INDEX<IdxT>};

          u_node.set(&u_value, &u_index);
          vss_node.set(&vss_value, &vss_index);

          auto data                      = makeData();
          data.parameters[Params::T6]    = test.T6;
          data.parameters[Params::Ks]    = test.Ks;
          data.parameters[Params::Lsmin] = test.Lsmin;
          data.parameters[Params::Lsmax] = test.Lsmax;

          PhasorDynamics::Stabilizer::Ieeest<ScalarT, IdxT> model(data);
          model.getSignals().template attachSignalNode<PhasorDynamics::Stabilizer::IeeestExternalVariables::U>(&u_node);
          model.getSignals().template assignSignalNode<PhasorDynamics::Stabilizer::IeeestInternalVariables::VSS>(&vss_node);

          model.allocate();
          success *= (model.verify() == 0);
          model.initialize();

          success *= vss_node.linked();
          success *= (vss_node.getVariableIndex() == 11);
          success *= isEqual(model.y()[10], test.raw_v7, tol_);
          success *= isEqual(model.y()[11], test.expected_vss, loose_tol);
          success *= isEqual(vss_node.read(), test.expected_vss, loose_tol);
        }

        return success.report(__func__);
      }

      TestOutcome residual()
      {
        TestStatus success = true;

        using Params = PhasorDynamics::Stabilizer::IeeestParameters;

        struct ResidualCase
        {
          const char*                 name;
          std::function<void(DataT&)> edit;
          const std::vector<ScalarT>  expected;
        };

        const std::vector<ResidualCase> cases = {
            {"baseline",
             [](DataT&) {},
             {0.19, 0.28, 0.37, 0.0878, 0.25, 0.24, -0.05, -0.42, -0.25, -0.31, 5.75, 0.0}},
            {"a4_zero",
             [](DataT& data)
             {
               data.parameters[Params::A4] = static_cast<RealT>(0.0);
             },
             {0.19, 0.28, 0.37, 0.227, 0.25, 0.24, -0.05, -0.42, -0.25, -0.31, 5.75, 0.0}},
            {"a3_a4_zero",
             [](DataT& data)
             {
               data.parameters[Params::A3] = static_cast<RealT>(0.0);
               data.parameters[Params::A4] = static_cast<RealT>(0.0);
             },
             {0.19, 0.28, 0.40, 0.32, 0.25, 0.24, -0.05, -0.42, -0.25, -0.31, 5.75, 0.0}},
            {"time_zero",
             [](DataT& data)
             {
               data.parameters[Params::T2] = static_cast<RealT>(0.0);
               data.parameters[Params::T4] = static_cast<RealT>(0.0);
               data.parameters[Params::T6] = static_cast<RealT>(0.0);
             },
             {0.19, 0.28, 0.37, 0.0878, 0.30, 0.30, 0.30, -0.42, -0.10, -0.10, 9.95, 0.0}},
        };

        const auto loose_tol = static_cast<RealT>(1.0e-4);
        for (const auto& test : cases)
        {
          PhasorDynamics::SignalNode<ScalarT, IdxT> u_node;
          PhasorDynamics::SignalNode<ScalarT, IdxT> vss_node;
          ScalarT                                   u_value{0.5};
          IdxT                                      u_index{12};
          ScalarT                                   vss_value{0.0};
          IdxT                                      vss_index{INVALID_INDEX<IdxT>};

          u_node.set(&u_value, &u_index);
          vss_node.set(&vss_value, &vss_index);

          auto data = makeData();
          test.edit(data);

          PhasorDynamics::Stabilizer::Ieeest<ScalarT, IdxT> model(data);
          model.getSignals().template attachSignalNode<PhasorDynamics::Stabilizer::IeeestExternalVariables::U>(&u_node);
          model.getSignals().template assignSignalNode<PhasorDynamics::Stabilizer::IeeestInternalVariables::VSS>(&vss_node);

          model.allocate();
          model.initialize();

          model.y()[0]  = 0.1;
          model.y()[1]  = 0.2;
          model.y()[2]  = 0.3;
          model.y()[3]  = 0.4;
          model.y()[4]  = 0.5;
          model.y()[5]  = 0.6;
          model.y()[6]  = 0.7;
          model.y()[7]  = 0.8;
          model.y()[8]  = 0.9;
          model.y()[9]  = 1.0;
          model.y()[10] = 0.05;
          model.y()[11] = 0.05;

          model.yp()[0] = 0.01;
          model.yp()[1] = 0.02;
          model.yp()[2] = 0.03;
          model.yp()[3] = 0.04;
          model.yp()[4] = 0.05;
          model.yp()[5] = 0.06;
          model.yp()[6] = 0.07;

          model.evaluateResidual();

          for (size_t i = 0; i < test.expected.size(); ++i)
          {
            auto test_tol = (i == 11) ? loose_tol : tol_;
            if (!isEqual(model.getResidual()[i], test.expected[i], test_tol))
            {
              std::cout << "Incorrect residual for " << test.name
                        << " row " << i << ": "
                        << std::setprecision(15) << model.getResidual()[i]
                        << " != " << test.expected[i] << "\n";
              success = false;
            }
          }
        }

        return success.report(__func__);
      }

      TestOutcome verify()
      {
        TestStatus success = true;
        using Params       = PhasorDynamics::Stabilizer::IeeestParameters;

        {
          PhasorDynamics::Stabilizer::Ieeest<ScalarT, IdxT> model(makeData());
          model.allocate();
          success *= (model.verify() != 0);
        }

        {
          PhasorDynamics::SignalNode<ScalarT, IdxT>         u_node;
          PhasorDynamics::Stabilizer::Ieeest<ScalarT, IdxT> model(makeData());
          model.getSignals().template attachSignalNode<PhasorDynamics::Stabilizer::IeeestExternalVariables::U>(&u_node);
          model.allocate();
          success *= (model.verify() != 0);
        }

        {
          PhasorDynamics::SignalNode<ScalarT, IdxT> u_node;
          ScalarT                                   u_value{0.0};
          IdxT                                      u_index{12};
          u_node.set(&u_value, &u_index);

          auto data                   = makeData();
          data.parameters[Params::A1] = static_cast<RealT>(1.0);
          data.parameters[Params::A2] = static_cast<RealT>(0.0);
          data.parameters[Params::A3] = static_cast<RealT>(0.0);
          data.parameters[Params::A4] = static_cast<RealT>(0.0);

          PhasorDynamics::Stabilizer::Ieeest<ScalarT, IdxT> model(data);
          model.getSignals().template attachSignalNode<PhasorDynamics::Stabilizer::IeeestExternalVariables::U>(&u_node);
          model.allocate();
          success *= (model.verify() != 0);
        }

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      TestOutcome jacobian()
      {
        TestStatus success = true;
        using DepVar       = DependencyTracking::Variable;

        std::vector<DependencyTracking::Variable::DependencyMap> dependency_tracking_jacobian;

        {
          PhasorDynamics::SignalNode<DepVar, IdxT> u_node;
          PhasorDynamics::SignalNode<DepVar, IdxT> vss_node;
          DepVar                                   u_value{0.5};
          IdxT                                     u_index{12};
          DepVar                                   vss_value{0.0};
          IdxT                                     vss_index{INVALID_INDEX<IdxT>};

          u_node.set(&u_value, &u_index);
          vss_node.set(&vss_value, &vss_index);

          PhasorDynamics::Stabilizer::Ieeest<DepVar, IdxT> model(makeData());
          model.getSignals().template attachSignalNode<PhasorDynamics::Stabilizer::IeeestExternalVariables::U>(&u_node);
          model.getSignals().template assignSignalNode<PhasorDynamics::Stabilizer::IeeestInternalVariables::VSS>(&vss_node);

          model.allocate();
          model.initialize();

          for (size_t i = 0; i < model.size(); ++i)
          {
            model.y()[i].setVariableNumber(i);
          }
          u_value.setVariableNumber(model.size());
          u_value.setValue(0.5);

          model.y()[0].setValue(0.1);
          model.y()[1].setValue(0.2);
          model.y()[2].setValue(0.3);
          model.y()[3].setValue(0.4);
          model.y()[4].setValue(0.5);
          model.y()[5].setValue(0.6);
          model.y()[6].setValue(0.7);
          model.y()[7].setValue(0.8);
          model.y()[8].setValue(0.9);
          model.y()[9].setValue(1.0);
          model.y()[10].setValue(0.05);
          model.y()[11].setValue(0.05);

          model.yp()[0].setValue(0.01);
          model.yp()[1].setValue(0.02);
          model.yp()[2].setValue(0.03);
          model.yp()[3].setValue(0.04);
          model.yp()[4].setValue(0.05);
          model.yp()[5].setValue(0.06);
          model.yp()[6].setValue(0.07);

          model.evaluateResidual();
          std::vector<DepVar> residual_y = model.getResidual();

          model.initialize();
          for (size_t i = 0; i < model.size(); ++i)
          {
            model.y()[i] = model.y()[i].getValue();
            model.yp()[i].setVariableNumber(i);
          }
          u_value = 0.5;

          model.y()[0].setValue(0.1);
          model.y()[1].setValue(0.2);
          model.y()[2].setValue(0.3);
          model.y()[3].setValue(0.4);
          model.y()[4].setValue(0.5);
          model.y()[5].setValue(0.6);
          model.y()[6].setValue(0.7);
          model.y()[7].setValue(0.8);
          model.y()[8].setValue(0.9);
          model.y()[9].setValue(1.0);
          model.y()[10].setValue(0.05);
          model.y()[11].setValue(0.05);

          model.yp()[0].setValue(0.01);
          model.yp()[1].setValue(0.02);
          model.yp()[2].setValue(0.03);
          model.yp()[3].setValue(0.04);
          model.yp()[4].setValue(0.05);
          model.yp()[5].setValue(0.06);
          model.yp()[6].setValue(0.07);

          model.evaluateResidual();
          std::vector<DepVar> residual_yp = model.getResidual();

          dependency_tracking_jacobian.resize(residual_y.size());
          for (IdxT i = 0; i < residual_y.size(); ++i)
          {
            auto dependency_y  = residual_y[i].getDependencies();
            auto dependency_yp = residual_yp[i].getDependencies();

            for (const auto& pair_y : dependency_y)
            {
              auto it_yp = dependency_yp.find(pair_y.first);
              if (it_yp != dependency_yp.end())
              {
                dependency_tracking_jacobian[i].insert(std::make_pair(pair_y.first, pair_y.second + it_yp->second));
              }
              else
              {
                dependency_tracking_jacobian[i].insert(std::make_pair(pair_y.first, pair_y.second));
              }
            }

            for (const auto& pair_yp : dependency_yp)
            {
              if (!dependency_y.contains(pair_yp.first))
              {
                dependency_tracking_jacobian[i].insert(std::make_pair(pair_yp.first, pair_yp.second));
              }
            }
          }
        }

        std::vector<DependencyTracking::Variable::DependencyMap> enzyme_jacobian;

        {
          PhasorDynamics::SignalNode<ScalarT, IdxT> u_node;
          PhasorDynamics::SignalNode<ScalarT, IdxT> vss_node;
          ScalarT                                   u_value{0.5};
          IdxT                                      u_index{12};
          ScalarT                                   vss_value{0.0};
          IdxT                                      vss_index{INVALID_INDEX<IdxT>};

          u_node.set(&u_value, &u_index);
          vss_node.set(&vss_value, &vss_index);

          PhasorDynamics::Stabilizer::Ieeest<ScalarT, IdxT> model(makeData());
          model.getSignals().template attachSignalNode<PhasorDynamics::Stabilizer::IeeestExternalVariables::U>(&u_node);
          model.getSignals().template assignSignalNode<PhasorDynamics::Stabilizer::IeeestInternalVariables::VSS>(&vss_node);

          model.allocate();
          model.initialize();

          model.y()[0]  = 0.1;
          model.y()[1]  = 0.2;
          model.y()[2]  = 0.3;
          model.y()[3]  = 0.4;
          model.y()[4]  = 0.5;
          model.y()[5]  = 0.6;
          model.y()[6]  = 0.7;
          model.y()[7]  = 0.8;
          model.y()[8]  = 0.9;
          model.y()[9]  = 1.0;
          model.y()[10] = 0.05;
          model.y()[11] = 0.05;

          model.yp()[0] = 0.01;
          model.yp()[1] = 0.02;
          model.yp()[2] = 0.03;
          model.yp()[3] = 0.04;
          model.yp()[4] = 0.05;
          model.yp()[5] = 0.06;
          model.yp()[6] = 0.07;

          model.updateTime(0.0, 1.0);
          model.evaluateResidual();
          model.evaluateJacobian();

          auto model_jacobian = model.getJacobian();
          model_jacobian.deduplicate();
          enzyme_jacobian = GridKit::Testing::MapFromCOO(model_jacobian);
        }

        for (size_t i = 0; i < dependency_tracking_jacobian.size(); ++i)
        {
          success *= GridKit::Testing::isEqual(dependency_tracking_jacobian[i], enzyme_jacobian[i], tol_);
        }

        return success.report(__func__);
      }
#endif

    private:
      static constexpr ScalarT tol_ = 10 * std::numeric_limits<ScalarT>::epsilon();

      auto makeData() -> DataT
      {
        using Params = PhasorDynamics::Stabilizer::IeeestParameters;

        DataT data;
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
    }; // class StabilizerIeeestTests

  } // namespace Testing
} // namespace GridKit
