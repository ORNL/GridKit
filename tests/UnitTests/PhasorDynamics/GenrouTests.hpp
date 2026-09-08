#include <iomanip>
#include <iostream>
#include <sstream>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENROU/Genrou.hpp>
#include <GridKit/Model/VariableMonitorController.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Testing/Tokenizer.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {

    template <class ScalarT, typename IdxT>
    class GenrouTests
    {
    private:
      using RealT                   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;
      using GenrouDataT             = PhasorDynamics::GenrouData<RealT, IdxT>;
      static constexpr ScalarT tol_ = 10 * std::numeric_limits<ScalarT>::epsilon();

      static GenrouDataT makeGenrouData()
      {
        using Parameter = typename GenrouDataT::Parameters;
        using Buses     = typename GenrouDataT::Buses;

        GenrouDataT data;
        data.device_class                 = "Genrou";
        data.disambiguation_string        = "1";
        data.buses[Buses::bus]            = 1;
        data.parameters[Parameter::p0]    = RealT{1.0};
        data.parameters[Parameter::q0]    = RealT{0.05013};
        data.parameters[Parameter::H]     = RealT{3.0};
        data.parameters[Parameter::D]     = RealT{0.0};
        data.parameters[Parameter::Ra]    = RealT{0.0};
        data.parameters[Parameter::Tdop]  = RealT{7.0};
        data.parameters[Parameter::Tdopp] = RealT{0.04};
        data.parameters[Parameter::Tqop]  = RealT{0.75};
        data.parameters[Parameter::Tqopp] = RealT{0.05};
        data.parameters[Parameter::Xd]    = RealT{2.1};
        data.parameters[Parameter::Xdp]   = RealT{0.2};
        data.parameters[Parameter::Xdpp]  = RealT{0.18};
        data.parameters[Parameter::Xq]    = RealT{0.5};
        data.parameters[Parameter::Xqp]   = RealT{0.5};
        data.parameters[Parameter::Xqpp]  = RealT{0.18};
        data.parameters[Parameter::Xl]    = RealT{0.15};
        data.parameters[Parameter::S10]   = RealT{0.0};
        data.parameters[Parameter::S12]   = RealT{0.0};

        return data;
      }

    public:
      GenrouTests()  = default;
      ~GenrouTests() = default;

      TestOutcome constructor()
      {
        TestStatus success = true;

        auto* bus = new PhasorDynamics::Bus<ScalarT, IdxT>(1.0, 0.0);

        PhasorDynamics::Component<ScalarT, IdxT>* machine =
            new PhasorDynamics::Genrou<ScalarT, IdxT>(bus);

        success *= (machine != nullptr);

        if (machine)
        {
          delete machine;
        }
        delete bus;

        return success.report(__func__);
      }

      /**
       * @brief Checks residual evaluation.
       *
       * The test instantiates and initializes Genrou model. Properly
       * initialized model should have residual equal to zero within machine
       * precision.
       *
       * @return TestOutcome - wheter test was successful
       */
      TestOutcome residual()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT>    bus(1.0, 0.0);
        PhasorDynamics::Genrou<ScalarT, IdxT> gen(&bus,
                                                  1,
                                                  0.05013,
                                                  3,
                                                  0,
                                                  0,
                                                  7,
                                                  0.04,
                                                  0.05,
                                                  0.75,
                                                  2.1,
                                                  0.2,
                                                  0.18,
                                                  0.5,
                                                  0.5,
                                                  0.18,
                                                  0.15,
                                                  0,
                                                  0);

        bus.allocate();
        bus.initialize();
        bus.evaluateResidual();

        gen.allocate();
        gen.initialize();
        gen.evaluateResidual();

        // Require results to be within machine precision
        auto tol = 10 * std::numeric_limits<RealT>::epsilon();

        const auto& f      = gen.getResidual();
        const auto* f_data = f.getData();
        for (std::size_t i = 0; i < f.getSize(); ++i)
        {
          if (!isEqual(f_data[i], 0.0, tol))
            success = false;
        }

        return success.report(__func__);
      }

      /**
       * @brief Checks monitored terminal current and power use system base.
       */
      TestOutcome monitor_system_base()
      {
        TestStatus success = true;

        using Parameter = typename GenrouDataT::Parameters;
        using Variable  = typename GenrouDataT::MonitorableVariables;

        auto data                       = makeGenrouData();
        data.parameters[Parameter::mva] = RealT{50.0};
        data.monitored_variables.insert(Variable::ir);
        data.monitored_variables.insert(Variable::p);

        PhasorDynamics::Bus<ScalarT, IdxT>    bus(1.0, 0.0);
        PhasorDynamics::Genrou<ScalarT, IdxT> gen(&bus, data);

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
        monitor.addSink({Model::VariableMonitorFormat::CSV}, os);
        monitor.print();

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

      // A test to verify that the hard coded answers match those given by the residual functions
      // Hard coded parameters, differential, and algebraic terms
      TestOutcome hard_coded_residual()
      {
        TestStatus success = true;

        // GenRou generator parameters
        RealT p0{1};
        RealT q0{.05013};
        RealT H{3};
        RealT D{.5};
        RealT Ra{.1};
        RealT Tdop{7};
        RealT Tdopp{.04};
        RealT Tqopp{.05};
        RealT Tqop{.75};
        RealT Xd{2.1};
        RealT Xdp{.2};
        RealT Xdpp{.5};
        RealT Xq{.18};
        RealT Xqp{.3};
        RealT Xqpp{.5};
        RealT Xl{.5};
        RealT S10{.1};
        RealT S12{.2};

        ScalarT Vr1{1.0}; ///< Bus real voltage
        ScalarT Vi1{0};   ///< Bus imaginary voltage

        PhasorDynamics::Bus<ScalarT, IdxT>    bus(Vr1, Vi1);
        PhasorDynamics::Genrou<ScalarT, IdxT> gen(&bus, p0, q0, H, D, Ra, Tdop, Tdopp, Tqopp, Tqop, Xd, Xdp, Xdpp, Xq, Xqp, Xqpp, Xl, S10, S12);

        // Answer key is available only in double precision.
        // Therefore, only double precision tests are done at this time.
        const std::vector<ScalarT> res_answer = {
            -2.0 * std::numbers::pi_v<RealT> * 60.0,
            -static_cast<ScalarT>(10.) / static_cast<ScalarT>(9.),
            -static_cast<ScalarT>(223.) / static_cast<ScalarT>(525.),
            -54.75,
            -9.6,
            static_cast<ScalarT>(892.) / static_cast<ScalarT>(375.),
            0.21,
            -0.07,
            -0.19223748416156686,
            2.0,
            1.4,
            0.31,
            2.211,
            0.85,
            1.2,
            -static_cast<ScalarT>(13.) / static_cast<ScalarT>(130.),
            -static_cast<ScalarT>(143.) / static_cast<ScalarT>(52.)};

        bus.allocate();
        bus.initialize();

        // Allocate but not initialize generator model
        gen.allocate();
        // TODO: Set pmech and efd. They are currently not set in this test, as we are not
        // calling gen.initialize(). The private members are initialized to 0 as a workaround,
        // but this needs to be better handled in the model implementation.

        // Set variable values matching the answer key
        auto* y  = gen.y().getData();
        auto* yp = gen.yp().getData();

        static constexpr auto pi = std::numbers::pi_v<RealT>;

        y[0]  = pi;  // delta
        y[1]  = 2.0; // omega
        y[2]  = 2.0; // Eqp
        y[3]  = .1;  // psidp
        y[4]  = .01; // psiqp
        y[5]  = .6;  // Edp
        y[6]  = .2;  // psiqp
        y[7]  = .03; // psidpp
        y[8]  = .01; // psipp
        y[9]  = 2.0; // ksat
        y[10] = .8;  // vd
        y[11] = .4;  // vq
        y[12] = 2.0; // telec
        y[13] = 1.1; // id
        y[14] = .3;  // iq
        y[15] = .9;  // ir
        y[16] = .25; // ii

        // Set derivative values matching the answer key
        yp[0] = 2.0 * pi * 60.0; // delta_dot
        yp[1] = -1.5;            // omega_dot
        yp[2] = 1.0;             // Eqp_dot
        yp[3] = 1.0;             // psidp_dot
        yp[4] = 1.0;             // psiqp_dot
        yp[5] = 1.0;             // Edp_dot

        gen.y().setDataUpdated();
        gen.yp().setDataUpdated();
        gen.evaluateResidual();
        auto&       residual      = gen.getResidual();
        const auto* residual_data = residual.getData();

        for (size_t i = 0; i < res_answer.size(); ++i)
        {
          if (!isEqual(residual_data[i], res_answer[i], tol_))
          {
            std::cout << "Incorrect result for residual " << i << ": "
                      << residual_data[i] << " != " << res_answer[i] << "\n";
            success = false;
            break;
          }
        }

        return success.report(__func__);
      }

      TestOutcome accessors()
      {
        TestStatus success = true;
        success.skipTest();

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

        // Jacobian via DependencyTracking
        std::vector<DependencyTracking::Variable::DependencyMap> dependency_tracking_jacobian = DependencyTrackingJacobian();

        // Jacobian via Enzyme
        std::vector<DependencyTracking::Variable::DependencyMap> enzyme_jacobian = EnzymeJacobian();

        /// Compare DependencyTracking dependencies to Enzyme's
        for (size_t i = 0; i < dependency_tracking_jacobian.size(); ++i)
        {
          success *= (GridKit::Testing::isEqual(dependency_tracking_jacobian[i], enzyme_jacobian[i], tol));
        }

        return success.report(__func__);
      }

    private:
      std::vector<DependencyTracking::Variable::DependencyMap> DependencyTrackingJacobian()
      {
        DependencyTracking::Variable                               Vr1{1.0}; ///< Bus real voltage
        DependencyTracking::Variable                               Vi1{0.0}; ///< Bus imaginary voltage
        PhasorDynamics::Bus<DependencyTracking::Variable, IdxT>    bus(Vr1, Vi1);
        PhasorDynamics::Genrou<DependencyTracking::Variable, IdxT> gen(&bus,
                                                                       1,
                                                                       0.05013,
                                                                       3,
                                                                       0,
                                                                       0,
                                                                       7,
                                                                       0.04,
                                                                       0.05,
                                                                       0.75,
                                                                       2.1,
                                                                       0.2,
                                                                       0.18,
                                                                       0.5,
                                                                       0.5,
                                                                       0.18,
                                                                       0.15,
                                                                       0,
                                                                       0);

        bus.allocate();
        gen.allocate();

        for (size_t i = 0; i < bus.size(); ++i)
        {
          bus.setVariableIndex(i, i + gen.size()); // Reset bus variable indices
          bus.setResidualIndex(i, i + gen.size()); // Reset bus residual indices
        }

        bus.initialize();
        gen.initialize();

        gen.updateTime(0.0, 1.0); // Set alpha to 1.0 to verify d/dy' term

        bus.evaluateResidual();
        gen.evaluateResidual();

        bus.evaluateJacobian();
        gen.evaluateJacobian();
        auto* model_jacobian = gen.getCsrJacobian();
        std::cout << "Sparse Csr Matrix: Genrou DependencyTracking Jacobian\n";
        model_jacobian->print();

        return GridKit::Testing::MapFromCsr(model_jacobian);      
      }

      std::vector<DependencyTracking::Variable::DependencyMap> EnzymeJacobian()
      {
        ScalarT                               Vr1{1.0}; ///< Bus real voltage
        ScalarT                               Vi1{0.0}; ///< Bus imaginary voltage
        PhasorDynamics::Bus<ScalarT, IdxT>    bus(Vr1, Vi1);
        PhasorDynamics::Genrou<ScalarT, IdxT> gen(&bus,
                                                  1,
                                                  0.05013,
                                                  3,
                                                  0,
                                                  0,
                                                  7,
                                                  0.04,
                                                  0.05,
                                                  0.75,
                                                  2.1,
                                                  0.2,
                                                  0.18,
                                                  0.5,
                                                  0.5,
                                                  0.18,
                                                  0.15,
                                                  0,
                                                  0);

        bus.allocate();
        gen.allocate();

        for (size_t i = 0; i < bus.size(); ++i)
        {
          bus.setVariableIndex(i, i + gen.size()); // Reset bus variable indices
          bus.setResidualIndex(i, i + gen.size()); // Reset bus residual indices
        }

        bus.initialize();
        gen.initialize();

        gen.updateTime(0.0, 1.0); // Set alpha to 1.0 to verify d/dy' term

        bus.evaluateResidual();
        gen.evaluateResidual();

        bus.evaluateJacobian();
        gen.evaluateJacobian();
        gen.constructCsr();
        auto* model_jacobian = gen.getCsrJacobian();
        std::cout << "Sparse Csr Matrix: Genrou Enzyme Jacobian\n";
        model_jacobian->print();

        return GridKit::Testing::MapFromCsr(model_jacobian);
      }
#endif
    }; // class GenrouTest

  } // namespace Testing
} // namespace GridKit
