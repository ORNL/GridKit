#include <iomanip>
#include <iostream>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusData.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <GridKit/Model/PhasorDynamics/Governor/Tgov1/Tgov1.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENROU/Genrou.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENROU/GenrouData.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {

    template <class ScalarT, typename IdxT>
    class GovernorTgov1Tests
    {
    private:
      using RealT                   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;
      using real_type               = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;
      static constexpr ScalarT tol_ = 10 * std::numeric_limits<ScalarT>::epsilon();

    public:
      GovernorTgov1Tests()  = default;
      ~GovernorTgov1Tests() = default;

      TestOutcome constructor()
      {
        TestStatus success = true;

        PhasorDynamics::Bus<ScalarT, IdxT>        bus(1.0, 0.0);
        PhasorDynamics::SignalNode<ScalarT, IdxT> pmech;
        PhasorDynamics::SignalNode<ScalarT, IdxT> omega;

        PhasorDynamics::Governor::Tgov1<ScalarT, IdxT>* gov =
            new PhasorDynamics::Governor::Tgov1<ScalarT, IdxT>(&pmech, &omega);

        success *= (gov != nullptr);

        if (gov)
        {
          delete gov;
        }

        return success.report(__func__);
      }

      /**
       * For Tgov:
       * This section is the updated residual, not residual initialization.
       */

      TestOutcome residual()
      {
        TestStatus success = true;

        using BusType     = PhasorDynamics::BusData<ScalarT, IdxT>::BusType;
        using GenrouDataT = PhasorDynamics::GenrouData<ScalarT, IdxT>;
        using Parameter   = typename GenrouDataT::Parameters;
        using Buses       = typename GenrouDataT::Buses;

        PhasorDynamics::BusData<ScalarT, IdxT> busdata;
        busdata.bus_id   = 0;
        busdata.bus_type = BusType::DEFAULT;
        busdata.Vr0      = 1.0;
        busdata.Vi0      = 0.0;

        PhasorDynamics::GenrouData<ScalarT, IdxT> gendata;
        gendata.buses[Buses::bus] = 0;

        gendata.parameters[Parameter::p0]    = 1.;
        gendata.parameters[Parameter::q0]    = 0.05013;
        gendata.parameters[Parameter::H]     = 3.;
        gendata.parameters[Parameter::D]     = 0.;
        gendata.parameters[Parameter::Ra]    = 0.;
        gendata.parameters[Parameter::Tdop]  = 7.;
        gendata.parameters[Parameter::Tdopp] = .04;
        gendata.parameters[Parameter::Tqopp] = .05;
        gendata.parameters[Parameter::Tqop]  = .75;
        gendata.parameters[Parameter::Xd]    = 2.1;
        gendata.parameters[Parameter::Xdp]   = 0.2;
        gendata.parameters[Parameter::Xdpp]  = 0.18;
        gendata.parameters[Parameter::Xq]    = 0.5;
        gendata.parameters[Parameter::Xqp]   = 0.5;
        gendata.parameters[Parameter::Xqpp]  = 0.18;
        gendata.parameters[Parameter::Xl]    = 0.15;
        gendata.parameters[Parameter::S10]   = 0.;
        gendata.parameters[Parameter::S12]   = 0.;

        PhasorDynamics::Bus<ScalarT, IdxT>             bus(busdata);
        PhasorDynamics::SignalNode<ScalarT, IdxT>      pmech;
        PhasorDynamics::SignalNode<ScalarT, IdxT>      omega;
        PhasorDynamics::Genrou<ScalarT, IdxT>          gen(&bus, &omega, &pmech, gendata);
        PhasorDynamics::Governor::Tgov1<ScalarT, IdxT> gov(&pmech, &omega);

        // Test answer keys
        const std::vector<ScalarT> res_answer = {static_cast<ScalarT>(2.0) / static_cast<ScalarT>(15.0),
                                                 -2.0,
                                                 static_cast<ScalarT>(1.0) / static_cast<ScalarT>(3.0)};

        bus.allocate();
        gen.allocate();
        gov.allocate();

        bus.initialize();
        gen.initialize();
        gov.initialize();

        bus.evaluateResidual();
        gen.evaluateResidual();
        gov.evaluateResidual();

        // Set variable values matching the answer key
        auto* y  = gov.y().getData();
        auto* yp = gov.yp().getData();

        y[0] = 1.0;                                                     // Ptx
        y[1] = 1.0;                                                     // Pv
        y[2] = static_cast<ScalarT>(10.0) / static_cast<ScalarT>(15.0); // Pmech

        // Set derivative values matching the answer key
        yp[0] = static_cast<ScalarT>(8.0) / static_cast<ScalarT>(15.0); // ptx_dot
        yp[1] = 2.0;                                                    // pv_dot

        gov.y().setDataUpdated();
        gov.yp().setDataUpdated();
        gov.evaluateResidual();
        auto&       residual      = gov.getResidual();
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

      /**
       * @brief Checks residual evaluation.
       *
       * The test instantiates and initializes Genrou model. Properly
       * initialized model should have residual equal to zero within machine
       * precision.
       *
       * @return TestOutcome - wheter test was successful
       *
       * (Verifies the residual evaluates to zero for the initial conditions)
       */
      TestOutcome zeroInitialResidual()
      {
        TestStatus success = true;

        using BusType     = typename PhasorDynamics::BusData<ScalarT, IdxT>::BusType;
        using GenrouDataT = PhasorDynamics::GenrouData<ScalarT, IdxT>;
        using Parameter   = typename GenrouDataT::Parameters;
        using Buses       = typename GenrouDataT::Buses;

        PhasorDynamics::BusData<ScalarT, IdxT> busdata;
        busdata.bus_id   = 0;
        busdata.bus_type = BusType::DEFAULT;
        busdata.Vr0      = 1.0;
        busdata.Vi0      = 0.0;

        PhasorDynamics::GenrouData<ScalarT, IdxT> gendata;
        gendata.buses[Buses::bus] = 0;

        gendata.parameters[Parameter::p0]    = 1.;
        gendata.parameters[Parameter::q0]    = 0.05013;
        gendata.parameters[Parameter::H]     = 3.;
        gendata.parameters[Parameter::D]     = 0.;
        gendata.parameters[Parameter::Ra]    = 0.;
        gendata.parameters[Parameter::Tdop]  = 7.;
        gendata.parameters[Parameter::Tdopp] = .04;
        gendata.parameters[Parameter::Tqopp] = .05;
        gendata.parameters[Parameter::Tqop]  = .75;
        gendata.parameters[Parameter::Xd]    = 2.1;
        gendata.parameters[Parameter::Xdp]   = 0.2;
        gendata.parameters[Parameter::Xdpp]  = 0.18;
        gendata.parameters[Parameter::Xq]    = 0.5;
        gendata.parameters[Parameter::Xqp]   = 0.5;
        gendata.parameters[Parameter::Xqpp]  = 0.18;
        gendata.parameters[Parameter::Xl]    = 0.15;
        gendata.parameters[Parameter::S10]   = 0.;
        gendata.parameters[Parameter::S12]   = 0.;

        PhasorDynamics::Bus<ScalarT, IdxT>             bus(busdata);
        PhasorDynamics::SignalNode<ScalarT, IdxT>      pmech;
        PhasorDynamics::SignalNode<ScalarT, IdxT>      omega;
        PhasorDynamics::Genrou<ScalarT, IdxT>          gen(&bus, &omega, &pmech, gendata);
        // Create governor to be tested
        PhasorDynamics::Governor::Tgov1<ScalarT, IdxT> gov(&pmech, &omega);

        bus.allocate();
        gov.allocate();
        gen.allocate();

        bus.initialize();
        gen.initialize();
        gov.initialize();

        bus.evaluateResidual();
        gen.evaluateResidual();
        gov.evaluateResidual();

        // Require results to be within machine precision
        auto tol = 10 * std::numeric_limits<RealT>::epsilon();

        const auto& f      = gov.getResidual();
        const auto* f_data = f.getData();
        for (std::size_t i = 0; i < f.getSize(); ++i)
        {
          if (!isEqual(f_data[i], 0.0, tol))
            success = false;
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
       * A test case to verify Jacobian values
       */
      TestOutcome jacobian()
      {
        TestStatus success = true;

        using BusType     = PhasorDynamics::BusData<ScalarT, IdxT>::BusType;
        using GenrouDataT = PhasorDynamics::GenrouData<ScalarT, IdxT>;
        using Parameter   = typename GenrouDataT::Parameters;
        using Buses       = typename GenrouDataT::Buses;

        PhasorDynamics::BusData<ScalarT, IdxT> busdata;
        busdata.bus_id   = 0;
        busdata.bus_type = BusType::DEFAULT;
        busdata.Vr0      = 1.0;
        busdata.Vi0      = 0.0;

        PhasorDynamics::GenrouData<ScalarT, IdxT> gendata;
        gendata.buses[Buses::bus] = 0;

        gendata.parameters[Parameter::p0]    = 1.;
        gendata.parameters[Parameter::q0]    = 0.05013;
        gendata.parameters[Parameter::H]     = 3.;
        gendata.parameters[Parameter::D]     = 0.;
        gendata.parameters[Parameter::Ra]    = 0.;
        gendata.parameters[Parameter::Tdop]  = 7.;
        gendata.parameters[Parameter::Tdopp] = .04;
        gendata.parameters[Parameter::Tqopp] = .05;
        gendata.parameters[Parameter::Tqop]  = .75;
        gendata.parameters[Parameter::Xd]    = 2.1;
        gendata.parameters[Parameter::Xdp]   = 0.2;
        gendata.parameters[Parameter::Xdpp]  = 0.18;
        gendata.parameters[Parameter::Xq]    = 0.5;
        gendata.parameters[Parameter::Xqp]   = 0.5;
        gendata.parameters[Parameter::Xqpp]  = 0.18;
        gendata.parameters[Parameter::Xl]    = 0.15;
        gendata.parameters[Parameter::S10]   = 0.;
        gendata.parameters[Parameter::S12]   = 0.;

        /// Jacobian via DependencyTracking
        std::vector<DependencyTracking::Variable::DependencyMap> dependency_tracking_jacobian = DependencyTrackingJacobian(busdata, gendata);

        /// Jacobian via Enzyme
        std::vector<DependencyTracking::Variable::DependencyMap> enzyme_jacobian = EnzymeJacobian(busdata, gendata);

        /// Compare DependencyTracking dependencies to Enzyme's
        auto tol = 10 * std::numeric_limits<RealT>::epsilon();
        for (size_t i = 0; i < dependency_tracking_jacobian.size(); ++i)
        {
          success *= (GridKit::Testing::isEqual(dependency_tracking_jacobian[i], enzyme_jacobian[i], tol));
        }
        return success.report(__func__);
      }

    private:
      std::vector<DependencyTracking::Variable::DependencyMap> DependencyTrackingJacobian(
          PhasorDynamics::BusData<ScalarT, IdxT>    busdata,
          PhasorDynamics::GenrouData<ScalarT, IdxT> gendata)
      {
        PhasorDynamics::Bus<DependencyTracking::Variable, IdxT>             bus(busdata);
        PhasorDynamics::SignalNode<DependencyTracking::Variable, IdxT>      pmech;
        PhasorDynamics::SignalNode<DependencyTracking::Variable, IdxT>      omega;
        PhasorDynamics::Genrou<DependencyTracking::Variable, IdxT>          gen(&bus, &omega, &pmech, gendata);
        // Create governor to be tested
        PhasorDynamics::Governor::Tgov1<DependencyTracking::Variable, IdxT> gov(&pmech, &omega);

        bus.allocate();
        gov.allocate();
        gen.allocate();

        gen.setVariableIndex(1, gov.size()); // Reset omega index

        bus.initialize();
        gen.initialize();
        gov.initialize();

        gov.updateTime(0.0, 1.0); // Set alpha to 1.0 to verify d/dy' term

        bus.evaluateResidual();
        gen.evaluateResidual();
        gov.evaluateResidual();

        gov.evaluateJacobian();
        auto* model_jacobian = gov.getCsrJacobian();
        std::cout << "Sparse Csr Matrix: Tgov1 DependencyTracking Jacobian\n";
        model_jacobian->print();

        return GridKit::Testing::MapFromCsr(model_jacobian);      
      }

      std::vector<DependencyTracking::Variable::DependencyMap> EnzymeJacobian(
          PhasorDynamics::BusData<ScalarT, IdxT>    busdata,
          PhasorDynamics::GenrouData<ScalarT, IdxT> gendata)
      {
        PhasorDynamics::Bus<ScalarT, IdxT>             bus(busdata);
        PhasorDynamics::SignalNode<ScalarT, IdxT>      pmech;
        PhasorDynamics::SignalNode<ScalarT, IdxT>      omega;
        PhasorDynamics::Genrou<ScalarT, IdxT>          gen(&bus, &omega, &pmech, gendata);
        // Create governor to be tested
        PhasorDynamics::Governor::Tgov1<ScalarT, IdxT> gov(&pmech, &omega);

        bus.allocate();
        gov.allocate();
        gen.allocate();

        gen.setVariableIndex(1, gov.size()); // Reset omega index

        bus.initialize();
        gen.initialize();
        gov.initialize();

        gov.updateTime(0.0, 1.0); // Set alpha to 1.0 to verify d/dy' term

        bus.evaluateResidual();
        gen.evaluateResidual();
        gov.evaluateResidual();

        gov.evaluateJacobian();
        gov.constructCsr();
        auto* model_jacobian = gov.getCsrJacobian();
        std::cout << "Sparse Csr Matrix: Tgov1 Enzyme Jacobian\n";
        model_jacobian->print();

        return GridKit::Testing::MapFromCsr(model_jacobian);
      }
#endif
    }; // class GovernorTgov1Tests

  } // namespace Testing
} // namespace GridKit
