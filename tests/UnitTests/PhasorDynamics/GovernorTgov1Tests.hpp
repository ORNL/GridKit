#include <iomanip>
#include <iostream>

#include <AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <Definitions.hpp>
#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/Bus/BusData.hpp>
#include <Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <Model/PhasorDynamics/Governor/Tgov1/Tgov1.hpp>
#include <Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <Model/PhasorDynamics/SynchronousMachine/GENROUwS/Genrou.hpp>
#include <Model/PhasorDynamics/SynchronousMachine/GENROUwS/GenrouData.hpp>
#include <Utilities/MapFromCOO.hpp>
#include <Utilities/TestHelpers.hpp>
#include <Utilities/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {

    template <class ScalarT, typename IdxT>
    class GovernorTgov1Tests
    {
    private:
      using real_type = typename PhasorDynamics::Component<ScalarT, IdxT>::real_type;

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

        using BusType          = PhasorDynamics::BusData<ScalarT, IdxT>::BusType;
        using GenrouParameters = PhasorDynamics::GenrouData<ScalarT, IdxT>::Parameters;
        using GenrouPorts      = PhasorDynamics::GenrouData<ScalarT, IdxT>::Ports;

        PhasorDynamics::BusData<ScalarT, IdxT> busdata;
        busdata.bus_id   = 0;
        busdata.bus_type = BusType::DEFAULT;
        busdata.Vr0      = 1.0;
        busdata.Vi0      = 0.0;

        PhasorDynamics::GenrouData<ScalarT, IdxT> gendata;
        gendata.ports[GenrouPorts::bus] = 0;

        gendata.parameters[GenrouParameters::p0]    = 1.;
        gendata.parameters[GenrouParameters::q0]    = 0.05013;
        gendata.parameters[GenrouParameters::H]     = 3.;
        gendata.parameters[GenrouParameters::D]     = 0.;
        gendata.parameters[GenrouParameters::Ra]    = 0.;
        gendata.parameters[GenrouParameters::Tdop]  = 7.;
        gendata.parameters[GenrouParameters::Tdopp] = .04;
        gendata.parameters[GenrouParameters::Tqopp] = .05;
        gendata.parameters[GenrouParameters::Tqop]  = .75;
        gendata.parameters[GenrouParameters::Xd]    = 2.1;
        gendata.parameters[GenrouParameters::Xdp]   = 0.2;
        gendata.parameters[GenrouParameters::Xdpp]  = 0.18;
        gendata.parameters[GenrouParameters::Xq]    = 0.5;
        gendata.parameters[GenrouParameters::Xqp]   = 0.5;
        gendata.parameters[GenrouParameters::Xqpp]  = 0.18;
        gendata.parameters[GenrouParameters::Xl]    = 0.15;
        gendata.parameters[GenrouParameters::S10]   = 0.;
        gendata.parameters[GenrouParameters::S12]   = 0.;

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
        auto tol = 10 * std::numeric_limits<real_type>::epsilon();

        const std::vector<ScalarT>& f = gov.getResidual();
        for (const auto& f_val : f)
        {
          if (!isEqual(f_val, 0.0, tol))
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

        using BusType          = PhasorDynamics::BusData<ScalarT, IdxT>::BusType;
        using GenrouParameters = PhasorDynamics::GenrouData<ScalarT, IdxT>::Parameters;
        using GenrouPorts      = PhasorDynamics::GenrouData<ScalarT, IdxT>::Ports;

        PhasorDynamics::BusData<ScalarT, IdxT> busdata;
        busdata.bus_id   = 0;
        busdata.bus_type = BusType::DEFAULT;
        busdata.Vr0      = 1.0;
        busdata.Vi0      = 0.0;

        PhasorDynamics::GenrouData<ScalarT, IdxT> gendata;
        gendata.ports[GenrouPorts::bus] = 0;

        gendata.parameters[GenrouParameters::p0]    = 1.;
        gendata.parameters[GenrouParameters::q0]    = 0.05013;
        gendata.parameters[GenrouParameters::H]     = 3.;
        gendata.parameters[GenrouParameters::D]     = 0.;
        gendata.parameters[GenrouParameters::Ra]    = 0.;
        gendata.parameters[GenrouParameters::Tdop]  = 7.;
        gendata.parameters[GenrouParameters::Tdopp] = .04;
        gendata.parameters[GenrouParameters::Tqopp] = .05;
        gendata.parameters[GenrouParameters::Tqop]  = .75;
        gendata.parameters[GenrouParameters::Xd]    = 2.1;
        gendata.parameters[GenrouParameters::Xdp]   = 0.2;
        gendata.parameters[GenrouParameters::Xdpp]  = 0.18;
        gendata.parameters[GenrouParameters::Xq]    = 0.5;
        gendata.parameters[GenrouParameters::Xqp]   = 0.5;
        gendata.parameters[GenrouParameters::Xqpp]  = 0.18;
        gendata.parameters[GenrouParameters::Xl]    = 0.15;
        gendata.parameters[GenrouParameters::S10]   = 0.;
        gendata.parameters[GenrouParameters::S12]   = 0.;

        /// Jacobian via DependencyTracking
        std::vector<DependencyTracking::Variable> dependency_tracking_residuals = DependencyTrackingJacobian(busdata, gendata);

        /// Jacobian via Enzyme
        std::vector<DependencyTracking::Variable::DependencyMap> enzyme_jacobian = EnzymeJacobian(busdata, gendata);

        /// Compare DependencyTracking dependencies to Enzyme's
        auto tol = 1000 * std::numeric_limits<real_type>::epsilon(); ///< Todo: Check why the results don't match to machine precision
        for (size_t i = 0; i < dependency_tracking_residuals.size(); ++i)
        {
          DependencyTracking::Variable                       res           = dependency_tracking_residuals[i];
          const DependencyTracking::Variable::DependencyMap& dependencies  = res.getDependencies();
          success                                                         *= (GridKit::Testing::isEqual(dependencies, enzyme_jacobian[i], tol));
        }

        return success.report(__func__);
      }

    private:
      std::vector<DependencyTracking::Variable> DependencyTrackingJacobian(
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

        bus.initialize();
        gen.initialize();
        gov.initialize();

        for (size_t i = 0; i < gov.size(); ++i)
        {
          gov.y()[i].setVariableNumber(i); ///< Independent variables
        }

        bus.evaluateResidual();
        gen.evaluateResidual();
        gov.evaluateResidual(); ///< Computes the residual and the Jacobian values by tracking
                                ///< the dependencies
        std::vector<DependencyTracking::Variable> residual = gov.getResidual();

        /// Print the dependencies
        for (size_t i = 0; i < residual.size(); ++i)
        {
          std::cout << i << "th residual: ";
          (residual[i]).print(std::cout);
          std::cout << "\n";
        }

        return residual;
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

        bus.initialize();
        gen.initialize();
        gov.initialize();

        bus.evaluateResidual();
        gen.evaluateResidual();
        gov.evaluateResidual();

        gov.evaluateJacobian();
        GridKit::LinearAlgebra::COO_Matrix<ScalarT, IdxT> model_jacobian = gov.getJacobian();
        model_jacobian.printMatrix("Model Jacobian");

        return GridKit::Testing::MapFromCOO(model_jacobian);
      }
#endif
    }; // class GovernorTgov1Tests

  } // namespace Testing
} // namespace GridKit
