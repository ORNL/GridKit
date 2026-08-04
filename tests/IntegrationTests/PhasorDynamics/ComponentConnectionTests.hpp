#pragma once

#include <limits>

#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/ESDC1A/Esdc1a.hpp>
#include <GridKit/Model/PhasorDynamics/Exciter/ESDC1A/Esdc1aData.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SynchronousMachine/GENROU/Genrou.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    /// Verify that a pair of components wired through a shared signal node
    /// resolves that node and agrees on its value across system assembly.
    /// These cases exercise component-to-component connections only; the
    /// solver-driven cases live in @ref PDIntegrationTests.
    template <typename scalar_type, typename index_type>
    class ComponentConnectionTests
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      ComponentConnectionTests()  = default;
      ~ComponentConnectionTests() = default;

      // Initialization and evaluateResidual() reassociate the same
      // expressions, so a steady state lands within a few ulps of exact zero
      // rather than on it. The cases here are exact to within one ulp; 100
      // eps is the shared margin used across the phasor-dynamics suites.
      static constexpr RealT kTol =
          static_cast<RealT>(100.0) * std::numeric_limits<RealT>::epsilon();

      /// GENROU seeds the shared field-voltage node before ESDC1A consumes
      /// and preserves that value during system initialization.
      TestOutcome genrouEsdc1a()
      {
        using MachineExternal = PhasorDynamics::GenrouExternalVariables;
        using ExciterInternal = PhasorDynamics::Exciter::Esdc1aInternalVariables;
        using ExciterParams   = PhasorDynamics::Exciter::Esdc1aParameters;

        TestStatus success = true;

        PhasorDynamics::SystemModel<ScalarT, IdxT> system;
        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus(
            static_cast<ScalarT>(1.0),
            static_cast<ScalarT>(0.0));
        PhasorDynamics::SignalNode<ScalarT, IdxT> efd;
        PhasorDynamics::Genrou<ScalarT, IdxT>     machine(&bus);

        PhasorDynamics::Exciter::Esdc1aData<RealT, IdxT> exciter_data;
        exciter_data.parameters[ExciterParams::Tr] = static_cast<RealT>(0.02);
        exciter_data.parameters[ExciterParams::Tb] = static_cast<RealT>(0.5);

        PhasorDynamics::Exciter::Esdc1a<ScalarT, IdxT> exciter(&bus, exciter_data);

        machine.getSignals().template attachSignalNode<MachineExternal::EFD>(&efd);
        exciter.getSignals().template assignSignalNode<ExciterInternal::EFD>(&efd);

        system.addBus(&bus);
        system.addComponent(&machine);
        system.addComponent(&exciter);

        success *= system.allocate() == 0;
        success *= efd.linked();
        success *= system.initialize() == 0;
        success *= system.evaluateResidual() == 0;
        success *= isEqual(efd.read(), static_cast<ScalarT>(1.0), kTol);

        const auto* residual = exciter.getResidual().getData();
        for (IdxT row = 0; row < exciter.size(); ++row)
        {
          success *= isEqual(residual[row], static_cast<ScalarT>(0.0), kTol);
        }

        return success.report(__func__);
      }
    };

  } // namespace Testing
} // namespace GridKit
