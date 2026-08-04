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
    /// Connection tests for pairs of components that share a signal node.
    /// Each case checks that the node links both components and that they
    /// agree on its value after initialization. Solver-driven cases live
    /// in @ref PDIntegrationTests.
    template <typename scalar_type, typename index_type>
    class ComponentConnectionTests
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;

      ComponentConnectionTests()  = default;
      ~ComponentConnectionTests() = default;

      // The tolerance only absorbs floating-point roundoff.
      static constexpr RealT kTol = std::numeric_limits<RealT>::epsilon();

      /// GENROU initializes first and writes the field voltage it needs to
      /// the shared node. ESDC1A then initializes around that value and
      /// must leave it unchanged at a steady state.
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

        // At zero power the required field voltage equals the terminal voltage.
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
