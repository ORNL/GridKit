#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class BusTests
    {
    private:
      using RealT = typename PhasorDynamics::Bus<ScalarT, IdxT>::RealT;

    public:
      BusTests()  = default;
      ~BusTests() = default;

      /// Constructor, allocation, and initialization checks
      TestOutcome constructor()
      {
        TestStatus success = true;

        ScalarT Vr{1.0};
        ScalarT Vi{2.0};

        PhasorDynamics::BusBase<ScalarT, IdxT>* bus = nullptr;

        // Create an infinite bus
        bus      = new PhasorDynamics::BusInfinite<ScalarT, IdxT>();
        success *= isEqual(bus->Vr(), 0.0);
        success *= isEqual(bus->Vi(), 0.0);
        delete bus;

        bus      = new PhasorDynamics::BusInfinite<ScalarT, IdxT>(Vr, Vi);
        success *= isEqual(bus->Vr(), Vr);
        success *= isEqual(bus->Vi(), Vi);
        delete bus;

        // Create an default bus
        bus = new PhasorDynamics::Bus<ScalarT, IdxT>();
        bus->allocate();
        bus->initialize();
        success *= isEqual(bus->Vr(), 0.0);
        success *= isEqual(bus->Vi(), 0.0);
        delete bus;

        bus = new PhasorDynamics::Bus<ScalarT, IdxT>(Vr, Vi);
        bus->allocate();
        bus->initialize();
        success *= isEqual(bus->Vr(), Vr);
        success *= isEqual(bus->Vi(), Vi);
        delete bus;

        bus = nullptr;

        return success.report(__func__);
      }

      /// Accessor method tests
      TestOutcome residual()
      {
        TestStatus success = true;

        ScalarT Vr{1.0};
        ScalarT Vi{2.0};
        ScalarT Ir{1.0};
        ScalarT Ii{2.0};

        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus_inf;
        bus_inf.Ir()  = Ir;
        success      *= isEqual(bus_inf.Ir(), Ir);
        bus_inf.Ii()  = Ii;
        success      *= isEqual(bus_inf.Ii(), Ii);

        bus_inf.evaluateResidual();
        success *= isEqual(bus_inf.Ir(), 0.0);
        success *= isEqual(bus_inf.Ii(), 0.0);

        PhasorDynamics::Bus<ScalarT, IdxT> bus(Vr, Vi);
        bus.allocate();
        bus.initialize();
        success *= isEqual(bus.Vr(), Vr);
        success *= isEqual(bus.Vi(), Vi);

        bus.Ir()  = Ir;
        success  *= isEqual(bus.Ir(), Ir);
        bus.Ii()  = Ii;
        success  *= isEqual(bus.Ii(), Ii);
        bus.getResidual().setDataUpdated();

        bus.evaluateResidual();
        success *= isEqual(bus.Ir(), 0.0);
        success *= isEqual(bus.Ii(), 0.0);

        return success.report(__func__);
      }

      /// Fault current applied and cleared through setFault
      TestOutcome fault()
      {
        TestStatus success = true;

        ScalarT Vr{1.0};
        ScalarT Vi{2.0};
        RealT   R{1.0};
        RealT   X{2.0};

        PhasorDynamics::Bus<ScalarT, IdxT> bus(Vr, Vi);
        bus.allocate();
        bus.initialize();

        // Fault current I = -V / (R + jX)
        const std::complex<RealT> current = -std::complex<RealT>(Vr, Vi) / std::complex<RealT>(R, X);

        success *= bus.setFault(true, R, X) == 0;
        bus.evaluateResidual();
        success *= isEqual(bus.Ir(), current.real());
        success *= isEqual(bus.Ii(), current.imag());

        success *= bus.setFault(false, R, X) == 0;
        bus.evaluateResidual();
        success *= isEqual(bus.Ir(), 0.0);
        success *= isEqual(bus.Ii(), 0.0);

        // A zero fault impedance is rejected
        success *= bus.setFault(true, 0.0, 0.0) != 0;

        // An infinite bus cannot be faulted
        PhasorDynamics::BusInfinite<ScalarT, IdxT> bus_inf;
        success *= bus_inf.setFault(true, R, X) != 0;

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// Fault Jacobian matches DependencyTracking and keeps its pattern when cleared
      TestOutcome jacobian()
      {
        TestStatus success = true;

        RealT R{1.0};
        RealT X{2.0};

        // Jacobian via DependencyTracking, which numbers y entries 2 * index
        DependencyTracking::Variable Vr{1.0};
        DependencyTracking::Variable Vi{2.0};

        PhasorDynamics::Bus<DependencyTracking::Variable, IdxT> dependency_tracking_bus(Vr, Vi);
        dependency_tracking_bus.allocate();
        dependency_tracking_bus.initialize();
        dependency_tracking_bus.setFault(true, R, X);
        dependency_tracking_bus.evaluateResidual();

        const auto* f    = dependency_tracking_bus.getResidual().getData();
        const auto  size = static_cast<size_t>(dependency_tracking_bus.size());

        std::vector<DependencyTracking::Variable::DependencyMap> dependency_tracking_jacobian(size);
        for (size_t row = 0; row < size; ++row)
        {
          for (const auto& [number, value] : f[row].getDependencies())
          {
            dependency_tracking_jacobian[row][number / 2] = value;
          }
        }

        // Jacobian from the bus
        PhasorDynamics::Bus<ScalarT, IdxT> bus(1.0, 2.0);
        bus.allocate();
        bus.initialize();
        bus.setFault(true, R, X);
        bus.evaluateResidual();
        bus.evaluateJacobian();

        auto*        jacobian = bus.getCooJacobian();
        const IdxT*  rows     = jacobian->getRowData();
        const IdxT*  cols     = jacobian->getColData();
        const RealT* vals     = jacobian->getValues();

        std::vector<DependencyTracking::Variable::DependencyMap> bus_jacobian(size);
        for (IdxT i = 0; i < jacobian->getNnz(); ++i)
        {
          bus_jacobian[static_cast<size_t>(rows[i])][static_cast<size_t>(cols[i])] = vals[i];
        }

        for (size_t row = 0; row < size; ++row)
        {
          success *= isEqual(dependency_tracking_jacobian[row], bus_jacobian[row]);
        }

        // Clearing the fault keeps the same entries with zero values
        const IdxT nnz = jacobian->getNnz();
        bus.setFault(false, R, X);
        bus.evaluateResidual();
        bus.evaluateJacobian();
        success *= bus.nnz() == nnz;
        for (IdxT i = 0; i < nnz; ++i)
        {
          success *= isEqual(vals[i], 0.0);
        }

        return success.report(__func__);
      }
#endif
    };

  } // namespace Testing
} // namespace GridKit
