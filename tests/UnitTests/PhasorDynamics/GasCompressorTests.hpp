/**
 * @file GasCompressorTests.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Tests for the gas compressor model.
 */
#include <limits>
#include <type_traits>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Constants.hpp>
#include <GridKit/Definitions.hpp>
#include <GridKit/Model/PhasorDynamics/Gas/GasBus/GasBus.hpp>
#include <GridKit/Model/PhasorDynamics/Gas/GasBusInfinite/GasBusInfinite.hpp>
#include <GridKit/Model/PhasorDynamics/Gas/GasCompressor/GasCompressor.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class GasCompressorTests
    {
    private:
      using RealT              = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;
      using GasBusDataT        = PhasorDynamics::Gas::GasBusData<RealT, IdxT>;
      using GasCompressorDataT = PhasorDynamics::Gas::GasCompressorData<RealT, IdxT>;
      using GasBusT            = PhasorDynamics::Gas::GasBus<ScalarT, IdxT>;
      using GasBusInfiniteT    = PhasorDynamics::Gas::GasBusInfinite<ScalarT, IdxT>;
      using GasCompressorT     = PhasorDynamics::Gas::GasCompressor<ScalarT, IdxT>;
      using DependencyMap      = DependencyTracking::Variable::DependencyMap;

      static constexpr RealT tol_   = 10 * std::numeric_limits<RealT>::epsilon();
      static constexpr RealT alpha_ = 2.0;

      static GasBusDataT busData(RealT pressure)
      {
        GasBusDataT data;
        data.p0 = pressure;
        return data;
      }

      /// A compressor with ratio 1.5, between buses at 4 and 6 MPa
      static GasCompressorDataT compressorData()
      {
        GasCompressorDataT data;
        data.device_class                                                    = "GasCompressor";
        data.disambiguation_string                                           = "C";
        data.parameters[PhasorDynamics::Gas::GasCompressorParameters::ratio] = RealT{1.5};
        return data;
      }

      /// Allocate and initialize, numbering the bus variables after the compressor's
      template <typename T>
      static int connect(PhasorDynamics::Gas::GasBusBase<T, IdxT>&    bus1,
                         PhasorDynamics::Gas::GasBusBase<T, IdxT>&    bus2,
                         PhasorDynamics::Gas::GasCompressor<T, IdxT>& compressor)
      {
        bus1.allocate();
        bus2.allocate();
        compressor.allocate();

        IdxT next = compressor.size();
        for (auto* bus : {&bus1, &bus2})
        {
          if (bus->size() > 0)
          {
            bus->setVariableIndex(0, next);
            bus->setResidualIndex(0, next);
            ++next;
          }
        }

        bus1.initialize();
        bus2.initialize();
        return compressor.initialize();
      }

      static void assign(double& variable, double value)
      {
        variable = value;
      }

      static void assign(DependencyTracking::Variable& variable, double value)
      {
        variable.setValue(value);
      }

      /// Bus pressures rising at 2000 and 2500 Pa/s, with 10 kg/s through the compressor
      template <typename T>
      static void setState(PhasorDynamics::Gas::GasBus<T, IdxT>&        bus1,
                           PhasorDynamics::Gas::GasBus<T, IdxT>&        bus2,
                           PhasorDynamics::Gas::GasCompressor<T, IdxT>& compressor)
      {
        assign(bus1.yp().getData()[0], 2000.0);
        assign(bus2.yp().getData()[0], 2500.0);
        assign(compressor.y().getData()[0], 10.0);
      }

    public:
      GasCompressorTests()  = default;
      ~GasCompressorTests() = default;

      /// The mass flow is the compressor's one variable, and it is algebraic
      TestOutcome constructor()
      {
        TestStatus success = true;

        GasBusT        bus1(busData(4.0e6));
        GasBusT        bus2(busData(6.0e6));
        GasCompressorT compressor(&bus1, &bus2, compressorData());
        compressor.allocate();
        compressor.tagDifferentiable();

        success *= (compressor.size() == 1);
        success *= !compressor.tag()[0];

        return success.report(__func__);
      }

      TestOutcome parameterValidation()
      {
        TestStatus success = true;

        GasBusT         bus1(busData(4.0e6));
        GasBusT         bus2(busData(6.0e6));
        GasBusInfiniteT supply1(busData(4.0e6));
        GasBusInfiniteT supply2(busData(6.0e6));

        success *= GasCompressorT(&bus1, &bus2, compressorData()).verify().passed();
        success *= GasCompressorT(&supply1, &bus2, compressorData()).verify().passed();
        success *= !GasCompressorT(nullptr, &bus2, compressorData()).verify().passed();
        success *= !GasCompressorT(&supply1, &supply2, compressorData()).verify().passed();

        auto data                                                             = compressorData();
        data.parameters[PhasorDynamics::Gas::GasCompressorParameters::ratio]  = RealT{0.0};
        success                                                              *= !GasCompressorT(&bus1, &bus2, data).verify().passed();

        return success.report(__func__);
      }

      /// Initialization requires the bus pressures to satisfy the ratio
      TestOutcome initialization()
      {
        TestStatus success = true;

        GasBusT        bus1(busData(4.0e6));
        GasBusT        bus2(busData(6.0e6));
        GasCompressorT compressor(&bus1, &bus2, compressorData());
        success *= (connect(bus1, bus2, compressor) == 0);
        success *= isEqual(compressor.y().getData()[0], ScalarT{0.0});

        GasBusT        bus3(busData(5.0e6));
        GasCompressorT mismatched(&bus1, &bus3, compressorData());
        success *= (connect(bus1, bus3, mismatched) != 0);

        return success.report(__func__);
      }

      /// Residual away from steady flow, derived by hand
      TestOutcome residual()
      {
        TestStatus success = true;

        GasBusT        bus1(busData(4.0e6));
        GasBusT        bus2(busData(6.0e6));
        GasCompressorT compressor(&bus1, &bus2, compressorData());
        connect(bus1, bus2, compressor);
        setState(bus1, bus2, compressor);

        bus1.evaluateResidual();
        bus2.evaluateResidual();
        compressor.evaluateResidual();

        // 1.5 (2000) - 2500
        success *= isEqual(compressor.getResidual().getData()[0], ScalarT{500.0}, tol_);
        success *= isEqual(bus1.q(), ScalarT{-10.0}, tol_);
        success *= isEqual(bus2.q(), ScalarT{10.0}, tol_);

        return success.report(__func__);
      }

#ifdef GRIDKIT_ENABLE_ENZYME
      /// The hand-written Jacobian matches DependencyTracking, bus rows included
      TestOutcome jacobian()
      {
        TestStatus success = true;

        const auto reference    = jacobianRows<DependencyTracking::Variable>();
        const auto hand_written = jacobianRows<ScalarT>();

        success *= (reference.size() == hand_written.size());
        for (size_t i = 0; i < reference.size(); ++i)
        {
          success *= isEqual(reference[i], hand_written[i], tol_);
        }

        return success.report(__func__);
      }

      /// A fixed-pressure suction bus has no row and no column in the Jacobian
      TestOutcome fixedPressureBus()
      {
        TestStatus success = true;

        GasBusInfiniteT supply(busData(4.0e6));
        GasBusT         bus2(busData(6.0e6));
        GasCompressorT  compressor(&supply, &bus2, compressorData());
        connect(supply, bus2, compressor);

        compressor.updateTime(0.0, alpha_);
        bus2.evaluateResidual();
        compressor.evaluateResidual();
        compressor.evaluateJacobian();

        // The ratio row keeps the bus-2 column, and bus 2 keeps the flow column
        success *= (compressor.getCooJacobian()->getNnz() == 2);

        return success.report(__func__);
      }

    private:
      /// Jacobian rows of the compressor followed by the rows of buses 1 and 2
      template <typename T>
      std::vector<DependencyMap> jacobianRows()
      {
        PhasorDynamics::Gas::GasBus<T, IdxT>        bus1(busData(4.0e6));
        PhasorDynamics::Gas::GasBus<T, IdxT>        bus2(busData(6.0e6));
        PhasorDynamics::Gas::GasCompressor<T, IdxT> compressor(&bus1, &bus2, compressorData());
        connect(bus1, bus2, compressor);
        setState(bus1, bus2, compressor);

        compressor.updateTime(0.0, alpha_);
        bus1.evaluateResidual();
        bus2.evaluateResidual();
        compressor.evaluateResidual();
        compressor.evaluateJacobian();

        if constexpr (std::is_same_v<T, DependencyTracking::Variable>)
        {
          auto rows = MapFromCsr(compressor.getCsrJacobian());
          rows.push_back(withDerivatives(bus1.q().getDependencies()));
          rows.push_back(withDerivatives(bus2.q().getDependencies()));
          return rows;
        }
        else
        {
          auto*                      jacobian = compressor.getCooJacobian();
          std::vector<DependencyMap> rows(static_cast<size_t>(compressor.size()) + 2);
          for (IdxT i = 0; i < jacobian->getNnz(); ++i)
          {
            rows[jacobian->getRowData()[i]][jacobian->getColData()[i]] += jacobian->getValues()[i];
          }
          return rows;
        }
      }

      /// Fold derivative dependencies (odd numbers) into their variables, scaled by alpha
      static DependencyMap withDerivatives(const DependencyMap& dependencies)
      {
        DependencyMap row;
        for (const auto& [number, value] : dependencies)
        {
          RealT scale = ONE<RealT>;
          if (number % 2 == 1)
          {
            scale = alpha_;
          }
          row[number / 2] += scale * value;
        }
        return row;
      }
#endif
    };
  } // namespace Testing
} // namespace GridKit
