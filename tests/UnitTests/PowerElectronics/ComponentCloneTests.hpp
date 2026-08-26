
#include <GridKit/Model/PowerElectronics/Bus/MicrogridBus.hpp>
#include <GridKit/Model/PowerElectronics/Bus/SignalNode.hpp>
#include <GridKit/Model/PowerElectronics/DistributedGenerator/DistributedGenerator.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridBusDQ/MicrogridBusDQ.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridLine/MicrogridLine.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridLoad/MicrogridLoad.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  template <typename ComponentT, typename ScalarT, typename IdxT>
  bool verifyComponentClone(ComponentT& component)
  {
    using RealT = typename CircuitComponent<ScalarT, IdxT>::RealT;

    bool success = true;

    auto* clone = dynamic_cast<ComponentT*>(component.clone());

    if (clone == nullptr)
    {
      return false;
    }

    /**************************************************************************
     * Verify the Clone Initially Matches the Original
     **************************************************************************/

    success &= clone != &component;
    success &= clone->size() == component.size();
    success &= clone->nnz() == component.nnz();
    success &= clone->getInternalSize() == component.getInternalSize();
    success &= clone->getExternSize() == component.getExternSize();
    success &= clone->getExternIndices() == component.getExternIndices();
    success &= clone->getIDcomponent() == component.getIDcomponent();

    /**************************************************************************
     * Verify Connection Independence
     **************************************************************************/

    for (IdxT i = 0; i < component.size(); ++i)
    {
      const IdxT connection = component.getNodeConnection(i);

      success &= clone->getNodeConnection(i) == connection;

      clone->setConnectionNodes(i, connection + 1);

      success &= component.getNodeConnection(i) == connection;
      success &= clone->getNodeConnection(i) == connection + 1;

      clone->setConnectionNodes(i, connection);
    }

    /**************************************************************************
     * Verify State Independence
     **************************************************************************/

    auto checkVectorIndependence = [&success](auto& original, auto& copy)
    {
      success &= original.getSize() == copy.getSize();

      if (original.getSize() == 0)
      {
        return;
      }

      auto* original_data = original.getData();
      auto* copy_data     = copy.getData();

      success &= original_data != copy_data;

      const auto original_value = original_data[0];

      copy_data[0] = original_value + 1.0;

      success &= original_data[0] == original_value;
      success &= copy_data[0] == original_value + 1.0;

      copy_data[0] = original_value;
    };

    checkVectorIndependence(component.y(), clone->y());
    checkVectorIndependence(component.yp(), clone->yp());
    checkVectorIndependence(component.getResidual(), clone->getResidual());
    checkVectorIndependence(component.absoluteTolerance(), clone->absoluteTolerance());
    checkVectorIndependence(component.param(), clone->param());

    /**************************************************************************
     * Verify Jacobian Independence
     **************************************************************************/

    if (component.nnz() > 0)
    {
      auto* original_rows   = component.jacobianCooRows();
      auto* original_cols   = component.jacobianCooCols();
      auto* original_values = component.jacobianCooValues();

      auto* clone_rows   = clone->jacobianCooRows();
      auto* clone_cols   = clone->jacobianCooCols();
      auto* clone_values = clone->jacobianCooValues();

      success &= clone_rows != original_rows;
      success &= clone_cols != original_cols;
      success &= clone_values != original_values;

      for (IdxT i = 0; i < component.nnz(); ++i)
      {
        success &= clone_rows[i] == original_rows[i];
        success &= clone_cols[i] == original_cols[i];
        success &= clone_values[i] == original_values[i];
      }

      const IdxT  original_row   = original_rows[0];
      const IdxT  original_col   = original_cols[0];
      const RealT original_value = original_values[0];

      clone_rows[0]   = original_row + 1;
      clone_cols[0]   = original_col + 1;
      clone_values[0] = original_value + 1.0;

      success &= original_rows[0] == original_row;
      success &= original_cols[0] == original_col;
      success &= original_values[0] == original_value;

      clone_rows[0]   = original_row;
      clone_cols[0]   = original_col;
      clone_values[0] = original_value;
    }

    delete clone;

    return success;
  }

  namespace Testing
  {
    template <typename ScalarT, typename IdxT>
    class CircuitComponentCloneTests
    {
      using SignalNode          = PowerElectronics::SignalNode<ScalarT, IdxT>;
      using Bus                 = PowerElectronics::MicrogridBus<ScalarT, IdxT>;
      using BusDQ               = MicrogridBusDQ<ScalarT, IdxT>;
      using Generator           = DistributedGenerator<ScalarT, IdxT>;
      using GeneratorParameters = DistributedGeneratorParameters<ScalarT, IdxT>;
      using Line                = MicrogridLine<ScalarT, IdxT>;
      using Load                = MicrogridLoad<ScalarT, IdxT>;

    public:
      CircuitComponentCloneTests()
      {
        /**************************************************************************
         * Construct Network Nodes
         **************************************************************************/

        signal_.allocate();

        bus1_.allocate();
        bus2_.allocate();

        /**************************************************************************
         * Distributed Generator Parameters
         **************************************************************************/

        generator_parameters_.wb_  = 2.0 * M_PI * 50.0;
        generator_parameters_.wc_  = 31.41;
        generator_parameters_.mp_  = 9.4e-5;
        generator_parameters_.Vn_  = 380.0;
        generator_parameters_.nq_  = 1.3e-3;
        generator_parameters_.F_   = 0.75;
        generator_parameters_.Kiv_ = 420.0;
        generator_parameters_.Kpv_ = 0.1;
        generator_parameters_.Kic_ = 2.0e4;
        generator_parameters_.Kpc_ = 15.0;
        generator_parameters_.Cf_  = 5.0e-5;
        generator_parameters_.rLf_ = 0.1;
        generator_parameters_.Lf_  = 1.35e-3;
        generator_parameters_.rLc_ = 0.03;
        generator_parameters_.Lc_  = 0.35e-3;

        /**************************************************************************
         * Construct Components
         **************************************************************************/

        generator_ = new Generator(1, generator_parameters_, true, &signal_, &bus1_);

        line_ = new Line(2, 0.23, 0.1 / (2.0 * M_PI * 50.0), &signal_, &bus1_, &bus2_);

        load_ = new Load(3, 3.0, 2.0 / (2.0 * M_PI * 50.0), &signal_, &bus1_);

        bus_dq_ = new BusDQ(4, 1.0e4, &bus1_);

        /**************************************************************************
         * Allocate Components
         **************************************************************************/

        generator_->allocate();
        line_->allocate();
        load_->allocate();
        bus_dq_->allocate();
      }

      ~CircuitComponentCloneTests()
      {
        delete generator_;
        delete line_;
        delete load_;
        delete bus_dq_;
      }

      TestOutcome distributedGeneratorClone()
      {
        TestStatus success = true;

        success *= verifyComponentClone<Generator, ScalarT, IdxT>(*generator_);

        return success.report(__func__);
      }

      TestOutcome microgridLineClone()
      {
        TestStatus success = true;

        success *= verifyComponentClone<Line, ScalarT, IdxT>(*line_);

        return success.report(__func__);
      }

      TestOutcome microgridLoadClone()
      {
        TestStatus success = true;

        success *= verifyComponentClone<Load, ScalarT, IdxT>(*load_);

        return success.report(__func__);
      }

      TestOutcome microgridBusDQClone()
      {
        TestStatus success = true;

        success *= verifyComponentClone<BusDQ, ScalarT, IdxT>(*bus_dq_);

        return success.report(__func__);
      }

    private:
      SignalNode signal_;

      Bus bus1_;
      Bus bus2_;

      GeneratorParameters generator_parameters_;

      Generator* generator_{nullptr};
      Line*      line_{nullptr};
      Load*      load_{nullptr};
      BusDQ*     bus_dq_{nullptr};
    };
  } // namespace Testing
} // namespace GridKit
