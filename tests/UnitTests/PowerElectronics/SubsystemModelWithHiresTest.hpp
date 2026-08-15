#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include <GridKit/Model/PowerElectronics/Bus/MicrogridBus.hpp>
#include <GridKit/Model/PowerElectronics/CircuitComponent.hpp>
#include <GridKit/Model/PowerElectronics/NodeBase.hpp>
#include <GridKit/Model/PowerElectronics/PartitionInterface/BusPartitionInterface.hpp>
#include <GridKit/Model/PowerElectronics/SubsystemModel.hpp>
#include <GridKit/Model/PowerElectronics/SystemModelPowerElectronics.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{

  /*!
   * @brief Hires Component 1 class.
   *
   */
  template <class ScalarT, typename IdxT>
  class HiresComponent1 : public CircuitComponent<ScalarT, IdxT>
  {
    using RealT = typename CircuitComponent<ScalarT, IdxT>::RealT;
    using NodeT = typename PowerElectronics::NodeBase<ScalarT, IdxT>;

    using CircuitComponent<ScalarT, IdxT>::size_;
    using CircuitComponent<ScalarT, IdxT>::nnz_;
    using CircuitComponent<ScalarT, IdxT>::alpha_;
    using CircuitComponent<ScalarT, IdxT>::y_ext_;
    using CircuitComponent<ScalarT, IdxT>::y_int_;
    using CircuitComponent<ScalarT, IdxT>::yp_ext_;
    using CircuitComponent<ScalarT, IdxT>::yp_int_;
    using CircuitComponent<ScalarT, IdxT>::abs_tol_;
    using CircuitComponent<ScalarT, IdxT>::f_ext_;
    using CircuitComponent<ScalarT, IdxT>::f_int_;
    using CircuitComponent<ScalarT, IdxT>::idc_;

    using CircuitComponent<ScalarT, IdxT>::extern_indices_;
    using CircuitComponent<ScalarT, IdxT>::n_extern_;
    using CircuitComponent<ScalarT, IdxT>::n_intern_;

  public:
    HiresComponent1(NodeT* bus, IdxT id)
      : node_ref_(bus)
    {
      size_           = 5;
      n_intern_       = 3;
      n_extern_       = 2;
      extern_indices_ = {0, 1};
      idc_            = id;
      nnz_            = 12;
    }

    ~HiresComponent1()
    {
    }

    int allocate() final
    {
      CircuitComponent<ScalarT, IdxT>::allocate();

      this->setExternalConnectionNodes(0, node_ref_->getNodeConnection(0));
      this->setExternalConnectionNodes(1, node_ref_->getNodeConnection(1));

      return 0;
    }

    int initialize() final
    {
      return 0;
    }

    int tagDifferentiable() final
    {
      return 0;
    }

    int evaluateInternalResidual() final
    {
      // Internals
      f_int_[0] = -yp_int_[0] - 1.71 * y_int_[0] + 0.43 * y_int_[1] + 8.32 * y_int_[2] + 0.0007;
      f_int_[1] = -yp_int_[1] + 1.71 * y_int_[0] - 8.75 * y_int_[1];
      f_int_[2] = -yp_int_[2] - 10.03 * y_int_[2] + 0.43 * *y_ext_[0] + 0.035 * *y_ext_[1];

      return 0;
    }

    int evaluateExternalResidual()
    {
      // outputs
      *f_ext_[0] += 8.32 * y_int_[1] + 1.71 * y_int_[2] - 0.1 * *y_ext_[0];
      *f_ext_[1] += -0.7 * *y_ext_[1];

      return 0;
    }

    int evaluateJacobian() final
    {

      this->zeroJacMatrix();

      // Internal Jacobian Entries
      std::vector<IdxT>  row = {2, 2, 2, 3, 3, 4, 4, 4};
      std::vector<IdxT>  col = {2, 3, 4, 2, 3, 4, 0, 1};
      std::vector<RealT> val = {-1.71 - alpha_, 0.43, 8.32, 1.71, -8.75 - alpha_, -10.03 - alpha_, 0.43, 0.035};

      this->setJacValues(row, col, val);

      // External Jacobian Entries
      row = {0, 0, 0, 1};
      col = {3, 4, 0, 1};
      val = {8.32, 1.71, -0.1, -0.7};

      this->setJacValues(row, col, val);

      return 0;
    }

    int evaluateIntegrand() final
    {
      return 0;
    }

    int initializeAdjoint() final
    {
      return 0;
    }

    int evaluateAdjointResidual() final
    {
      return 0;
    }

    int evaluateAdjointIntegrand() final
    {
      return 0;
    }

    /**
     * @brief Compute the absolute tolerance for each variable in the model
     */
    int setAbsoluteTolerance(RealT rel_tol) final
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return 0;
    }

  private:
    NodeT* node_ref_;
  };

  /*!
   * @brief Hires Bus Component (Component 2).
   *
   */
  template <class ScalarT, typename IdxT>
  class HiresBus : public CircuitComponent<ScalarT, IdxT>
  {

    using RealT = typename CircuitComponent<ScalarT, IdxT>::RealT;
    using NodeT = typename PowerElectronics::NodeBase<ScalarT, IdxT>;

    using CircuitComponent<ScalarT, IdxT>::size_;
    using CircuitComponent<ScalarT, IdxT>::nnz_;
    using CircuitComponent<ScalarT, IdxT>::alpha_;
    using CircuitComponent<ScalarT, IdxT>::y_ext_;
    using CircuitComponent<ScalarT, IdxT>::y_int_;
    using CircuitComponent<ScalarT, IdxT>::yp_ext_;
    using CircuitComponent<ScalarT, IdxT>::yp_int_;
    using CircuitComponent<ScalarT, IdxT>::abs_tol_;
    using CircuitComponent<ScalarT, IdxT>::f_ext_;
    using CircuitComponent<ScalarT, IdxT>::f_int_;
    using CircuitComponent<ScalarT, IdxT>::idc_;

    using CircuitComponent<ScalarT, IdxT>::extern_indices_;
    using CircuitComponent<ScalarT, IdxT>::n_extern_;
    using CircuitComponent<ScalarT, IdxT>::n_intern_;

  public:
    HiresBus(NodeT* bus, IdxT id)
      : node_ref_(bus)
    {
      size_           = 2;
      n_intern_       = 0;
      n_extern_       = 2;
      extern_indices_ = {0, 1};
      idc_            = id;
      nnz_            = 2;
    }

    ~HiresBus()
    {
    }

    int allocate() final
    {
      CircuitComponent<ScalarT, IdxT>::allocate();

      this->setExternalConnectionNodes(0, node_ref_->getNodeConnection(0));
      this->setExternalConnectionNodes(1, node_ref_->getNodeConnection(1));

      return 0;
    }

    int initialize() final
    {
      return 0;
    }

    int tagDifferentiable() final
    {
      return 0;
    }

    int evaluateInternalResidual() final
    {
      return 0;
    }

    int evaluateExternalResidual() final
    {
      *f_ext_[0] += -*yp_ext_[0] - *y_ext_[0];
      *f_ext_[1] += -*yp_ext_[1] - *y_ext_[1];

      return 0;
    }

    int evaluateJacobian() final
    {
      this->zeroJacMatrix();

      std::vector<IdxT>  row = {0, 1};
      std::vector<IdxT>  col = {0, 1};
      std::vector<RealT> val = {-alpha_ - 1.0, -alpha_ - 1.0};

      this->setJacValues(row, col, val);

      return 0;
    }

    int evaluateIntegrand() final
    {
      return 0;
    }

    int initializeAdjoint() final
    {
      return 0;
    }

    int evaluateAdjointResidual() final
    {
      return 0;
    }

    int evaluateAdjointIntegrand() final
    {
      return 0;
    }

    /**
     * @brief Compute the absolute tolerance for each variable in the model
     */
    int setAbsoluteTolerance(RealT rel_tol) final
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return 0;
    }

  private:
    NodeT* node_ref_;
  };

  /*!
   * @brief Hires Component 3 class.
   *
   */
  template <class ScalarT, typename IdxT>
  class HiresComponent3 : public CircuitComponent<ScalarT, IdxT>
  {
    using RealT = typename CircuitComponent<ScalarT, IdxT>::RealT;
    using NodeT = typename PowerElectronics::NodeBase<ScalarT, IdxT>;

    using CircuitComponent<ScalarT, IdxT>::size_;
    using CircuitComponent<ScalarT, IdxT>::nnz_;
    using CircuitComponent<ScalarT, IdxT>::alpha_;
    using CircuitComponent<ScalarT, IdxT>::y_ext_;
    using CircuitComponent<ScalarT, IdxT>::y_int_;
    using CircuitComponent<ScalarT, IdxT>::yp_ext_;
    using CircuitComponent<ScalarT, IdxT>::yp_int_;
    using CircuitComponent<ScalarT, IdxT>::abs_tol_;
    using CircuitComponent<ScalarT, IdxT>::f_ext_;
    using CircuitComponent<ScalarT, IdxT>::f_int_;
    using CircuitComponent<ScalarT, IdxT>::idc_;

    using CircuitComponent<ScalarT, IdxT>::extern_indices_;
    using CircuitComponent<ScalarT, IdxT>::n_extern_;
    using CircuitComponent<ScalarT, IdxT>::n_intern_;

  public:
    HiresComponent3(NodeT* bus, IdxT id)
      : node_ref_(bus)
    {
      size_           = 5;
      n_intern_       = 3;
      n_extern_       = 2;
      extern_indices_ = {0, 1};
      idc_            = id;
      nnz_            = 15;
    }

    ~HiresComponent3()
    {
    }

    int allocate()
    {
      CircuitComponent<ScalarT, IdxT>::allocate();

      this->setExternalConnectionNodes(0, node_ref_->getNodeConnection(0));
      this->setExternalConnectionNodes(1, node_ref_->getNodeConnection(1));

      return 0;
    }

    int initialize()
    {
      return 0;
    }

    int tagDifferentiable()
    {
      return 0;
    }

    int evaluateInternalResidual()
    {

      // Internals
      f_int_[0] = -yp_int_[0] - 280 * y_int_[0] * y_int_[2] + 0.69 * *y_ext_[0] + 1.71 * *y_ext_[1] - 0.43 * y_int_[0] + 0.69 * y_int_[1];
      f_int_[1] = -yp_int_[1] + 280 * y_int_[0] * y_int_[2] - 1.81 * y_int_[1];
      f_int_[2] = -yp_int_[2] - 280 * y_int_[0] * y_int_[2] + 1.81 * y_int_[1];

      return 0;
    }

    int evaluateExternalResidual()
    {
      // Externals
      *f_ext_[0] += -0.02 * *y_ext_[0];
      *f_ext_[1] += -0.045 * *y_ext_[1] + 0.43 * y_int_[0] + 0.43 * y_int_[1];

      return 0;
    }

    int evaluateJacobian()
    {
      this->zeroJacMatrix();

      // Internal Jacobian Entries [row 1]
      std::vector<IdxT>  row = {2, 2, 2, 2, 2};
      std::vector<IdxT>  col = {2, 3, 4, 0, 1};
      std::vector<RealT> val = {-280 * y_int_[2] - 0.43 - alpha_, 0.69, -280 * y_int_[0], 0.69, 1.71};

      this->setJacValues(row, col, val);

      // Internal Jacobian Entries [row 2]
      row = {3, 3, 3};
      col = {2, 3, 4};
      val = {280 * y_int_[2], -1.81 - alpha_, 280 * y_int_[0]};

      this->setJacValues(row, col, val);

      // Internal Jacobian Entries [row 3]
      row = {4, 4, 4};
      col = {2, 3, 4};
      val = {-280 * y_int_[2], 1.81, -280 * y_int_[0] - alpha_};

      this->setJacValues(row, col, val);

      // External Jacobian Entries
      row = {0, 1, 1, 1};
      col = {0, 2, 3, 1};
      val = {-0.02, 0.43, 0.43, -0.045};

      this->setJacValues(row, col, val);

      return 0;
    }

    int evaluateIntegrand()
    {
      return 0;
    }

    int initializeAdjoint()
    {
      return 0;
    }

    int evaluateAdjointResidual()
    {
      return 0;
    }

    int evaluateAdjointIntegrand()
    {
      return 0;
    }

    /**
     * @brief Compute the absolute tolerance for each variable in the model
     *
     */
    int setAbsoluteTolerance(RealT rel_tol)
    {
      abs_tol_.setToConst(static_cast<ScalarT>(rel_tol));
      return 0;
    }

  private:
    NodeT* node_ref_;
  };

  namespace Testing
  {

    /**
     * This example assembles the HIRES problem and partitions it to demonstrate
     * that the partitioned system can reconstruct the full system.
     *
     * The HIRES problem consists of eight ODE equations and is divided into three
     * components. HiresComponent1 contains equations 1, 2, and 3,
     * HiresComponent3 contains equations 6, 7, and 8, and HiresBus contains
     * equations 4 and 5. HiresBus plays a role similar to MicrogridBusDQ, where
     * other components connected to the bus contribute terms to its equations.
     *
     * In the equations below, terms enclosed in parentheses ( ) are contributions
     * from HiresComponent1, terms enclosed in square brackets [ ] are contributions
     * from HiresComponent3, and the remaining terms in equations 4 and 5 belong
     * to HiresBus.
     *\f[
     * f_1 = dy_1/dt + 1.71y_1 - 0.43y_2 - 8.32y_3 - 0.0007
     *
     * f_2 = dy_2/dt - 1.71y_1 + 8.75y_2
     *
     * f_3 = dy_3/dt + 10.03y_3 - 0.43y_4 - 0.035y_5
     *
     * f_4 = dy_4/dt + y_4 + (0.1y_4 - 8.32y_2 - 1.71y_3) + [0.02y_4]
     *
     * f_5 = dy_5/dt + y_5 + \left(0.7y_5 \right) + \left[0.045y_5 - 0.43y_6 - 0.43y_7 \right]
     *
     * f_6 = dy_6/dt - 280y_6y_8 + 0.69y_4 + 1.71y_5 - 0.43y_6 + 0.69y_7
     *
     * f_7 = dy_7/dt + 280y_6y_8 - 1.81y_7
     *
     * f_8 = dy_8/dt - 280y_6y_8 + 1.81y_7
     *\f]
     *
     * The assembled system has the following structure:
     *
     * (HiresComponent1) -------- (HiresBus) -------- (HiresComponent3)
     *
     * The system is partitioned between HiresBus and HiresComponent3. A bus
     * partition interface is introduced in the partition containing HiresBus to
     * preserve the contribution of HiresComponent3 to the bus equations. The
     * residuals evaluated independently by the partitions are then printed with
     * the residual of the full system for eye-ball comparision.
     *
     * This example will also be useful later for testing the order of accuracy of co-simulation methods.
     */
    template <class ScalarT, typename IdxT>
    class SubsystemModelWithHiresTests
    {
      using RealT     = typename CircuitComponent<ScalarT, IdxT>::RealT;
      using Bus       = PowerElectronics::MicrogridBus<ScalarT, IdxT>;
      using Subsystem = SubsystemModel<ScalarT, IdxT>;
      using System    = PowerElectronicsModel<ScalarT, IdxT>;

    public:
      /**
       * @brief Construct the HIRES reference system and subsystem partitions.
       */
      SubsystemModelWithHiresTests()
        : system_(new System()),
          partition1_(new Subsystem()),
          partition2_(new Subsystem()),
          comp1_(new HiresComponent1<ScalarT, IdxT>(&bus_, 1)),
          bus1_(new HiresBus<ScalarT, IdxT>(&bus_, 2)),
          comp3_(new HiresComponent3<ScalarT, IdxT>(&bus_, 3))
      {

        y_  = {1, 2, 3, 4, 5, 6, 7, 8};
        yp_ = {1, 2, 3, 4, 5, 6, 7, 8};
        // ---------------------------------------------------------------------
        // Assemble and allocate the monolithic reference system
        // ---------------------------------------------------------------------

        system_->addComponent(comp1_);
        system_->addComponent(comp3_);
        system_->addComponent(bus1_);
        system_->addNode(&bus_);

        system_->allocate();

        distributeVariables(y_, yp_);
        system_->updateTime(1.0, 2.0);
        system_->evaluateResidual();

        // ---------------------------------------------------------------------
        // Construct the partition interface
        // ---------------------------------------------------------------------

        auto* comp3_copy = new HiresComponent3<ScalarT, IdxT>(*comp3_);
        bus_interface_   = new BusPartitionInterface<ScalarT, IdxT>(&bus_, comp3_copy, 4);

        bus_interface_->allocate();

        // ---------------------------------------------------------------------
        // Assemble the subsystem partitions
        // ---------------------------------------------------------------------

        partition1_->addComponent(comp1_);
        partition1_->addComponent(bus1_);
        partition1_->addInterface(bus_interface_);
        partition1_->addNode(&bus_);

        partition2_->addComponent(comp3_);

        partition1_->allocate();
        partition2_->allocate();
        partitions_ = {partition1_, partition2_};
      }

      ~SubsystemModelWithHiresTests()
      {
        delete partition1_;
        delete partition2_;
        delete system_;
      }

      void distributeVariables(const std::vector<ScalarT>& y, const std::vector<ScalarT>& yp)
      {
        auto* system_y  = system_->y().getData();
        auto* system_yp = system_->yp().getData();

        for (size_t i = 0; i < system_->size(); ++i)
        {
          system_y[i]  = y[i];
          system_yp[i] = yp[i];
        }

        system_->y().setDataUpdated();
        system_->yp().setDataUpdated();

        for (auto* partition : partitions_)
        {
          for (size_t i = 0; i < partition->getExternSize(); ++i)
          {
            const auto global_index = partition->getExternalDataIndices()[i];

            partition->getExternalDataY()[i]  = y[global_index];
            partition->getExternalDataYP()[i] = yp[global_index];
          }

          auto* partition_y  = partition->y().getData();
          auto* partition_yp = partition->yp().getData();

          for (size_t i = 0; i < partition->getInternalSize(); ++i)
          {
            const auto global_index = partition->getNodeConnection(i);

            partition_y[i]  = y[global_index];
            partition_yp[i] = yp[global_index];
          }

          partition->y().setDataUpdated();
          partition->yp().setDataUpdated();
        }
      }

      /**
       * @brief Verify that the partitioned HIRES residual matches the
       *        monolithic residual.
       */
      TestOutcome residual()
      {
        TestStatus success = true;

        // Distribute variables to all partitions
        distributeVariables(y_, yp_);

        std::vector<ScalarT> partition_residual(system_->size(), 0.0);

        for (auto* partition : partitions_)
        {
          partition->updateTime(1.0, 2.0);
          partition->evaluateResidual();

          auto* residual = partition->getResidual().getData();

          for (size_t i = 0; i < partition->getInternalSize(); ++i)
          {
            partition_residual[partition->getNodeConnection(i)] = residual[i];
          }

          partition->getResidual().setDataUpdated();
        }

        auto* reference_residual = system_->getResidual().getData();

        RealT max_error = 0.0;

        for (size_t i = 0; i < system_->size(); ++i)
        {
          double error = std::abs(partition_residual[i] - reference_residual[i]) / std::abs(reference_residual[i] + 1);
          max_error    = std::max(max_error, error);
        }

        std::cout << "max error " << max_error << std::endl;

        success *= max_error <= std::numeric_limits<RealT>::epsilon();

        return success.report(__func__);
      }

      /**
       * @brief Verify that each subsystem Jacobian matches the corresponding
       *        entries of the monolithic HIRES Jacobian.
       */
      TestOutcome jacobian()
      {
        TestStatus success = true;

        RealT alpha = 2.0;

        distributeVariables(y_, yp_);

        /*
         * GridKit stores the HIRES variables in component assembly order:
         *
         *   System index:  0   1   2   3   4   5   6   7
         *   HIRES index:   0   1   2   5   6   7   3   4
         *   Variable:      y1  y2  y3  y6  y7  y8  y4  y5
         */
        const std::array<IdxT, 8> sysmodel_to_hires = {
            0, 1, 2, 5, 6, 7, 3, 4};

        const std::array<IdxT, 8> hires_to_sysmodel = {
            0, 1, 2, 6, 7, 3, 4, 5};

        const RealT y6 = y_[hires_to_sysmodel[5]];
        const RealT y8 = y_[hires_to_sysmodel[7]];

        std::array<std::array<RealT, 8>, 8> reference_jac =
            {{{-alpha - 1.71, 0.43, 8.32, 0.0, 0.0, 0.0, 0.0, 0.0},
              {1.71, -alpha - 8.75, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
              {0.0, 0.0, -alpha - 10.03, 0.43, 0.035, 0.0, 0.0, 0.0},
              {0.0, 8.32, 1.71, -alpha - 1.12, 0.0, 0.0, 0.0, 0.0},
              {0.0, 0.0, 0.0, 0.0, -alpha - 1.745, 0.43, 0.43, 0.0},
              {0.0, 0.0, 0.0, 0.69, 1.71, -alpha - 280.0 * y8 - 0.43, 0.69, -280.0 * y6},
              {0.0, 0.0, 0.0, 0.0, 0.0, 280.0 * y8, -alpha - 1.81, 280.0 * y6},
              {0.0, 0.0, 0.0, 0.0, 0.0, -280.0 * y8, 1.81, -alpha - 280.0 * y6}}};

        for (auto* partition : partitions_)
        {
          partition->updateTime(1.0, alpha);
          partition->evaluateJacobian();

          auto* partition_jac = partition->getCsrJacobian();

          const auto* row_ptr = partition_jac->getRowData();
          const auto* cols    = partition_jac->getColData();
          const auto* vals    = partition_jac->getValues();

          const size_t n = partition->getInternalSize();

          std::vector<std::vector<RealT>> dense_jac(n, std::vector<RealT>(n, 0.0));

          // Convert the partition CSR Jacobian to a dense matrix.
          for (size_t row = 0; row < n; ++row)
          {
            for (IdxT k = row_ptr[row]; k < row_ptr[row + 1]; ++k)
            {
              dense_jac[row][cols[k]] = vals[k];
            }
          }

          // Compare with the corresponding entries of the full Jacobian.
          for (size_t row = 0; row < n; ++row)
          {
            const IdxT global_row = partition->getNodeConnection(row);
            const IdxT ref_row    = sysmodel_to_hires[global_row];

            for (size_t col = 0; col < n; ++col)
            {
              const IdxT global_col = partition->getNodeConnection(col);
              const IdxT ref_col    = sysmodel_to_hires[global_col];

              success *= std::abs(dense_jac[row][col] - reference_jac[ref_row][ref_col]) <= std::numeric_limits<RealT>::epsilon();
            }
          }
        }

        return success.report(__func__);
      }

    private:
      Bus bus_;

      System* system_;

      Subsystem* partition1_;
      Subsystem* partition2_;

      HiresComponent1<ScalarT, IdxT>* comp1_;
      HiresBus<ScalarT, IdxT>*        bus1_;
      HiresComponent3<ScalarT, IdxT>* comp3_;

      BusPartitionInterface<ScalarT, IdxT>* bus_interface_;

      std::vector<Subsystem*> partitions_;

      std::vector<ScalarT> y_;
      std::vector<ScalarT> yp_;
    };

  } // namespace Testing

} // namespace GridKit
