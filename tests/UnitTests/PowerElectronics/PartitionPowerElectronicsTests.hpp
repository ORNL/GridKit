
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <numbers>
#include <utility>
#include <vector>

#include <GridKit/Model/PowerElectronics/Bus/MicrogridBus.hpp>
#include <GridKit/Model/PowerElectronics/Bus/SignalNode.hpp>
#include <GridKit/Model/PowerElectronics/Component.hpp>
#include <GridKit/Model/PowerElectronics/DistributedGenerator/DistributedGenerator.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridBusDQ/MicrogridBusDQ.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridLine/MicrogridLine.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridLoad/MicrogridLoad.hpp>
#include <GridKit/Model/PowerElectronics/PartitionInterface/BusLinePartitionInterface.hpp>
#include <GridKit/Model/PowerElectronics/SubsystemModel.hpp>
#include <GridKit/Model/PowerElectronics/SystemModel.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

/*
 *
 * This network is partitioned into two subsystems, and the partitioned residuals and Jacobians
 * are compared against those of the monolithic system.
 *
 *                  Signal Node
 *                       |
 *                 (shared signal)
 *
 *       DG1                         DG2
 *        |                           |
 *        |                           |
 *       Bus1 -------- Line -------- Bus2
 *        |                           |
 *        |                           |
 *      Load1                       Load2
 *
 * The partition boundary is located between the transmission line and Bus2.
 *
 */
namespace GridKit
{
  namespace PowerElectronics
  {
    namespace Testing
    {
      template <class ScalarT, typename IdxT>
      class PartitionPowerElectronicsTests
      {
        using TestOutcome = GridKit::Testing::TestOutcome;
        using TestStatus  = GridKit::Testing::TestStatus;
        using RealT       = typename Component<ScalarT, IdxT>::RealT;
        using System      = SystemModel<ScalarT, IdxT>;
        using Subsystem   = SubsystemModel<ScalarT, IdxT>;

        using Bus       = MicrogridBus<ScalarT, IdxT>;
        using Signal    = SignalNode<ScalarT, IdxT>;
        using Generator = DistributedGenerator<ScalarT, IdxT>;
        using BusDQ     = MicrogridBusDQ<ScalarT, IdxT>;
        using Line      = MicrogridLine<ScalarT, IdxT>;
        using Load      = MicrogridLoad<ScalarT, IdxT>;

        using Parameters = DistributedGeneratorParameters<ScalarT, IdxT>;
        using Interface  = BusLinePartitionInterface<ScalarT, IdxT>;

      public:
        PartitionPowerElectronicsTests()
        {
          createComponents();
          assembleMonolithic();
          saveReferenceResidual();
          saveReferenceJacobian();
          assemblePartitions();
          distributeVariables();
        }

        /**
         * @brief Verify sequential partition residuals against
         *        the saved monolithic reference.
         */
        TestOutcome residual()
        {
          const auto partition_residual = evaluatePartitions(false);

          return compareResiduals(partition_residual, reference_residual_).report(__func__);
        }

        /**
         * @brief Verify that parallel partition evaluation produces
         *        the same result as sequential evaluation.
         */
        TestOutcome parallelResidual()
        {
          const auto sequential = evaluatePartitions(false);
          const auto parallel   = evaluatePartitions(true);

          TestStatus success = true;

          success *= compareResiduals(parallel, sequential);
          success *= compareResiduals(parallel, reference_residual_);

          return success.report(__func__);
        }

        /**
         * @brief Verify that the subsystem Jacobians match the corresponding
         *        blocks of the monolithic Jacobian.
         */
        TestOutcome jacobian()
        {
          TestStatus success = true;

          const RealT tolerance = 100 * std::numeric_limits<RealT>::epsilon();

          for (auto* partition : partitions_)
          {
            partition->updateTime(1.0, 2.0);
            partition->evaluateJacobian();

            const auto jac_map        = GridKit::Testing::MapFromCsr(partition->getCsrJacobian());
            const auto partition_size = partition->getInternalSize();

            for (size_t row = 0; row < partition_size; ++row)
            {
              const auto global_row = partition->getNodeConnection(row);

              for (size_t col = 0; col < partition_size; ++col)
              {
                const auto global_col = partition->getNodeConnection(col);

                const auto& partition_row = jac_map[row];
                const auto& reference_row = reference_jacobian_[global_row];

                const RealT partition_value = partition_row.contains(col) ? partition_row.at(col) : RealT{0};
                const RealT reference_value = reference_row.contains(global_col) ? reference_row.at(global_col) : RealT{0};

                success *= GridKit::Testing::isEqual(partition_value, reference_value, tolerance);
              }
            }
          }

          return success.report(__func__);
        }

      private:
        /**
         * @brief Construct the physical microgrid components.
         */
        void createComponents()
        {
          constexpr RealT pi = std::numbers::pi_v<RealT>;

          Parameters parms{};

          parms.wb_  = 2.0 * pi * 50.0;
          parms.wc_  = 31.41;
          parms.mp_  = 9.4e-5;
          parms.Vn_  = 380.0;
          parms.nq_  = 1.3e-3;
          parms.F_   = 0.75;
          parms.Kiv_ = 420.0;
          parms.Kpv_ = 0.1;
          parms.Kic_ = 2.0e4;
          parms.Kpc_ = 15.0;
          parms.Cf_  = 5.0e-5;
          parms.rLf_ = 0.1;
          parms.Lf_  = 1.35e-3;
          parms.rLc_ = 0.03;
          parms.Lc_  = 0.35e-3;

          dg1_ = new Generator(1, parms, true, &signal_, &bus1_);
          dg2_ = new Generator(2, parms, false, &signal_, &bus2_);

          line_ = new Line(3, 0.23, 0.1 / (2.0 * pi * 50.0), &signal_, &bus1_, &bus2_);

          load1_ = new Load(4, 3.0, 2.0 / (2.0 * pi * 50.0), &signal_, &bus1_);
          load2_ = new Load(5, 2.0, 1.0 / (2.0 * pi * 50.0), &signal_, &bus2_);

          busDQ1_ = new BusDQ(6, 1.0e4, &bus1_);
          busDQ2_ = new BusDQ(7, 1.0e4, &bus2_);
        }

        /**
         * @brief Assemble the monolithic reference system.
         */
        void assembleMonolithic()
        {
          system_ = std::make_unique<System>(true);

          system_->addComponent(dg1_);
          system_->addComponent(load1_);
          system_->addComponent(busDQ1_);
          system_->addComponent(line_);
          system_->addComponent(dg2_);
          system_->addComponent(load2_);
          system_->addComponent(busDQ2_);

          system_->addNode(&signal_);
          system_->addNode(&bus1_);
          system_->addNode(&bus2_);

          system_->allocate();

          y_.resize(system_->size());
          yp_.resize(system_->size());

          for (size_t i = 0; i < system_->size(); ++i)
          {
            const auto index = static_cast<RealT>(i + 1);

            y_[i]  = static_cast<ScalarT>(0.01 * index);
            yp_[i] = static_cast<ScalarT>(0.001 * index);
          }
        }

        /**
         * @brief Evaluate the monolithic residual exactly once,
         *        before constructing the partitions.
         */
        void saveReferenceResidual()
        {
          auto* system_y  = system_->y().getData();
          auto* system_yp = system_->yp().getData();

          for (size_t i = 0; i < system_->size(); ++i)
          {
            system_y[i]  = y_[i];
            system_yp[i] = yp_[i];
          }

          system_->y().setDataUpdated();
          system_->yp().setDataUpdated();

          system_->updateTime(1.0, 2.0);
          system_->evaluateResidual();

          const auto* residual = system_->getResidual().getData();

          reference_residual_.assign(residual, residual + system_->size());
        }

        /**
         * @brief Save the monolithic Jacobian before partitioning.
         */
        void saveReferenceJacobian()
        {
          system_->updateTime(1.0, 2.0);
          system_->evaluateJacobian();

          reference_jacobian_ = GridKit::Testing::MapFromCsr(system_->getCsrJacobian());
        }

        /**
         * @brief Assemble the two subsystem partitions.
         */
        void assemblePartitions()
        {
          partition1_ = std::make_unique<Subsystem>();
          partition2_ = std::make_unique<Subsystem>();

          // Boundary between the line and Bus2.
          interface_ = new Interface(&bus2_, line_, 8);
          interface_->allocate();

          // Partition 1: DG1, Load1, Bus1, Line
          partition1_->addComponent(dg1_);
          partition1_->addComponent(load1_);
          partition1_->addComponent(busDQ1_);
          partition1_->addComponent(line_);

          partition1_->addNode(&signal_);
          partition1_->addNode(&bus1_);

          // Partition 2: DG2, Load2, Bus2
          partition2_->addComponent(dg2_);
          partition2_->addComponent(load2_);
          partition2_->addComponent(busDQ2_);
          partition2_->addInterface(interface_);

          partition2_->addNode(&bus2_);

          partition1_->allocate();
          partition2_->allocate();

          partitions_ = {partition1_.get(), partition2_.get()};
        }

        /**
         * @brief Distribute global states to the partitions.
         */
        void distributeVariables()
        {
          for (auto* partition : partitions_)
          {
            auto forcing = [this, partition](ScalarT) -> void
            {
              auto* y_ext  = partition->getExternalDataY().getData();
              auto* yp_ext = partition->getExternalDataYP().getData();

              const auto& indices = partition->getExternalDataIndices();

              for (size_t i = 0; i < partition->getExternSize(); ++i)
              {
                y_ext[i]  = y_[indices[i]];
                yp_ext[i] = yp_[indices[i]];
              }
            };

            partition->setForcingFunction(std::move(forcing));

            auto* y  = partition->y().getData();
            auto* yp = partition->yp().getData();

            for (size_t i = 0; i < partition->getInternalSize(); ++i)
            {
              const auto index = partition->getNodeConnection(i);

              y[i]  = y_[index];
              yp[i] = yp_[index];
            }

            partition->y().setDataUpdated();
            partition->yp().setDataUpdated();
          }
        }

        /**
         * @brief Evaluate partitions and reconstruct the global residual.
         *
         * @param parallel Enable OpenMP parallel evaluation.
         */
        std::vector<ScalarT> evaluatePartitions([[maybe_unused]] bool parallel)
        {
          const int num_partitions = static_cast<int>(partitions_.size());

#ifdef _OPENMP
#pragma omp parallel for if (parallel) schedule(static)
#endif
          for (int p = 0; p < num_partitions; ++p)
          {
            auto* partition = partitions_[p];

            partition->updateTime(1.0, 2.0);
            partition->evaluateResidual();
          }

          std::vector<ScalarT> residual(reference_residual_.size(), 0.0);

          for (auto* partition : partitions_)
          {
            const auto* f_values = partition->getResidual().getData();

            for (size_t i = 0; i < partition->getInternalSize(); ++i)
            {
              const auto index = partition->getNodeConnection(i);
              residual[index]  = f_values[i];
            }
          }

          return residual;
        }

        /**
         * @brief Compare two residual vectors.
         */
        TestStatus compareResiduals(
            const std::vector<ScalarT>& actual,
            const std::vector<ScalarT>& expected) const
        {
          TestStatus success = true;

          const RealT tolerance = 100 * std::numeric_limits<RealT>::epsilon();

          for (size_t i = 0; i < expected.size(); ++i)
          {
            success *= GridKit::Testing::isEqual(actual[i], expected[i], tolerance);
          }

          return success;
        }

        // nodes
        Signal signal_;
        Bus    bus1_;
        Bus    bus2_;

        // components
        Generator* dg1_    = nullptr;
        Generator* dg2_    = nullptr;
        Line*      line_   = nullptr;
        Load*      load1_  = nullptr;
        Load*      load2_  = nullptr;
        BusDQ*     busDQ1_ = nullptr;
        BusDQ*     busDQ2_ = nullptr;

        // Partition interface
        Interface* interface_ = nullptr;

        // Reference system
        std::unique_ptr<System> system_;

        // Partitioned system
        std::unique_ptr<Subsystem> partition1_;
        std::unique_ptr<Subsystem> partition2_;
        std::vector<Subsystem*>    partitions_;

        // Test data
        std::vector<ScalarT> y_;
        std::vector<ScalarT> yp_;
        std::vector<ScalarT> reference_residual_;

        std::vector<DependencyTracking::Variable::DependencyMap> reference_jacobian_;
      };

    } // namespace Testing
  } // namespace PowerElectronics
} // namespace GridKit
