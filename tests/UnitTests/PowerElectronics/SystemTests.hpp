#pragma once

#include <format>
#include <iostream>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/PowerElectronics/DistributedGenerator/DistributedGenerator.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridBusDQ/MicrogridBusDQ.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridLine/MicrogridLine.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridLoad/MicrogridLoad.hpp>
#include <GridKit/Model/PowerElectronics/SystemModelPowerElectronics.hpp>
#include <GridKit/Utilities/TestHelpers.hpp>
#include <GridKit/Utilities/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {

    template <class ScalarT, typename IdxT>
    class SystemTests
    {
    public:
      SystemTests()  = default;
      ~SystemTests() = default;

      TestOutcome jacobianSparsity()
      {
        TestStatus success = true;

        double abs_tol         = 1.0e-8;
        double rel_tol         = 1.0e-8;
        size_t max_step_number = 3000;
        bool   use_jac         = true;

        // Create model
        auto* sysmodel = new GridKit::PowerElectronicsModel<DependencyTracking::Variable, size_t>(rel_tol, abs_tol, use_jac, max_step_number);

        // Modeled after the problem in the paper
        double RN = 1.0e4;

        // DG Params
        GridKit::DistributedGeneratorParameters<double, size_t> parms1;
        parms1.wb_  = 2.0 * M_PI * 50.0;
        parms1.wc_  = 31.41;
        parms1.mp_  = 9.4e-5;
        parms1.Vn_  = 380.0;
        parms1.nq_  = 1.3e-3;
        parms1.F_   = 0.75;
        parms1.Kiv_ = 420.0;
        parms1.Kpv_ = 0.1;
        parms1.Kic_ = 2.0e4;
        parms1.Kpc_ = 15.0;
        parms1.Cf_  = 5.0e-5;
        parms1.rLf_ = 0.1;
        parms1.Lf_  = 1.35e-3;
        parms1.rLc_ = 0.03;
        parms1.Lc_  = 0.35e-3;

        GridKit::DistributedGeneratorParameters<double, size_t> parms2;
        // Parameters from MATLAB Microgrid code for first DG
        parms2.wb_  = 2.0 * M_PI * 50.0;
        parms2.wc_  = 31.41;
        parms2.mp_  = 12.5e-5;
        parms2.Vn_  = 380.0;
        parms2.nq_  = 1.5e-3;
        parms2.F_   = 0.75;
        parms2.Kiv_ = 390.0;
        parms2.Kpv_ = 0.05;
        parms2.Kic_ = 16.0e3;
        parms2.Kpc_ = 10.5;
        parms2.Cf_  = 50.0e-6;
        parms2.rLf_ = 0.1;
        parms2.Lf_  = 1.35e-3;
        parms2.rLc_ = 0.03;
        parms2.Lc_  = 0.35e-3;

        // Line params
        double rline1 = 0.23;
        double Lline1 = 0.1 / (2.0 * M_PI * 50.0);

        double rline2 = 0.35;
        double Lline2 = 0.58 / (2.0 * M_PI * 50.0);

        double rline3 = 0.23;
        double Lline3 = 0.1 / (2.0 * M_PI * 50.0);

        // load parms
        double rload1 = 3.0;
        double Lload1 = 2.0 / (2.0 * M_PI * 50.0);

        double rload2 = 2.0;
        double Lload2 = 1.0 / (2.0 * M_PI * 50.0);

        // indexing sets
        size_t Nsize              = 2;
        //							DGs	+		- refframe	   Lines +				Loads
        size_t vec_size_internals = 13 * (2 * Nsize) - 1 + (2 + 4 * (Nsize - 1)) + 2 * Nsize;
        //							\omegaref + BusDQ
        size_t vec_size_externals = 1 + 2 * (2 * Nsize);
        size_t dqbus1             = vec_size_internals + 1;
        size_t dqbus2             = vec_size_internals + 3;
        size_t dqbus3             = vec_size_internals + 5;
        size_t dqbus4             = vec_size_internals + 7;

        size_t vec_size_total = vec_size_internals + vec_size_externals;

        size_t indexv = 0;

        // dg 1
        GridKit::DistributedGenerator<DependencyTracking::Variable, size_t>* dg1 = new GridKit::DistributedGenerator<DependencyTracking::Variable, size_t>(0, parms1, true);
        // ref motor
        dg1->setExternalConnectionNodes(0, vec_size_internals);
        // outputs
        dg1->setExternalConnectionNodes(1, dqbus1);
        dg1->setExternalConnectionNodes(2, dqbus1 + 1);
        //"grounding" of the difference
        dg1->setExternalConnectionNodes(3, static_cast<size_t>(-1));
        // internal connections
        for (size_t i = 0; i < 12; i++)
        {
          dg1->setExternalConnectionNodes(4 + i, indexv + i);
        }
        indexv += 12;
        sysmodel->addComponent(dg1);

        // dg 2
        GridKit::DistributedGenerator<DependencyTracking::Variable, size_t>* dg2 = new GridKit::DistributedGenerator<DependencyTracking::Variable, size_t>(1, parms1, false);
        // ref motor
        dg2->setExternalConnectionNodes(0, vec_size_internals);
        // outputs
        dg2->setExternalConnectionNodes(1, dqbus2);
        dg2->setExternalConnectionNodes(2, dqbus2 + 1);
        // internal connections
        for (size_t i = 0; i < 13; i++)
        {
          dg2->setExternalConnectionNodes(3 + i, indexv + i);
        }
        indexv += 13;
        sysmodel->addComponent(dg2);

        // dg 3
        GridKit::DistributedGenerator<DependencyTracking::Variable, size_t>* dg3 = new GridKit::DistributedGenerator<DependencyTracking::Variable, size_t>(2, parms2, false);
        // ref motor
        dg3->setExternalConnectionNodes(0, vec_size_internals);
        // outputs
        dg3->setExternalConnectionNodes(1, dqbus3);
        dg3->setExternalConnectionNodes(2, dqbus3 + 1);
        // internal connections
        for (size_t i = 0; i < 13; i++)
        {
          dg3->setExternalConnectionNodes(3 + i, indexv + i);
        }
        indexv += 13;
        sysmodel->addComponent(dg3);

        // dg 4
        GridKit::DistributedGenerator<DependencyTracking::Variable, size_t>* dg4 = new GridKit::DistributedGenerator<DependencyTracking::Variable, size_t>(3, parms2, false);
        // ref motor
        dg4->setExternalConnectionNodes(0, vec_size_internals);
        // outputs
        dg4->setExternalConnectionNodes(1, dqbus4);
        dg4->setExternalConnectionNodes(2, dqbus4 + 1);

        // internal connections
        for (size_t i = 0; i < 13; i++)
        {
          dg4->setExternalConnectionNodes(3 + i, indexv + i);
        }
        indexv += 13;
        sysmodel->addComponent(dg4);

        // Lines

        // line 1
        GridKit::MicrogridLine<DependencyTracking::Variable, size_t>* l1 = new GridKit::MicrogridLine<DependencyTracking::Variable, size_t>(4, rline1, Lline1);
        // ref motor
        l1->setExternalConnectionNodes(0, vec_size_internals);
        // input connections
        l1->setExternalConnectionNodes(1, dqbus1);
        l1->setExternalConnectionNodes(2, dqbus1 + 1);
        // output connections
        l1->setExternalConnectionNodes(3, dqbus2);
        l1->setExternalConnectionNodes(4, dqbus2 + 1);
        // internal connections
        for (size_t i = 0; i < 2; i++)
        {
          l1->setExternalConnectionNodes(5 + i, indexv + i);
        }
        indexv += 2;
        sysmodel->addComponent(l1);

        // line 2
        GridKit::MicrogridLine<DependencyTracking::Variable, size_t>* l2 = new GridKit::MicrogridLine<DependencyTracking::Variable, size_t>(5, rline2, Lline2);
        // ref motor
        l2->setExternalConnectionNodes(0, vec_size_internals);
        // input connections
        l2->setExternalConnectionNodes(1, dqbus2);
        l2->setExternalConnectionNodes(2, dqbus2 + 1);
        // output connections
        l2->setExternalConnectionNodes(3, dqbus3);
        l2->setExternalConnectionNodes(4, dqbus3 + 1);
        // internal connections
        for (size_t i = 0; i < 2; i++)
        {
          l2->setExternalConnectionNodes(5 + i, indexv + i);
        }
        indexv += 2;
        sysmodel->addComponent(l2);

        // line 3
        GridKit::MicrogridLine<DependencyTracking::Variable, size_t>* l3 = new GridKit::MicrogridLine<DependencyTracking::Variable, size_t>(6, rline3, Lline3);
        // ref motor
        l3->setExternalConnectionNodes(0, vec_size_internals);
        // input connections
        l3->setExternalConnectionNodes(1, dqbus3);
        l3->setExternalConnectionNodes(2, dqbus3 + 1);
        // output connections
        l3->setExternalConnectionNodes(3, dqbus4);
        l3->setExternalConnectionNodes(4, dqbus4 + 1);
        // internal connections
        for (size_t i = 0; i < 2; i++)
        {
          l3->setExternalConnectionNodes(5 + i, indexv + i);
        }
        indexv += 2;
        sysmodel->addComponent(l3);

        //  loads

        // load 1
        GridKit::MicrogridLoad<DependencyTracking::Variable, size_t>* load1 = new GridKit::MicrogridLoad<DependencyTracking::Variable, size_t>(7, rload1, Lload1);
        // ref motor
        load1->setExternalConnectionNodes(0, vec_size_internals);
        // input connections
        load1->setExternalConnectionNodes(1, dqbus1);
        load1->setExternalConnectionNodes(2, dqbus1 + 1);
        // internal connections
        for (size_t i = 0; i < 2; i++)
        {
          load1->setExternalConnectionNodes(3 + i, indexv + i);
        }
        indexv += 2;
        sysmodel->addComponent(load1);

        // load 2
        GridKit::MicrogridLoad<DependencyTracking::Variable, size_t>* load2 = new GridKit::MicrogridLoad<DependencyTracking::Variable, size_t>(8, rload2, Lload2);
        // ref motor
        load2->setExternalConnectionNodes(0, vec_size_internals);
        // input connections
        load2->setExternalConnectionNodes(1, dqbus3);
        load2->setExternalConnectionNodes(2, dqbus3 + 1);
        // internal connections
        for (size_t i = 0; i < 2; i++)
        {
          load2->setExternalConnectionNodes(3 + i, indexv + i);
        }
        indexv += 2;
        sysmodel->addComponent(load2);

        // Virtual PQ Buses
        GridKit::MicrogridBusDQ<DependencyTracking::Variable, size_t>* bus1 = new GridKit::MicrogridBusDQ<DependencyTracking::Variable, size_t>(9, RN);

        bus1->setExternalConnectionNodes(0, dqbus1);
        bus1->setExternalConnectionNodes(1, dqbus1 + 1);
        sysmodel->addComponent(bus1);

        GridKit::MicrogridBusDQ<DependencyTracking::Variable, size_t>* bus2 = new GridKit::MicrogridBusDQ<DependencyTracking::Variable, size_t>(10, RN);

        bus2->setExternalConnectionNodes(0, dqbus2);
        bus2->setExternalConnectionNodes(1, dqbus2 + 1);
        sysmodel->addComponent(bus2);

        GridKit::MicrogridBusDQ<DependencyTracking::Variable, size_t>* bus3 = new GridKit::MicrogridBusDQ<DependencyTracking::Variable, size_t>(11, RN);

        bus3->setExternalConnectionNodes(0, dqbus3);
        bus3->setExternalConnectionNodes(1, dqbus3 + 1);
        sysmodel->addComponent(bus3);

        GridKit::MicrogridBusDQ<DependencyTracking::Variable, size_t>* bus4 = new GridKit::MicrogridBusDQ<DependencyTracking::Variable, size_t>(12, RN);

        bus4->setExternalConnectionNodes(0, dqbus4);
        bus4->setExternalConnectionNodes(1, dqbus4 + 1);
        sysmodel->addComponent(bus4);

        sysmodel->allocate(vec_size_total);

        // Create Intial points for states
        for (size_t i = 0; i < vec_size_total; i++)
        {
          sysmodel->y()[i]  = 0.0;
          sysmodel->yp()[i] = 0.0;
          sysmodel->y()[i].setVariableNumber(i);
        }

        // Create Intial derivatives specifics generated in MATLAB
        // DGs 1
        sysmodel->yp()[2]      = parms1.Vn_;
        sysmodel->yp()[4]      = parms1.Kpv_ * parms1.Vn_;
        sysmodel->yp()[6]      = (parms1.Kpc_ * parms1.Kpv_ * parms1.Vn_) / parms1.Lf_;
        sysmodel->yp()[12 + 3] = parms1.Vn_;
        sysmodel->yp()[12 + 5] = parms1.Kpv_ * parms1.Vn_;
        sysmodel->yp()[12 + 7] = (parms1.Kpc_ * parms1.Kpv_ * parms1.Vn_) / parms1.Lf_;
        for (size_t i = 2; i < 4; i++)
        {
          sysmodel->yp()[13 * i - 1 + 3] = parms2.Vn_;
          sysmodel->yp()[13 * i - 1 + 5] = parms2.Kpv_ * parms2.Vn_;
          sysmodel->yp()[13 * i - 1 + 7] = (parms2.Kpc_ * parms2.Kpv_ * parms2.Vn_) / parms2.Lf_;
        }

        // since the intial P_com = 0
        sysmodel->y()[vec_size_internals] = parms1.wb_;
        sysmodel->y()[vec_size_internals].setVariableNumber(vec_size_internals);

        for (size_t i = 0; i < vec_size_total; i++)
        {
          sysmodel->yp()[i].setVariableNumber(i);
        }

        sysmodel->initialize();
        sysmodel->evaluateResidual();

        std::vector<DependencyTracking::Variable>& residuals = sysmodel->getResidual();

        size_t* row_indices = sysmodel->getJacRowIndices();
        size_t* col_indices = sysmodel->getJacColIndices();

        for (size_t row = 0; row < residuals.size(); row++)
        {
          bool row_correct = true;

          DependencyTracking::Variable& row_residual   = residuals[row];
          auto&                         dependency_map = row_residual.getDependencies();
          auto                          dependencies   = std::vector(dependency_map.begin(), dependency_map.end());

          // Verify both rows are the same size
          row_correct &= (row_indices[row + 1] - row_indices[row]) == dependencies.size();

          if (row_correct)
          {
            for (size_t col_idx = 0; col_idx < dependencies.size(); col_idx++)
            {
              row_correct &= col_indices[col_idx + row_indices[row]] == std::get<0>(dependencies[col_idx]);
            }
          }

          if (!row_correct)
          {
            std::cerr << "Row " << row << " is incorrect:\n";
            std::cerr << "            CSR Row: {";
            for (size_t col_idx = row_indices[row]; col_idx < row_indices[row + 1]; col_idx++)
            {
              std::cerr << std::format("{:3}{}", col_indices[col_idx], col_idx == row_indices[row + 1] - 1 ? "" : ", ");
            }
            std::cerr << "}\n"
                      << "    Dep Tracker Row: {";
            for (size_t col_idx = 0; col_idx < dependencies.size(); col_idx++)
            {
              std::cerr << std::format("{:3}{}", std::get<0>(dependencies[col_idx]), col_idx == dependencies.size() - 1 ? "" : ", ");
            }
            std::cerr << "}\n\n";
          }

          success *= row_correct;
        }

        delete sysmodel;

        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
