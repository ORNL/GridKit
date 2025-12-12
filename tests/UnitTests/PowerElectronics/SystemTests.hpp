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

      TestOutcome jacobianSparsity(size_t Nsize)
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

        // DG Params Vector
        // All DGs have the same set of parameters except for the first two.
        GridKit::DistributedGeneratorParameters<double, size_t> DG_parms1;
        DG_parms1.wb_  = 2.0 * M_PI * 50.0;
        DG_parms1.wc_  = 31.41;
        DG_parms1.mp_  = 9.4e-5;
        DG_parms1.Vn_  = 380.0;
        DG_parms1.nq_  = 1.3e-3;
        DG_parms1.F_   = 0.75;
        DG_parms1.Kiv_ = 420.0;
        DG_parms1.Kpv_ = 0.1;
        DG_parms1.Kic_ = 2.0e4;
        DG_parms1.Kpc_ = 15.0;
        DG_parms1.Cf_  = 5.0e-5;
        DG_parms1.rLf_ = 0.1;
        DG_parms1.Lf_  = 1.35e-3;
        DG_parms1.rLc_ = 0.03;
        DG_parms1.Lc_  = 0.35e-3;

        GridKit::DistributedGeneratorParameters<double, size_t> DG_parms2;
        DG_parms2.wb_  = 2.0 * M_PI * 50.0;
        DG_parms2.wc_  = 31.41;
        DG_parms2.mp_  = 12.5e-5;
        DG_parms2.Vn_  = 380.0;
        DG_parms2.nq_  = 1.5e-3;
        DG_parms2.F_   = 0.75;
        DG_parms2.Kiv_ = 390.0;
        DG_parms2.Kpv_ = 0.05;
        DG_parms2.Kic_ = 16.0e3;
        DG_parms2.Kpc_ = 10.5;
        DG_parms2.Cf_  = 50.0e-6;
        DG_parms2.rLf_ = 0.1;
        DG_parms2.Lf_  = 1.35e-3;
        DG_parms2.rLc_ = 0.03;
        DG_parms2.Lc_  = 0.35e-3;

        std::vector<GridKit::DistributedGeneratorParameters<double, size_t>> DGParams_list(2 * Nsize, DG_parms2);

        DGParams_list[0] = DG_parms1;
        DGParams_list[1] = DG_parms1;

        // line vector params
        // Every odd line has the same parameters and every even line has the same parameters
        double              rline1 = 0.23;
        double              Lline1 = 0.1 / (2.0 * M_PI * 50.0);
        double              rline2 = 0.35;
        double              Lline2 = 0.58 / (2.0 * M_PI * 50.0);
        std::vector<double> rline_list(2 * Nsize - 1, 0.0);
        std::vector<double> Lline_list(2 * Nsize - 1, 0.0);
        for (size_t i = 0; i < rline_list.size(); i++)
        {
          rline_list[i] = (i % 2) ? rline2 : rline1;
          Lline_list[i] = (i % 2) ? Lline2 : Lline1;
        }

        // load parms
        // Only the first load has the same paramaters.
        double rload1 = 3.0;
        double Lload1 = 2.0 / (2.0 * M_PI * 50.0);
        double rload2 = 2.0;
        double Lload2 = 1.0 / (2.0 * M_PI * 50.0);

        std::vector<double> rload_list(Nsize, rload2);
        std::vector<double> Lload_list(Nsize, Lload2);
        rload_list[0] = rload1;
        Lload_list[0] = Lload1;

        size_t num_dgs   = 2 * Nsize;
        size_t num_lines = 2 * Nsize - 1;
        size_t num_loads = Nsize;
        size_t num_buses = 2 * Nsize;

        size_t vec_size_internals =
            //                                                                                  - refframe
            num_dgs * DistributedGenerator<DependencyTracking::Variable, size_t>::NUM_INTERNALS - 1
            + num_lines * MicrogridLine<DependencyTracking::Variable, size_t>::NUM_INTERNALS
            + num_loads * MicrogridLoad<DependencyTracking::Variable, size_t>::NUM_INTERNALS
            + num_buses * MicrogridBusDQ<DependencyTracking::Variable, size_t>::NUM_INTERNALS;

        size_t vec_size_externals =
            //                                                                              + refframe
            num_buses * MicrogridBusDQ<DependencyTracking::Variable, size_t>::NUM_EXTERNALS + 1;

        std::vector<size_t> vdqbus_index(2 * Nsize, 0);
        vdqbus_index[0] = vec_size_internals + 1;
        for (size_t i = 1; i < vdqbus_index.size(); i++)
        {
          vdqbus_index[i] = vdqbus_index[i - 1] + 2;
        }

        // Total size of the vector setup
        size_t vec_size_total = vec_size_internals + vec_size_externals;

        // Create the reference DG
        auto* dg_ref = new DistributedGenerator<DependencyTracking::Variable, size_t>(0,
                                                                                      DGParams_list[0],
                                                                                      true);
        // ref motor
        dg_ref->setExternalConnectionNodes(0, vec_size_internals);
        // outputs
        dg_ref->setExternalConnectionNodes(1, vdqbus_index[0]);
        dg_ref->setExternalConnectionNodes(2, vdqbus_index[0] + 1);
        //"grounding" of the difference
        dg_ref->setExternalConnectionNodes(3, static_cast<size_t>(-1));
        // internal connections
        for (size_t i = 0; i < dg_ref->NUM_INTERNALS - 1; i++)
        {
          dg_ref->setExternalConnectionNodes(4 + i, i);
        }
        sysmodel->addComponent(dg_ref);

        // Keep track of models and index location
        size_t indexv   = dg_ref->NUM_INTERNALS - 1;
        size_t model_id = 1;
        // Add all other DGs
        for (size_t i = 1; i < 2 * Nsize; i++)
        {
          // current DG to add
          auto* dg = new DistributedGenerator<DependencyTracking::Variable, size_t>(model_id++,
                                                                                    DGParams_list[i],
                                                                                    false);
          // ref motor
          dg->setExternalConnectionNodes(0, vec_size_internals);
          // outputs
          dg->setExternalConnectionNodes(1, vdqbus_index[i]);
          dg->setExternalConnectionNodes(2, vdqbus_index[i] + 1);
          // internal connections
          for (size_t j = 0; j < dg_ref->NUM_INTERNALS; j++)
          {
            dg->setExternalConnectionNodes(3 + j, indexv + j);
          }
          indexv += dg_ref->NUM_INTERNALS;
          sysmodel->addComponent(dg);
        }

        // Load all the Line compoenents
        for (size_t i = 0; i < 2 * Nsize - 1; i++)
        {
          // line
          auto* line_model = new MicrogridLine<DependencyTracking::Variable, size_t>(model_id++,
                                                                                     rline_list[i],
                                                                                     Lline_list[i]);
          // ref motor
          line_model->setExternalConnectionNodes(0, vec_size_internals);
          // input connections
          line_model->setExternalConnectionNodes(1, vdqbus_index[i]);
          line_model->setExternalConnectionNodes(2, vdqbus_index[i] + 1);
          // output connections
          line_model->setExternalConnectionNodes(3, vdqbus_index[i + 1]);
          line_model->setExternalConnectionNodes(4, vdqbus_index[i + 1] + 1);
          // internal connections
          for (size_t j = 0; j < 2; j++)
          {
            line_model->setExternalConnectionNodes(5 + j, indexv + j);
          }
          indexv += 2;
          sysmodel->addComponent(line_model);
        }

        //  Load all the Load components
        for (size_t i = 0; i < Nsize; i++)
        {
          auto* load_model = new MicrogridLoad<DependencyTracking::Variable, size_t>(model_id++,
                                                                                     rload_list[i],
                                                                                     Lload_list[i]);
          // ref motor
          load_model->setExternalConnectionNodes(0, vec_size_internals);
          // input connections
          load_model->setExternalConnectionNodes(1, vdqbus_index[2 * i]);
          load_model->setExternalConnectionNodes(2, vdqbus_index[2 * i] + 1);
          // internal connections
          for (size_t j = 0; j < 2; j++)
          {
            load_model->setExternalConnectionNodes(3 + j, indexv + j);
          }
          indexv += 2;
          sysmodel->addComponent(load_model);
        }

        // Add all the microgrid Virtual DQ Buses
        for (size_t i = 0; i < 2 * Nsize; i++)
        {
          auto* virDQbus_model = new MicrogridBusDQ<DependencyTracking::Variable, size_t>(model_id++, RN);

          virDQbus_model->setExternalConnectionNodes(0, vdqbus_index[i]);
          virDQbus_model->setExternalConnectionNodes(1, vdqbus_index[i] + 1);
          sysmodel->addComponent(virDQbus_model);
        }

        sysmodel->allocate(vec_size_total);

        // Create Intial points for states
        for (size_t i = 0; i < vec_size_total; i++)
        {
          sysmodel->y()[i]  = 0.0;
          sysmodel->yp()[i] = 0.0;
          sysmodel->y()[i].setVariableNumber(i);
        }

        // Create Intial derivatives specifics generated in MATLAB
        for (size_t i = 0; i < 2 * Nsize; i++)
        {
          sysmodel->yp()[dg_ref->NUM_INTERNALS * i - 1 + 3] = DGParams_list[i].Vn_;
          sysmodel->yp()[dg_ref->NUM_INTERNALS * i - 1 + 5] = DGParams_list[i].Kpv_ * DGParams_list[i].Vn_;
          sysmodel->yp()[dg_ref->NUM_INTERNALS * i - 1 + 7] = (DGParams_list[i].Kpc_ * DGParams_list[i].Kpv_ * DGParams_list[i].Vn_) / DGParams_list[i].Lf_;
        }

        // since the intial P_com = 0, the set the intial vector to the reference frame
        sysmodel->y()[vec_size_internals] = DG_parms1.wb_;
        sysmodel->y()[vec_size_internals].setVariableNumber(vec_size_internals);

        for (size_t i = 0; i < vec_size_total; i++)
        {
          sysmodel->yp()[i].setVariableNumber(i);
        }

        sysmodel->initialize();
        sysmodel->evaluateResidual();

        std::vector<DependencyTracking::Variable>& residuals = sysmodel->getResidual();

        size_t* row_indices = sysmodel->getCsrJac().getRowData();
        size_t* col_indices = sysmodel->getCsrJac().getColData();

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

        std::string test_name = __func__ + std::format(" (Nsize = {})", Nsize);
        return success.report(test_name.c_str());
      }
    };
  } // namespace Testing
} // namespace GridKit
