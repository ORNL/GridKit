

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <filesystem>

#include <Utilities/Testing.hpp>
#include "SolutionKeys.hpp"

#include <ComponentLib/PowerElectronicsComponents/DistributedGenerator/DistributedGenerator.hpp>
#include <ComponentLib/PowerElectronicsComponents/MicrogridLine/MicrogridLine.hpp>
#include <ComponentLib/PowerElectronicsComponents/MicrogridLoad/MicrogridLoad.hpp>
#include <ComponentLib/PowerElectronicsComponents/MicrogridBusDQ/MicrogridBusDQ.hpp>

#include <PowerElectronicsModel.hpp>
#include <Solver/Dynamic/Ida.hpp>
#include <Solver/Dynamic/DynamicSolver.hpp>


static int test(size_t Nsize);

int main(int argc, char const *argv[])
{
    int retval = 0;
    retval += test(2);
    retval += test(4);
    retval += test(8);
    if (retval > 0)
    {
        std::cout << "Some tests fail!!\n";
    }
    else
    {
        std::cout << "All tests pass!!\n";
    }
    return retval;
}

/**
 * @brief Tests network of distributed generators.
 * 
 * @param Nsize - The number of DG line load cobinations to generate for scale
 * @return int returns 0 if successful, >0 otherwise
 */
int test(size_t Nsize)
{
    double abstol = 1.0e-8;
    double reltol = 1.0e-8;
    size_t max_step_amount = 3000;
    bool usejac = true;

    //TODO:setup as named parameters
    //Create circuit model
    ModelLib::PowerElectronicsModel<double, size_t>* sysmodel = new ModelLib::PowerElectronicsModel<double, size_t>(reltol, abstol, usejac, max_step_amount);

    const std::vector<double>* true_vec = &answer_key_N8;

    switch (Nsize)
    {
        case 2:
            true_vec = &answer_key_N2;
            break;
        case 4:
            true_vec = &answer_key_N4;
            break;
        case 8:
            // true_vec = &answer_key_N8;
            break;
        default:
            std::cout << "No reference solution for Nsize = " << Nsize << ".\n";
            std::cout << "Using default Nsize = 8.\n";
    }

    
    // Modeled after the problem in the paper
    // Every Bus has the same virtual resistance. This is due to the numerical stability as mentioned in the paper.
    double RN = 1.0e4;

    // DG Params Vector
    // All DGs have the same set of parameters except for the first two.
    ModelLib::DistributedGeneratorParameters<double, size_t> DG_parms1;
    DG_parms1.wb_ = 2.0*M_PI*50.0;
    DG_parms1.wc_ = 31.41;
    DG_parms1.mp_ = 9.4e-5;
    DG_parms1.Vn_ = 380.0;
    DG_parms1.nq_ = 1.3e-3;
    DG_parms1.F_ = 0.75;
    DG_parms1.Kiv_ = 420.0;
    DG_parms1.Kpv_ = 0.1;
    DG_parms1.Kic_ = 2.0e4;
    DG_parms1.Kpc_ = 15.0;
    DG_parms1.Cf_ = 5.0e-5;
    DG_parms1.rLf_ = 0.1;
    DG_parms1.Lf_ = 1.35e-3;
    DG_parms1.rLc_ = 0.03;
    DG_parms1.Lc_ = 0.35e-3;

    ModelLib::DistributedGeneratorParameters<double, size_t> DG_parms2;
    DG_parms2.wb_ = 2.0*M_PI*50.0;
    DG_parms2.wc_ = 31.41;
    DG_parms2.mp_ = 12.5e-5;
    DG_parms2.Vn_ = 380.0;
    DG_parms2.nq_ = 1.5e-3;
    DG_parms2.F_ = 0.75;
    DG_parms2.Kiv_ = 390.0;
    DG_parms2.Kpv_ = 0.05;
    DG_parms2.Kic_ = 16.0e3;
    DG_parms2.Kpc_ = 10.5;
    DG_parms2.Cf_ = 50.0e-6;
    DG_parms2.rLf_ = 0.1;
    DG_parms2.Lf_ = 1.35e-3;
    DG_parms2.rLc_ = 0.03;
    DG_parms2.Lc_ = 0.35e-3;

    std::vector<ModelLib::DistributedGeneratorParameters<double, size_t>> DGParams_list(2*Nsize, DG_parms2);

    DGParams_list[0] = DG_parms1;
    DGParams_list[1] = DG_parms1;

    //line vector params
    //Every odd line has the same parameters and every even line has the same parameters
    double rline1 = 0.23;
    double Lline1 = 0.1 / (2.0 * M_PI * 50.0);
    double rline2 = 0.35;
    double Lline2 = 0.58 / (2.0 * M_PI * 50.0);
    std::vector<double> rline_list(2*Nsize-1, 0.0);
    std::vector<double> Lline_list(2*Nsize-1, 0.0);
    for (size_t i = 0; i < rline_list.size(); i++)
    {
        rline_list[i] = (i % 2) ? rline2 : rline1;
        Lline_list[i] = (i % 2) ? Lline2 : Lline1;
    }
    

    //load parms
    //Only the first load has the same paramaters.
    double rload1 = 3.0;
    double Lload1 = 2.0 / (2.0 * M_PI * 50.0);
    double rload2 = 2.0;
    double Lload2 = 1.0 / (2.0 * M_PI * 50.0);

    std::vector<double> rload_list(Nsize, rload2);
    std::vector<double> Lload_list(Nsize, Lload2);
    rload_list[0] = rload1;
    Lload_list[0] = Lload1;

    //							DGs	+		- refframe	   Lines +				Loads
    size_t vec_size_internals = 13*(2*Nsize) - 1 + (2 + 4*(Nsize - 1)) + 2*Nsize;
    //							\omegaref + BusDQ
    size_t vec_size_externals = 1 +	2*(2*Nsize);

    std::vector<size_t> vdqbus_index(2*Nsize,0);
    vdqbus_index[0] = vec_size_internals + 1;
    for (size_t i = 1; i < vdqbus_index.size(); i++)
    {
        vdqbus_index[i] = vdqbus_index[i-1] + 2;
    }

    //Total size of the vector setup
    size_t vec_size_total = vec_size_internals + vec_size_externals;


    //Create the reference DG
    ModelLib::DistributedGenerator<double, size_t>* dg_ref = new ModelLib::DistributedGenerator<double, size_t>(0, DGParams_list[0], true);
    //ref motor
    dg_ref->setExternalConnectionNodes(0, vec_size_internals);
    //outputs
    dg_ref->setExternalConnectionNodes(1, vdqbus_index[0]);
    dg_ref->setExternalConnectionNodes(2, vdqbus_index[0] + 1);
    //"grounding" of the difference
    dg_ref->setExternalConnectionNodes(3, -1);
    //internal connections
    for (size_t i = 0; i < 12; i++)
    {
        dg_ref->setExternalConnectionNodes(4 + i, i);
    }
    sysmodel->addComponent(dg_ref);

    //Keep track of models and index location
    size_t indexv = 12;
    size_t model_id = 1;
    //Add all other DGs
    for (size_t i = 1; i < 2*Nsize; i++)
    {
        //current DG to add
        ModelLib::DistributedGenerator<double, size_t>* dg = new ModelLib::DistributedGenerator<double, size_t>(model_id++, DGParams_list[i], false);
        //ref motor
        dg->setExternalConnectionNodes(0,vec_size_internals);
        //outputs
        dg->setExternalConnectionNodes(1,vdqbus_index[i]);
        dg->setExternalConnectionNodes(2,vdqbus_index[i] + 1);
        //internal connections
        for (size_t j = 0; j < 13; j++)
        {
            dg->setExternalConnectionNodes(3 + j,indexv + j);
        }
        indexv += 13;
        sysmodel->addComponent(dg);
    }

    // Load all the Line compoenents
    for (size_t i = 0; i < 2*Nsize - 1; i++)
    {
        //line
        ModelLib::MicrogridLine<double, size_t>* line_model = new ModelLib::MicrogridLine<double, size_t>(model_id++, rline_list[i], Lline_list[i]);
        //ref motor
        line_model->setExternalConnectionNodes(0, vec_size_internals);
        //input connections
        line_model->setExternalConnectionNodes(1, vdqbus_index[i]);
        line_model->setExternalConnectionNodes(2, vdqbus_index[i] + 1);
        //output connections
        line_model->setExternalConnectionNodes(3, vdqbus_index[i+1]);
        line_model->setExternalConnectionNodes(4, vdqbus_index[i+1] + 1);
        //internal connections
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
        ModelLib::MicrogridLoad<double, size_t>* load_model = new ModelLib::MicrogridLoad<double, size_t>(model_id++, rload_list[i], Lload_list[i]);
        //ref motor
        load_model->setExternalConnectionNodes(0, vec_size_internals);
        //input connections
        load_model->setExternalConnectionNodes(1, vdqbus_index[2*i]);
        load_model->setExternalConnectionNodes(2, vdqbus_index[2*i] + 1);
        //internal connections
        for (size_t j = 0; j < 2; j++)
        {
            load_model->setExternalConnectionNodes(3 + j, indexv + j);
        }
        indexv += 2;
        sysmodel->addComponent(load_model);
    }

    //Add all the microgrid Virtual DQ Buses
    for (size_t i = 0; i < 2*Nsize; i++)
    {
        ModelLib::MicrogridBusDQ<double, size_t>* virDQbus_model = new ModelLib::MicrogridBusDQ<double, size_t>(
        model_id++, RN);

        virDQbus_model->setExternalConnectionNodes(0, vdqbus_index[i]);
        virDQbus_model->setExternalConnectionNodes(1, vdqbus_index[i] + 1);
        sysmodel->addComponent(virDQbus_model);
    }
    
    //allocate all the intial conditions
    sysmodel->allocate(vec_size_total);

    std::cout << sysmodel->y().size() << std::endl;
    std::cout << vec_size_internals << ", " << vec_size_externals << "\n";

    //Create Intial points for states. Every state is to specified to the zero intially
    for (size_t i = 0; i < vec_size_total; i++)
    {
        sysmodel->y()[i] = 0.0;
        sysmodel->yp()[i] = 0.0;
    }

    // Create Intial derivatives specifics generated in MATLAB
    for (size_t i = 0; i < 2*Nsize; i++)
    {
        sysmodel->yp()[13*i - 1 + 3] = DGParams_list[i].Vn_;
        sysmodel->yp()[13*i - 1 + 5] = DGParams_list[i].Kpv_ * DGParams_list[i].Vn_;
        sysmodel->yp()[13*i - 1 + 7] = (DGParams_list[i].Kpc_ * DGParams_list[i].Kpv_ * DGParams_list[i].Vn_) / DGParams_list[i].Lf_;
    }

    //since the intial P_com = 0, the set the intial vector to the reference frame
    sysmodel->y()[vec_size_internals] = DG_parms1.wb_;

    sysmodel->initialize();
    sysmodel->evaluateResidual();

    std::vector<double>& fres = sysmodel->getResidual(); 
    std::cout << "Verify Intial Resisdual is Zero: {\n";
    for (size_t i = 0; i < fres.size(); i++)
    {
        printf("%lu : %.16e \n", i, fres[i]);
    }
    std::cout << "}\n";

    sysmodel->updateTime(0.0, 1.0e-8);
    sysmodel->evaluateJacobian();
    std::cout << "Intial Jacobian with alpha:\n";


    //Create numerical integrator and configure it for the generator model
    AnalysisManager::Sundials::Ida<double, size_t>* idas = new AnalysisManager::Sundials::Ida<double, size_t>(sysmodel);

    double t_init  = 0.0;
    double t_final = 1.0;

    // setup simulation
    idas->configureSimulation();
    idas->getDefaultInitialCondition();
    idas->initializeSimulation(t_init);

    idas->runSimulation(t_final);

    std::vector<double>& yfinial = sysmodel->y(); 

    std::cout << "Final Vector y\n";
    for (size_t i = 0; i < yfinial.size(); i++)
    {
        printf("%lu: % 2.16e\n", i, yfinial[i]);
    }

    bool test_pass = true;
    double tol = 3e-3;

    // check relative error
    std::cout << "Test the Relative Error\n";
    for (size_t i = 0; i < true_vec->size(); i++)
    {
        test_pass *= GridKit::Testing::isEqual(yfinial[i], true_vec->at(i), tol);
        printf("%lu: % 2.16e\n", i, abs(true_vec->at(i) - yfinial[i]) / abs(true_vec->at(i)));
    }
    
    delete idas;
    delete sysmodel; 

    if(test_pass)
    {
        std::cout << "Test with Nsize = " << Nsize << " passes!\n";
        return 0;
    }
    else
    {
        std::cout << "Test with Nsize = " << Nsize << " fails!\n";
        return 1;
    }
}
