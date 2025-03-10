#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <filesystem>
#include<vector>
using namespace std;
//#include <Solver/Dynamic/Ida.hpp>
//#include <Solver/Dynamic/DynamicSolver.hpp>

#define _DEL1    0
#define _W1      1
#define _E_qp1   2
#define _PSI_dp1 3
#define _PSI_qp1 4
#define _E_dp1   5
#define _I_Gr1   6
#define _I_Gi1   7
#define _I_r21   8
#define _I_i21   9
#define _I_r31   10
#define _I_i31   11
#define _I_Lr1   12
#define _I_Li1   13
#define _I_Lr2   14
#define _I_Li2   15
#define _I_r12   16
#define _I_i12   17
#define _I_r32   18
#define _I_i32   19
#define _I_Gr3   20
#define _I_Gi3   21
#define _I_r13   22
#define _I_i13   23
#define _I_r23   24
#define _I_i23   25
#define _V_r1    26
#define _V_i1    27
#define _V_r2    28
#define _V_i2    29
#define _V_r3    30
#define _V_i3    31
#define _I_d1    32
#define _I_q1    33
#define _DEL2    34
#define _W2      35
#define _E_qp2   36
#define _PSI_dp2 37
#define _PSI_qp2 38
#define _E_dp2   39
#define _I_d2    40
#define _I_q2    41


struct generator_params 
{
    double w_0;
    double H;
    double D;
    double T_dop;
    double T_dopp;
    double T_qopp;
    double T_qop;
    double X_d1;
    double X_d2;
    double X_d3;
    double X_d4;
    double X_d5;
    double X_q1;
    double X_q2;
    double X_q3;
    double X_q4;
    double X_q5;
    double X_d;
    double X_dp;
    double X_dpp;
    double X_q;
    double X_qp;
    double X_qpp;
    double X_L;
    double X_qd;
    double R_a;
    double S_A;
    double S_B;
};

int get_psi_qpp(double &psi_qpp,double psi_qp, double X_q4, double E_dp, double X_q5)
{
    psi_qpp = -psi_qp*X_q4 - E_dp*X_q5;
    return 0;
}

int get_psi_dpp(double &psi_dpp, double psi_dp, double X_d4, double E_qp, double X_d5)
{
    psi_dpp = psi_dp*X_d4 + E_qp*X_d5;
    return 0;
}


int diff_gen(std::vector<double> Y, std::vector<double> &Yp, std::vector<double> Z, generator_params P, double P_mech, double E_fd)
{   
    double delta = Y[_DEL1];
    double w = Y[_W1];
    double E_qp = Y[_E_qp1];
    double psi_dp = Y[_PSI_dp1];
    double psi_qp = Y[_PSI_qp1];
    double E_dp = Y[_E_dp1];

    double I_d = Z[0];
    double I_q = Z[1];

    double psi_qpp;
    double psi_dpp;
    double psi_pp;
    double k_sat;
    double T_elec;

    get_psi_qpp(psi_qpp, psi_qp, P.X_q4, E_dp, P.X_q5);
    get_psi_dpp(psi_dpp, psi_dp, P.X_d4, E_qp, P.X_d5);
    psi_pp = sqrt(psi_qpp*psi_qpp + psi_dpp*psi_dpp);
    k_sat = P.S_B*(psi_pp - P.S_A)*(psi_pp - P.S_A);
    T_elec = (psi_dpp - I_d*P.X_dpp)*I_q - (psi_qpp - I_q*P.X_dpp)*I_d;

    Yp[0] = w*P.w_0;
    Yp[1] = (1/(2*P.H))*((P_mech - P.D*w)/(1+w) - T_elec);
    Yp[2] = (1/P.T_dop)*(E_fd  - (E_qp + P.X_d1*(I_d + P.X_d3*(E_qp - psi_dp - P.X_d2*I_d)) + psi_dpp*k_sat));
    Yp[3] = (1/P.T_dopp)*(E_qp - psi_dp - P.X_d2*I_d);
    Yp[4] = (1/P.T_qopp)*(E_dp - psi_qp + P.X_q2*I_q);
    Yp[5] = (1/P.T_qop)*(-E_dp + P.X_qd*psi_qpp*k_sat + P.X_q1*(I_q - P.X_q3*(E_dp + I_q*P.X_q2 - psi_qp)));

    return 0;

}

int get_algebraic_res(std::vector<double> &res, std::vector<double> Y, std::vector<double> Z, std::vector<double> V, std::vector<double> I, generator_params P, double R_a)
{
    double delta = Y[_DEL1];
    double w = Y[_W1];
    double E_qp = Y[_E_qp1];
    double psi_dp = Y[_PSI_dp1];
    double psi_qp = Y[_PSI_qp1];
    double E_dp = Y[_E_dp1];

    double I_d = Z[0];
    double I_q = Z[1];

    double I_r = I[0];
    double I_i = I[1];

    double V_r = V[0];
    double V_i = V[1];

    double psi_qpp;
    double psi_dpp;

    get_psi_qpp(psi_qpp, psi_qp, P.X_q4, E_dp, P.X_q5);
    get_psi_dpp(psi_dpp, psi_dp, P.X_d4, E_qp, P.X_d5);

    res[0] = I_d - I_r*sin(delta) + I_i*cos(delta);
    res[1] = I_q - I_r*cos(delta) - I_i*sin(delta);
    res[2] = -psi_qpp*(1 + w) - V_r*sin(delta) + V_i*cos(delta) - I_d*R_a + I_q*P.X_qpp;
    res[3] =  psi_dpp*(1 + w) - V_r*cos(delta) - V_i*sin(delta) - I_d*P.X_qpp - I_q*R_a;

    return 0;
}


int buseq(std::vector<double> &I, std::vector<double> &r)
{
    r[0] = 0.0;
    r[1] = 0.0;

    for (size_t i = 0; i < I.size(); i++)
    {
        if(i%2 == 0){
            r[0] += I[i];
        }
        else{
            r[1] += I[i];
        }
    }

    return 0;
    
}

/**
 * @brief Calculates the current given the voltages
 * I[0]: I_r1
 * I[1]: I_i1
 * I[2]: I_r2
 * I[3]: I_i2
 * 
 * V[0]: V_r1
 * V[1]: V_i1
 * V[2]: V_r2
 * V[3]: V_i2
 */
int brancheq(std::vector<double> &res, std::vector<double> &I, std::vector<double> &V, double g, double b, double G, double B)
{
    res[0] = -I[0] -(g + G/2)*V[0] + (b + B/2)*V[1] + g*V[2] - b*V[3];
    res[1] = -I[1] -(b + B/2)*V[0] - (g + G/2)*V[1] + b*V[2] + g*V[3];
    res[2] = -I[2] + g*V[0] - b*V[1] - (g + G/2)*V[2] + (b + B/2)*V[3];
    res[3] = -I[3] + b*V[0] + g*V[1] - (b + B/2)*V[2] - (g + G/2)*V[3];

    return 0;
}

/**
 * @brief calculates the current for the load given the voltages and Resistances 
 * I[0]: I_r1
 * I[1]: I_i1
 * 
 * V[0]: V_r1
 * V[1]: V_i1
 */
int loadeq(std::vector<double> &I, std::vector<double> &V, std::vector<double> &res, double R1, double R2)
{
    res[0] = I[0] - (V[0]*R1 + V[1]*R2)/(R1*R1 + R2*R2);
    res[1] = I[1] - (V[1]*R1 - V[0]*R2)/(R1*R1 + R2*R2);

    return 0;
}

int copy_to_res_vector(std::vector<double> &gen_res, std::vector<double> &towrite, int &at)
{
    for(int i=0; i < towrite.size(); i++){
        gen_res[at] = towrite[i];
        at++;
    }
    return 0;
}


int diff_res(double t, std::vector<double> &Y, std::vector<double> &Yp, generator_params params, std::vector<double> &gen_res )
{

    int residual_index = 0;

    std::vector<double> sub_Y(Y.begin(), Y.begin() + 6);
    std::vector<double> sub_Yp(Yp.begin(), Yp.begin() + 6);
    std::vector<double> Z(2);
    std::vector<double> diff_Yp(6);
    Z[0] = Y[_I_d1];
    Z[1] = Y[_I_q1]; 

    
    double P_mech = 3.0;
    double E_fd = 1.0;
    
    diff_gen(sub_Y, diff_Yp, Z, params, P_mech, E_fd);

    for(int i=0; i < diff_Yp.size(); i++)
    {
        gen_res[i] = diff_Yp[i] - Yp[i];
        residual_index++;
    }


    std::vector<double> algeb_res(4);
    std::vector<double> V = {Y[_V_r1], Y[_V_i1]};
    std::vector<double> I = {Y[_I_Gr1], Y[_I_Gi1]};
    double R_a = 2.0;

    get_algebraic_res(algeb_res, sub_Y, Z, V, I, params, R_a);
    copy_to_res_vector(gen_res, algeb_res, residual_index);

    

    //Calculate residuals from bus equations and
    //append to the global residual vector gen_res

    std::vector<double> bus_res(2);

    std::vector<double> I_bus_eqn(Y.begin() + 6, Y.begin() + 14);
    buseq(I_bus_eqn, bus_res);
    copy_to_res_vector(gen_res, bus_res, residual_index);


    I_bus_eqn.assign(Y.begin() + 14, Y.begin() + 20);
    buseq(I_bus_eqn, bus_res);
    copy_to_res_vector(gen_res, bus_res, residual_index);


    I_bus_eqn.assign(Y.begin() + 20, Y.begin() + 26);
    buseq(I_bus_eqn, bus_res);
    copy_to_res_vector(gen_res, bus_res, residual_index);

    

    //Calculate residuals from branch equations and
    //append to the global residual vector gen_res
    std::vector<double> branch_res(4);
    double b = 1.0;
    double g = 1.0;
    double B = 2.0;
    double G = 2.0;

    std::vector<double> I_branch = {Y[_I_r21], Y[_I_i21], Y[_I_r12], Y[_I_i12]};
    std::vector<double> V_branch(Y.begin()+26, Y.begin() + 30);
    brancheq(branch_res, I_branch, V_branch, g, b, G, B);
    copy_to_res_vector(gen_res, branch_res, residual_index);

    
    I_branch = {Y[_I_r32], Y[_I_i32], Y[_I_r23], Y[_I_i23]};
    V_branch.assign(Y.begin() + 28, Y.begin() + 32);
    brancheq(branch_res, I_branch, V_branch, g, b, G, B);
    copy_to_res_vector(gen_res, branch_res, residual_index);

    
    I_branch = {Y[_I_r31], Y[_I_i31], Y[_I_r13], Y[_I_i13]};
    V_branch = {Y[_V_r1], Y[_V_i1], Y[_V_r3], Y[_V_i3]};
    brancheq(branch_res, I_branch, V_branch, g, b, G, B);
    copy_to_res_vector(gen_res, branch_res, residual_index);

    
    //Calculate residuals from load equations and
    //append to the global residual vector gen_res
    std::vector<double> load_res(2);
    double R1 = 1.0;
    double R2 = 1.0;

    std::vector<double> I_load = {Y[_I_Lr1], Y[_I_Li1]};
    std::vector<double> V_load = {Y[_V_r1], Y[_V_i1]};
    loadeq(I_load, V_load, load_res, R1, R2);
    copy_to_res_vector(gen_res, load_res, residual_index);


    I_load = {Y[_I_Lr2], Y[_I_Li2]};
    V_load = {Y[_V_r2], Y[_V_i2]};
    loadeq(I_load, V_load, load_res, R1, R2);
    copy_to_res_vector(gen_res, load_res, residual_index);


    sub_Y.assign(Y.begin() + 34, Y.begin() + 40);
    sub_Yp.assign(Yp.begin() + 34, Yp.begin() + 40);
    Z[0] = Y[_I_d2];
    Z[1] = Y[_I_q2]; 

    
    P_mech = 3.0;
    E_fd = 1.0;
    
    diff_gen(sub_Y, diff_Yp, Z, params, P_mech, E_fd);

    for(int i=0; i < diff_Yp.size(); i++)
    {
        gen_res[i] = diff_Yp[i] - Yp[i];
        residual_index++;
    }


    V = {Y[_V_r3], Y[_V_i3]};
    I = {Y[_I_Gr3], Y[_I_Gi3]};
    R_a = 2.0;

    get_algebraic_res(algeb_res, sub_Y, Z, V, I, params, R_a);
    copy_to_res_vector(gen_res, algeb_res, residual_index);



    for(int i=0; i < gen_res.size(); i++)
    {
        std::cout << "index:  " << i << " - value: " << gen_res[i] << std::endl;
    }

    return 0;

}




int main(int argc, char const *argv[])
{
    //These are for branch 1
    double b1 = 0.0;
    double g1 = 0.0;
    double B1 = 0.0;
    double G1 = 0.0;

    //These are for branch 2
    double b2 = 0.0;
    double g2 = 0.0;
    double B2 = 0.0;
    double G2 = 0.0;

    //These are for branch 3
    double b3 = 0.0;
    double g3 = 0.0;
    double B3 = 0.0;
    double G3 = 0.0;

    generator_params params;

    params.w_0 = 1;
    params.H = 1;
    params.D = 1;
    params.T_dop = 1;
    params.T_dopp = 1;
    params.T_qopp = 1;
    params.T_qop = 1;
    params.X_d1 = 1;
    params.X_d2 = 1;
    params.X_d3 = 1;
    params.X_d4 = 1;
    params.X_d5 = 1;
    params.X_q1 = 1;
    params.X_q2 = 1;
    params.X_q3 = 1;
    params.X_q4 = 1;
    params.X_q5 = 1;
    params.X_d = 1;
    params.X_dp = 1;
    params.X_dpp = 1;
    params.X_q = 1;
    params.X_qp = 1;
    params.X_qpp = 1;
    params.X_L = 1;
    params.X_qd = 1;
    params.R_a = 1;
    params.S_A = 1;
    params.S_B = 1;

    std::vector<double> gen_res(42);
    std::vector<double> Y = {1, 1, 2, 2, 2, 1, 1, 2, -6, -5, 4, 1, -2,
         -2, 0.5, 4.5, -1, 15, 0, -18, 1, 1, 2, 2, -3, -2, 2, 1, 1, 2, 
        3, 4, 1 , 2, 1, 1, 2, 2, 2, 1, 1, 2};

    std::vector<double> Yp(42);
    Yp[0]=1, Yp[1]=-5, Yp[2]=-65, Yp[3]=-1, Yp[4]=1, Yp[5]=-48;
    Yp[34]=1, Yp[35]=-5, Yp[36]=-65, Yp[37]=-1, Yp[38]=1, Yp[39]=-48;

    diff_res(1,  Y, Yp, params, gen_res );


    // vector<double> Y = {1, 1, 2, 2, 2, 1};
    // vector<double> Z = {1, 2};
    // vector<double> Yp(6);
    // double P_mech = 3.0;
    // double E_fd = 1.0;

    // diff_gen(Y, Yp, Z, params, P_mech, E_fd);
    

    // for(auto val: Yp)
    // {
    //     std::cout << val << " ";
    // }

    // std::cout << "\n--------------------Algebraics----------------------" << std::endl;

    // std::vector<double> algeb_res(4);
    // std::vector<double> V = {2, 1};
    // std::vector<double> I = {1, 2};
    // double R_a = 2;
    // get_algebraic_res(algeb_res, Y, Z, V, I, params, R_a);

    // for(auto val: algeb_res)
    // {
    //     std::cout << val << " ";
    // }

    // std::cout << "\n--------------------Buses----------------------" << std::endl;

    // std::vector<double> bus_res(2);
    // I = {3, 4, -6, -5, 4, 1, -2, -2};
    // buseq(I, bus_res);
    // for(auto val: bus_res)
    // {
    //     std::cout << val << " ";
    // }
    // std::cout <<std::endl;

    // I = {0.5, 4.5, -1, 15, 0, -18};
    // buseq(I, bus_res);
    // for(auto val: bus_res)
    // {
    //     std::cout << val << " ";
    // }

    // std::cout <<std::endl;

    // I = {1, 1, 2, 2, -3, -2};
    // buseq(I, bus_res);
    // for(auto val: bus_res)
    // {
    //     std::cout << val << " ";
    // }

    // std::cout << "\n--------------------Branches----------------------" << std::endl;
    // double g, b, B, G;
    // g = 1, b=1, G=2, B=2;

    // std::vector<double> branch_res(4);
    // I = {-6, -5, -1, 15};
    // V = {2, 1, 1, 2};
    // brancheq(branch_res, I, V, g, b, B, G);

    // for(auto val: branch_res)
    // {
    //     std::cout << val << " ";
    // }
    // std::cout <<std::endl;

    // I = {0, -18, -2, -2};
    // V = {1, 2, 3, 4};
    // brancheq(branch_res, I, V, g, b, B, G);
    
    // for(auto val: branch_res)
    // {
    //     std::cout << val << " ";
    // }
    // std::cout <<std::endl;

    // I = {4, 1, 2, 2};
    // V = {2, 1, 3, 4};
    // brancheq(branch_res, I, V, g, b, B, G);
    
    // for(auto val: branch_res)
    // {
    //     std::cout << val << " ";
    // }
    // std::cout <<std::endl;

    // std::cout << "\n--------------------Branches----------------------" << std::endl;
    // //int loadeq(std::vector<double> &I, std::vector<double> &V, std::vector<double> &res, double R1, double R2)

    // std::vector<double> load_res(2);
    // I = {-2, -2};
    // V = {2, 1};
    // double R1 = 1;
    // double R2 = 1;

    // loadeq(I, V, load_res, R1, R2);

    // for(auto val: load_res)
    // {
    //     std::cout << val << " ";
    // }
    // std::cout <<std::endl;


    // I = {0.5, 4.5};
    // V = {1, 2};

    // loadeq(I, V, load_res, R1, R2);

    // for(auto val: load_res)
    // {
    //     std::cout << val << " ";
    // }
    // std::cout <<std::endl;


    
    return 0;
}
