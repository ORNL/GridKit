#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <filesystem>

#include <Solver/Dynamic/Ida.hpp>
#include <Solver/Dynamic/DynamicSolver.hpp>


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
    double delta = Y[0];
    double w = Y[1];
    double E_qp = Y[2];
    double psi_dp = Y[3];
    double psi_qp = Y[4];
    double E_dp = Y[5];
    double I_d = Z[0];
    double I_q = Z[1];

    double psi_qpp;
    double psi_dpp;
    double psi_pp;
    double k_sat;
    double T_elec;

    get_psi_qpp(&psi_qpp, psi_qp, P.X_q4, E_dp, P.X_q5);
    get_psi_dpp(&psi_dpp, psi_dp, P.X_d4, E_qp, P.X_d5);
    psi_pp = sqrt(psi_qpp*psi_qpp + psi_dpp*psi_dpp);
    k_sat = P.S_B*(psi_pp - P.S_A)*(psi_pp - P.S_A);
    T_elec = (psi_dpp - I_d*P.X_dpp)*I_q - (psi_qpp - I_q*P.X_dpp)*I_d;
    
    Yp[0] = w*P.w_0;
    Yp[1] = (1/(2*P.H))*((P_mech - P.D*w)/(1+w) - T_elec);
    Yp[2] = (1/P.T_dop)*(E_fd  - (E_qp + P.X_d1*(I_d + P.X_d3*(E_qp - psi_dp - P.X_d2*I_d)) + psi_dpp*k_sat));
    Yp[3] = (1/P.T_dopp)*(E_qp - psi_dp - P.X_d2*I_d);
    Yp[4] = (1/P.T_qopp)*(E_dp - psi_qp + X_q2*I_q);
    Yp[5] = (1/P.T_qop)*(-E_dp + X_qd*psi_qpp*k_sat + P.X_q1*(I_q - P.X_q3*(E_dp + I_q*P.X_q2 - psi_qp)))

    return 0;

}

int get_algebraic_res(sdt::vector<double> &res, std::vector<double> Y, std::vector<double> Z, std::vector<double> V, std::vector<double> I, double R_a, double X_qpp)
{
    double delta = Y[0];
    double w = Y[1];
    double E_qp = Y[2];
    double psi_dp = Y[3];
    double psi_qp = Y[4];
    double E_dp = Y[5];

    double I_d = Z[0];
    double I_q = Z[1];

    double I_r = I[0];
    double I_i = I[1];

    double V_r = V[0];
    double V_i = V[1];

    double psi_qpp;
    double psi_dpp;

    get_psi_qpp(&psi_qpp, psi_qp, P.X_q4, E_dp, P.X_q5);
    get_psi_dpp(&psi_dpp, psi_dp, P.X_d4, E_qp, P.X_d5);

    res[0] = I_d - I_r*sin(delta) + I_i*cos(delta);
    res[1] = I_q - I_r*cos(delta) - I_i*sin(delta);
    res[2] = -psi_qpp*(1 + w) - V_r*sin(delta) + V_i*cos(delta) - I_d*R_a + I_q*X_qpp;
    res[3] =  psi_dpp*(1 + w) - V_r*cos(delta) - V_i*sin(delta) - I_d*X_qpp - I_q*R_a;
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
int brancheq(std::vector<double> &res, std::vector<double> &I, std::vector<double> &V, double g, double b, double B, double G)
{
    res[0] = -I[0] -(g + G/2)*V[0] + (b + B/2)*V[1] + g*V[2] - b*V[3];
    res[1] = -I[1] -(b + B/2)*V[0] - (g + G/2)*V[1] + b*V[2] + g*V[3];
    res[2] = -I[2] = g*V[0] - b*V[1] - (g + G/2)*V[2] + (b + B/2)*V[3];
    res[3] = -I[3] = b*V[0] + g*V[1] - (b + B/2)*V[2] - (g + G/2)*V[3];

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
    res[1] = I[1] - (V[1]*R1 + V[0]*R2)/(R1*R1 + R2*R2);

    return 0;
}

int loadeq_res()

std::vector<double> diff_res(double t, std::vector<double> &Y, std::vector<double> &Yp)
{
    std::vector<double> gen_res;

    std::vector<double> sub_Y = (Y.begin(), Y.begin() + 6);
    std::vector<double> sub_Yp = (Yp.begin(), Yp.begin() + 6);
    std::vector<double> Z = (Y.begin() + 32, Yp.begin() + 34);
    std::vector<double> diff_Yp(6);
    generator_params params;
    diff_gen(sub_Y, &diff_Yp, Z, params, 1, 1);

    for(int i=0; i < sub_Yp.size(); i++)
    {
        gen_res.push_back(sub_Yp[i] - Yp[i] );
    }

    std::vector<double> algeb_res(4);
    std::vector<double> V = {Y[34], Y[35]};
    std::vector<double> I = {Y[6], Y[7]};
    get_algebraic_res(&algeb_res, sub_Y, Z, V, I, 1, 1);

    for(int i=0; i < algeb_res.size(); i++)
    {
        gen_res.push_back(algeb_res[i]);
    }

    //Calculate residuals from bus equations and
    //append to the global residual vector gen_res
    int buseq(std::vector<double> &I, std::vector<double> &r)

    std::vector<double> bus_res(2);

    std::vector<double> I_bus_eqn = (Y.begin() + 6, Y.begin() + 14);
    buseq(I_bus_eqn, &bus_res);
    gen_res.insert(gen_res.end(), bus_res.begin(), bus_res.end());

    I_bus_eqn = (Y.begin() + 14, Y.begin() + 20);
    buseq(I_bus_eqn, &bus_res);
    gen_res.insert(gen_res.end(), bus_res.begin(), bus_res.end());

    I_bus_eqn = (Y.begin() + 20, Y.begin() + 26);
    buseq(I_bus_eqn, &bus_res);
    gen_res.insert(gen_res.end(), bus_res.begin(), bus_res.end());

    //Calculate residuals from branch equations and
    //append to the global residual vector gen_res
    std::vector<double> branch_res(4);

    std::vector<double> I_branch = {Y[8], Y[9], Y[16], Y[17]};
    std::vector<double> V_branch = {Y[26], Y[27], Y[28], Y[29]};
    brancheq(&branch_res, I_branch, V_branch, 1, 1, 1, 1);
    gen_res.insert(gen_res.end(), branch_res.begin(), branch_res.end());
    
    I_branch = {Y[18], Y[19], Y[24], Y[25]};
    V_branch = {Y[28], Y[29], Y[30], Y[31]};
    brancheq(&branch_res, I_branch, V_branch, 1, 1, 1, 1);
    gen_res.insert(gen_res.end(), branch_res.begin(), branch_res.end());

    
    I_branch = {Y[10], Y[11], Y[22], Y[23]};
    V_branch = {Y[26], Y[27], Y[30], Y[31]};
    brancheq(&branch_res, I_branch, V_branch, 1, 1, 1, 1);
    gen_res.insert(gen_res.end(), branch_res.begin(), branch_res.end());

    //Calculate residuals from load equations and
    //append to the global residual vector gen_res
    int loadeq(std::vector<double> &I, std::vector<double> &V, std::vector<double> &res, double R1, double R2)

    std::vector<double> load_res(2);
    std::vector<double> I_load = {Y[12], Y[13]};
    std::vector<double> V_load = {Y[26], Y[27]};
    loadeq(I_load, V_load, load_res, 1, 1);
    gen_res.insert(gen_res.end(), load_res.begin(), load_res.end());

    return gen_res;

}

int get_params()


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

    return 0;
}
