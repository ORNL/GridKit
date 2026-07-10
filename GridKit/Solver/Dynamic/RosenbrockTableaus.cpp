#include "Rosenbrock.hpp"

namespace Integrator
{
  template <class ScalarT, typename IdxT>
  Rosenbrock<ScalarT, IdxT>::Tableau Rosenbrock<ScalarT, IdxT>::Tableau::linImplicitEuler()
  {
    constexpr size_t num_stages = 1;

    Tableau re = {
        .num_stages_ = num_stages,
        .gamma_      = 1.0,
        .alpha_sum_  = std::make_unique_for_overwrite<RealT[]>(num_stages),
        .gamma_sum_  = std::make_unique_for_overwrite<RealT[]>(num_stages),
        .m_          = std::make_unique_for_overwrite<RealT[]>(num_stages),
        .e_          = {},
        .A_          = std::make_unique_for_overwrite<RealT[]>(num_stages * num_stages),
        .C_          = std::make_unique_for_overwrite<RealT[]>(num_stages * num_stages),
        .H_          = {},
        .order_      = 1,
        .is_krylov_  = true,
        .is_w_       = false,
        .is_dae_     = true,
    };

    re.alpha_sum_[0] = 0.0;

    re.gamma_sum_[0] = 1.0;

    re.m_[0] = 1.0;

    re.A_[0] = 0.0;

    re.C_[0] = 2.0;

    return re;
  }

  /**
   * @brief A 5th order Rosenbrock method suitable for solving parabolic PDEs. Created in \cite steinebach2023construction.
   *
   *
   */
  template <class ScalarT, typename IdxT>
  Rosenbrock<ScalarT, IdxT>::Tableau Rosenbrock<ScalarT, IdxT>::Tableau::rodas5p()
  {
    constexpr size_t num_stages = 8;

    Tableau re = {
        .num_stages_ = num_stages,
        .gamma_      = 0.21193756319429014,
        .alpha_sum_  = std::make_unique_for_overwrite<RealT[]>(num_stages),
        .gamma_sum_  = std::make_unique_for_overwrite<RealT[]>(num_stages),
        .m_          = std::make_unique_for_overwrite<RealT[]>(num_stages),
        .e_          = std::make_unique_for_overwrite<RealT[]>(num_stages),
        .A_          = std::make_unique_for_overwrite<RealT[]>(num_stages * num_stages),
        .C_          = std::make_unique_for_overwrite<RealT[]>(num_stages * num_stages),
        .H_          = std::make_unique_for_overwrite<RealT[]>(num_stages * 3),
        .order_      = 5,
        .is_krylov_  = false,
        .is_w_       = false,
        .is_dae_     = true,
    };

    re.alpha_sum_[0] = 0.0;
    re.alpha_sum_[1] = 0.6358126895828704;
    re.alpha_sum_[2] = 0.4095798393397535;
    re.alpha_sum_[3] = 0.9769306725060716;
    re.alpha_sum_[4] = 0.4288403609558664;
    re.alpha_sum_[5] = 1.0;
    re.alpha_sum_[6] = 1.0;
    re.alpha_sum_[7] = 1.0;

    re.gamma_sum_[0] = 0.21193756319429014;
    re.gamma_sum_[1] = -0.42387512638858027;
    re.gamma_sum_[2] = -0.3384627126235924;
    re.gamma_sum_[3] = 1.8046452872882734;
    re.gamma_sum_[4] = 2.325825639765069;
    re.gamma_sum_[5] = 0.0;
    re.gamma_sum_[6] = 0.0;
    re.gamma_sum_[7] = 0.0;

    re.m_[0] = -7.502846399306121;
    re.m_[1] = 2.561846144803919;
    re.m_[2] = -11.627539656261098;
    re.m_[3] = -0.18268767659942256;
    re.m_[4] = 0.030198172008377946;
    re.m_[5] = 1.0;
    re.m_[6] = 1.0;
    re.m_[7] = 1.0;

    re.e_[0] = 0.0;
    re.e_[1] = 0.0;
    re.e_[2] = 0.0;
    re.e_[3] = 0.0;
    re.e_[4] = 0.0;
    re.e_[5] = 0.0;
    re.e_[6] = 0.0;
    re.e_[7] = 1.0;

    re.A_[1 * 8 + 0] = 3.0;

    re.A_[2 * 8 + 0] = 2.849394379747939;
    re.A_[2 * 8 + 1] = 0.45842242204463923;

    re.A_[3 * 8 + 0] = -6.954028509809101;
    re.A_[3 * 8 + 1] = 2.489845061869568;
    re.A_[3 * 8 + 2] = -10.358996098473584;

    re.A_[4 * 8 + 0] = 2.8029986275628964;
    re.A_[4 * 8 + 1] = 0.5072464736228206;
    re.A_[4 * 8 + 2] = -0.3988312541770524;
    re.A_[4 * 8 + 3] = -0.04721187230404641;

    re.A_[5 * 8 + 0] = -7.502846399306121;
    re.A_[5 * 8 + 1] = 2.561846144803919;
    re.A_[5 * 8 + 2] = -11.627539656261098;
    re.A_[5 * 8 + 3] = -0.18268767659942256;
    re.A_[5 * 8 + 4] = 0.030198172008377946;

    re.A_[6 * 8 + 0] = -7.502846399306121;
    re.A_[6 * 8 + 1] = 2.561846144803919;
    re.A_[6 * 8 + 2] = -11.627539656261098;
    re.A_[6 * 8 + 3] = -0.18268767659942256;
    re.A_[6 * 8 + 4] = 0.030198172008377946;
    re.A_[6 * 8 + 5] = 1.0;

    re.A_[7 * 8 + 0] = -7.502846399306121;
    re.A_[7 * 8 + 1] = 2.561846144803919;
    re.A_[7 * 8 + 2] = -11.627539656261098;
    re.A_[7 * 8 + 3] = -0.18268767659942256;
    re.A_[7 * 8 + 4] = 0.030198172008377946;
    re.A_[7 * 8 + 5] = 1.0;
    re.A_[7 * 8 + 6] = 1.0;

    re.C_[1 * 8 + 0] = -14.155112264123755;

    re.C_[2 * 8 + 0] = -17.97296035885952;
    re.C_[2 * 8 + 1] = -2.859693295451294;

    re.C_[3 * 8 + 0] = 147.12150275711716;
    re.C_[3 * 8 + 1] = -1.41221402718213;
    re.C_[3 * 8 + 2] = 71.68940251302358;

    re.C_[4 * 8 + 0] = 165.43517024871676;
    re.C_[4 * 8 + 1] = -0.4592823456491126;
    re.C_[4 * 8 + 2] = 42.90938336958603;
    re.C_[4 * 8 + 3] = -5.961986721573306;

    re.C_[5 * 8 + 0] = 24.854864614690072;
    re.C_[5 * 8 + 1] = -3.0009227002832186;
    re.C_[5 * 8 + 2] = 47.4931110020768;
    re.C_[5 * 8 + 3] = 5.5814197821558125;
    re.C_[5 * 8 + 4] = -0.6610691825249471;

    re.C_[6 * 8 + 0] = 30.91273214028599;
    re.C_[6 * 8 + 1] = -3.1208243349937974;
    re.C_[6 * 8 + 2] = 77.79954646070892;
    re.C_[6 * 8 + 3] = 34.28646028294783;
    re.C_[6 * 8 + 4] = -19.097331116725623;
    re.C_[6 * 8 + 5] = -28.087943162872662;

    re.C_[7 * 8 + 0] = 37.80277123390563;
    re.C_[7 * 8 + 1] = -3.2571969029072276;
    re.C_[7 * 8 + 2] = 112.26918849496327;
    re.C_[7 * 8 + 3] = 66.9347231244047;
    re.C_[7 * 8 + 4] = -40.06618937091002;
    re.C_[7 * 8 + 5] = -54.66780262877968;
    re.C_[7 * 8 + 6] = -9.48861652309627;

    re.H_[0 * 8 + 0] = 25.948786856663858;
    re.H_[0 * 8 + 1] = -2.5579724845846235;
    re.H_[0 * 8 + 2] = 10.433815404888879;
    re.H_[0 * 8 + 3] = -2.3679251022685204;
    re.H_[0 * 8 + 4] = 0.524948541321073;
    re.H_[0 * 8 + 5] = 1.1241088310450404;
    re.H_[0 * 8 + 6] = 0.4272876194431874;
    re.H_[0 * 8 + 7] = -0.17202221070155493;

    re.H_[1 * 8 + 0] = -9.91568850695171;
    re.H_[1 * 8 + 1] = -0.9689944594115154;
    re.H_[1 * 8 + 2] = 3.0438037242978453;
    re.H_[1 * 8 + 3] = -24.495224566215796;
    re.H_[1 * 8 + 4] = 20.176138334709044;
    re.H_[1 * 8 + 5] = 15.98066361424651;
    re.H_[1 * 8 + 6] = -6.789040303419874;
    re.H_[1 * 8 + 7] = -6.710236069923372;

    re.H_[2 * 8 + 0] = 11.419903575922262;
    re.H_[2 * 8 + 1] = 2.8879645146136994;
    re.H_[2 * 8 + 2] = 72.92137995996029;
    re.H_[2 * 8 + 3] = 80.12511834622643;
    re.H_[2 * 8 + 4] = -52.072871366152654;
    re.H_[2 * 8 + 5] = -59.78993625266729;
    re.H_[2 * 8 + 6] = -0.15582684282751913;
    re.H_[2 * 8 + 7] = 4.883087185713722;

    return re;
  }

  template class Rosenbrock<double, int>;
} // namespace Integrator
