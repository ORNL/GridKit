#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#include <GridKit/Model/EMT/Constants.hpp>
#include <GridKit/Model/EMT/Parameters/Effects/SkinEffect/SkinEffect.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Conductor/Conductor.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace
{
  using scalar_type = double;
  using index_type  = size_t;

  constexpr scalar_type pi  = GridKit::EMT::Constants::pi<scalar_type>();
  constexpr scalar_type mu0 = GridKit::EMT::Constants::mu0<scalar_type>();

  bool close(scalar_type value, scalar_type ref, scalar_type tol = 1.0e-10)
  {
    return std::abs(value - ref) <= tol * (1.0 + std::abs(ref));
  }

  GridKit::EMT::Parameters::ConductorData<scalar_type, index_type>
  makeConductorData(std::vector<scalar_type> radius,
                    std::vector<scalar_type> sigma,
                    std::vector<scalar_type> mu)
  {
    GridKit::EMT::Parameters::ConductorData<scalar_type, index_type> data;
    data.K            = static_cast<index_type>(radius.size());
    data.radius       = std::move(radius);
    data.inner_radius = std::vector<scalar_type>(static_cast<size_t>(data.K), 0.0);
    data.sigma        = std::move(sigma);
    data.mu           = std::move(mu);
    return data;
  }

  struct SkinLayout
  {
    SkinLayout(index_type conductor_count, index_type retained_roots)
      : K(conductor_count),
        S(retained_roots),
        q(0),
        g(K * S),
        h(g + K),
        w(h + K),
        r(w + K),
        l(r + K),
        n(l + K)
    {
    }

    const index_type K;
    const index_type S;
    const index_type q;
    const index_type g;
    const index_type h;
    const index_type w;
    const index_type r;
    const index_type l;
    const index_type n;
  };

  template <typename ModelT>
  bool zeroResidual(ModelT& skin, const std::vector<scalar_type>& f, index_type base)
  {
    skin.evaluateResidual();
    for (index_type i = 0; i < skin.size(); ++i)
    {
      if (!close(f[base + i], 0.0, 1.0e-8))
      {
        return false;
      }
    }
    return true;
  }

  GridKit::Testing::TestOutcome skin_effect_fixed_root_closed_form()
  {
    using GridKit::Testing::TestStatus;
    using namespace GridKit::EMT::Parameters;

    const index_type  base  = 3;
    const scalar_type r     = 0.01;
    const scalar_type sigma = 3.7e7;
    const scalar_type omega = 2.0 * pi * 60.0;

    Conductor<scalar_type, index_type>  conductor(makeConductorData({r}, {sigma}, {mu0}));
    SkinEffect<scalar_type, index_type> skin(conductor);
    const SkinLayout                    L(1, SkinEffect<scalar_type, index_type>::rootCount());

    std::vector<scalar_type> y(base + skin.size(), 0.0);
    std::vector<scalar_type> yp(base + skin.size(), 0.0);
    std::vector<scalar_type> f(base + skin.size(), 0.0);
    std::vector<index_type>  index(base + skin.size(), 0);
    std::iota(index.begin(), index.end(), index_type{0});

    skin.bind(y.data(), yp.data(), f.data(), index.data(), base);
    skin.setCoordinate(omega, 1.0);
    skin.initialize();

    const scalar_type lb = mu0 / (4.0 * pi);
    scalar_type       g  = 0.0;
    scalar_type       h  = 0.0;
    for (index_type k = 0; k < L.S; ++k)
    {
      const scalar_type xi  = (static_cast<scalar_type>(k) + 0.75) * pi;
      const scalar_type rb  = xi * xi / (4.0 * pi * sigma * r * r);
      const scalar_type q   = rb * rb + omega * omega * lb * lb;
      g                    += rb / q;
      h                    += lb / q;
    }
    const scalar_type w_ref = g * g + omega * omega * h * h;
    const scalar_type r_ref = g / w_ref;
    const scalar_type l_ref = h / w_ref;
    const auto        rskin = skin.rskin();
    const auto        lskin = skin.lskin();

    TestStatus success  = true;
    success            *= (skin.size() == L.n);
    success            *= (rskin.offset == base + L.r);
    success            *= (rskin.rows == 1);
    success            *= (rskin.cols == 1);
    success            *= (lskin.offset == base + L.l);
    success            *= (lskin.rows == 1);
    success            *= (lskin.cols == 1);
    success            *= close(y[rskin.offset], r_ref, 1.0e-12);
    success            *= close(y[lskin.offset], l_ref, 1.0e-12);
    success            *= zeroResidual(skin, f, base);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome skin_effect_residual_form_after_state_perturbation()
  {
    using GridKit::Testing::TestStatus;
    using namespace GridKit::EMT::Parameters;

    const index_type  base  = 2;
    const scalar_type omega = 2.0 * pi * 400.0;
    const SkinLayout  L(1, SkinEffect<scalar_type, index_type>::rootCount());

    Conductor<scalar_type, index_type>  conductor(makeConductorData({0.012}, {5.8e7}, {mu0}));
    SkinEffect<scalar_type, index_type> skin(conductor);

    std::vector<scalar_type> y(base + skin.size(), 0.0);
    std::vector<scalar_type> yp(base + skin.size(), 0.0);
    std::vector<scalar_type> f(base + skin.size(), 0.0);
    std::vector<index_type>  index(base + skin.size(), 0);
    std::iota(index.begin(), index.end(), index_type{0});

    skin.bind(y.data(), yp.data(), f.data(), index.data(), base);
    skin.setCoordinate(omega, 1.0);
    skin.initialize();

    const std::vector<scalar_type> y0 = y;
    const scalar_type              dq = 1.0e-6;
    const scalar_type              dr = 2.0e-6;
    const scalar_type              dl = 3.0e-9;

    TestStatus success  = true;
    success            *= zeroResidual(skin, f, base);

    y              = y0;
    y[base + L.q] += dq;
    skin.evaluateResidual();
    success *= close(f[base + L.q], -dq, 1.0e-12);

    y              = y0;
    y[base + L.r] += dr;
    skin.evaluateResidual();
    success *= close(f[base + L.r], -y0[base + L.w] * dr, 1.0e-12);

    y              = y0;
    y[base + L.l] += dl;
    skin.evaluateResidual();
    success *= close(f[base + L.l], -y0[base + L.w] * dl, 1.0e-12);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome skin_effect_multiconductor_signal_contract()
  {
    using GridKit::Testing::TestStatus;
    using namespace GridKit::EMT::Parameters;

    const index_type  base  = 5;
    const scalar_type omega = 2.0 * pi * 60.0;

    Conductor<scalar_type, index_type> conductor(
        makeConductorData({0.01, 0.02}, {3.7e7, 5.8e7}, {mu0, 2.0 * mu0}));
    SkinEffect<scalar_type, index_type> skin(conductor);
    const SkinLayout                    L(2, SkinEffect<scalar_type, index_type>::rootCount());

    std::vector<scalar_type> y(base + skin.size(), 0.0);
    std::vector<scalar_type> yp(base + skin.size(), 0.0);
    std::vector<scalar_type> f(base + skin.size(), 0.0);
    std::vector<index_type>  index(base + skin.size(), 0);
    std::iota(index.begin(), index.end(), index_type{0});

    skin.bind(y.data(), yp.data(), f.data(), index.data(), base);
    skin.setCoordinate(omega, 1.0);
    skin.initialize();

    const auto rskin = skin.rskin();
    const auto lskin = skin.lskin();

    TestStatus success  = true;
    success            *= (skin.size() == L.n);
    success            *= (skin.inputs().empty());
    success            *= (rskin.offset == base + L.r);
    success            *= (lskin.offset == rskin.offset + rskin.size());
    success            *= (rskin.rows == 2);
    success            *= (rskin.cols == 1);
    success            *= (lskin.rows == 2);
    success            *= (lskin.cols == 1);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome skin_effect_low_frequency_dc_limit()
  {
    using GridKit::Testing::TestStatus;
    using namespace GridKit::EMT::Parameters;

    const index_type  base  = 0;
    const scalar_type r     = 0.01;
    const scalar_type sigma = 3.7e7;
    const scalar_type omega = 2.0 * pi * 1.0e-3;

    Conductor<scalar_type, index_type>  conductor(makeConductorData({r}, {sigma}, {mu0}));
    SkinEffect<scalar_type, index_type> skin(conductor);

    std::vector<scalar_type> y(base + skin.size(), 0.0);
    std::vector<scalar_type> yp(base + skin.size(), 0.0);
    std::vector<scalar_type> f(base + skin.size(), 0.0);
    std::vector<index_type>  index(base + skin.size(), 0);
    std::iota(index.begin(), index.end(), index_type{0});

    skin.bind(y.data(), yp.data(), f.data(), index.data(), base);
    skin.setCoordinate(omega, 1.0);
    skin.initialize();

    const scalar_type r_dc  = 1.0 / (pi * sigma * r * r);
    const scalar_type l_max = mu0 / (4.0 * pi);
    const auto        rskin = skin.rskin();
    const auto        lskin = skin.lskin();

    TestStatus success  = true;
    success            *= close(y[rskin.offset], r_dc, 1.0e-2);
    success            *= (0.0 < y[lskin.offset]);
    success            *= (y[lskin.offset] < l_max);
    success            *= zeroResidual(skin, f, base);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome skin_effect_frequency_trend()
  {
    using GridKit::Testing::TestStatus;
    using namespace GridKit::EMT::Parameters;

    const index_type  base       = 0;
    const scalar_type omega_low  = 2.0 * pi * 60.0;
    const scalar_type omega_high = 2.0 * pi * 10000.0;

    Conductor<scalar_type, index_type>  conductor(makeConductorData({0.01}, {3.7e7}, {mu0}));
    SkinEffect<scalar_type, index_type> skin(conductor);

    std::vector<scalar_type> y(base + skin.size(), 0.0);
    std::vector<scalar_type> yp(base + skin.size(), 0.0);
    std::vector<scalar_type> f(base + skin.size(), 0.0);
    std::vector<index_type>  index(base + skin.size(), 0);
    std::iota(index.begin(), index.end(), index_type{0});

    skin.bind(y.data(), yp.data(), f.data(), index.data(), base);

    skin.setCoordinate(omega_low, 1.0);
    skin.initialize();
    const scalar_type r_low        = y[skin.rskin().offset];
    const scalar_type l_low        = y[skin.lskin().offset];
    const bool        low_residual = zeroResidual(skin, f, base);

    skin.setCoordinate(omega_high, 1.0);
    skin.initialize();
    const scalar_type r_high = y[skin.rskin().offset];
    const scalar_type l_high = y[skin.lskin().offset];

    TestStatus success  = true;
    success            *= low_residual;
    success            *= zeroResidual(skin, f, base);
    success            *= (r_high > 3.0 * r_low);
    success            *= (l_high < 0.5 * l_low);
    success            *= (r_low > 0.0);
    success            *= (r_high > 0.0);
    success            *= (l_low > 0.0);
    success            *= (l_high > 0.0);

    return success.report(__func__);
  }
} // namespace

int main()
{
  GridKit::Testing::TestingResults result;
  result += skin_effect_fixed_root_closed_form();
  result += skin_effect_residual_form_after_state_perturbation();
  result += skin_effect_multiconductor_signal_contract();
  result += skin_effect_low_frequency_dc_limit();
  result += skin_effect_frequency_trend();
  return result.summary();
}
