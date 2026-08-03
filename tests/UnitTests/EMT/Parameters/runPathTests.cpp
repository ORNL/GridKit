#include <cmath>
#include <vector>

#include <GridKit/Model/EMT/Constants.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Conductor/Conductor.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Path/Path.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Tower/Tower.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace
{
  using scalar_type = double;
  using index_type  = size_t;
  using PathData    = GridKit::EMT::Parameters::PathData<scalar_type, index_type>;
  using PathModel   = GridKit::EMT::Parameters::Path<scalar_type, index_type>;
  using PathPoint   = GridKit::EMT::Parameters::PathPoint<scalar_type>;
  using Conductor   = GridKit::EMT::Parameters::Conductor<scalar_type, index_type>;
  using Tower       = GridKit::EMT::Parameters::Tower<scalar_type, index_type>;

  constexpr scalar_type pi = GridKit::EMT::Constants::pi<scalar_type>();

  bool close(scalar_type value, scalar_type ref, scalar_type tol = 1.0e-12)
  {
    return std::abs(value - ref) <= tol * (1.0 + std::abs(ref));
  }

  template <typename FuncT>
  bool throws(FuncT&& func)
  {
    try
    {
      func();
    }
    catch (const std::exception&)
    {
      return true;
    }
    return false;
  }

  GridKit::EMT::Parameters::TowerData<scalar_type, index_type> makeTowerData()
  {
    GridKit::EMT::Parameters::TowerData<scalar_type, index_type> data;
    data.K        = 2;
    data.span     = 100.0;
    data.height   = {20.0, 20.0};
    data.position = {-5.0, 5.0};
    return data;
  }

  GridKit::EMT::Parameters::ConductorData<scalar_type, index_type> makeConductorData()
  {
    GridKit::EMT::Parameters::ConductorData<scalar_type, index_type> data;
    data.K            = 2;
    data.radius       = {0.01, 0.01};
    data.inner_radius = {0.0, 0.0};
    data.sigma        = {1.0, 1.0};
    data.mu           = {1.0, 1.0};
    data.weight       = {10.0, 10.0};
    data.phase        = {"a", "b"};
    return data;
  }

  struct PathFixture
  {
    PathFixture()
      : conductor(makeConductorData()),
        tower(makeTowerData(), conductor)
    {
      conductor.initialize();
      tower.initialize();
    }

    Conductor conductor;
    Tower     tower;
  };

  GridKit::Testing::TestOutcome path_explicit_length()
  {
    using GridKit::Testing::TestStatus;

    PathFixture fixture;
    PathData    data;
    data.length = 1234.5;

    PathModel path(data, fixture.tower);
    path.initialize();

    TestStatus success  = true;
    success            *= (path.size() == 0);
    success            *= path.inputs().empty();
    success            *= close(path.length(), 1234.5);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome path_gis_equator_length()
  {
    using GridKit::Testing::TestStatus;

    PathData data;
    data.path = {PathPoint{0.0, 0.0}, PathPoint{0.0, 1.0}};

    PathFixture fixture;
    PathModel   path(data, fixture.tower);
    path.initialize();

    const scalar_type ref =
        PathModel::earthRadius() * pi / 180.0 * fixture.tower.saggedToSpanRatio();

    TestStatus success  = true;
    success            *= close(path.length(), ref, 1.0e-12);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome path_gis_multi_interval_dateline()
  {
    using GridKit::Testing::TestStatus;

    PathData data;
    data.path = {PathPoint{0.0, 179.0}, PathPoint{0.0, -179.0}, PathPoint{1.0, -179.0}};

    PathFixture fixture;
    PathModel   path(data, fixture.tower);
    path.initialize();

    const scalar_type ref =
        3.0 * PathModel::earthRadius() * pi / 180.0 * fixture.tower.saggedToSpanRatio();

    TestStatus success  = true;
    success            *= close(path.length(), ref, 1.0e-12);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome path_explicit_length_overrides_gis()
  {
    using GridKit::Testing::TestStatus;

    PathData data;
    data.length = 50.0;
    data.path   = {PathPoint{0.0, 0.0}, PathPoint{0.0, 1.0}};

    PathFixture fixture;
    PathModel   path(data, fixture.tower);
    path.initialize();

    TestStatus success  = true;
    success            *= close(path.length(), 50.0);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome path_zero_state_behavior()
  {
    using GridKit::Testing::TestStatus;

    PathData data;
    data.length = 100.0;
    PathFixture fixture;
    PathModel   path(data, fixture.tower);
    path.initialize();

    std::vector<bool> tag;

    TestStatus success  = true;
    success            *= (path.size() == 0);
    success            *= path.inputs().empty();
    path.tagDifferentiable(tag);
    success *= tag.empty();
    success *= (path.evaluateResidual() == 0);
    success *= (path.evaluateJacobian() == 0);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome path_invalid_data()
  {
    using GridKit::Testing::TestStatus;

    TestStatus success = true;

    success *= throws([&]
                      {
                        PathFixture fixture;
                        PathData data;
                        data.length = -1.0;
                        PathModel path(data, fixture.tower);
                        path.initialize(); });
    success *= throws([&]
                      {
                        PathFixture fixture;
                        PathData data;
                        PathModel path(data, fixture.tower);
                        path.initialize(); });
    success *= throws([&]
                      {
                        PathFixture fixture;
                        PathData data;
                        data.path = {PathPoint{0.0, 0.0}};
                        PathModel path(data, fixture.tower);
                        path.initialize(); });
    success *= throws([&]
                      {
                        PathFixture fixture;
                        PathData data;
                        data.path = {PathPoint{91.0, 0.0}, PathPoint{0.0, 0.0}};
                        PathModel path(data, fixture.tower);
                        path.initialize(); });
    success *= throws([&]
                      {
                        PathFixture fixture;
                        PathData data;
                        data.path = {PathPoint{0.0, -181.0}, PathPoint{0.0, 0.0}};
                        PathModel path(data, fixture.tower);
                        path.initialize(); });
    success *= throws([&]
                      {
                        PathFixture fixture;
                        PathData data;
                        data.path = {PathPoint{0.0, 0.0}, PathPoint{0.0, 0.0}};
                        PathModel path(data, fixture.tower);
                        path.initialize(); });

    return success.report(__func__);
  }
} // namespace

int main()
{
  GridKit::Testing::TestingResults result;
  result += path_explicit_length();
  result += path_gis_equator_length();
  result += path_gis_multi_interval_dateline();
  result += path_explicit_length_overrides_gis();
  result += path_zero_state_behavior();
  result += path_invalid_data();
  return result.summary();
}
