#include <cmath>
#include <optional>
#include <vector>

#include <GridKit/Model/EMT/Parameters/Geometry/Conductor/Conductor.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Tower/Tower.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace
{
  using scalar_type = double;
  using index_type  = size_t;
  using Conductor   = GridKit::EMT::Parameters::Conductor<scalar_type, index_type>;
  using ConductorData =
      GridKit::EMT::Parameters::ConductorData<scalar_type, index_type>;
  using Tower     = GridKit::EMT::Parameters::Tower<scalar_type, index_type>;
  using TowerData = GridKit::EMT::Parameters::TowerData<scalar_type, index_type>;

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

  TowerData makeTowerData()
  {
    TowerData data;
    data.K        = 2;
    data.span     = 100.0;
    data.height   = {20.0, 24.0};
    data.position = {-5.0, 5.0};
    return data;
  }

  ConductorData makeConductorData()
  {
    ConductorData data;
    data.K            = 2;
    data.radius       = {0.01, 0.01};
    data.inner_radius = {0.0, 0.0};
    data.sigma        = {1.0, 1.0};
    data.mu           = {1.0, 1.0};
    data.weight       = {10.0, 12.0};
    data.phase        = {"a", "b"};
    return data;
  }

  GridKit::Testing::TestOutcome tower_geometry_outputs()
  {
    using GridKit::Testing::TestStatus;

    Conductor conductor(makeConductorData());
    conductor.initialize();
    Tower tower(makeTowerData(), conductor);
    tower.initialize();

    const scalar_type direct_distance = std::sqrt(10.0 * 10.0 + 4.0 * 4.0);
    const scalar_type image_distance  = std::sqrt(10.0 * 10.0 + 44.0 * 44.0);

    TestStatus success  = true;
    success            *= (tower.size() == 0);
    success            *= tower.inputs().empty();
    success            *= close(tower.span(), 100.0);
    success            *= close(tower.height(1), 24.0);
    success            *= close(tower.position(0), -5.0);
    success            *= close(tower.horizontalDistance(0, 1), 10.0);
    success            *= close(tower.distance(0, 1), direct_distance);
    success            *= close(tower.imageDistance(0, 1), image_distance);
    success            *= close(tower.spanLength(0), 100.0);
    success            *= close(tower.spanLength(1), 100.0);
    success            *= close(tower.sag(0), 0.0);
    success            *= close(tower.minimumHeight(0), 20.0);
    success            *= close(tower.height(0), 20.0);
    success            *= close(tower.saggedToSpanRatio(), 1.0);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome tower_optional_tension_outputs()
  {
    using GridKit::Testing::TestStatus;

    Conductor conductor(makeConductorData());
    conductor.initialize();
    auto tower_data    = makeTowerData();
    tower_data.tension = {2000.0, 1800.0};
    Tower tower(tower_data, conductor);
    tower.initialize();

    const scalar_type a0        = 2000.0 / 10.0;
    const scalar_type a1        = 1800.0 / 12.0;
    const scalar_type eta0      = 100.0 / (2.0 * a0);
    const scalar_type eta1      = 100.0 / (2.0 * a1);
    const scalar_type sag0      = a0 * (std::cosh(eta0) - 1.0);
    const scalar_type sag1      = a1 * (std::cosh(eta1) - 1.0);
    const scalar_type length0   = 2.0 * a0 * std::sinh(eta0);
    const scalar_type length1   = 2.0 * a1 * std::sinh(eta1);
    const scalar_type eff0      = 20.0 - sag0 + (2.0 * a0 * a0 / 100.0) * std::sinh(eta0) - a0;
    const scalar_type eff1      = 24.0 - sag1 + (2.0 * a1 * a1 / 100.0) * std::sinh(eta1) - a1;
    const scalar_type distance  = std::sqrt(10.0 * 10.0 + (eff0 - eff1) * (eff0 - eff1));
    const scalar_type ratio_ref = 0.5 * (length0 + length1) / 100.0;

    TestStatus success  = true;
    success            *= close(tower.sag(0), sag0);
    success            *= close(tower.spanLength(0), length0);
    success            *= close(tower.spanLength(1), length1);
    success            *= close(tower.minimumHeight(0), 20.0 - sag0);
    success            *= close(tower.height(0), eff0);
    success            *= close(tower.height(1), eff1);
    success            *= close(tower.distance(0, 1), distance);
    success            *= close(tower.saggedToSpanRatio(), ratio_ref);
    success            *= (tower.saggedToSpanRatio() > 1.0);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome tower_mixed_tension_outputs()
  {
    using GridKit::Testing::TestStatus;

    Conductor conductor(makeConductorData());
    conductor.initialize();
    auto tower_data    = makeTowerData();
    tower_data.tension = {2000.0, std::nullopt};
    Tower tower(tower_data, conductor);
    tower.initialize();

    const scalar_type a0        = 2000.0 / 10.0;
    const scalar_type eta0      = 100.0 / (2.0 * a0);
    const scalar_type sag0      = a0 * (std::cosh(eta0) - 1.0);
    const scalar_type length0   = 2.0 * a0 * std::sinh(eta0);
    const scalar_type eff0      = 20.0 - sag0 + (2.0 * a0 * a0 / 100.0) * std::sinh(eta0) - a0;
    const scalar_type ratio_ref = 0.5 * (length0 / 100.0 + 1.0);

    TestStatus success  = true;
    success            *= close(tower.sag(0), sag0);
    success            *= close(tower.sag(1), 0.0);
    success            *= close(tower.spanLength(0), length0);
    success            *= close(tower.spanLength(1), 100.0);
    success            *= close(tower.minimumHeight(0), 20.0 - sag0);
    success            *= close(tower.minimumHeight(1), 24.0);
    success            *= close(tower.height(0), eff0);
    success            *= close(tower.height(1), 24.0);
    success            *= close(tower.saggedToSpanRatio(), ratio_ref);
    success            *= (tower.saggedToSpanRatio() > 1.0);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome tower_invalid_data()
  {
    using GridKit::Testing::TestStatus;

    TestStatus success = true;

    success *= throws([&]
                      {
                        Conductor conductor(makeConductorData());
                        auto tower_data = makeTowerData();
                        tower_data.span = 0.0;
                        Tower tower(tower_data, conductor);
                        tower.initialize(); });
    success *= throws([&]
                      {
                        Conductor conductor(makeConductorData());
                        auto tower_data    = makeTowerData();
                        tower_data.height = {0.0, 24.0};
                        Tower tower(tower_data, conductor);
                        tower.initialize(); });
    success *= throws([&]
                      {
                        Conductor conductor(makeConductorData());
                        auto tower_data    = makeTowerData();
                        tower_data.tension = {2000.0};
                        Tower tower(tower_data, conductor);
                        tower.initialize(); });
    success *= throws([&]
                      {
                        Conductor conductor(makeConductorData());
                        auto tower_data    = makeTowerData();
                        tower_data.tension = {2000.0, 0.0};
                        Tower tower(tower_data, conductor);
                        tower.initialize(); });

    return success.report(__func__);
  }
} // namespace

int main()
{
  GridKit::Testing::TestingResults result;
  result += tower_geometry_outputs();
  result += tower_optional_tension_outputs();
  result += tower_mixed_tension_outputs();
  result += tower_invalid_data();
  return result.summary();
}
