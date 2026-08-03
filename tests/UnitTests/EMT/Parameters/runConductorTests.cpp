#include <cstddef>
#include <exception>

#include <GridKit/Model/EMT/Parameters/Geometry/Conductor/Conductor.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace
{
  using scalar_type   = double;
  using index_type    = size_t;
  using Conductor     = GridKit::EMT::Parameters::Conductor<scalar_type, index_type>;
  using ConductorData = GridKit::EMT::Parameters::ConductorData<scalar_type, index_type>;

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

  ConductorData makeData()
  {
    ConductorData data;
    data.K            = 2;
    data.radius       = {0.01, 0.02};
    data.inner_radius = {0.0, 0.005};
    data.sigma        = {3.7e7, 5.8e7};
    data.mu           = {1.0, 1.0};
    data.weight       = {11.0, 12.0};
    data.phase        = {"a", "b"};
    return data;
  }

  GridKit::Testing::TestOutcome conductor_valid_data()
  {
    using GridKit::Testing::TestStatus;

    Conductor conductor(makeData());

    TestStatus success  = true;
    success            *= (conductor.initialize() == 0);
    success            *= (conductor.size() == 0);
    success            *= conductor.inputs().empty();
    success            *= (conductor.conductorCount() == 2);
    success            *= (conductor.phase(0) == "a");

    auto no_phase = makeData();
    no_phase.phase.clear();
    Conductor conductor_without_phase(no_phase);
    success *= (conductor_without_phase.initialize() == 0);

    return success.report(__func__);
  }

  GridKit::Testing::TestOutcome conductor_invalid_data()
  {
    using GridKit::Testing::TestStatus;

    TestStatus success = true;

    success *= throws([&]
                      {
                        auto data = makeData();
                        data.K    = 0;
                        Conductor conductor(data);
                        conductor.initialize(); });
    success *= throws([&]
                      {
                        auto data = makeData();
                        data.radius.pop_back();
                        Conductor conductor(data);
                        conductor.initialize(); });
    success *= throws([&]
                      {
                        auto data             = makeData();
                        data.inner_radius[0] = data.radius[0];
                        Conductor conductor(data);
                        conductor.initialize(); });
    success *= throws([&]
                      {
                        auto data    = makeData();
                        data.sigma[0] = 0.0;
                        Conductor conductor(data);
                        conductor.initialize(); });
    success *= throws([&]
                      {
                        auto data     = makeData();
                        data.weight[0] = -1.0;
                        Conductor conductor(data);
                        conductor.initialize(); });
    success *= throws([&]
                      {
                        auto data    = makeData();
                        data.phase[0] = "x";
                        Conductor conductor(data);
                        conductor.initialize(); });

    return success.report(__func__);
  }
} // namespace

int main()
{
  GridKit::Testing::TestingResults result;
  result += conductor_valid_data();
  result += conductor_invalid_data();
  return result.summary();
}
