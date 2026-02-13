/**
 * @file ThreeBusBasicJson.cpp
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Example running a 3-bus system
 *
 * Simulates a 3-bus system with two Genrou 6th order generator models and
 * compares results with data generated for the same system by Poweworld.
 *
 */
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <vector>

#include <GridKit/Model/PhasorDynamics/ComponentLibrary.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Utilities/CliOptions/CliOptions.hpp>
#include <GridKit/Utilities/Testing.hpp>

#define ERROR_TOL 1.0e-4

using scalar_type = double;
using real_type   = double;
using index_type  = size_t;

template <typename T = std::string>
class Tokenizer
{
public:
  Tokenizer() = delete;

  explicit Tokenizer(const std::string& in, char delimiter = ' ')
  {
    std::istringstream iss(in);
    for (std::string item; std::getline(iss, item, delimiter);)
    {
      std::istringstream(item) >> tokens_.emplace_back();
    }
  }

  const std::vector<T>& operator()() const
  {
    return tokens_;
  }

private:
  std::vector<T> tokens_;
};

enum class Norm
{
  L1,
  L2,
  LInf
};

template <std::size_t N>
struct TimeDataGroup
{
  real_type                t;
  std::array<real_type, N> data;

  TimeDataGroup(const std::vector<real_type>& v) : t{v[0]}
  {
    assert(v.size() == N + 1);
    std::copy(next(begin(v)), end(v), begin(data));
  }

  template <typename... TArgs>
  TimeDataGroup(real_type t, TArgs&&... args) : t{t}, data{args...}
  {
    static_assert(sizeof...(args) == N);
  }

private:
  friend TimeDataGroup operator-(const TimeDataGroup& a, const TimeDataGroup& b)
  {
    assert(GridKit::Testing::isEqual(a.t, b.t, ERROR_TOL));
    TimeDataGroup ret(a);
    for (std::size_t i = 0; i < N; ++i)
    {
      ret.data[i] -= b.data[i];
    }
    return ret;
  }

  friend real_type l1Norm(const TimeDataGroup& a)
  {
    real_type mx = 0.0;
    for (auto v : a.data)
    {
      mx = std::max(mx, std::abs(v));
    }
    return mx;
  }

  friend real_type l2Norm(const TimeDataGroup& a)
  {
    real_type ret = 0.0;
    for (auto v : a.data)
    {
      ret += v * v;
    }
    return std::sqrt(ret);
  }

  friend real_type lInfNorm(const TimeDataGroup& a)
  {
    real_type ret = 0.0;
    for (auto v : a.data)
    {
      ret += std::abs(v);
    }
    return ret;
  }

  friend real_type errorNorm(
      const TimeDataGroup& a, const TimeDataGroup& b, Norm norm = Norm::L1)
  {
    if (norm == Norm::L2)
    {
      return l2Norm(a - b);
    }

    if (norm == Norm::LInf)
    {
      return lInfNorm(a - b);
    }

    return l1Norm(a - b);
  }
};

int main(int argc, const char* argv[])
{
  using namespace GridKit::PhasorDynamics;
  using namespace AnalysisManager::Sundials;
  using namespace GridKit::Utilities;

  auto error_allowed = static_cast<real_type>(1e-4);

  //
  // Input file
  //

  std::filesystem::path input_file;
  if (argc < 2)
  {
    if (std::filesystem::exists("ThreeBusBasic.json"))
    {
      input_file = std::filesystem::current_path() / "ThreeBusBasic.json";
    }
    else
    {
      std::cout << "\n"
                   "ERROR: No input file found or provided.\n"
                   "\n"
                   "Usage:\n"
                   "       ThreeBusBasicJson <json-input-file>\n"
                   "\n"
                   "Please provide a JSON input file as a positional command-line \n"
                   "argument.\n"
                   "\n"
                   "By default this example will look for \"ThreeBusBasic.json\" in the \n"
                   "current working directory and use that if found.\n"
                   "\n";
      exit(1);
    }
  }
  else
  {
    input_file = argv[1];
  }

  std::cout << "Example: ThreeBusBasicJson\n";
  std::cout << "Input file: " << input_file << '\n';

  //
  // Create model data
  //

  auto data = parseSystemModelData(input_file);

  //
  // Instantiate system
  //

  SystemModel<scalar_type, index_type> sys(data);
  sys.allocate();

  // Get access to fault 0
  auto* fault = sys.getBusFault(0);

  real_type dt = 1.0 / 4.0 / 60.0;

  // Set up simulation
  Ida<scalar_type, index_type> ida(&sys);
  ida.configureSimulation();

  // Run simulation, output each `dt` interval
  scalar_type start = static_cast<scalar_type>(clock());
  ida.initializeSimulation(0.0, false);

  // Run for 1s
  int nout = static_cast<int>(std::round((1.0 - 0.0) / dt));
  ida.runSimulation(1.0, nout);

  // Introduce fault to ground and run for 0.1s
  fault->setStatus(true);
  ida.initializeSimulation(1.0, false);
  nout = static_cast<int>(std::round((1.1 - 1.0) / dt));
  ida.runSimulation(1.1, nout);

  // Clear fault and run until t = 10s.
  fault->setStatus(false);
  ida.initializeSimulation(1.1, false);
  nout = static_cast<int>(std::round((10.0 - 1.1) / dt));
  ida.runSimulation(10.0, nout);
  double stop = static_cast<double>(clock());

  sys.stopMonitor();

  /* Check worst-case error */
  real_type max_error      = 0;
  real_type max_error_time = 0;

  std::ifstream ifs("mon.csv");
  std::ifstream ifs_ref("ThreeBusBasic.ref.csv");

  std::string line;
  std::string line_ref;
  std::getline(ifs, line);
  std::getline(ifs_ref, line_ref);
  for (; ifs && ifs_ref; std::getline(ifs, line), std::getline(ifs_ref, line_ref))
  {
    TimeDataGroup<5> ref(Tokenizer<real_type>(line_ref, ',')());
    TimeDataGroup<5> grp(Tokenizer<real_type>(line, ',')());

    auto err = errorNorm(grp, ref);

    if (err > max_error)
    {
      max_error      = err;
      max_error_time = grp.t;
    }
  }

  std::cout
      << "Max error " << max_error
      << " at time t = " << max_error_time << "\n";
  std::cout << "\n\nComplete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";

  return max_error < error_allowed ? 0 : 1;
}
