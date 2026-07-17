/**
 * @file ThreeBusClassical.cpp
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @brief Example running a 3-bus system with Zip load model
 *
 * Simulates a 3-bus system with classical generators and ZIP loads.
 *
 */
#include <ctime>
#include <filesystem>
#include <fstream>
#include <vector>

#include <GridKit/Model/PhasorDynamics/ComponentLibrary.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/Testing.hpp>

#include "ThreeBusZipLoad.hpp"

using scalar_type = double;
using real_type   = double;
using index_type  = size_t;

struct OutputData
{
  real_type t;
  real_type gen2speed;
  real_type gen3speed;
  real_type v2mag;
  real_type v3mag;

  OutputData& operator-=(const OutputData& other)
  {
    assert(GridKit::Testing::isEqual(t, other.t, reference_tol));
    gen2speed -= other.gen2speed;
    gen3speed -= other.gen3speed;
    v2mag     -= other.v2mag;
    v3mag     -= other.v3mag;
    return *this;
  }

  real_type norm() const
  {
    return std::max({
        std::abs(gen2speed),
        std::abs(gen3speed),
        std::abs(v2mag),
        std::abs(v3mag),
    });
  }
};

const OutputData operator-(const OutputData& lhs, const OutputData& rhs)
{
  return OutputData(lhs) -= rhs;
}

std::ostream& operator<<(std::ostream& out, const OutputData& data)
{
  out << data.t << ","
      << data.gen2speed << ","
      << data.gen3speed << ","
      << data.v2mag << ","
      << data.v3mag;
  return out;
}

int main(int argc, const char* argv[])
{
  using namespace GridKit::PhasorDynamics;
  using namespace AnalysisManager::Sundials;

  auto error_allowed = static_cast<real_type>(0.003);

  //
  // Input file
  //

  std::filesystem::path input_file;
  if (argc < 2)
  {
    if (std::filesystem::exists("ThreeBusZipLoad.json"))
    {
      input_file = std::filesystem::current_path() / "ThreeBusZipLoad.json";
    }
    else
    {
      std::cout << "\n"
                   "ERROR: No input file found or provided.\n"
                   "\n"
                   "Usage:\n"
                   "       ThreeBusZipLoad <json-input-file>\n"
                   "\n"
                   "Please provide a JSON input file as a positional command-line \n"
                   "argument.\n"
                   "\n"
                   "By default this example will look for \"ThreeBusZipLoad.json\" in the \n"
                   "current working directory and use that if found.\n"
                   "\n";
      exit(1);
    }
  }
  else
  {
    input_file = argv[1];
  }

  std::cout << "Example: ThreeBusZipLoadJson\n";
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

  std::vector<OutputData> output;

  auto output_cb = [&](real_type t)
  {
    const auto* y_val = sys.y().getData();

    output.push_back(OutputData{t,
                                1 + static_cast<real_type>(y_val[7]),
                                1 + static_cast<real_type>(y_val[12]),
                                std::hypot(static_cast<real_type>(y_val[0]), static_cast<real_type>(y_val[1])),
                                std::hypot(static_cast<real_type>(y_val[2]), static_cast<real_type>(y_val[3]))});
  };

  // Set up simulation
  Ida<scalar_type, index_type> ida(&sys);
  ida.configureSimulation();

  // Run simulation, output each `dt` interval
  real_type start = static_cast<real_type>(clock());
  ida.initializeSimulation(0.0, false);

  // Run for 1s
  ida.runSimulation(1.0, dt, output_cb);

  // Introduce fault to ground and run for 0.1s
  fault->setStatus(true);
  ida.initializeSimulation(1.0, false);
  ida.runSimulation(1.1, dt, output_cb);

  // Clear fault and run until t = 10s.
  fault->setStatus(false);
  ida.initializeSimulation(1.1, false);
  ida.runSimulation(10.0, dt, output_cb);
  real_type stop = static_cast<real_type>(clock());

  /* Check worst-case error */
  real_type worst_error      = 0;
  real_type worst_error_time = 0;

  std::ostream nullout(nullptr);
  // std::ostream& out = nullout;

  // Uncomment code below to print output to a file:
  std::ofstream fileout;
  fileout.open("Example_ThreeBus_ZipLoad_results.csv");
  std::ostream& out = fileout;

  out << "Time,gen2speed,gen3speed,v2mag,v3mag\n";
  out << 0. << "," << 1. << "," << 1. << "," << 1. << "," << 1. << "\n";

  for (index_type i = 0; i < output.size(); ++i)
  {
    OutputData ref{reference_solution[i + 1][0],
                   reference_solution[i + 1][1],
                   reference_solution[i + 1][2],
                   reference_solution[i + 1][4],
                   reference_solution[i + 1][5]};
    OutputData out_data = output[i];

    out << out_data << '\n';

    real_type err = (out_data - ref).norm();
    if (err > worst_error)
    {
      worst_error      = err;
      worst_error_time = out_data.t;
    }
  }
  fileout.close();

  std::cout << "Max error " << worst_error
            << " at time t = " << worst_error_time << "\n";
  std::cout << "\n\nComplete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n\n";

  int success = 0;
  if (worst_error < error_allowed)
  {
    std::cout << "Test result: PASS\n\n";
  }
  else
  {
    std::cout << "Test result: FAIL\n\n";
    success = 1;
  }

  return success;
}
