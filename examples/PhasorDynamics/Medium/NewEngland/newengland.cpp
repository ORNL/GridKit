/**
 * @file newengland.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Example running the New England IEEE 39-Bus Case
 *
 * Simulates a New England IEEE 39-Bus with Genrou 6th order generator model and
 * compares results with data generated for the same system by Poweworld.
 *
 */
#include <ctime>
#include <filesystem>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/ComponentLibrary.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Utilities/Testing.hpp>

#include "newengland.hpp"

int main(int argc, const char* argv[])
{
  using namespace GridKit::PhasorDynamics;
  using namespace AnalysisManager::Sundials;

  using scalar_type = double;
  // using real_type   = double;
  using index_type  = size_t;

  // Read Input JSON File
  std::filesystem::path input_file;
  if (argc < 2)
  {
    if (std::filesystem::exists("newengland.json"))
    {
      input_file = std::filesystem::current_path() / "newengland.json";
    }
    else
    {
      std::cout << "\n"
                   "ERROR: No input file found or provided.\n"
                   "\n"
                   "Usage:\n"
                   "       newengland <json-input-file>\n"
                   "\n"
                   "Please provide a JSON input file as a positional command-line \n"
                   "argument.\n"
                   "\n"
                   "By default this example will look for \"newengland.json\" in the \n"
                   "current working directory and use that if found.\n"
                   "\n";
      exit(1);
    }
  }
  else
  {
    input_file = argv[1];
  }

  std::cout << "Example: newengland\n";
  std::cout << "Input file: " << input_file << '\n';

  // Parse Model Data
  auto data = parseSystemModelData(input_file);

  //  Instantiate System Model
  SystemModel<scalar_type, index_type> sys(data);
  sys.allocate();


  int status = 0;
  return status;
}
