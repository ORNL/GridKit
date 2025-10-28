/**
 * @file hawaiiJson.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Example running a 2-bus system
 *
 * Simulates a 2-bus system with Genrou 6th order generator model and
 * compares results with data generated for the same system by Poweworld.
 *
 */
#include <ctime>
#include <filesystem>
#include <iostream>

#include "hawaii.hpp"
#include <Model/PhasorDynamics/ComponentLibrary.hpp>
#include <Model/PhasorDynamics/SystemModel.hpp>
#include <Model/PhasorDynamics/SystemModelData.hpp>
#include <Model/PhasorDynamics/SystemModelDataJSONParser.hpp>
#include <Solver/Dynamic/Ida.hpp>
#include <Utilities/Testing.hpp>
#include <nlohmann/json.hpp>

int main(int argc, const char* argv[])
{
    using namespace GridKit::PhasorDynamics;
    using namespace AnalysisManager::Sundials;

    using scalar_type = double;
    using real_type   = double;
    using index_type  = size_t;

    //
    // Input file
    //

    std::filesystem::path input_file;
    if (argc < 2)
    {
    if (std::filesystem::exists("hawaii.json"))
    {
        input_file = std::filesystem::current_path() / "hawaii.json";
    }
    else
    {
        std::cout << "\n"
                    "ERROR: No input file found or provided.\n"
                    "\n"
                    "Usage:\n"
                    "       hawaiiJson <json-input-file>\n"
                    "\n"
                    "Please provide a JSON input file as a positional command-line \n"
                    "argument.\n"
                    "\n"
                    "By default this example will look for \"hawaii.json\" in the \n"
                    "current working directory and use that if found.\n"
                    "\n";
        exit(1);
    }
    }
    else
    {
    input_file = argv[1];
    }

    std::cout << "Example: hawaiiJson\n";
    std::cout << "Input file: " << input_file << '\n';

    //
    // Create model data
    //

    SystemModelData<scalar_type, index_type> data(json::parse(std::ifstream(input_file)));

    //
    // Instantiate system model
    //

    SystemModel<scalar_type, index_type> sys(data);
    sys.allocate();

    int status=0;
    return status;
}
