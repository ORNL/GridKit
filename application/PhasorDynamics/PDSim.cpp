#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

#include "CTGStudy.hpp"

using Log = GridKit::Utilities::Logger;

using namespace GridKit::PhasorDynamics;

int main(int argc, const char* argv[])
{
  // Study file
  if (argc < 2)
  {
    Log::error() << "No input file provided" << std::endl;
    std::cout << "\n"
                 "Usage:\n"
                 "       pdsim <json-input-file>\n"
                 "\n"
                 "Please provide a json input file for the study to run.\n"
                 "\n";
    exit(1);
  }

  auto study = CTGStudy(parseStudyData(argv[1]));

  study.run();

  std::cout << "\n\nComplete in " << study.getRunTime() << " seconds\n";

  return study.checkStatus().get();
}
