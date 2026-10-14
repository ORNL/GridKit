#include <filesystem>
#include <string>
#include <vector>

#include <GridKit/Utilities/CliArgs/CliArgs.hpp>

#include "SystemJacobianTests.hpp"

int main(int argc, const char* argv[])
{
  using Log = GridKit::Utilities::Logger;

  // Print timings
  Log::raiseVerbosity(Log::SUMMARY);

  GridKit::Utilities::CliArgs args{{.name     = {"--cases", "-c"},
                                    .help     = "Case file, or directory searched recursively for *.case.json",
                                    .required = true}};
  args.parseArgs(argc, argv);

  const std::filesystem::path path = args["cases"]();

  std::vector<std::filesystem::path> case_files;
  if (std::filesystem::is_directory(path))
  {
    for (const auto& entry : std::filesystem::recursive_directory_iterator(path))
    {
      if (entry.is_regular_file() && entry.path().filename().string().ends_with(".case.json"))
      {
        case_files.push_back(entry.path());
      }
    }
  }
  else
  {
    case_files.push_back(path);
  }

  if (case_files.empty())
  {
    Log::error() << "No case files found\n";
    return 1;
  }

  GridKit::Testing::TestingResults                      result;
  GridKit::Testing::SystemJacobianTests<double, size_t> test;

  for (const auto& case_file : case_files)
  {
    result += test.jacobian(case_file);
  }

  test.printTimings();

  return result.summary();
}
