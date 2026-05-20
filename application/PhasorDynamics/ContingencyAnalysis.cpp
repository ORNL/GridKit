#include <future>
#include <thread>
#include <tuple>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

#include "CTGStudy.hpp"

using Clock = std::chrono::high_resolution_clock;
using Log   = GridKit::Utilities::Logger;

using namespace GridKit::PhasorDynamics;
using namespace GridKit::Testing;

TestStatus runStudy(std::size_t faultId, StudyData studyData)
{
  // Change id in schedule to current fault id
  for (auto& event : studyData.events)
  {
    event.element_id = faultId;
  }

  try
  {
    auto study = CTGStudy(studyData);
    study.run();
    return study.checkStatus(CTGStudy::Print{false});
  }
  catch (...)
  {
    Log::warning() << "exception caught at fault id: " << faultId << std::endl;
    return {false};
  }
}

using Dur = std::chrono::duration<double>;

void runStudySerial(const StudyData& studyData, std::vector<TestStatus>& statVec)
{
  for (std::size_t i = 0; i < studyData.model_data.bus_fault.size(); ++i)
  {
    auto stat  = runStudy(i, studyData);
    statVec[i] = stat;
  }
}

void runStudyAsync(const StudyData& studyData, std::vector<TestStatus>& statVec)
{
  auto nFaults = studyData.model_data.bus_fault.size();

  std::vector<std::future<TestStatus>> futures;
  futures.reserve(nFaults);
  for (std::size_t i = 0; i < nFaults; ++i)
  {
    futures.emplace_back(std::async(std::launch::async, runStudy, i, studyData));
  }

  for (std::size_t i = 0; i < nFaults; ++i)
  {
    auto stat  = futures[i].get();
    statVec[i] = stat;
  }
}

#ifdef _OPENMP
void runStudyOpenMP(const StudyData& studyData, std::vector<TestStatus>& statVec)
{
  auto nFaults = studyData.model_data.bus_fault.size();
#pragma omp parallel for
  for (std::size_t i = 0; i < nFaults; ++i)
  {
    auto stat  = runStudy(i, studyData);
    statVec[i] = stat;
  }
}
#endif

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

  auto studyData = parseStudyData(argv[1]);

  const auto start = Clock::now();

  auto faults  = studyData.model_data.bus_fault;
  auto statVec = std::vector<TestStatus>(faults.size(), true);

  // runStudySerial(studyData, statVec);
#ifdef _OPENMP
  std::cout << "omp\n";
  runStudyOpenMP(studyData, statVec);
#else
  std::cout << "async\n";
  runStudyAsync(studyData, statVec);
#endif

  const auto stop = Clock::now();
  const auto dur  = std::chrono::duration<double>(stop - start);
  std::cout << "\n\nComplete in " << dur << "\n";

  TestStatus status;
  for (std::size_t i = 0; i < statVec.size(); ++i)
  {
    status *= statVec[i];
    if (!statVec[i])
    {
      std::cout << "Study failed for fault: "
                << faults[i].disambiguation_string << '\n';
    }
  }

  return status.get();
}
