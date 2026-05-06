#include <omp.h>

#include <future>
#include <thread>
#include <tuple>

#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

#include "CTGStudy.hpp"

using Log = GridKit::Utilities::Logger;

using namespace GridKit::PhasorDynamics;
using namespace GridKit::Testing;

using StudyResult = std::pair<TestStatus, CTGStudy::DurationT>;
using Clock       = std::chrono::high_resolution_clock;

StudyResult runStudy(std::size_t faultId, StudyData studyData)
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
    return std::make_pair(study.checkStatus(), study.getRunTime());
  }
  catch (...)
  {
    std::cout << "exception caught (" << faultId << ")" << std::endl;
    return {false, {}};
  }
}

using Dur = std::chrono::duration<double>;

std::tuple<Dur, Dur> runStudySerial(
    const StudyData& studyData, std::vector<TestStatus>& statVec)
{
  Dur maxDur{};
  Dur avgDur{};
  for (std::size_t i = 0; i < studyData.model_data.bus_fault.size(); ++i)
  {
    auto [stat, dur] = runStudy(i, studyData);

    statVec[i]  = stat;
    avgDur     += dur;
    maxDur      = (maxDur < dur) ? dur : maxDur;
  }
  avgDur /= static_cast<double>(studyData.model_data.bus_fault.size());

  return std::make_tuple(maxDur, avgDur);
}

std::tuple<Dur, Dur> runStudyAsync(
    const StudyData& studyData, std::vector<TestStatus>& statVec)
{
  auto nFaults = studyData.model_data.bus_fault.size();

  std::vector<std::future<StudyResult>> futures;
  futures.reserve(nFaults);
  for (std::size_t i = 0; i < nFaults; ++i)
  {
    futures.emplace_back(std::async(std::launch::async, runStudy, i, studyData));
  }

  Dur maxDur{};
  Dur avgDur{};
  for (std::size_t i = 0; i < nFaults; ++i)
  {
    auto [stat, dur] = futures[i].get();

    statVec[i]  = stat;
    avgDur     += dur;
    maxDur      = (maxDur < dur) ? dur : maxDur;
  }
  avgDur /= static_cast<double>(nFaults);

  return std::make_tuple(maxDur, avgDur);
}

std::tuple<Dur, Dur> runStudyOpenMP(
    const StudyData& studyData, std::vector<TestStatus>& statVec)
{
  Dur::rep maxTicks = 0;
  Dur::rep avgTicks = 0;
  auto     nFaults  = studyData.model_data.bus_fault.size();
#pragma omp parallel for reduction(+ : avgTicks) reduction(max : maxTicks)
  for (std::size_t i = 0; i < nFaults; ++i)
  {
    auto [stat, dur] = runStudy(i, studyData);

    statVec[i] = stat;

    auto ticks  = dur.count();
    avgTicks   += ticks;
    maxTicks    = (maxTicks < ticks) ? ticks : maxTicks;
  }
  avgTicks /= static_cast<double>(nFaults);

  return std::make_tuple(Dur(maxTicks), Dur(avgTicks));
}

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

  auto [maxDur, avgDur] = runStudySerial(studyData, statVec);
  // auto [maxDur, avgDur] = runStudyAsync(studyData, statVec);
  // auto [maxDur, avgDur] = runStudyOpenMP(studyData, statVec);

  const auto stop = Clock::now();
  const auto dur  = std::chrono::duration<double>(stop - start);
  std::cout << "\n\nComplete in " << dur << "\n";
  std::cout << "Max study time: " << maxDur << "\n";
  std::cout << "Avg study time: " << avgDur << "\n";

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
