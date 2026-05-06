#pragma once

#include <chrono>
#include <filesystem>
#include <fstream>
#include <vector>

#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace fs = ::std::filesystem;

    /**
     * @brief Describes an event that is used to modify the simulation at the
     * given time point
     */
    struct SystemEvent
    {
      /// Type of event determines action performed
      enum class Type
      {
        FAULT_ON,
        FAULT_OFF
      };

      /// Time event takes place
      double      time;
      /// Event type
      Type        type;
      /// ID of element used in event (e.g., bus fault id)
      std::size_t element_id;
    };

    /**
     * @brief Data defined in JSON file for parameterized study
     */
    struct StudyData
    {
      /// path to system model JSON file
      fs::path                 system_model_file;
      /// time step size
      double                   dt;
      /// max time
      double                   tmax;
      /// set of system events
      std::vector<SystemEvent> events;
      /// path to output file
      fs::path                 output_file;
      /// path to reference file for validation
      fs::path                 reference_file;
      /// Error tolerance (between output file and reference file)
      double                   error_tol;
      /// Instance of model data
      SystemModelData<>        model_data;
    };

    /**
     * @brief Wrapper function to parse `StudyData` from JSON and perform
     * follow-up configuration
     */
    StudyData parseStudyData(const fs::path& file_path);

    class CTGStudy
    {
    public:
      using ScalarT      = double;
      using IdxT         = std::size_t;
      using RealT        = double;
      using DurationT    = std::chrono::duration<RealT>;
      using SystemModelT = SystemModel<ScalarT, IdxT>;
      using IdaT         = AnalysisManager::Sundials::Ida<ScalarT, IdxT>;

      CTGStudy(const StudyData& data);

      void run();

      GridKit::Testing::TestStatus checkStatus(bool print = true) const;

      DurationT getRunTime() const
      {
        return dur_;
      }

    private:
      StudyData    study_data_;
      SystemModelT sys_;
      IdaT         ida_;
      DurationT    dur_{};
    };

  } // namespace PhasorDynamics
} // namespace GridKit
