#pragma once

#include <filesystem>

#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Utilities/TestHelpers.hpp>
#include <GridKit/Utilities/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    class SystemModelTests
    {
    public:
      SystemModelTests()
      {
      }

      virtual ~SystemModelTests()
      {
      }

      /**
       * @brief Test for exception thrown when signals are incorrectly
       * configured
       */
      TestOutcome signalError()
      {
        using namespace std::filesystem;
        using namespace GridKit::PhasorDynamics;
        auto                        input_file = current_path() / "ThreeBusBasicBad.json";
        auto                        data       = parseSystemModelData(input_file);
        SystemModel<double, size_t> sys(data);

        TestStatus status{true};
        status *= throws<std::runtime_error>(
            [&]()
            { sys.allocate(); });

        return status.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
