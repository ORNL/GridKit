#pragma once

#include <filesystem>
#include <fstream>
#include <limits>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SignalSource/SampledSignalSource/SampledSignalSource.hpp>
#include <GridKit/Model/PhasorDynamics/SignalSource/SampledSignalSource/SampledSignalSourceData.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class SampledSignalSourceTests
    {
    private:
      using RealT        = typename PhasorDynamics::Component<ScalarT, IdxT>::RealT;
      using Source       = PhasorDynamics::SignalSource::SampledSignalSource<ScalarT, IdxT>;
      using SourceData   = typename Source::ModelDataT;
      using SourcePorts  = typename SourceData::Ports;
      using SourceParams = typename SourceData::Parameters;

    public:
      TestOutcome inlineSource()
      {
        TestStatus success = true;

        SourceData data;
        data.device_class                     = "SampledSignalSource";
        data.disambiguation_string            = "SRC";
        data.parameters[SourceParams::scale]  = static_cast<RealT>(2.0);
        data.parameters[SourceParams::offset] = static_cast<RealT>(-1.0);
        data.ports[SourcePorts::output]       = 1;
        data.source_type                      = "samples";
        data.samples                          = {{0.0, 0.0}, {0.1, 1.0}, {0.2, 0.0}};

        Source                                    source(data);
        PhasorDynamics::SignalNode<ScalarT, IdxT> output;

        source.getSignals().template assignSignalNode<PhasorDynamics::SignalSource::SampledSignalSourceInternalVariables::OUT>(&output);
        source.allocate();
        success *= (source.verify() == 0);

        source.updateTime(0.05, 0.0);
        success *= isEqual(source.valueAt(0.05), static_cast<RealT>(0.0));
        success *= isEqual(static_cast<RealT>(output.read()), static_cast<RealT>(0.0));

        return success.report(__func__);
      }

      TestOutcome csvSource()
      {
        TestStatus success = true;

        const auto file_path = std::filesystem::current_path() / "sampled_signal_source_test.csv";
        {
          std::ofstream file(file_path);
          file << "t,u\n";
          file << "0.0,0.0\n";
          file << "0.5,1.0\n";
          file << "1.0,0.0\n";
        }

        SourceData data;
        data.device_class               = "SampledSignalSource";
        data.disambiguation_string      = "CSV";
        data.ports[SourcePorts::output] = 1;
        data.source_type                = "csv";
        data.file                       = file_path;
        data.time_column                = "t";
        data.value_column               = "u";

        Source                                    source(data);
        PhasorDynamics::SignalNode<ScalarT, IdxT> output;
        source.getSignals().template assignSignalNode<PhasorDynamics::SignalSource::SampledSignalSourceInternalVariables::OUT>(&output);
        source.allocate();
        source.updateTime(0.25, 0.0);

        success *= (source.verify() == 0);
        success *= isEqual(source.valueAt(0.25), static_cast<RealT>(0.5));
        success *= isEqual(static_cast<RealT>(output.read()), static_cast<RealT>(0.5));

        std::filesystem::remove(file_path);

        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
