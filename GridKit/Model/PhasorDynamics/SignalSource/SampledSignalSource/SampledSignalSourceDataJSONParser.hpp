#pragma once

#include <sstream>
#include <stdexcept>

#include <nlohmann/json.hpp>

#include <GridKit/Model/PhasorDynamics/ComponentDataJSONParser.hpp>
#include <GridKit/Model/PhasorDynamics/SignalSource/SampledSignalSource/SampledSignalSourceData.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace SignalSource
    {
      using json = nlohmann::json;
      using Log  = ::GridKit::Utilities::Logger;

      template <typename RealT, typename IdxT>
      void from_json(const json& j, SampledSignalSourceData<RealT, IdxT>& data)
      {
        using DataT = SampledSignalSourceData<RealT, IdxT>;
        auto& base  = static_cast<ComponentData<RealT,
                                                IdxT,
                                                typename DataT::Parameters,
                                                typename DataT::Ports,
                                                typename DataT::MonitorableVariables>&>(data);
        GridKit::PhasorDynamics::from_json(j, base);

        std::stringstream error_context;
        error_context << "\n\tSee the \"SampledSignalSource\" device with \"id\": \""
                      << data.disambiguation_string
                      << "\" in the \"devices\" list of your JSON file.";

        auto fail = [&](const std::string& message)
        {
          std::stringstream ss;
          ss << "\n\t" << message << error_context.str();
          Log::error() << ss.str() << std::endl;
          throw std::runtime_error(ss.str());
        };

        if (!j.contains("source"))
        {
          fail("SampledSignalSource requires a \"source\" object.");
        }

        const auto& source = j.at("source");
        source.at("type").get_to(data.source_type);

        if (data.source_type == "csv")
        {
          data.file         = source.at("file").get<std::string>();
          data.time_column  = source.value("time_column", std::string{"t"});
          data.value_column = source.value("value_column", std::string{"u"});
        }
        else if (data.source_type == "samples")
        {
          if (!source.contains("samples") || !source.at("samples").is_array())
          {
            fail("Inline SampledSignalSource data requires a \"samples\" array.");
          }

          for (const auto& sample : source.at("samples"))
          {
            if (!sample.is_array() || sample.size() != 2)
            {
              fail("Each inline sample must be a two-entry [time, value] array.");
            }

            data.samples.emplace_back(sample.at(0).template get<RealT>(),
                                      sample.at(1).template get<RealT>());
          }
        }
        else
        {
          fail("SampledSignalSource source \"type\" must be \"csv\" or \"samples\".");
        }
      }
    } // namespace SignalSource
  } // namespace PhasorDynamics
} // namespace GridKit
