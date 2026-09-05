#pragma once

#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

#include <GridKit/Model/EMT/SystemModel.hpp>

#include "AnalysisUtilities.hpp"

namespace GridKit::EMT
{
  /** Optional complete DAE recording, including states without model monitors. */
  template <typename ScalarT, typename IdxT>
  class StateMonitor
  {
    using System = SystemModel<ScalarT, IdxT>;

  public:
    StateMonitor(System& system, const StudyData& study)
      : system_(system)
    {
      if (study.state_output_file.empty())
      {
        return;
      }
      if (!study.output_file.empty()
          && (fs::absolute(study.state_output_file).lexically_normal()
                  == fs::absolute(study.output_file).lexically_normal()
              || (fs::exists(study.state_output_file) && fs::exists(study.output_file)
                  && fs::equivalent(study.state_output_file, study.output_file))))
      {
        throw std::invalid_argument("State and variable monitors require distinct output files");
      }
      output_.open(study.state_output_file);
      if (!output_)
      {
        throw std::runtime_error("Cannot open state output file: " + study.state_output_file.string());
      }
      output_.exceptions(std::ios::badbit | std::ios::failbit);
      output_ << std::setprecision(17);

      std::vector<std::string> names(system_.size());
      json                     metadata;
      metadata["variables"] = json::array();
      auto describe         = [&](auto&& self, const json& body, const std::string& prefix) -> void
      {
        for (const auto& device : body.at("devices"))
        {
          const auto  path      = prefix + device.at("id").get<std::string>();
          const auto& component = system_.component(path);
          if (dynamic_cast<const Container<ScalarT, IdxT>*>(&component) != nullptr)
          {
            self(self, device, path + ".");
            continue;
          }
          const auto& indices = component.getVariableIndices();
          for (std::size_t local = 0; local < indices.size(); ++local)
          {
            const auto global = indices[local];
            names.at(global)  = path + "[" + std::to_string(local) + "]";
            metadata["variables"].push_back({{"index", global}, {"component", path}, {"class", device.at("class")}, {"local_index", local}, {"differential", system_.tag().at(global)}});
          }
        }
      };
      describe(describe, json::parse(openFile(study.system_model_file)), "");
      metadata["description"] = "CSV columns are time, all y, then all yp in global index order. Local indices and units follow each component's internal variable definitions. Algebraic yp values are solver interpolants, not independent states.";
      std::ofstream schema(study.state_output_file.string() + ".json");
      schema.exceptions(std::ios::badbit | std::ios::failbit);
      schema << metadata.dump(2) << '\n';
      output_ << "time";
      for (const auto& kind : {"y:", "yp:"})
      {
        for (const auto& name : names)
        {
          output_ << ",\"" << kind;
          for (const auto character : name)
          {
            if (character == '"')
            {
              output_ << '"';
            }
            output_ << character;
          }
          output_ << '"';
        }
      }
      output_ << '\n';
    }

    void write(double time)
    {
      if (!output_.is_open())
      {
        return;
      }
      output_ << time;
      for (const auto* values : {system_.y().getData(), system_.yp().getData()})
      {
        for (IdxT i = 0; i < system_.size(); ++i)
        {
          output_ << ',' << values[i];
        }
      }
      output_ << '\n';
    }

  private:
    System&       system_;
    std::ofstream output_;
  };
} // namespace GridKit::EMT
