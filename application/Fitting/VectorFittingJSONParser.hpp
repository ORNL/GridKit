/**
 * @file VectorFittingJSONParser.hpp
 *
 * @brief JSON parsing for the VectorFitting application specification.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#pragma once

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#include <nlohmann/json.hpp>

#include "VectorFittingData.hpp"

namespace GridKit
{
  namespace Optimization
  {
    namespace Application
    {
      namespace Detail
      {
        namespace fs = std::filesystem;

        using json = ::nlohmann::json;

        inline std::ifstream openJsonFile(const fs::path& file_path)
        {
          if (!fs::exists(file_path))
          {
            throw std::runtime_error("VectorFitting solver file not found: " + file_path.string());
          }

          auto stream = std::ifstream(file_path);
          if (!stream)
          {
            throw std::runtime_error(
                "Failed to open VectorFitting solver file: " + file_path.string());
          }
          return stream;
        }

        inline void rejectUnsupportedVocabulary(const json& j)
        {
          for (const auto* key : {"frequency_hz", "hz", "weights_column"})
          {
            if (j.contains(key))
            {
              throw std::runtime_error(
                  "VectorFitting solver JSON does not support \"" + std::string(key) + "\"");
            }
          }
        }

        inline RationalTerms parseTerms(const std::string& name)
        {
          if (name == "none")
          {
            return RationalTerms::NONE;
          }
          if (name == "constant")
          {
            return RationalTerms::CONSTANT;
          }
          if (name == "linear")
          {
            return RationalTerms::LINEAR;
          }
          throw std::runtime_error("Invalid VectorFitting terms: " + name);
        }

        inline Weighting parseWeighting(const std::string& name)
        {
          if (name == "uniform")
          {
            return Weighting::UNIFORM;
          }
          if (name == "inverse_magnitude")
          {
            return Weighting::INVERSE_MAGNITUDE;
          }
          if (name == "inverse_squared_magnitude")
          {
            return Weighting::INVERSE_SQUARED_MAGNITUDE;
          }
          throw std::runtime_error("Invalid VectorFitting weighting: " + name);
        }
      } // namespace Detail

      /**
       * @brief Parse a VectorFitting solver specification file.
       *
       * The specification names the sampled-response CSV, the response
       * shape, the fitting request, and the output paths.
       */
      inline VectorFittingData
      parseVectorFittingData(const std::filesystem::path& solver_file)
      {
        using json = Detail::json;

        const auto j = json::parse(Detail::openJsonFile(solver_file));
        Detail::rejectUnsupportedVocabulary(j);

        VectorFittingData data;

        data.samples =
            j.at("samples").template get<std::filesystem::path>();
        if (data.samples.empty())
        {
          throw std::runtime_error(
              "VectorFitting solver JSON requires nonempty \"samples\"");
        }

        data.rows = j.value("rows", data.rows);
        data.cols = j.value("cols", data.cols);
        if (data.rows < 1 || data.cols < 1)
        {
          throw std::runtime_error(
              "VectorFitting requires positive \"rows\" and \"cols\"");
        }

        if (j.contains("fit"))
        {
          const auto& fit = j.at("fit");
          Detail::rejectUnsupportedVocabulary(fit);

          data.fit.poles     = fit.value("poles", data.fit.poles);
          data.fit.min_poles = fit.value("min_poles", data.fit.min_poles);
          data.fit.max_poles = fit.value("max_poles", data.fit.max_poles);
          data.fit.target_rel_rms =
              fit.value("target_rel_rms", data.fit.target_rel_rms);
          data.fit.terms = Detail::parseTerms(
              fit.value("terms", std::string{"constant"}));
          data.fit.weighting = Detail::parseWeighting(
              fit.value("weighting", std::string{"uniform"}));

          if (data.fit.poles < 1)
          {
            throw std::runtime_error(
                "fit: \"poles\" must be positive");
          }
          if (data.fit.max_poles != 0 && (data.fit.min_poles < 1 || data.fit.max_poles < data.fit.min_poles))
          {
            throw std::runtime_error(
                "fit: pole bounds must satisfy 1 <= min <= max");
          }
          if (!(data.fit.target_rel_rms > 0.0))
          {
            throw std::runtime_error(
                "fit: \"target_rel_rms\" must be positive");
          }
        }

        const auto& output = j.at("output");
        Detail::rejectUnsupportedVocabulary(output);
        data.output_model =
            output.at("model").template get<std::filesystem::path>();
        if (data.output_model.empty())
        {
          throw std::runtime_error(
              "VectorFitting output requires nonempty \"model\"");
        }
        if (output.contains("state_space"))
        {
          data.output_state_space =
              output.at("state_space")
                  .template get<std::filesystem::path>();
        }

        const auto base_path = solver_file.parent_path();
        if (!data.samples.is_absolute())
        {
          data.samples = base_path / data.samples;
        }
        if (!data.output_model.is_absolute())
        {
          data.output_model = base_path / data.output_model;
        }
        if (!data.output_state_space.empty() && !data.output_state_space.is_absolute())
        {
          data.output_state_space = base_path / data.output_state_space;
        }

        return data;
      }
    } // namespace Application
  } // namespace Optimization
} // namespace GridKit
