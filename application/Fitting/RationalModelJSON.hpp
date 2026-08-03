/**
 * @file RationalModelJSON.hpp
 *
 * @brief JSON serialization of fitted rational models, shared by the
 *        fitting applications.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#pragma once

#include <nlohmann/json.hpp>

#include <GridKit/Solver/Optimization/Rational/RationalModel.hpp>

namespace GridKit
{
  namespace Optimization
  {
    namespace Application
    {
      /**
       * @brief Serialize a rational model with the VectorFit operator
       *        vocabulary: D, E, poles as [real, imag] pairs, residues
       *        pole-major over row-major channels.
       */
      template <typename scalar_type, typename index_type>
      nlohmann::json
      modelToJson(const RationalModel<scalar_type, index_type>& model)
      {
        nlohmann::json j;
        j["rows"] = model.rows;
        j["cols"] = model.cols;

        if (!model.d.empty())
        {
          j["D"] = model.d;
        }
        if (!model.e.empty())
        {
          j["E"] = model.e;
        }

        auto poles = nlohmann::json::array();
        for (const auto& pole : model.poles)
        {
          poles.push_back({pole.real(), pole.imag()});
        }
        j["poles"] = poles;

        const auto channels = static_cast<size_t>(model.rows) * static_cast<size_t>(model.cols);

        auto residues = nlohmann::json::array();
        for (size_t q = 0; q < model.poles.size(); ++q)
        {
          auto per_channel = nlohmann::json::array();
          for (size_t ch = 0; ch < channels; ++ch)
          {
            const auto value = model.residues[q * channels + ch];
            per_channel.push_back({value.real(), value.imag()});
          }
          residues.push_back(per_channel);
        }
        j["residues"] = residues;
        return j;
      }
    } // namespace Application
  } // namespace Optimization
} // namespace GridKit
