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
       *        vocabulary and shapes: D and E as rows x cols real
       *        matrices, poles as [real, imag] pairs, and residues as
       *        one rows x cols matrix of [real, imag] entries per pole.
       */
      template <typename scalar_type, typename index_type>
      nlohmann::json
      modelToJson(const RationalModel<scalar_type, index_type>& model)
      {
        const auto rows = static_cast<size_t>(model.rows);
        const auto cols = static_cast<size_t>(model.cols);

        const auto real_matrix =
            [rows, cols](const std::vector<typename RationalModel<
                             scalar_type,
                             index_type>::RealT>& flat)
        {
          auto matrix = nlohmann::json::array();
          for (size_t row = 0; row < rows; ++row)
          {
            auto values = nlohmann::json::array();
            for (size_t col = 0; col < cols; ++col)
            {
              values.push_back(flat[row * cols + col]);
            }
            matrix.push_back(values);
          }
          return matrix;
        };

        nlohmann::json j;
        j["rows"] = model.rows;
        j["cols"] = model.cols;

        if (!model.d.empty())
        {
          j["D"] = real_matrix(model.d);
        }
        if (!model.e.empty())
        {
          j["E"] = real_matrix(model.e);
        }

        auto poles = nlohmann::json::array();
        for (const auto& pole : model.poles)
        {
          poles.push_back({pole.real(), pole.imag()});
        }
        j["poles"] = poles;

        auto residues = nlohmann::json::array();
        for (size_t q = 0; q < model.poles.size(); ++q)
        {
          auto matrix = nlohmann::json::array();
          for (size_t row = 0; row < rows; ++row)
          {
            auto values = nlohmann::json::array();
            for (size_t col = 0; col < cols; ++col)
            {
              const auto value =
                  model.residues[(q * rows + row) * cols + col];
              values.push_back({value.real(), value.imag()});
            }
            matrix.push_back(values);
          }
          residues.push_back(matrix);
        }
        j["residues"] = residues;
        return j;
      }
    } // namespace Application
  } // namespace Optimization
} // namespace GridKit
