#pragma once

#include <array>

#include "VectorFit.hpp"

namespace GridKit::Examples::EMT
{
  struct LineData
  {
    size_t              N, K;
    std::vector<size_t> conductors; // Zero-based internally.
    VectorFit           yc, input, output;
    std::vector<double> delays;

    static LineData parse(const Json& value, const std::filesystem::path& directory)
    {
      if (value.is_string())
      {
        const auto path = directory / value.get<std::string>();
        return parse(readJson(path), path.parent_path());
      }
      fields(value, {"N", "K", "conductors", "Yc", "H"});
      LineData data;
      data.N              = dimension(value.at("N"));
      data.K              = dimension(value.at("K"));
      const auto& mapping = value.at("conductors");
      if (!mapping.is_array() || mapping.size() != data.K)
      {
        throw std::invalid_argument("Line conductor mapping must have K entries");
      }
      std::vector<bool> present(data.N, false);
      for (const auto& c : mapping)
      {
        const size_t phase = dimension(c) - 1;
        if (phase >= data.N)
        {
          throw std::invalid_argument("Line conductor phase out of range");
        }
        data.conductors.push_back(phase);
        present[phase] = true;
      }
      if (std::find(present.begin(), present.end(), false) != present.end())
      {
        throw std::invalid_argument("Line conductor mapping must cover every phase");
      }
      data.yc          = VectorFit::parse(value.at("Yc"), directory);
      Json h           = value.at("H");
      auto h_directory = directory;
      if (h.is_string())
      {
        const auto path = directory / h.get<std::string>();
        h               = readJson(path);
        h_directory     = path.parent_path();
      }
      if (h.contains("fit"))
      {
        fields(h, {"fit", "delay"});
        data.input  = VectorFit::identity(data.K);
        data.output = VectorFit::parse(h.at("fit"), h_directory);
        data.delays.assign(data.K, real(h.at("delay")));
      }
      else
      {
        fields(h, {"input", "delays", "output"});
        data.input  = VectorFit::parse(h.at("input"), h_directory);
        data.output = VectorFit::parse(h.at("output"), h_directory);
        if (!h.at("delays").is_array())
        {
          throw std::invalid_argument("Line delays must be an array");
        }
        for (const auto& delay : h.at("delays"))
        {
          data.delays.push_back(real(delay));
        }
      }
      for (const auto* fit : {&data.yc, &data.input, &data.output})
      {
        if (fit->rows != data.K || fit->cols != data.K)
        {
          throw std::invalid_argument("Line fits must be K by K");
        }
      }
      if (data.delays.size() != data.K)
      {
        throw std::invalid_argument("Line requires one delay per mode");
      }
      for (double delay : data.delays)
      {
        if (delay < 0)
        {
          throw std::invalid_argument("Line delays must be nonnegative");
        }
      }
      return data;
    }
  };

  /// Reciprocal line with independent terminal and directional states.
  /// All currents i[terminal] are positive from the line into the circuit.
  class LineDistributed
  {
  public:
    explicit LineDistributed(LineData coefficients)
      : data(std::move(coefficients))
    {
    }

    LineData                               data;
    std::array<size_t, 2>                  reflected, modal;
    std::array<std::vector<Expression>, 2> characteristic, incident, current;

    void allocate(LinearDAE& dae, double voltage_scale, double current_scale)
    {
      for (size_t end = 0; end < 2; ++end)
      {
        reflected[end]     = dae.variables(data.K, false, current_scale);
        modal[end]         = dae.variables(data.K, false, current_scale);
        yc_state_[end]     = data.yc.allocate(dae, voltage_scale);
        input_state_[end]  = data.input.allocate(dae, current_scale);
        output_state_[end] = data.output.allocate(dae, current_scale);
      }
    }

    void stamp(LinearDAE& dae, std::array<size_t, 2> voltage, size_t delayed)
    {
      for (size_t end = 0; end < 2; ++end)
      {
        characteristic[end] = data.yc.stamp(dae, yc_state_[end], voltage[end], data.conductors);
        // delayed[(1-end)*K] carries waves launched at the opposite end.
        incident[end]       = data.output.stamp(dae, output_state_[end], delayed + (1 - end) * data.K);
        const auto outgoing = data.input.stamp(dae, input_state_[end], reflected[end]);
        current[end].resize(data.N);
        for (size_t k = 0; k < data.K; ++k)
        {
          dae.add(reflected[end] + k, reflected[end] + k, -1);
          dae.add(reflected[end] + k, characteristic[end][k], 2);
          dae.add(reflected[end] + k, incident[end][k], -1);
          dae.add(modal[end] + k, modal[end] + k, 1);
          dae.add(modal[end] + k, outgoing[k], -1);
          auto& injection = current[end][data.conductors[k]];
          injection.insert(injection.end(), incident[end][k].begin(), incident[end][k].end());
          for (const auto& [index, coefficient] : characteristic[end][k])
          {
            injection.push_back({index, -coefficient});
          }
        }
      }
    }

    /// Diagnostic name of a local variable owned by the line, or an empty string.
    std::string name(size_t local) const
    {
      for (size_t end = 0; end < 2; ++end)
      {
        const std::string e = std::to_string(end + 1);
        const struct
        {
          const char* label;
          size_t      first, count;
        } blocks[] = {{"ref", reflected[end], data.K},
                      {"modal", modal[end], data.K},
                      {"yc_state", yc_state_[end], data.K * data.yc.poles.size()},
                      {"hin_state", input_state_[end], data.K * data.input.poles.size()},
                      {"hout_state", output_state_[end], data.K * data.output.poles.size()}};
        for (const auto& block : blocks)
        {
          if (local >= block.first && local < block.first + block.count)
          {
            const size_t offset = local - block.first;
            return std::string(block.label) + e + "_p" + std::to_string(offset / data.K) + "_k" + std::to_string(offset % data.K + 1);
          }
        }
      }
      return "";
    }

  private:
    std::array<size_t, 2> yc_state_, input_state_, output_state_;
  };
} // namespace GridKit::Examples::EMT
