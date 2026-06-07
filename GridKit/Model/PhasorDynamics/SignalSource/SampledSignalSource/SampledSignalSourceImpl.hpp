#pragma once

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <variant>

#include <GridKit/Model/PhasorDynamics/SignalSource/SampledSignalSource/SampledSignalSource.hpp>
#include <GridKit/Model/PhasorDynamics/SignalSource/SampledSignalSource/SampledSignalSourceData.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace SignalSource
    {
      using Log = ::GridKit::Utilities::Logger;

      namespace Detail
      {
        inline std::string trim(const std::string& text)
        {
          const auto begin = text.find_first_not_of(" \t\r\n");
          if (begin == std::string::npos)
          {
            return {};
          }
          const auto end = text.find_last_not_of(" \t\r\n");
          return text.substr(begin, end - begin + 1);
        }

        inline std::vector<std::string> splitCsvLine(const std::string& line)
        {
          std::vector<std::string> values;
          std::stringstream        ss(line);
          std::string              value;
          while (std::getline(ss, value, ','))
          {
            values.push_back(trim(value));
          }
          return values;
        }
      } // namespace Detail

      template <typename scalar_type, typename index_type>
      SampledSignalSource<scalar_type, index_type>::SampledSignalSource()
        : monitor_(std::make_unique<MonitorT>())
      {
      }

      template <typename scalar_type, typename index_type>
      SampledSignalSource<scalar_type, index_type>::SampledSignalSource(const ModelDataT& data)
        : monitor_(std::make_unique<MonitorT>(data))
      {
        initializeParameters(data);
        loadSamples(data);
        initializeMonitor();
      }

      template <typename scalar_type, typename index_type>
      SampledSignalSource<scalar_type, index_type>::~SampledSignalSource()
      {
      }

      template <typename scalar_type, typename index_type>
      int SampledSignalSource<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int SampledSignalSource<scalar_type, index_type>::allocate()
      {
        size_ = 0;
        nnz_  = 0;
        y_.clear();
        yp_.clear();
        f_.clear();
        tag_.clear();
        variable_indices_.clear();
        residual_indices_.clear();

        if (signals_.template isAssigned<SampledSignalSourceInternalVariables::OUT>())
        {
          signals_.template getSignalNode<SampledSignalSourceInternalVariables::OUT>()->set(
              &value_,
              &invalid_index_);
        }

        return 0;
      }

      template <typename scalar_type, typename index_type>
      int SampledSignalSource<scalar_type, index_type>::verify() const
      {
        int ret = configuration_errors_;

        if (!signals_.template isAssigned<SampledSignalSourceInternalVariables::OUT>())
        {
          Log::error() << "SampledSignalSource: required output signal is not assigned\n";
          ret += 1;
        }
        else if (!signals_.template getSignalNode<SampledSignalSourceInternalVariables::OUT>()->linked())
        {
          Log::error() << "SampledSignalSource: output signal is assigned but not linked\n";
          ret += 1;
        }

        if (samples_.empty())
        {
          Log::error() << "SampledSignalSource: at least one sample is required\n";
          ret += 1;
        }

        return ret;
      }

      template <typename scalar_type, typename index_type>
      int SampledSignalSource<scalar_type, index_type>::initialize()
      {
        updateTime(time_, 0.0);
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int SampledSignalSource<scalar_type, index_type>::tagDifferentiable()
      {
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int SampledSignalSource<scalar_type, index_type>::evaluateResidual()
      {
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int SampledSignalSource<scalar_type, index_type>::evaluateJacobian()
      {
        return 0;
      }

      template <typename scalar_type, typename index_type>
      void SampledSignalSource<scalar_type, index_type>::updateTime(RealT t, RealT a)
      {
        Component<ScalarT, IdxT>::updateTime(t, a);
        value_ = ScalarT{valueAt(t)};
      }

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* SampledSignalSource<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      template <typename scalar_type, typename index_type>
      auto SampledSignalSource<scalar_type, index_type>::sample(RealT t) const -> RealT
      {
        if (samples_.empty())
        {
          return 0.0;
        }

        if (t <= samples_.front().first)
        {
          return samples_.front().second;
        }
        if (t >= samples_.back().first)
        {
          return samples_.back().second;
        }

        auto upper = std::lower_bound(samples_.begin(),
                                      samples_.end(),
                                      t,
                                      [](const auto& sample, RealT time)
                                      {
                                        return sample.first < time;
                                      });

        if (upper != samples_.end() && upper->first == t)
        {
          return upper->second;
        }

        auto       lower = std::prev(upper);
        const auto span  = upper->first - lower->first;
        const auto theta = (t - lower->first) / span;
        return lower->second + theta * (upper->second - lower->second);
      }

      template <typename scalar_type, typename index_type>
      auto SampledSignalSource<scalar_type, index_type>::valueAt(RealT t) const -> RealT
      {
        return scale_ * sample(t) + offset_;
      }

      template <typename scalar_type, typename index_type>
      void SampledSignalSource<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        auto read_real = [](const auto& value)
        {
          return std::visit([](const auto& stored) -> RealT
                            { return static_cast<RealT>(stored); },
                            value);
        };

        if (data.parameters.contains(ModelDataT::Parameters::scale))
        {
          scale_ = read_real(data.parameters.at(ModelDataT::Parameters::scale));
        }
        if (data.parameters.contains(ModelDataT::Parameters::offset))
        {
          offset_ = read_real(data.parameters.at(ModelDataT::Parameters::offset));
        }
      }

      template <typename scalar_type, typename index_type>
      void SampledSignalSource<scalar_type, index_type>::loadSamples(const ModelDataT& data)
      {
        if (data.source_type == "csv")
        {
          loadCsvSamples(data);
        }
        else
        {
          samples_ = data.samples;
        }
        validateSamples();
      }

      template <typename scalar_type, typename index_type>
      void SampledSignalSource<scalar_type, index_type>::loadCsvSamples(const ModelDataT& data)
      {
        std::ifstream file(data.file);
        if (!file)
        {
          Log::error() << "SampledSignalSource: could not open CSV file: " << data.file << "\n";
          configuration_errors_ += 1;
          return;
        }

        std::string line;
        if (!std::getline(file, line))
        {
          Log::error() << "SampledSignalSource: CSV file is empty: " << data.file << "\n";
          configuration_errors_ += 1;
          return;
        }

        const auto header   = Detail::splitCsvLine(line);
        auto       time_it  = std::find(header.begin(), header.end(), data.time_column);
        auto       value_it = std::find(header.begin(), header.end(), data.value_column);
        if (time_it == header.end() || value_it == header.end())
        {
          Log::error() << "SampledSignalSource: CSV columns \"" << data.time_column
                       << "\" and \"" << data.value_column << "\" were not both found in "
                       << data.file << "\n";
          configuration_errors_ += 1;
          return;
        }

        const auto time_column  = static_cast<size_t>(std::distance(header.begin(), time_it));
        const auto value_column = static_cast<size_t>(std::distance(header.begin(), value_it));

        while (std::getline(file, line))
        {
          if (Detail::trim(line).empty())
          {
            continue;
          }

          const auto values = Detail::splitCsvLine(line);
          if (time_column >= values.size() || value_column >= values.size())
          {
            Log::error() << "SampledSignalSource: malformed CSV row in " << data.file
                         << ": " << line << "\n";
            configuration_errors_ += 1;
            continue;
          }

          try
          {
            samples_.emplace_back(static_cast<RealT>(std::stod(values[time_column])),
                                  static_cast<RealT>(std::stod(values[value_column])));
          }
          catch (const std::exception&)
          {
            Log::error() << "SampledSignalSource: nonnumeric CSV row in " << data.file
                         << ": " << line << "\n";
            configuration_errors_ += 1;
          }
        }
      }

      template <typename scalar_type, typename index_type>
      void SampledSignalSource<scalar_type, index_type>::validateSamples()
      {
        if (samples_.empty())
        {
          return;
        }

        for (size_t i = 1; i < samples_.size(); ++i)
        {
          if (samples_[i].first <= samples_[i - 1].first)
          {
            Log::error() << "SampledSignalSource: sample times must be strictly increasing\n";
            configuration_errors_ += 1;
            return;
          }
        }
      }

      template <typename scalar_type, typename index_type>
      void SampledSignalSource<scalar_type, index_type>::initializeMonitor()
      {
        using Variable = typename ModelDataT::MonitorableVariables;
        monitor_->set(Variable::out, [this]
                      { return value_; });
      }
    } // namespace SignalSource
  } // namespace PhasorDynamics
} // namespace GridKit
