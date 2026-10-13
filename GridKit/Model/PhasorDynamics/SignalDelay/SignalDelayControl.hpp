#pragma once

#include <algorithm>
#include <array>
#include <memory>
#include <optional>
#include <sstream>
#include <vector>

#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>

namespace GridKit::PhasorDynamics
{
  /// IDAS binding for source histories. Bindings and sample buffers are allocated
  /// once; the model and its signal topology must outlive this control.
  template <typename scalar_type, typename index_type>
  class SignalHistoryControl
  {
    using System     = SystemModel<scalar_type, index_type>;
    using Solver     = AnalysisManager::Sundials::Ida<scalar_type, index_type>;
    using RealT      = typename System::RealT;
    using IdxT       = index_type;
    using History    = SignalHistory<RealT>;
    using Source     = typename System::SourceHistory;
    using Connection = typename System::DelayedConnection;

    struct SourceBinding
    {
      IdxT    id;
      Source* source;
    };

    struct TapBinding
    {
      const Connection*    connection;
      History*             history;
      std::array<RealT, 6> values{};
    };

    enum class Check
    {
      SOURCE_RECONSTRUCTION,
      DELAYED_OUTPUT
    };

    struct Assessment
    {
      RealT               error{0};
      IdxT                source{};
      std::optional<IdxT> signal;
      Check               check{Check::SOURCE_RECONSTRUCTION};
    };

  public:
    SignalHistoryControl(System& system, SignalHistoryOptions<RealT> options)
      : system_(system)
    {
      for (auto& [id, source] : system_.signal_histories_)
      {
        source.history.setOptions(options);
        sources_.push_back({id, &source});
      }
      for (const auto& connection : system_.delayed_connections_)
      {
        if (connection.component)
        {
          taps_.push_back({&connection, &system_.signal_histories_.at(connection.source_id).history, {}});
        }
      }
    }

    void initialize(RealT t, bool consistent)
    {
      for (auto& binding : sources_)
      {
        auto& source = *binding.source;
        if (consistent)
        {
          source.history.initialize(t, static_cast<RealT>(source.source->read()));
        }
        else
        {
          source.history.begin(t);
        }
      }
      if (!consistent)
      {
        for (const auto& tap : taps_)
        {
          if (!tap.history->initialized())
          {
            // The actual initial time can differ from configureSimulation's zero.
            tap.connection->component->initialize();
          }
        }
        system_.y().setDataUpdated();
        system_.yp().setDataUpdated();
      }
      state_sources_ = std::all_of(sources_.begin(), sources_.end(), [](const auto& binding)
                                   { return binding.source->source->getVariableIndex() != INVALID_INDEX<IdxT>; });
    }

    typename Solver::StepAssessment assess(RealT begin, RealT end, const typename Solver::DenseOutput& dense)
    {
      assessment_ = {};
      begin_      = begin;
      end_        = end;
      IdxT source_id{};
      try
      {
        unsigned order = 4;
        for (auto& binding : sources_)
        {
          auto& source = *binding.source;
          source_id    = binding.id;
          source.history.prepare(end, static_cast<RealT>(source.source->read()));
          order = std::min(order, source.history.errorOrder());
        }
        // Six samples reconstruct IDAS's degree <= 5 polynomial. A seventh checks
        // prescribed sources at the midpoint. Computed sources retain the full
        // model-state fallback; state variables need only y interpolation.
        for (size_t i = 0; i <= 6; ++i)
        {
          const RealT            t = i == 6 ? begin + (end - begin) / 2
                                            : begin + (end - begin) * History::dense_nodes[i];
          std::span<const RealT> values;
          if (state_sources_)
          {
            values = dense.values(t);
          }
          else
          {
            dense(t);
          }
          for (const auto& binding : sources_)
          {
            const auto& source = *binding.source;
            const RealT value  = state_sources_ ? values[static_cast<size_t>(source.source->getVariableIndex())]
                                                : static_cast<RealT>(source.source->read());
            source_id          = binding.id;
            const RealT error  = source.history.estimate(t, value);
            if (error > assessment_.error)
            {
              assessment_ = {error, binding.id, {}, Check::SOURCE_RECONSTRUCTION};
            }
          }
          if (assessment_.error > 1)
          {
            return {assessment_.error, 1}; // Startup and corners need conservative rejection.
          }
          if (i < 6)
          {
            for (auto& tap : taps_)
            {
              const auto* output = tap.connection->output;
              tap.values[i]      = state_sources_ ? values[static_cast<size_t>(output->getVariableIndex())]
                                                  : static_cast<RealT>(output->read());
            }
          }
        }
        for (const auto& tap : taps_)
        {
          const auto& connection = *tap.connection;
          source_id              = connection.source_id;
          const RealT error      = tap.history->estimateDelayed(begin, end, connection.delay, tap.values);
          if (error > assessment_.error)
          {
            assessment_ = {error, connection.source_id, connection.output->signalId(), Check::DELAYED_OUTPUT};
          }
        }
        return {assessment_.error, assessment_.error > 1 ? 1 : order};
      }
      catch (const std::exception& error)
      {
        std::ostringstream message;
        message << "Signal history check failed for source=" << source_id
                << " on [" << begin << ", " << end << "]: " << error.what();
        throw std::runtime_error(message.str());
      }
    }

    std::string describe() const
    {
      std::ostringstream message;
      message << (assessment_.check == Check::SOURCE_RECONSTRUCTION ? "source reconstruction" : "delayed output")
              << ", source=" << assessment_.source
              << ", interval=[" << begin_ << ", " << end_ << "]";
      if (assessment_.signal)
      {
        message << ", signal=" << *assessment_.signal;
      }
      return message.str();
    }

    void accept()
    {
      for (auto& binding : sources_)
      {
        binding.source->history.accept();
      }
    }

    void finish()
    {
      for (auto& binding : sources_)
      {
        binding.source->history.maintain(binding.source->delays);
      }
    }

    void discard()
    {
      for (auto& binding : sources_)
      {
        binding.source->history.discard();
      }
    }

  private:
    System&                    system_;
    std::vector<SourceBinding> sources_;
    std::vector<TapBinding>    taps_;
    Assessment                 assessment_;
    RealT                      begin_{0};
    RealT                      end_{0};
    bool                       state_sources_{false};
  };

  /// Attach before configuring the forward solver. The system must outlive it.
  template <typename scalar_type, typename index_type>
  void configureSignalDelays(SystemModel<scalar_type, index_type>&                                      system,
                             AnalysisManager::Sundials::Ida<scalar_type, index_type>&                   solver,
                             SignalHistoryOptions<typename SystemModel<scalar_type, index_type>::RealT> options)
  {
    if (!system.hasSignalDelays())
    {
      return;
    }
    using Control = SignalHistoryControl<scalar_type, index_type>;
    using RealT   = typename SystemModel<scalar_type, index_type>::RealT;
    using Solver  = AnalysisManager::Sundials::Ida<scalar_type, index_type>;
    auto control  = std::make_shared<Control>(system, options);
    solver.setStepControl({[control](RealT t, bool consistent)
                           { control->initialize(t, consistent); },
                           [control](RealT begin, RealT end, const typename Solver::DenseOutput& dense)
                           { return control->assess(begin, end, dense); },
                           [control]()
                           { control->accept(); },
                           [control]()
                           { control->finish(); },
                           [control]()
                           { control->discard(); },
                           4,
                           [control]()
                           { return control->describe(); }});
  }
} // namespace GridKit::PhasorDynamics
