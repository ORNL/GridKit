#pragma once

#include <numbers>

#include "LineDistributed.hpp"

namespace GridKit::Examples::EMT
{
  inline std::array<double, 3> phaseValues(const Json& j)
  {
    if (!j.is_array() || j.size() != 3)
    {
      throw std::invalid_argument("Expected three phase values");
    }
    return {real(j[0]), real(j[1]), real(j[2])};
  }

  inline double transition(double t, double start, double duration)
  {
    const double x = std::clamp((t - start) / duration, 0.0, 1.0);
    return x * x * x * (10 + x * (-15 + 6 * x));
  }

  struct Event
  {
    std::string         type;
    double              on, off, rise, fall, value;
    std::array<bool, 3> phases{false, false, false};

    double envelope(double t) const
    {
      return transition(t, on, rise) - transition(t, off, fall);
    }
  };

  struct CircuitData
  {
    LineData              line;
    std::array<double, 3> rms, phase, rs, resistance, inductance;
    double                frequency, t0, tf, dt, rtol, history_atol, hmax;
    bool                  compress, steady;
    std::vector<Event>    events;

    static CircuitData parse(const Json& j, const std::filesystem::path& directory)
    {
      fields(j, {"name", "line", "source", "load", "initialization", "events", "solver"});
      CircuitData d;
      d.line = LineData::parse(j.at("line"), directory);
      if (d.line.N != 3)
      {
        throw std::invalid_argument("The demonstration circuit requires N = 3");
      }
      const auto& source = j.at("source");
      fields(source, {"rms", "phase", "frequency", "resistance"});
      d.rms            = phaseValues(source.at("rms"));
      d.phase          = phaseValues(source.at("phase"));
      d.rs             = phaseValues(source.at("resistance"));
      d.frequency      = real(source.at("frequency"));
      const auto& load = j.at("load");
      fields(load, {"resistance", "inductance"});
      d.resistance                     = phaseValues(load.at("resistance"));
      d.inductance                     = phaseValues(load.at("inductance"));
      const std::string initialization = j.at("initialization");
      if (initialization != "zero" && initialization != "sinusoidal")
      {
        throw std::invalid_argument("Initialization must be zero or sinusoidal");
      }
      d.steady           = initialization == "sinusoidal";
      const auto& solver = j.at("solver");
      fields(solver, {"t0", "tf", "dt", "rel_tol", "history_abs_tol", "h_max", "compress"});
      d.t0           = real(solver.at("t0"));
      d.tf           = real(solver.at("tf"));
      d.dt           = real(solver.at("dt"));
      d.rtol         = real(solver.at("rel_tol"));
      d.history_atol = real(solver.at("history_abs_tol"));
      d.hmax         = real(solver.at("h_max"));
      d.compress     = solver.at("compress").get<bool>();
      if (!(d.frequency > 0 && d.tf > d.t0 && d.dt > 0 && d.rtol > 0 && d.rtol < 1 && d.history_atol > 0 && d.hmax >= 0))
      {
        throw std::invalid_argument("Invalid source frequency or solver settings");
      }
      for (size_t p = 0; p < 3; ++p)
      {
        if (!(d.rms[p] >= 0 && d.rs[p] > 0 && d.resistance[p] > 0 && d.inductance[p] >= 0))
        {
          throw std::invalid_argument("Require RMS >= 0, Rs/Rload > 0 and Lload >= 0");
        }
      }
      if (!j.at("events").is_array())
      {
        throw std::invalid_argument("Events must be an array");
      }
      for (const auto& item : j.at("events"))
      {
        fields(item, {"type", "on", "off", "rise", "fall", "value", "phases"});
        Event event;
        event.type  = item.at("type");
        event.on    = real(item.at("on"));
        event.off   = real(item.at("off"));
        event.rise  = real(item.at("rise"));
        event.fall  = real(item.at("fall"));
        event.value = real(item.at("value"));
        if (!(event.on >= d.t0 && event.rise > 0 && event.fall > 0 && event.off >= event.on + event.rise))
        {
          throw std::invalid_argument("Events require positive transition times and off >= on + rise");
        }
        if (event.type != "source" && event.type != "fault")
        {
          throw std::invalid_argument("Event type must be source or fault");
        }
        if (event.type == "fault" && event.value < 0)
        {
          throw std::invalid_argument("Fault conductance must be nonnegative");
        }
        if (!item.at("phases").is_array() || item.at("phases").empty())
        {
          throw std::invalid_argument("Event phases must be a nonempty array");
        }
        for (const auto& phase : item.at("phases"))
        {
          const size_t p = dimension(phase) - 1;
          if (p >= 3 || event.phases[p])
          {
            throw std::invalid_argument("Invalid or duplicate event phase");
          }
          event.phases[p] = true;
        }
        d.events.push_back(event);
      }
      return d;
    }

    double omega() const
    {
      return 2 * std::numbers::pi * frequency;
    }

    double source(double t, size_t phase_index) const
    {
      double amplitude = steady ? 1 : 0;
      for (const auto& event : events)
      {
        if (event.type == "source" && event.phases[phase_index])
        {
          amplitude += event.value * event.envelope(t);
        }
      }
      return amplitude * std::sqrt(2.0) * rms[phase_index] * std::cos(omega() * t + phase[phase_index]);
    }

    double fault(double t, size_t phase_index) const
    {
      double conductance = 0;
      for (const auto& event : events)
      {
        if (event.type == "fault" && event.phases[phase_index])
        {
          conductance += event.value * event.envelope(t);
        }
      }
      return conductance;
    }
  };
} // namespace GridKit::Examples::EMT
