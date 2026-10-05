#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/BusFault/BusFault.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelDataJSONParser.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace
{
  using RealT = double;
  using IdxT  = std::size_t;
  using namespace GridKit::PhasorDynamics;
  using SystemT = SystemModel<RealT, IdxT>;
  using DataT   = SystemModelData<RealT, IdxT>;
  using IdaT    = AnalysisManager::Sundials::Ida<RealT, IdxT>;

  struct Error
  {
    RealT voltage{};
    RealT state{};
    RealT speed{};
    RealT angle{};
  };

  void require(bool condition, const char* message)
  {
    if (!condition)
    {
      throw std::runtime_error(message);
    }
  }

  DataT readCase(const json& document)
  {
    auto data = document.get<DataT>();
    for (const auto& bus : data.bus)
    {
      require(bus.bus_type == BusData<RealT, IdxT>::BusType::DEFAULT,
              "Partition test requires ordinary buses");
    }
    return data;
  }

  void checkPartition(const json& a, const json& b, const json& reference)
  {
    auto records = [](const json& objects, const char* key)
    {
      auto result = json::object();
      for (const auto& object : objects)
      {
        if (object.at("class") == "BusToSignalAdapter")
        {
          continue;
        }
        const auto id = object.at(key).dump();
        require(!result.contains(id), "Duplicate case record");
        result[id] = object;
      }
      return result;
    };
    for (const auto* group : {"buses", "devices"})
    {
      const bool  buses    = std::string(group) == "buses";
      const auto* key      = buses ? "number" : "id";
      auto        combined = records(a.at(group), key);
      const auto  other    = records(b.at(group), key);
      for (const auto& [id, record] : other.items())
      {
        require(!combined.contains(id) || (buses && combined.at(id) == record),
                "Partitions duplicate a device or disagree on a bus");
        combined[id] = record;
      }
      require(combined == records(reference.at(group), key),
              "Partitions do not reconstruct the intact case");
    }
  }

  void bindCurrent(SystemT& system, const DataT& data, RealT& ir, RealT& ii, IdxT& index)
  {
    const auto& inputs = data.adapter.at(0).signal_inputs;
    system.getSignalNode(inputs.at(BusToSignalAdapterSignalInputs::ir))->link(&ir, &index);
    system.getSignalNode(inputs.at(BusToSignalAdapterSignalInputs::ii))->link(&ii, &index);
  }

  std::vector<RealT> differentialStates(SystemT& system)
  {
    std::vector<RealT> states;
    for (IdxT i = 0; i < system.size(); ++i)
    {
      if (system.tag()[i])
      {
        states.push_back(system.y().getData()[i]);
      }
    }
    return states;
  }

  RealT maxDifference(const std::vector<RealT>& a, const std::vector<RealT>& b)
  {
    require(a.size() == b.size(), "State counts differ");
    RealT error = 0.0;
    for (IdxT i = 0; i < a.size(); ++i)
    {
      const RealT difference = std::abs(a[i] - b[i]);
      require(std::isfinite(difference), "Nonfinite solution");
      error = std::max(error, difference);
    }
    return error;
  }

  void compareBuses(SystemT& system, const DataT& data, SystemT& reference, Error& error)
  {
    system.evaluateResidual();
    for (const auto& bus : data.bus)
    {
      auto* actual   = system.getBus(bus.bus_id);
      auto* expected = reference.getBus(bus.bus_id);
      require(std::hypot(actual->Ir(), actual->Ii()) < 1e-8, "Partition violates KCL");
      const RealT difference = std::hypot(actual->Vr() - expected->Vr(), actual->Vi() - expected->Vi());
      require(std::isfinite(difference), "Nonfinite bus voltage");
      error.voltage = std::max(error.voltage, difference);
    }
  }

  void compareStates(SystemT& a, SystemT& b, SystemT& reference, Error& error)
  {
    auto       states   = differentialStates(a);
    const auto states_b = differentialStates(b);
    states.insert(states.end(), states_b.begin(), states_b.end());
    const auto expected = differentialStates(reference);
    error.state         = std::max(error.state, maxDifference(states, expected));
    // These fixtures contain only GENROU differential states, in verified ID order.
    for (IdxT i = 0; i < states.size(); i += 6)
    {
      error.speed = std::max(error.speed, std::abs(states[i + 1] - expected[i + 1]));
      error.angle = std::max(error.angle, std::abs((states[i] - states[0]) - (expected[i] - expected[0])));
    }
  }

  void replay(IdaT& ida, RealT t0, RealT t1)
  {
    ida.getSavedInitialCondition();
    ida.initializeSimulation(t0);
    if (t1 > t0)
    {
      ida.runSimulation(t1);
    }
  }

  // Match the two terminal voltages by solving for real and imaginary current.
  template <typename Residual>
  bool solveInterface(Residual&& residual, RealT& ir, RealT& ii)
  {
    for (IdxT iteration = 0; iteration < 20; ++iteration)
    {
      RealT vr, vi;
      residual(ir, ii, vr, vi);
      if (!std::isfinite(vr) || !std::isfinite(vi))
      {
        return false;
      }
      if (std::hypot(vr, vi) < 1e-8)
      {
        return true;
      }

      const RealT h = 1e-5 * std::max(RealT{1}, std::hypot(ir, ii));
      RealT       vr_r, vi_r, vr_i, vi_i;
      residual(ir + h, ii, vr_r, vi_r);
      residual(ir, ii + h, vr_i, vi_i);
      const RealT dvr_dir = (vr_r - vr) / h, dvi_dir = (vi_r - vi) / h;
      const RealT dvr_dii = (vr_i - vr) / h, dvi_dii = (vi_i - vi) / h;
      const RealT det = dvr_dir * dvi_dii - dvr_dii * dvi_dir;
      if (!std::isfinite(det) || std::abs(det) <= 1e-12 * std::hypot(dvr_dir, dvi_dir) * std::hypot(dvr_dii, dvi_dii))
      {
        return false;
      }
      ir -= (dvi_dii * vr - dvr_dii * vi) / det;
      ii -= (dvr_dir * vi - dvi_dir * vr) / det;
    }
    return false;
  }

  // Save once per accepted starting point; all trials and retries restore this state.
  template <typename Couple>
  RealT advance(IdaT& a, IdaT& b, Couple&& couple, RealT t0, RealT t1, RealT& ir, RealT& ii)
  {
    a.saveInitialCondition();
    b.saveInitialCondition();
    const RealT initial_ir = ir, initial_ii = ii;
    for (IdxT retry = 0; retry < 9 && t1 > t0; ++retry)
    {
      try
      {
        if (couple(t0, t1))
        {
          return t1;
        }
      }
      catch (const AnalysisManager::Sundials::SundialsException&)
      {
        // A failed solver trial leaves the saved starting point intact.
      }
      ir = initial_ir;
      ii = initial_ii;
      t1 = std::midpoint(t0, t1);
    }
    require(couple(t0, t0), "Could not restore the accepted coupling state");
    throw std::runtime_error("Partition coupling failed after step reduction");
  }

  Error runPartitioned(const DataT& a_data, const DataT& b_data, const DataT& reference_data, RealT dt, RealT tf)
  {
    RealT   ir_a = 0.0, ii_a = 0.0, ir_b = 0.0, ii_b = 0.0;
    IdxT    index = GridKit::INVALID_INDEX<IdxT>;
    SystemT a(a_data), b(b_data), reference(reference_data);
    bindCurrent(a, a_data, ir_a, ii_a, index);
    bindCurrent(b, b_data, ir_b, ii_b, index);
    require(a.allocate() == 0 && b.allocate() == 0 && reference.allocate() == 0, "Model allocation failed");
    require(a.hasJacobian() && b.hasJacobian() && reference.hasJacobian(),
            "Partitioned integration requires model Jacobians");
    IdaT a_ida(&a), b_ida(&b), reference_ida(&reference);
    for (auto* ida : {&a_ida, &b_ida, &reference_ida})
    {
      ida->setTolerance(1e-10, 1e-12);
      ida->configureSimulation();
    }
    require(differentialStates(a).size() == 6 * a_data.genrou.size()
                && differentialStates(b).size() == 6 * b_data.genrou.size()
                && differentialStates(reference).size() == 6 * reference_data.genrou.size(),
            "Expected only GENROU differential states");

    auto* bus_a = a.getBus(a_data.adapter.front().buses.at(BusToSignalAdapterBuses::bus));
    auto* bus_b = b.getBus(b_data.adapter.front().buses.at(BusToSignalAdapterBuses::bus));
    a.evaluateResidual();
    RealT ir = -bus_a->Ir(), ii = -bus_a->Ii();
    auto  couple = [&](RealT t0, RealT t1)
    {
      auto mismatch = [&](RealT trial_ir, RealT trial_ii, RealT& vr, RealT& vi)
      {
        // Equal and opposite injections, on the common system base.
        ir_a = trial_ir;
        ii_a = trial_ii;
        ir_b = -trial_ir;
        ii_b = -trial_ii;
        replay(a_ida, t0, t1);
        replay(b_ida, t0, t1);
        vr = bus_a->Vr() - bus_b->Vr();
        vi = bus_a->Vi() - bus_b->Vi();
      };
      return solveInterface(mismatch, ir, ii);
    };

    Error error;
    auto  compare = [&]()
    {
      compareBuses(a, a_data, reference, error);
      compareBuses(b, b_data, reference, error);
      compareStates(a, b, reference, error);
    };
    auto synchronize = [&](RealT t)
    {
      a_ida.saveInitialCondition();
      b_ida.saveInitialCondition();
      require(couple(t, t), "Interface initialization did not converge");
    };
    auto run_interval = [&](RealT t0, RealT end)
    {
      RealT previous = t0;
      while (previous < end)
      {
        RealT t = std::min(previous + dt, end);
        if (end - t < 16 * std::numeric_limits<RealT>::epsilon() * std::max(RealT{1}, std::abs(end)))
        {
          t = end;
        }
        t = advance(a_ida, b_ida, couple, previous, t, ir, ii);
        reference_ida.runSimulation(t);
        compare();
        previous = t;
      }
    };
    SystemT& faulted   = a_data.bus_fault.empty() ? b : a;
    auto     set_fault = [&](bool active, RealT t)
    {
      const auto before_a = differentialStates(a), before_b = differentialStates(b);
      faulted.getBusFault(0)->setStatus(active);
      reference.getBusFault(0)->setStatus(active);
      synchronize(t);
      reference_ida.initializeSimulation(t);
      require(maxDifference(before_a, differentialStates(a)) < 1e-12
                  && maxDifference(before_b, differentialStates(b)) < 1e-12,
              "Fault transition changed differential states");
      compare();
    };

    synchronize(0.0);
    reference_ida.initializeSimulation(0.0);
    compare();
    require(error.voltage < 1e-8 && error.state < 1e-12, "Partitioned initialization differs from the intact circuit");
    for (IdxT i = 0; i < reference.size(); ++i)
    {
      if (reference.tag()[i])
      {
        require(std::abs(reference.yp().getData()[i]) < 1e-8, "Initial operating point is not stationary");
      }
    }
    const auto initial = differentialStates(reference);
    run_interval(0.0, 1.0);
    require(maxDifference(initial, differentialStates(reference)) < 1e-8, "Intact operating point drifted");
    set_fault(true, 1.0);
    // A rejected completed trial must reproduce a clean half-step from the checkpoint.
    a_ida.saveInitialCondition();
    b_ida.saveInitialCondition();
    const RealT half_step = 1.0 + 0.5 * dt;
    require(couple(1.0, half_step), "Clean retry reference failed");
    const auto expected_a = differentialStates(a), expected_b = differentialStates(b);
    require(couple(1.0, 1.0), "Retry reference restoration failed");
    bool reject = true;
    auto trial  = [&](RealT t0, RealT t1)
    {
      const bool converged = couple(t0, t1);
      if (reject)
      {
        reject = false;
        return false;
      }
      return converged;
    };
    const RealT accepted = advance(a_ida, b_ida, trial, 1.0, 1.0 + dt, ir, ii);
    require(std::abs(accepted - half_step) < 1e-14
                && maxDifference(expected_a, differentialStates(a)) < 1e-10
                && maxDifference(expected_b, differentialStates(b)) < 1e-10,
            "Rejected trial changed the accepted starting state");
    require(couple(1.0, 1.0), "Retry test restoration failed");
    run_interval(1.0, 1.1);
    set_fault(false, 1.1);
    run_interval(1.1, tf);
    std::cout << "h = " << dt << ": max voltage error = " << error.voltage
              << ", max state error = " << error.state
              << ", speed error = " << error.speed << ", relative angle error = " << error.angle << '\n';
    return error;
  }
} // namespace

int main(int argc, char** argv)
{
  GridKit::Testing::TestStatus success = true;
  try
  {
    require(argc >= 4 && argc <= 6, "Expected partition A, partition B, intact case paths, and optional end time and communication step");
    const auto a_json         = json::parse(std::ifstream(argv[1]));
    const auto b_json         = json::parse(std::ifstream(argv[2]));
    const auto reference_json = json::parse(std::ifstream(argv[3]));
    checkPartition(a_json, b_json, reference_json);
    const auto  a = readCase(a_json), b = readCase(b_json), reference = readCase(reference_json);
    const RealT tf = argc >= 5 ? std::stod(argv[4]) : 2.0;
    const RealT dt = argc == 6 ? std::stod(argv[5]) : 1.0 / 240.0;
    require(std::isfinite(dt) && dt > 0.0, "Communication step must be positive and finite");
    require(std::isfinite(tf) && tf > 1.1, "End time must follow fault clearing");
    require(a.adapter.size() == 1 && b.adapter.size() == 1 && reference.adapter.empty(),
            "Expected one interface per partition and none in the intact case");
    const auto boundary = a.adapter.front().buses.at(BusToSignalAdapterBuses::bus);
    require(boundary == b.adapter.front().buses.at(BusToSignalAdapterBuses::bus), "Interface bus IDs differ");
    for (const auto& bus_a : a.bus)
    {
      for (const auto& bus_b : b.bus)
      {
        require(bus_a.bus_id != bus_b.bus_id || bus_a.bus_id == boundary,
                "Only the interface bus may be shared");
      }
    }
    require(!a.genrou.empty() && !b.genrou.empty(), "Each partition requires dynamic generation");
    require(a.bus_fault.size() + b.bus_fault.size() == 1 && reference.bus_fault.size() == 1,
            "Expected one physical fault");
    IdxT generator = 0;
    for (const auto* data : {&a, &b})
    {
      for (const auto& gen : data->genrou)
      {
        require(gen.disambiguation_string == reference.genrou.at(generator++).disambiguation_string,
                "Partition generator order differs from the intact case");
      }
    }
    require(generator == reference.genrou.size(), "Generator counts differ");
    require(a.freq_base == reference.freq_base && b.freq_base == reference.freq_base
                && a.va_base == reference.va_base && b.va_base == reference.va_base,
            "Partition bases differ");
    const auto coarse  = runPartitioned(a, b, reference, dt, tf);
    const auto medium  = runPartitioned(a, b, reference, 0.5 * dt, tf);
    const auto fine    = runPartitioned(a, b, reference, 0.25 * dt, tf);
    success           *= medium.voltage < 0.7 * coarse.voltage && fine.voltage < 0.7 * medium.voltage;
    success           *= medium.state < 0.7 * coarse.state && fine.state < 0.7 * medium.state;
    success           *= medium.speed < 0.7 * coarse.speed && fine.speed < 0.7 * medium.speed;
    success           *= medium.angle < 0.7 * coarse.angle && fine.angle < 0.7 * medium.angle;
    success           *= fine.voltage < 1e-3 && fine.state < 1e-3 && fine.angle < 1e-3 && fine.speed < 1e-5;
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    success = false;
  }
  GridKit::Testing::TestingResults result;
  result += success.report("partitionedIntegration");
  return result.summary();
}
