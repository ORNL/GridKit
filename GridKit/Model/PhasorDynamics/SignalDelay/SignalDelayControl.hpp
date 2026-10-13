#pragma once

#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>

namespace GridKit::PhasorDynamics
{
  /// Attach history control before configuring the forward solver.
  /// The system must outlive the solver and its callbacks.
  template <class ScalarT, class IdxT>
  void configureSignalDelays(SystemModel<ScalarT, IdxT>&                                      system,
                             AnalysisManager::Sundials::Ida<ScalarT, IdxT>&                   solver,
                             SignalHistoryOptions<typename SystemModel<ScalarT, IdxT>::RealT> options)
  {
    if (!system.hasSignalDelays())
    {
      return;
    }
    using RealT   = typename SystemModel<ScalarT, IdxT>::RealT;
    using SolverT = AnalysisManager::Sundials::Ida<ScalarT, IdxT>;
    system.setSignalHistoryOptions(options);
    solver.setStepControl({[&system](RealT t, bool consistent)
                           { system.initializeSignalHistory(t, consistent); },
                           [&system](RealT begin, RealT end, const typename SolverT::DenseOutput& dense)
                           { return system.assessSignalHistory(begin, end, dense); },
                           [&system]()
                           { system.acceptSignalHistory(); },
                           [&system]()
                           { system.finishSignalHistory(); },
                           [&system]()
                           { system.discardSignalHistory(); },
                           4}); // Cubic reconstruction has fourth-order local error.
  }
} // namespace GridKit::PhasorDynamics
