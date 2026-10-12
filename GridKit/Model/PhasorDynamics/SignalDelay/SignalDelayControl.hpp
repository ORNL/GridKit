#pragma once

#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>

namespace GridKit::PhasorDynamics
{
  /// Attach history control before configuring the forward solver.
  /// The system must outlive the solver and its callbacks.
  template <class ScalarT, class IdxT>
  void configureSignalDelays(SystemModel<ScalarT, IdxT>&                    system,
                             AnalysisManager::Sundials::Ida<ScalarT, IdxT>& solver,
                             typename SystemModel<ScalarT, IdxT>::RealT     relative,
                             typename SystemModel<ScalarT, IdxT>::RealT     absolute,
                             bool                                           compress = true)
  {
    if (!system.hasSignalDelays())
    {
      return;
    }
    using RealT   = typename SystemModel<ScalarT, IdxT>::RealT;
    using SolverT = AnalysisManager::Sundials::Ida<ScalarT, IdxT>;
    system.setSignalHistoryTolerance(relative, absolute, compress);
    solver.setStepControl({[&system](RealT t, bool consistent)
                           { system.initializeSignalHistory(t, consistent); },
                           [&system](RealT begin, RealT end, const typename SolverT::DenseOutput& dense)
                           { return system.assessSignalHistory(begin, end, dense); },
                           [&system]()
                           { system.acceptSignalHistory(); },
                           [&system]()
                           { system.finishSignalHistory(); },
                           [&system]()
                           { system.discardSignalHistory(); }});
  }
} // namespace GridKit::PhasorDynamics
