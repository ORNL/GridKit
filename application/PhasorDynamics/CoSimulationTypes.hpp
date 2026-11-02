#pragma once

#include <GridKit/Definitions.hpp>

#ifdef GRIDKIT_ENABLE_ZMQ

#include <functional>
#include <optional>

#include <GridKit/CoSim/Message.hpp>

#include "AnalysisUtilities.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Alias for the co-simulation message type used in the applications.
    using CoSimulationMessage = CoSim::MessageData<double>;

    using NextCoupledStudy     = std::function<std::optional<StudyData>()>;
    using CompleteCoupledStudy = std::function<void(Testing::TestStatus)>;
  } // namespace PhasorDynamics
} // namespace GridKit
#endif
