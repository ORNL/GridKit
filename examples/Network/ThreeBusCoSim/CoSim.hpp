#pragma once

namespace GridKit
{
  namespace CoSim
  {
    /*
     * The following is meant to be used to communicate driver information over
     * the network. It is crude, but it at least gives language in the code that
     * expresses the meaning of the intent.
     *
     * - STEP is used to initiate (and continue) the co-simulation (take a step)
     * - END is used to end the co-simulation
     *
     * Both of these are sent by the client to keep the server in step.
     */

    using StatusT = int;

    inline constexpr StatusT STEP = 1;
    inline constexpr StatusT END  = 0;
  } // namespace CoSim
} // namespace GridKit
