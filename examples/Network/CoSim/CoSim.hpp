#pragma once

namespace GridKit
{
  namespace CoSim
  {
    /*
     * The following is meant to be used to communicate simulation driver
     * information over the network. It is crude, but it at least gives language
     * in the code that expresses the intent.
     *
     * - INIT is used to initialize the co-simulation
     * - STEP is used to initiate (and continue) the co-simulation (take a step)
     * - END is used to end the co-simulation
     * - FAIL is used to signal an early termination due to a failure
     */
    enum class Status
    {
      INIT,
      STEP,
      END,
      FAIL
    };

    /**
     * Integer representation of Status
     */
    using StatusRep = std::underlying_type_t<Status>;

    /**
     * @brief Write a Status to a stream as an integer
     */
    std::ostream& operator<<(std::ostream& os, Status stat)
    {
      os << static_cast<StatusRep>(stat);
      return os;
    }

    /**
     * @brief Read a Status from a stream, in which it is represented as an
     * integer
     */
    std::istream& operator>>(std::istream& is, Status& stat)
    {
      StatusRep rep;
      is >> rep;
      stat = static_cast<Status>(rep);
      return is;
    }
  } // namespace CoSim
} // namespace GridKit
