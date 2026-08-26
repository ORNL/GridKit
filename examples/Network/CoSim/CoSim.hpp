#pragma once

namespace GridKit
{
  namespace CoSim
  {
    enum class Status
    {
      INIT,
      STEP,
      END,
      FAIL
    };

    std::ostream& operator<<(std::ostream& os, Status stat)
    {
      os << static_cast<std::underlying_type_t<Status>>(stat);
      return os;
    }

    std::istream& operator>>(std::istream& is, Status& stat)
    {
      std::underlying_type_t<Status> rep;
      is >> rep;
      stat = static_cast<Status>(rep);
      return is;
    }
  } // namespace CoSim
} // namespace GridKit
