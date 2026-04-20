/**
 * @file BranchBase.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Abstract base class for phasor-dynamics branch models.
 */

#pragma once

#include <vector>

#include <GridKit/Constants.hpp>
#include <GridKit/Model/PhasorDynamics/Component.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <class ScalarT, typename IdxT>
    class BusBase;

    /**
     * @brief Abstract base class for phasor-dynamics branch models.
     *
     * Provides the common structural contract for branches: `BusBase*` ports,
     * a branch ID, and a `verifyPorts()` hook for concrete subclasses to
     * assert per-port width compatibility.
     */
    template <class ScalarT, typename IdxT>
    class BranchBase : public Component<ScalarT, IdxT>
    {
    public:
      using bus_type = BusBase<ScalarT, IdxT>;

      BranchBase()          = default;
      virtual ~BranchBase() = default;

      bus_type* port(IdxT k)
      {
        return ports_[static_cast<size_t>(k)];
      }

      const bus_type* port(IdxT k) const
      {
        return ports_[static_cast<size_t>(k)];
      }

      IdxT portCount() const
      {
        return static_cast<IdxT>(ports_.size());
      }

      virtual int setBranchID(IdxT) = 0;

      virtual const IdxT branchID() const
      {
        return branch_id_;
      }

      /// Called from allocate() to assert per-port width matches the concrete
      /// branch's physics. Throws std::runtime_error on mismatch.
      virtual int verifyPorts() const = 0;

    protected:
      std::vector<bus_type*> ports_;
      IdxT                   branch_id_{INVALID_INDEX<IdxT>};
    };

  } // namespace PhasorDynamics
} // namespace GridKit
