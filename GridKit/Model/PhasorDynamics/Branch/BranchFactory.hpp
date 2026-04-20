/**
 * @file BranchFactory.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Factory that constructs concrete BranchBase subclasses from BranchData.
 */

#pragma once

#include <stdexcept>
#include <string>

#include <GridKit/Model/PhasorDynamics/Branch/Branch.hpp>
#include <GridKit/Model/PhasorDynamics/Branch/BranchData.hpp>
#include <GridKit/Model/PhasorDynamics/BranchBase.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Constructs a concrete `BranchBase` from `BranchData`, dispatching
     * on `BranchData::BranchType`.
     *
     * Structural bus/branch compatibility is not checked here; the concrete
     * branch's `verifyPorts()` (run from `allocate()`) is responsible.
     *
     * Throws `std::runtime_error` if `branch_type` is unrecognized so
     * misconfiguration fails at construction rather than as a `nullptr`
     * dereference later.
     */
    template <typename ScalarT = double, typename IdxT = int>
    class BranchFactory
    {
    public:
      using RealT       = typename Model::Evaluator<ScalarT, IdxT>::RealT;
      using BranchData  = GridKit::PhasorDynamics::BranchData<RealT, IdxT>;
      using BranchTypeT = typename GridKit::PhasorDynamics::BranchData<RealT, IdxT>::BranchType;
      using bus_type    = BusBase<ScalarT, IdxT>;
      using Log         = ::GridKit::Utilities::Logger;

      BranchFactory() = delete;

      /**
       * @brief Construct a concrete branch selected by `data.branch_type`.
       *
       * @param data branch data; `branch_type` selects the concrete type.
       * @param bus1 bus attached to port 0.
       * @param bus2 bus attached to port 1.
       * @return owning pointer to a newly-allocated concrete branch.
       * @throws std::runtime_error if `branch_type` is unrecognized.
       */
      static BranchBase<ScalarT, IdxT>* create(const BranchData& data,
                                               bus_type*         bus1,
                                               bus_type*         bus2)
      {
        switch (data.branch_type)
        {
        case BranchTypeT::LINE:
          return new Branch<ScalarT, IdxT>(bus1, bus2, data);
        default:
          Log::error() << "BranchFactory: unrecognized branch_type "
                       << static_cast<int>(data.branch_type) << std::endl;
          throw std::runtime_error(
              "BranchFactory::create: unrecognized branch_type "
              + std::to_string(static_cast<int>(data.branch_type)));
        }
      }
    };
  } // namespace PhasorDynamics
} // namespace GridKit
