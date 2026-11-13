/**
 * @file Bus.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of an EMT bus.
 */

#pragma once

#include <cstddef>
#include <string_view>

#include <GridKit/Model/EMT/Bus/BusData.hpp>
#include <GridKit/Model/EMT/Component.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `Bus`, each of size N
    enum class BusInternalVariables : size_t
    {
      v, ///< \f$\mathbf{v}\f$ Bus voltage [V]
    };

    /// N-phase bus, whose residuals collect the injected currents
    template <typename scalar_type, typename index_type>
    class Bus final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;

    public:
      using ModelDataT         = BusData<RealT, IdxT>;
      using InternalVariablesT = BusInternalVariables;
      using PortsT             = Utilities::EmptyEnum;

      static constexpr std::string_view CLASS_NAME = "Bus";

      explicit Bus(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        setAbsoluteTolerance(RealT rel_tol) override;

    private:
      IdxT  N_{0};
      RealT kv_{ZERO<RealT>};
    };
  } // namespace EMT
} // namespace GridKit
