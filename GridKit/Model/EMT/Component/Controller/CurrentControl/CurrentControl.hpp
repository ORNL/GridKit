/**
 * @file CurrentControl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of an EMT current control.
 */

#pragma once

#include <array>
#include <cstddef>
#include <string_view>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/CurrentControl/CurrentControlData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `CurrentControl`, each of size 2
    enum class CurrentControlInternalVariables : size_t
    {
      xdq, ///< \f$\mathbf{x}_\mathrm{dq}\f$ Integrator states [V]
      mdq, ///< \f$\mathbf{m}_\mathrm{dq}\f$ Modulation signals
    };

    /// Current control in dq coordinates with a limited modulation magnitude
    template <typename scalar_type, typename index_type>
    class CurrentControl final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;

    public:
      using ModelDataT         = CurrentControlData<RealT, IdxT>;
      using InternalVariablesT = CurrentControlInternalVariables;
      using TerminalsT         = Utilities::EmptyEnum;
      using InputsT            = CurrentControlInputs;
      using OutputsT           = Utilities::EmptyEnum;

      static constexpr std::string_view CLASS_NAME = "CurrentControl";

      static constexpr std::array<IdxT, Utilities::enum_size<CurrentControlInternalVariables>()> NUM_INTERNAL_VARIABLES{2, 2};
      static constexpr std::array<IdxT, Utilities::enum_size<TerminalsT>()>                      NUM_TERMINAL_VARIABLES{};
      static constexpr std::array<IdxT, Utilities::enum_size<CurrentControlInputs>()>            NUM_INPUT_VARIABLES{1, 1, 2, 2, 1, 2, 1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<OutputsT>()>                        NUM_OUTPUT_VARIABLES{};

      explicit CurrentControl(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        setAbsoluteTolerance(RealT rel_tol) override;
      int                        tagDifferentiable() override;

      int evaluateResidual() override;
      int evaluateJacobian() override;

      template <typename T>
      [[gnu::always_inline]] void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

    private:
      // Local indices of the internal variables, in CurrentControlInternalVariables order
      enum : size_t
      {
        XD,
        XQ,
        MD,
        MQ,
      };

      // Offsets of the external values, in CurrentControlInputs order
      enum : size_t
      {
        VPOS,
        VNEG,
        VD,
        VQ,
        ID,
        IQ,
        OMEGA,
        AXIS_D,
        AXIS_Q,
        IP_CMD,
        IQ_CMD,
      };

      RealT V_{ZERO<RealT>};
      RealT mva_{ZERO<RealT>};
      RealT Kp_{ZERO<RealT>};
      RealT Ki_{ZERO<RealT>};
      RealT L_{ZERO<RealT>};
      RealT C_{ZERO<RealT>};
      RealT Mmax_{ONE<RealT>};

      RealT current_peak_{ZERO<RealT>}; ///< Peak phase current at the rating [A]
    };
  } // namespace EMT
} // namespace GridKit
