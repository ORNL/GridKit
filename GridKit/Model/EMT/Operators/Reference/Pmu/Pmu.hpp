/**
 * @file Pmu.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of an EMT phasor measurement unit.
 */

#pragma once

#include <array>
#include <cstddef>
#include <string_view>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Operators/Reference/Pll/Pll.hpp>
#include <GridKit/Model/EMT/Operators/Reference/Pmu/PmuData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `Pmu`, each of size 1
    enum class PmuInternalVariables : size_t
    {
      vms, ///< \f$v_\mathrm{ms}\f$ Mean-square voltage [V^2]
      p,   ///< \f$p\f$ Active power [W]
      q,   ///< \f$q\f$ Reactive power [var]
    };

    /// Outputs of a `Pmu`, each of size 1
    enum class PmuOutputs : size_t
    {
      vrms,  ///< \f$v_\mathrm{rms}\f$ Line-to-line RMS voltage [V]
      theta, ///< \f$\theta\f$ Angle [rad]
      omega, ///< \f$\omega\f$ Frequency [rad/s]
    };

    /// Submodels of a `Pmu`
    enum class PmuSubmodels : size_t
    {
      pll, ///< Phase-locked loop on the measured voltage
    };

    /// Phasor measurement unit on a three-phase bus
    template <typename scalar_type, typename index_type>
    class Pmu final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;
      using PllT       = Pll<scalar_type, index_type>;

    public:
      using ModelDataT         = PmuData<RealT, IdxT>;
      using InternalVariablesT = PmuInternalVariables;
      using TerminalsT         = Utilities::EmptyEnum;
      using InputsT            = PmuInputs;
      using OutputsT           = PmuOutputs;
      using SubmodelsT         = PmuSubmodels;

      static constexpr std::string_view CLASS_NAME = "Pmu";

      static constexpr std::array<IdxT, Utilities::enum_size<PmuInternalVariables>()> NUM_INTERNAL_VARIABLES{1, 1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<TerminalsT>()>           NUM_TERMINAL_VARIABLES{};
      static constexpr std::array<IdxT, Utilities::enum_size<PmuInputs>()>            NUM_INPUT_VARIABLES{3, 3, 2};
      static constexpr std::array<IdxT, Utilities::enum_size<PmuOutputs>()>           NUM_OUTPUT_VARIABLES{1, 1, 1};

      explicit Pmu(const ModelDataT& data);

      PllT& submodel(SubmodelsT symbol);

      Model::ConfigurationChecks verify() const override;
      int                        setAbsoluteTolerance(RealT rel_tol) override;
      int                        tagDifferentiable() override;

      int evaluateResidual() override;
      int evaluateOutputs() override;
      int evaluateJacobian() override;

      template <typename T>
      [[gnu::always_inline]] void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

      template <typename T>
      [[gnu::always_inline]] void output(OutputsT output, const T* y, const T* yp, const T* w, const T* wp, T* f) const;

    private:
      // Offsets of the external values, in PmuInputs order
      enum : size_t
      {
        V   = 0,
        I   = 3,
        PLL = 6,
      };

      const PllT* pll_{nullptr};
    };
  } // namespace EMT
} // namespace GridKit
