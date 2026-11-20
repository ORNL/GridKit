/**
 * @file Ieeest.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of an EMT IEEE power system stabilizer.
 */

#pragma once

#include <array>
#include <cstddef>
#include <string_view>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/Ieeest/IeeestData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of an `Ieeest`, each of size 1
    enum class IeeestInternalVariables : size_t
    {
      X1,  ///< \f$x_1\f$ Notch-filter signal state
      X2,  ///< \f$x_2\f$ First derivative of the filtered signal
      X3,  ///< \f$x_3\f$ Second derivative of the filtered signal
      X4,  ///< \f$x_4\f$ Third derivative of the filtered signal
      X5,  ///< \f$x_5\f$ Lead-lag 1 state
      X6,  ///< \f$x_6\f$ Lead-lag 2 state
      X7,  ///< \f$x_7\f$ Washout state
      V4,  ///< \f$v_4\f$ Notch-filter output
      V5,  ///< \f$v_5\f$ Lead-lag 1 output
      V6,  ///< \f$v_6\f$ Lead-lag 2 output
      V7,  ///< \f$v_7\f$ Stabilizer signal before the limits
      VSS, ///< \f$V_{\mathrm{ss}}\f$ Limited stabilizer signal [p.u.]
    };

    /// IEEE power system stabilizer
    template <typename scalar_type, typename index_type>
    class Ieeest final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;

    public:
      using ModelDataT         = IeeestData<RealT, IdxT>;
      using InternalVariablesT = IeeestInternalVariables;
      using TerminalsT         = Utilities::EmptyEnum;
      using InputsT            = IeeestInputs;
      using OutputsT           = Utilities::EmptyEnum;

      static constexpr std::string_view CLASS_NAME = "Ieeest";

      static constexpr std::array<IdxT, Utilities::enum_size<IeeestInternalVariables>()> NUM_INTERNAL_VARIABLES{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<TerminalsT>()>              NUM_TERMINAL_VARIABLES{};
      static constexpr std::array<IdxT, Utilities::enum_size<IeeestInputs>()>            NUM_INPUT_VARIABLES{1};
      static constexpr std::array<IdxT, Utilities::enum_size<OutputsT>()>                NUM_OUTPUT_VARIABLES{};

      explicit Ieeest(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        tagDifferentiable() override;

      int evaluateResidual() override;
      int evaluateJacobian() override;

      template <typename T>
      [[gnu::always_inline]] void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

    private:
      // Local indices of the internal variables, in IeeestInternalVariables order
      enum : size_t
      {
        X1,
        X2,
        X3,
        X4,
        X5,
        X6,
        X7,
        V4,
        V5,
        V6,
        V7,
        VSS,
      };

      // Offsets of the external values, in IeeestInputs order
      enum : size_t
      {
        U,
      };

      RealT A5_{ZERO<RealT>};
      RealT A6_{ZERO<RealT>};
      RealT T1_{ZERO<RealT>};
      RealT T2_{ZERO<RealT>};
      RealT T3_{ZERO<RealT>};
      RealT T4_{ZERO<RealT>};
      RealT T5_{ZERO<RealT>};
      RealT T6_{ZERO<RealT>};
      RealT Ks_{ZERO<RealT>};
      RealT Lsmin_{ZERO<RealT>};
      RealT Lsmax_{ZERO<RealT>};
      RealT Vcl_{ZERO<RealT>};
      RealT Vcu_{ZERO<RealT>};
      RealT Tdelay_{ZERO<RealT>};

      RealT a1_{ZERO<RealT>}; ///< Coefficients of the notch-filter denominator, whose constant term is one
      RealT a2_{ZERO<RealT>};
      RealT a3_{ZERO<RealT>};
      RealT a4_{ZERO<RealT>};

      RealT use_notch_{ZERO<RealT>};     ///< One with a notch filter
      RealT use_4th_order_{ZERO<RealT>}; ///< One for each order of the notch filter
      RealT use_3rd_order_{ZERO<RealT>};
      RealT use_2nd_order_{ZERO<RealT>};
      RealT inv_a4_{ZERO<RealT>}; ///< Reciprocal of the leading coefficient of each order, zero for the others
      RealT inv_a3_{ZERO<RealT>};
      RealT inv_a2_{ZERO<RealT>};
      RealT use_T2_{ZERO<RealT>}; ///< One for each block with a denominator time constant
      RealT use_T4_{ZERO<RealT>};
      RealT use_T6_{ZERO<RealT>};
    };
  } // namespace EMT
} // namespace GridKit
