/**
 * @file Regca.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of an EMT renewable energy generator/converter model.
 */

#pragma once

#include <array>
#include <cstddef>
#include <numbers>
#include <string_view>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Converter/Regca/RegcaData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `Regca`
    enum class RegcaInternalVariables : size_t
    {
      vmeas, ///< Filtered terminal voltage [p.u.]
      ip,    ///< \f$i_p\f$ Active current [p.u.]
      iq,    ///< \f$i_q\f$ Reactive current [p.u.]
      i,     ///< \f$\mathbf{i}\f$ Currents into the grid, size 3 [A]
    };

    /// Converter as a current source that follows limited current orders
    template <typename scalar_type, typename index_type>
    class Regca final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;

    public:
      using ModelDataT         = RegcaData<RealT, IdxT>;
      using InternalVariablesT = RegcaInternalVariables;
      using TerminalsT         = RegcaTerminals;
      using InputsT            = RegcaInputs;
      using OutputsT           = Utilities::EmptyEnum;

      static constexpr std::string_view CLASS_NAME = "Regca";

      static constexpr size_t PHASES = 3;

      static constexpr std::array<IdxT, Utilities::enum_size<RegcaInternalVariables>()> NUM_INTERNAL_VARIABLES{1, 1, 1, PHASES};
      static constexpr std::array<IdxT, Utilities::enum_size<RegcaTerminals>()>         NUM_TERMINAL_VARIABLES{PHASES};
      static constexpr std::array<IdxT, Utilities::enum_size<RegcaInputs>()>            NUM_INPUT_VARIABLES{1, 1, 1, 1};
      static constexpr std::array<IdxT, Utilities::enum_size<OutputsT>()>               NUM_OUTPUT_VARIABLES{};

      explicit Regca(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        initialize() override;
      int                        setAbsoluteTolerance(RealT rel_tol) override;
      int                        tagDifferentiable() override;

      IdxT sizeRoots() override
      {
        return lvpl_ ? 1 : 0;
      }

      int  evaluateRootsAt(RealT, const scalar_type* y, RealT* roots) override;
      bool applyRootEvents(RealT, std::span<const int> crossed) override;

      int evaluateResidual() override;
      int evaluateJacobian() override;

      template <typename T>
      [[gnu::always_inline]] void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

      template <typename T>
      [[gnu::always_inline]] void externalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

    private:
      // Offsets of the internal variables, in RegcaInternalVariables order
      enum : size_t
      {
        VM = 0,
        IP,
        IQ,
        I,
      };

      // Offsets of the external values, the terminal and then RegcaInputs in order
      enum : size_t
      {
        BUS = 0,
        VT  = PHASES,
        THETA,
        IP_CMD,
        IQ_CMD,
      };

      /// Angle of each phase from the angle of the terminal voltage [rad]
      static constexpr std::array<RealT, PHASES> GAMMA{ZERO<RealT>,
                                                       -TWO<RealT>* std::numbers::pi_v<RealT> / THREE<RealT>,
                                                       TWO<RealT>*  std::numbers::pi_v<RealT> / THREE<RealT>};

      RealT V_{ZERO<RealT>};
      RealT mva_{ZERO<RealT>};
      RealT Tg_{ZERO<RealT>};
      RealT TM_{ZERO<RealT>};
      RealT Rqmax_{ZERO<RealT>};
      RealT Rqmin_{ZERO<RealT>};
      RealT Rpmax_{ZERO<RealT>};
      RealT IL1_{ZERO<RealT>};
      RealT VL0_{ZERO<RealT>};
      RealT VL1_{ZERO<RealT>};
      RealT VA0_{ZERO<RealT>};
      RealT VA1_{ZERO<RealT>};
      RealT Vhvmax_{ZERO<RealT>};
      RealT Khv_;
      RealT Qmin_;
      bool  lvpl_;

      RealT rq_max_on_{ZERO<RealT>};    ///< One with an upper rate limit of the reactive current
      RealT rq_min_on_{ZERO<RealT>};    ///< One with a lower rate limit of the reactive current
      RealT current_peak_{ZERO<RealT>}; ///< Peak phase current at the rating [A]
    };
  } // namespace EMT
} // namespace GridKit
