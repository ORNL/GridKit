/**
 * @file TransformerThreeWinding.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of an EMT three-winding transformer.
 */

#pragma once

#include <array>
#include <cstddef>
#include <string_view>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/TransformerThreeWinding/TransformerThreeWindingData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `TransformerThreeWinding`, each of size 3 and in per unit
    enum class TransformerThreeWindingInternalVariables : size_t
    {
      i13,  ///< \f$\mathbf{i}_{13}\f$ Leakage current from node 1 to node 3
      i23,  ///< \f$\mathbf{i}_{23}\f$ Leakage current from node 2 to node 3
      psi1, ///< \f$\boldsymbol{\psi}_1\f$ Node 1 magnetizing flux linkage
      psi2, ///< \f$\boldsymbol{\psi}_2\f$ Node 2 magnetizing flux linkage
      psi3, ///< \f$\boldsymbol{\psi}_3\f$ Node 3 magnetizing flux linkage
    };

    /// Outputs of a `TransformerThreeWinding`, each of size 3 and in per unit
    enum class TransformerThreeWindingOutputs : size_t
    {
      e1,  ///< \f$\mathbf{e}_1\f$ Node 1 voltage
      e2,  ///< \f$\mathbf{e}_2\f$ Node 2 voltage
      e3,  ///< \f$\mathbf{e}_3\f$ Node 3 voltage
      im1, ///< \f$\mathbf{i}_\mathrm{m1}\f$ Node 1 magnetizing current
      im2, ///< \f$\mathbf{i}_\mathrm{m2}\f$ Node 2 magnetizing current
      im3, ///< \f$\mathbf{i}_\mathrm{m3}\f$ Node 3 magnetizing current
      iw1, ///< \f$\mathbf{i}_\mathrm{w1}\f$ Winding 1 current
      iw2, ///< \f$\mathbf{i}_\mathrm{w2}\f$ Winding 2 current
      iw3, ///< \f$\mathbf{i}_\mathrm{w3}\f$ Winding 3 current
    };

    /// Three-phase bank of three-winding transformers with coupled leakage and a saturable core
    template <typename scalar_type, typename index_type>
    class TransformerThreeWinding final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;

      static constexpr size_t PHASES{3};

    public:
      using ModelDataT         = TransformerThreeWindingData<RealT, IdxT>;
      using InternalVariablesT = TransformerThreeWindingInternalVariables;
      using TerminalsT         = TransformerThreeWindingTerminals;
      using InputsT            = Utilities::EmptyEnum;
      using OutputsT           = TransformerThreeWindingOutputs;

      static constexpr std::string_view CLASS_NAME = "TransformerThreeWinding";

      static constexpr std::array<IdxT, Utilities::enum_size<TransformerThreeWindingInternalVariables>()> NUM_INTERNAL_VARIABLES{PHASES, PHASES, PHASES, PHASES, PHASES};
      static constexpr std::array<IdxT, Utilities::enum_size<TransformerThreeWindingTerminals>()>         NUM_TERMINAL_VARIABLES{PHASES, PHASES, PHASES};
      static constexpr std::array<IdxT, Utilities::enum_size<InputsT>()>                                  NUM_INPUT_VARIABLES{};
      static constexpr std::array<IdxT, Utilities::enum_size<TransformerThreeWindingOutputs>()>           NUM_OUTPUT_VARIABLES{PHASES, PHASES, PHASES, PHASES, PHASES, PHASES, PHASES, PHASES, PHASES};

      explicit TransformerThreeWinding(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        tagDifferentiable() override;

      int evaluateResidual() override;
      int evaluateOutputs() override;
      int evaluateJacobian() override;

      template <typename T>
      [[gnu::always_inline]] void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

      template <typename T>
      [[gnu::always_inline]] void externalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

      template <typename T>
      [[gnu::always_inline]] void output(OutputsT output, const T* y, const T* yp, const T* w, const T* wp, T* f) const;

    private:
      using MapT = std::array<RealT, PHASES * PHASES>;

      // Local indices of the internal variables, in TransformerThreeWindingInternalVariables order
      enum : size_t
      {
        I13  = 0,
        I23  = PHASES,
        PSI1 = 2 * PHASES,
        PSI2 = 3 * PHASES,
        PSI3 = 4 * PHASES,
      };

      // Offsets of the terminals, in TransformerThreeWindingTerminals order
      enum : size_t
      {
        BUS1 = 0,
        BUS2 = PHASES,
        BUS3 = 2 * PHASES,
      };

      /// One winding, with its connection map in row-major order
      struct Winding
      {
        MapT  P{};
        RealT Vpk{ONE<RealT>};
        RealT Ipk{ZERO<RealT>};
        RealT turns{ONE<RealT>};

        /// Per-unit voltage of phase `n` of the winding from the voltages `bus` of its terminal
        template <typename T>
        [[gnu::always_inline]] T voltage(const T* bus, size_t n) const
        {
          return (P[n] * bus[0] + P[PHASES + n] * bus[1] + P[2 * PHASES + n] * bus[2]) / Vpk;
        }

        /// Write the current injections of the terminal to `f` from the winding currents `iw`
        template <typename T>
        [[gnu::always_inline]] void inject(const T* iw, T* f) const
        {
#pragma clang loop unroll(full)
          for (size_t m = 0; m < PHASES; ++m)
          {
            f[m] = -Ipk * (P[PHASES * m] * iw[0] + P[PHASES * m + 1] * iw[1] + P[PHASES * m + 2] * iw[2]);
          }
        }
      };

      /// Node voltage, magnetizing current and winding current of one phase of a winding
      template <typename T>
      struct WindingValues
      {
        T e;
        T im;
        T iw;
      };

      template <typename T>
      [[gnu::always_inline]] T magnetizing(T psi) const;

      template <typename T>
      [[gnu::always_inline]] WindingValues<T> windingValues(const Winding& winding, T u, T psi, T series) const;

      static bool  connectionMap(const Matrix<RealT>& matrix, MapT& map);
      static void  shiftPhase(MapT& map, RealT angle);
      static RealT windingVoltage(const MapT& map, RealT V);

      IdxT          N_{PHASES};
      RealT         S_{ZERO<RealT>};
      RealT         V1_{ZERO<RealT>};
      RealT         V2_{ZERO<RealT>};
      RealT         V3_{ZERO<RealT>};
      RealT         f_b_{ZERO<RealT>};
      Matrix<RealT> P1_;
      Matrix<RealT> P2_;
      Matrix<RealT> P3_;
      RealT         tap_{ONE<RealT>};
      RealT         phase_{ZERO<RealT>};
      RealT         R12_{ZERO<RealT>};
      RealT         R23_{ZERO<RealT>};
      RealT         R31_{ZERO<RealT>};
      RealT         X12_{ZERO<RealT>};
      RealT         X23_{ZERO<RealT>};
      RealT         X31_{ZERO<RealT>};
      RealT         I0_{ZERO<RealT>};
      RealT         P0_{ZERO<RealT>};
      RealT         knee_{ZERO<RealT>};
      RealT         Lsat_{ZERO<RealT>};

      Winding w1_;
      Winding w2_;
      Winding w3_;
      RealT   omega_b_{ZERO<RealT>};
      RealT   R3_{ZERO<RealT>};   ///< Resistance that the two leakage loops share
      RealT   X3_{ZERO<RealT>};   ///< Reactance that the two leakage loops share
      RealT   beta_{ZERO<RealT>}; ///< Share of the magnetizing branch at each node
      RealT   Gc_{ZERO<RealT>};
      RealT   G_{ZERO<RealT>}; ///< Share of the core-loss conductance at each node
      RealT   inv_Lm_{ZERO<RealT>};
      RealT   slope_{ZERO<RealT>}; ///< 1 / Lsat - 1 / Lm
    };
  } // namespace EMT
} // namespace GridKit
