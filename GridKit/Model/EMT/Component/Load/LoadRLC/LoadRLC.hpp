/**
 * @file LoadRLC.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of an EMT parallel RLC load.
 */

#pragma once

#include <cstddef>
#include <optional>
#include <string_view>
#include <vector>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Load/LoadRLC/LoadRLCData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `LoadRLC`
    enum class LoadRLCInternalVariables : size_t
    {
      i_l, ///< \f$i_{l,n}\f$ Inductor current injection of each phase with an inductor, in phase order [A]
    };

    /// Outputs of a `LoadRLC`, each of size N
    enum class LoadRLCOutputs : size_t
    {
      i, ///< \f$\mathbf{i}\f$ Current injection into the bus [A]
    };

    /// N-phase load with a resistor, an inductor and a capacitor in parallel on each phase
    template <typename scalar_type, typename index_type>
    class LoadRLC final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;

    public:
      using ModelDataT         = LoadRLCData<RealT, IdxT>;
      using InternalVariablesT = LoadRLCInternalVariables;
      using TerminalsT         = LoadRLCTerminals;
      using InputsT            = Utilities::EmptyEnum;
      using OutputsT           = LoadRLCOutputs;

      static constexpr std::string_view CLASS_NAME = "LoadRLC";

      explicit LoadRLC(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        setAbsoluteTolerance(RealT rel_tol) override;
      int                        tagDifferentiable() override;

      int evaluateOutputs() override;

      template <typename T>
      void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

      template <typename T>
      void externalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

      template <typename T>
      void output(OutputsT output, const T* y, const T* yp, const T* w, const T* wp, T* f) const;

    private:
      using PhaseValues = typename ModelDataT::PhaseValues;

      /// Branches of one phase, with zero for an absent branch
      struct Phase
      {
        RealT G{ZERO<RealT>};
        RealT L{ZERO<RealT>};
        RealT C{ZERO<RealT>};
        IdxT  i_l{0}; ///< Local index of the inductor current, for a phase with an inductor

        template <typename T>
        T current(const T* y, const T& v, const T& vp) const;
      };

      IdxT                       N_{0};
      std::optional<PhaseValues> R_;
      std::optional<PhaseValues> L_;
      std::optional<PhaseValues> C_;
      std::vector<Phase>         phases_;
    };
  } // namespace EMT
} // namespace GridKit
