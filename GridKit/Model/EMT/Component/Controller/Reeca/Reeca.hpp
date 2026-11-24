/**
 * @file Reeca.hpp
 * @brief WECC type A renewable electrical controller.
 */

#pragma once

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Component/Controller/Reeca/ReecaData.hpp>

namespace GridKit
{
  namespace EMT
  {
    enum class ReecaInternalVariables : size_t
    {
      VTERM,
      VMEAS,
      PMEAS,
      XPIQ,
      XPIV,
      QV,
      PREF_RATE,
      PORD,
      IQV,
      IPMAX,
      IQCMD,
      IPCMD,
      IQMAX,
      ILCAP,
      VSAFE,
      EQ,
      EPIV,
    };

    /// WECC type A renewable electrical controller
    template <typename ScalarT, typename IdxT>
    class Reeca final : public Component<ScalarT, IdxT>
    {
      using ComponentT = Component<ScalarT, IdxT>;
      using RealT      = typename ComponentT::RealT;

    public:
      using ModelDataT         = ReecaData<RealT, IdxT>;
      using InternalVariablesT = ReecaInternalVariables;
      using TerminalsT         = Utilities::EmptyEnum;
      using InputsT            = ReecaInputs;
      using OutputsT           = Utilities::EmptyEnum;

      static constexpr std::string_view                                             CLASS_NAME = "Reeca";
      static constexpr std::array<IdxT, Utilities::enum_size<InternalVariablesT>()> NUM_INTERNAL_VARIABLES{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
      static constexpr std::array<IdxT, 0>                                          NUM_TERMINAL_VARIABLES{};
      static constexpr std::array<IdxT, Utilities::enum_size<InputsT>()>            NUM_INPUT_VARIABLES{1, 1, 1, 1, 1, 1, 1};
      static constexpr std::array<IdxT, 0>                                          NUM_OUTPUT_VARIABLES{};

      explicit Reeca(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        tagDifferentiable() override;
      int                        evaluateResidual() override;
      int                        evaluateJacobian() override;

      template <typename T>
      [[gnu::always_inline]] void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;
      int                         initialize() override;

      IdxT sizeRoots() override
      {
        return 4;
      }

      int  evaluateRootsAt(RealT t, const ScalarT* y, RealT* roots) override;
      bool applyRootEvents(RealT t, std::span<const int> crossed) override;

    private:
      enum : size_t
      {
        VTERM,
        VMEAS,
        PMEAS,
        XPIQ,
        XPIV,
        QV,
        PREF_RATE,
        PORD,
        IQV,
        IPMAX,
        IQCMD,
        IPCMD,
        IQMAX,
        ILCAP,
        VSAFE,
        EQ,
        EPIV,
      };

      enum : size_t
      {
        VT_IN     = 0,
        P_IN      = 1,
        Q_IN      = 2,
        QEXT_IN   = 3,
        PFAREF_IN = 4,
        PREF_IN   = 5,
        OMEGA_IN  = 6,
      };

      // Slopes of the three voltage-dependent current limit segments.
      std::array<RealT, 3> slope_q_{}, slope_p_{};
      ModelDataT           data_;
      void                 voltageDip(bool dip, RealT t);
      RealT                reference_voltage_{0};
      RealT                dip_{0};
      RealT                hold_p_{0}, hold_q_{0};
      RealT                release_p_{0}, release_q_{0};
      RealT                held_ipmax_{0};
    };
  } // namespace EMT
} // namespace GridKit
