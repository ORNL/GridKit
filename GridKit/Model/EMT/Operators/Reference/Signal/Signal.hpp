/**
 * @file Signal.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of an EMT signal.
 */

#pragma once

#include <cstddef>
#include <string_view>
#include <vector>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Operators/Reference/Signal/SignalData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `Signal`, each of size N
    enum class SignalInternalVariables : size_t
    {
      s, ///< \f$\mathbf{s}\f$ Signal
    };

    /// N signals, which are held or follow the output added to their residuals
    template <typename scalar_type, typename index_type>
    class Signal final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;

    public:
      using ModelDataT         = SignalData<RealT, IdxT>;
      using InternalVariablesT = SignalInternalVariables;
      using TerminalsT         = Utilities::EmptyEnum;
      using InputsT            = Utilities::EmptyEnum;
      using OutputsT           = Utilities::EmptyEnum;

      static constexpr std::string_view CLASS_NAME = "Signal";

      explicit Signal(const ModelDataT& data);

      Model::ConfigurationChecks verify() const override;
      int                        tagDifferentiable() override;

      template <typename T>
      void internalResidual(const T* y, const T* yp, const T* w, const T* wp, T* f) const;

    private:
      IdxT N_{0};
      bool held_{false};
    };
  } // namespace EMT
} // namespace GridKit
