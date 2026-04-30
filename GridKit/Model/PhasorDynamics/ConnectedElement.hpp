#pragma once

#include <memory>

#include <GridKit/Model/PhasorDynamics/ComponentSignals.hpp>
#include <GridKit/Model/PhasorDynamics/GridElement.hpp>
#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename ElementP>
    struct ConnectedElementTraits
    {
    };

    /**
     * @brief Model base class for all system constituents
     */
    template <typename ElementP>
    class ConnectedElement : public ConnectedElementTraits<ElementP>::InterfaceT
    {
    public:
      using Traits     = ConnectedElementTraits<ElementP>;
      using Superclass = typename Traits::InterfaceT;
      using ScalarT    = typename Traits::ScalarT;
      using IdxT       = typename Traits::IdxT;
      using RealT      = typename Model::Evaluator<ScalarT, IdxT>::RealT;
      using MatrixT    = typename Model::Evaluator<ScalarT, IdxT>::MatrixT;
      using ModelDataT = typename Traits::ModelDataT;
      using MonitorT   = Model::VariableMonitor<ElementP>;

      using InternalVariablesT = typename Traits::InternalVariablesT;
      using ExternalVariablesT = typename Traits::ExternalVariablesT;

      ConnectedElement() = default;

      ConnectedElement(const ModelDataT& data);

      virtual ~ConnectedElement();

      auto getSignals()
          -> ComponentSignals<ScalarT,
                              IdxT,
                              InternalVariablesT,
                              ExternalVariablesT>&
      {
        return signals_;
      }

      const Model::VariableMonitorBase* getMonitor() const override;

    protected:
      /// Component signal extension
      ComponentSignals<ScalarT, IdxT, InternalVariablesT, ExternalVariablesT> signals_;

      /// Variable monitor
      std::unique_ptr<MonitorT> monitor_;

      using NotImplementedError = GridKit::Utilities::NotImplementedError;

    public:
      // TODO: evaluate how this complies with xSDK guidelines

      [[noreturn]] IdxT sizeQuadrature() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] IdxT sizeParams() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] std::vector<ScalarT>& yB() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const std::vector<ScalarT>& yB() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] std::vector<ScalarT>& ypB() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const std::vector<ScalarT>& ypB() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] std::vector<ScalarT>& param() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const std::vector<ScalarT>& param() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] std::vector<ScalarT>& param_up() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const std::vector<ScalarT>& param_up() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] std::vector<ScalarT>& param_lo() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const std::vector<ScalarT>& param_lo() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] int evaluateIntegrand() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] int initializeAdjoint() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] int evaluateAdjointResidual() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] int evaluateAdjointIntegrand() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] std::vector<ScalarT>& getIntegrand() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const std::vector<ScalarT>& getIntegrand() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] std::vector<ScalarT>& getAdjointResidual() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const std::vector<ScalarT>& getAdjointResidual() const override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] std::vector<ScalarT>& getAdjointIntegrand() override
      {
        throw NotImplementedError(__func__);
      }

      [[noreturn]] const std::vector<ScalarT>& getAdjointIntegrand() const override
      {
        throw NotImplementedError(__func__);
      }
    };

  } // namespace PhasorDynamics
} // namespace GridKit
