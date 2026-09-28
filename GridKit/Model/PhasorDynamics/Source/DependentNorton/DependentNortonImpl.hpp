/**
 * @file DependentNortonImpl.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Definition of the Norton source model.
 */

#pragma once

#include <cmath>
#include <variant>

#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
#include <GridKit/Model/PhasorDynamics/Source/DependentNorton/DependentNorton.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Source
    {
      /**
       * @brief Construct an unconfigured Norton source.
       */
      template <typename scalar_type, typename index_type>
      DependentNorton<scalar_type, index_type>::DependentNorton(BusT* bus)
        : bus_(bus)
      {
        size_ = 0;
      }

      /**
       * @brief Construct a Norton source from model data.
       *
       * @param[in] bus Terminal bus the source injects into.
       * @param[in] data Parameters and monitored-variable selections.
       */
      template <typename scalar_type, typename index_type>
      DependentNorton<scalar_type, index_type>::DependentNorton(BusT* bus, const ModelDataT& data)
        : bus_(bus),
          monitor_(std::make_unique<MonitorT>(data))
      {
        size_ = 0;
        initializeParameters(data);
        initializeMonitor();
      }

      template <typename scalar_type, typename index_type>
      DependentNorton<scalar_type, index_type>::~DependentNorton()
      {
      }

      /**
       * @brief Set the component identifier assigned by the system model.
       */
      template <typename scalar_type, typename index_type>
      int DependentNorton<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
      {
        gridkit_component_id_ = component_id;
        return 0;
      }

      /**
       * @brief Allocate the voltage, source-current, and current-contribution buffers.
       */
      template <typename scalar_type, typename index_type>
      int DependentNorton<scalar_type, index_type>::allocate()
      {
        if (!allocated_)
        {
          this->allocateVectors(size_);
        }

        wb_.resize(2);
        wb_.setToZero();
        h_.resize(2);
        h_.setToZero();

        auto signal_size = Utilities::enum_size<DependentNortonExternalVariables>();
        ws_.resize(static_cast<IdxT>(signal_size));
        ws_.setToZero();
        ws_indices_.assign(signal_size, INVALID_INDEX<IdxT>);

        allocated_ = true;
        return 0;
      }

      /**
       * @brief Validate the admittance, terminal connection, and required inputs.
       */
      template <typename scalar_type, typename index_type>
      int DependentNorton<scalar_type, index_type>::verify() const
      {
        int ret = static_cast<int>(parameter_error_count_);

        auto check = [&](bool condition, const char* message)
        {
          if (!condition)
          {
            Log::error() << "DependentNorton: " << message << '\n';
            ret += 1;
          }
        };

        check(bus_ != nullptr, "bus pointer is null");
        check(std::isfinite(G_), "G must be finite");
        check(std::isfinite(B_), "B must be finite");

        const auto inr_port = ports_.in.template port<DependentNortonSignalInputs::inr>();
        const auto ini_port = ports_.in.template port<DependentNortonSignalInputs::ini>();
        check(inr_port && inr_port.linked(), "inr requires a linked source");
        check(ini_port && ini_port.linked(), "ini requires a linked source");

        return ret;
      }

      /**
       * @brief Verify the prescribed inputs; there are no internal variables.
       */
      template <typename scalar_type, typename index_type>
      int DependentNorton<scalar_type, index_type>::initialize()
      {
        return verify();
      }

      template <typename scalar_type, typename index_type>
      int DependentNorton<scalar_type, index_type>::tagDifferentiable()
      {
        return 0;
      }

      template <typename scalar_type, typename index_type>
      int DependentNorton<scalar_type, index_type>::setAbsoluteTolerance([[maybe_unused]] RealT rel_tol)
      {
        return 0;
      }

      /**
       * @brief Read the voltage and current inputs and accumulate terminal current.
       *
       * @pre The terminal-bus residual has been zeroed for this evaluation.
       */
      template <typename scalar_type, typename index_type>
      int DependentNorton<scalar_type, index_type>::evaluateResidual()
      {
        const auto INR = static_cast<size_t>(DependentNortonExternalVariables::INR);
        const auto INI = static_cast<size_t>(DependentNortonExternalVariables::INI);

        const auto inr_port = ports_.in.template port<DependentNortonSignalInputs::inr>();
        const auto ini_port = ports_.in.template port<DependentNortonSignalInputs::ini>();

        auto* ws         = ws_.getData();
        ws[INR]          = inr_port.readSignal();
        ws[INI]          = ini_port.readSignal();
        ws_indices_[INR] = inr_port.signalVariableIndex();
        ws_indices_[INI] = ini_port.signalVariableIndex();

        auto* wb = wb_.getData();
        wb[0]    = Vr();
        wb[1]    = Vi();

        auto* h = h_.getData();
        evaluateBusResidual(y_.getData(), yp_.getData(), wb, ws, h);
        bus_->Ir() += h[0];
        bus_->Ii() += h[1];
        if (bus_->size() > 0)
        {
          bus_->getResidual().setDataUpdated();
        }

        return 0;
      }

      template <typename scalar_type, typename index_type>
      const Model::VariableMonitorBase* DependentNorton<scalar_type, index_type>::getMonitor() const
      {
        return monitor_.get();
      }

      /**
       * @brief Evaluate the terminal-current contribution on system base.
       *
       * @param[in] y Internal variables (unused).
       * @param[in] yp Internal variable derivatives (unused).
       * @param[in] wb Terminal voltage components.
       * @param[in] ws Norton source-current components.
       * @param[out] h Current injected into the terminal bus.
       */
      template <typename scalar_type, typename index_type>
      __attribute__((always_inline)) inline int DependentNorton<scalar_type, index_type>::evaluateBusResidual(
          [[maybe_unused]] const ScalarT* y,
          [[maybe_unused]] const ScalarT* yp,
          const ScalarT*                  wb,
          const ScalarT*                  ws,
          ScalarT*                        h)
      {
        const auto INR = static_cast<size_t>(DependentNortonExternalVariables::INR);
        const auto INI = static_cast<size_t>(DependentNortonExternalVariables::INI);

        const ScalarT Vr  = wb[0];
        const ScalarT Vi  = wb[1];
        const ScalarT INr = ws[INR];
        const ScalarT INi = ws[INI];

        h[0] = INr - G_ * Vr + B_ * Vi;
        h[1] = INi - B_ * Vr - G_ * Vi;

        return 0;
      }

      //
      //  Private methods
      //

      /**
       * @brief Read required parameters, recording errors for verify().
       */
      template <typename scalar_type, typename index_type>
      void DependentNorton<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
      {
        using Params = typename ModelDataT::Parameters;

        parameter_error_count_ = 0;

        auto load_required_real = [&](auto key, RealT& target, const char* name)
        {
          if (!data.parameters.contains(key))
          {
            Log::error() << "DependentNorton: missing required parameter '" << name << "'\n";
            ++parameter_error_count_;
            return;
          }

          const auto& value = data.parameters.at(key);
          if (const auto* real_value = std::get_if<RealT>(&value))
          {
            target = *real_value;
          }
          else if (const auto* index_value = std::get_if<IdxT>(&value))
          {
            target = static_cast<RealT>(*index_value);
          }
          else
          {
            Log::error() << "DependentNorton: parameter '" << name << "' must be numeric\n";
            ++parameter_error_count_;
          }
        };

        load_required_real(Params::G, G_, "G");
        load_required_real(Params::B, B_, "B");
      }

      /**
       * @brief Monitor terminal current and power from the current voltage and inputs.
       */
      template <typename scalar_type, typename index_type>
      void DependentNorton<scalar_type, index_type>::initializeMonitor()
      {
        using Variable = typename ModelDataT::MonitorableVariables;

        monitor_->set(Variable::ir, [this]
                      { return INr() - G_ * Vr() + B_ * Vi(); });
        monitor_->set(Variable::ii, [this]
                      { return INi() - B_ * Vr() - G_ * Vi(); });
        monitor_->set(Variable::p, [this]
                      { return Vr() * INr() + Vi() * INi() - G_ * (Vr() * Vr() + Vi() * Vi()); });
        monitor_->set(Variable::q, [this]
                      { return Vi() * INr() - Vr() * INi() + B_ * (Vr() * Vr() + Vi() * Vi()); });
      }

      /**
       * @brief Read the real terminal voltage owned by the connected bus.
       */
      template <typename scalar_type, typename index_type>
      typename DependentNorton<scalar_type, index_type>::ScalarT&
      DependentNorton<scalar_type, index_type>::Vr()
      {
        return bus_->Vr();
      }

      /**
       * @brief Read the imaginary terminal voltage owned by the connected bus.
       */
      template <typename scalar_type, typename index_type>
      typename DependentNorton<scalar_type, index_type>::ScalarT&
      DependentNorton<scalar_type, index_type>::Vi()
      {
        return bus_->Vi();
      }

      /**
       * @brief Read the real Norton source-current input.
       */
      template <typename scalar_type, typename index_type>
      typename DependentNorton<scalar_type, index_type>::ScalarT
      DependentNorton<scalar_type, index_type>::INr() const
      {
        return ports_.in.template port<DependentNortonSignalInputs::inr>().readSignal();
      }

      /**
       * @brief Read the imaginary Norton source-current input.
       */
      template <typename scalar_type, typename index_type>
      typename DependentNorton<scalar_type, index_type>::ScalarT
      DependentNorton<scalar_type, index_type>::INi() const
      {
        return ports_.in.template port<DependentNortonSignalInputs::ini>().readSignal();
      }
    } // namespace Source
  } // namespace PhasorDynamics
} // namespace GridKit
