/**
 * @file OverheadImpl.hpp
 *
 * @brief Overhead method definitions: layout, binding, and residual
 * and Jacobian assembly driven by IDA over angular frequency.
 *
 */

#pragma once

#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

#include <GridKit/Model/EMT/SystemAssembly.hpp>
#include <GridKit/Model/VariableMonitorController.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      namespace Detail
      {
        template <typename ScalarT, typename IdxT>
        class OverheadMonitor : public Model::VariableMonitorBase
        {
        public:
          using RealT    = typename ScalarTraits<ScalarT>::RealT;
          using SignalT  = Signal<IdxT>;
          using Variable = typename OverheadData<ScalarT, IdxT>::MonitorableVariables;

          OverheadMonitor(std::string label, const std::set<Variable>& variables)
            : label_(std::move(label))
          {
            for (const auto variable : variables)
            {
              variables_.push_back(variable);
            }
          }

          void setMatrix(Variable       variable,
                         SignalT        signal,
                         const ScalarT* values)
          {
            if (!selected(variable))
            {
              return;
            }

            const std::string prefix = variableName(variable);
            for (IdxT i = 0; i < signal.rows; ++i)
            {
              for (IdxT j = 0; j < signal.cols; ++j)
              {
                const IdxT offset = signal.offset + i * signal.cols + j;
                entries_.push_back(
                    {prefix + "_" + std::to_string(i) + "_" + std::to_string(j),
                     &values[static_cast<size_t>(offset)]});
              }
            }
          }

          void setComplexMatrix(Variable       variable,
                                SignalT        real_signal,
                                SignalT        imag_signal,
                                const ScalarT* values)
          {
            if (!selected(variable))
            {
              return;
            }

            const std::string prefix = variableName(variable);
            for (IdxT i = 0; i < real_signal.rows; ++i)
            {
              for (IdxT j = 0; j < real_signal.cols; ++j)
              {
                const IdxT real_offset = real_signal.offset + i * real_signal.cols + j;
                const IdxT imag_offset = imag_signal.offset + i * imag_signal.cols + j;
                const auto suffix      = std::to_string(i) + "_" + std::to_string(j);
                entries_.push_back(
                    {prefix + "_real_" + suffix,
                     &values[static_cast<size_t>(real_offset)]});
                entries_.push_back(
                    {prefix + "_imag_" + suffix,
                     &values[static_cast<size_t>(imag_offset)]});
              }
            }
          }

          void setVector(Variable variable, SignalT signal, const ScalarT* values)
          {
            if (!selected(variable))
            {
              return;
            }

            const std::string prefix = variableName(variable);
            for (IdxT i = 0; i < signal.rows; ++i)
            {
              entries_.push_back({prefix + "_" + std::to_string(i),
                                  &values[static_cast<size_t>(signal.offset + i)]});
            }
          }

          void setComplexVector(Variable       variable,
                                SignalT        real_signal,
                                SignalT        imag_signal,
                                const ScalarT* values)
          {
            if (!selected(variable))
            {
              return;
            }

            const std::string prefix = variableName(variable);
            for (IdxT i = 0; i < real_signal.rows; ++i)
            {
              entries_.push_back({prefix + "_real_" + std::to_string(i),
                                  &values[static_cast<size_t>(real_signal.offset + i)]});
              entries_.push_back({prefix + "_imag_" + std::to_string(i),
                                  &values[static_cast<size_t>(imag_signal.offset + i)]});
            }
          }

          bool empty() const override
          {
            return entries_.empty();
          }

          bool selectedAny(std::initializer_list<Variable> variables) const
          {
            for (const auto variable : variables)
            {
              if (selected(variable))
              {
                return true;
              }
            }
            return false;
          }

        private:
          using Csv  = VariableMonitorBase::Csv;
          using Json = VariableMonitorBase::Json;
          using Yaml = VariableMonitorBase::Yaml;

          struct Entry
          {
            std::string    name;
            const ScalarT* value{nullptr};
          };

          static std::string variableName(Variable variable)
          {
            switch (variable)
            {
            case Variable::R:
              return "R";
            case Variable::L:
              return "L";
            case Variable::G:
              return "G";
            case Variable::C:
              return "C";
            case Variable::Tv:
              return "Tv";
            case Variable::Ti:
              return "Ti";
            case Variable::Alpha:
              return "Alpha";
            case Variable::Beta:
              return "Beta";
            case Variable::Tau:
              return "Tau";
            case Variable::H:
              return "H";
            case Variable::Yc:
              return "Yc";
            case Variable::Zc:
              return "Zc";
            }
            return {};
          }

          bool selected(Variable variable) const
          {
            for (const auto selected_variable : variables_)
            {
              if (selected_variable == variable)
              {
                return true;
              }
            }
            return false;
          }

          std::string value(const Entry& entry) const
          {
            return Model::VariableMonitorDetail::formatReal(static_cast<RealT>(*entry.value));
          }

          void appendHeader(std::string& out, Csv csv) const override
          {
            for (const auto& entry : entries_)
            {
              out += csv.delim;
              out += label_;
              out += '_';
              out += entry.name;
            }
          }

          void append(std::string& out, Csv csv) const override
          {
            for (const auto& entry : entries_)
            {
              out += csv.delim;
              out += value(entry);
            }
          }

          void append(std::string& out, Json) const override
          {
            out += "    \"" + label_ + "\": {\n";
            for (size_t i = 0; i < entries_.size(); ++i)
            {
              const auto& entry  = entries_[i];
              out               += "      \"" + entry.name + "\": " + value(entry);
              out               += (i + 1 == entries_.size()) ? "\n" : ",\n";
            }
            out += "    }";
          }

          void append(std::string& out, Yaml) const override
          {
            out += "    " + label_ + ":\n";
            for (const auto& entry : entries_)
            {
              out += "      " + entry.name + ": " + value(entry) + "\n";
            }
          }

          std::string           label_;
          std::vector<Variable> variables_;
          std::vector<Entry>    entries_;
        };
      } // namespace Detail

      template <typename ScalarT, typename IdxT>
      Overhead<ScalarT, IdxT>::Overhead(const OverheadData<ScalarT, IdxT>& data)
        : conductor_(data.conductor),
          tower_(data.tower, conductor_),
          path_(data.path, tower_),
          geometric_inductance_(tower_, conductor_),
          skin_effect_(conductor_),
          carson_(tower_, data.carson),
          series_impedance_(skin_effect_, geometric_inductance_, carson_),
          shunt_potential_(tower_, conductor_),
          shunt_admittance_(shunt_potential_),
          reduction_map_(
              ReductionMap<ScalarT, IdxT>::fromConductors(data.conductor)),
          series_reduction_(
              reduction_map_.identity()
                  ? nullptr
                  : std::make_unique<SeriesReduction<ScalarT, IdxT>>(
                        series_impedance_, reduction_map_)),
          shunt_reduction_(
              reduction_map_.identity()
                  ? nullptr
                  : std::make_unique<ShuntReduction<ScalarT, IdxT>>(
                        shunt_admittance_, reduction_map_)),
          gamma_(modalSeries(), modalShunt()),
          yc_(modalSeries(), modalShunt()),
          zc_(modalSeries(), modalShunt()),
          monitor_(std::make_unique<Detail::OverheadMonitor<ScalarT, IdxT>>(
              "Overhead", data.monitored_variables)),
          monitor_controller_(std::make_unique<MonitorControllerT>("omega", omega_))
      {
        elements_ = {&conductor_,
                     &tower_,
                     &path_,
                     &geometric_inductance_,
                     &skin_effect_,
                     &carson_,
                     &series_impedance_,
                     &shunt_potential_,
                     &shunt_admittance_};
        if (series_reduction_)
        {
          elements_.push_back(series_reduction_.get());
          elements_.push_back(shunt_reduction_.get());
        }
        elements_.insert(elements_.end(), {&gamma_, &yc_, &zc_});
        for (const auto& sink : data.monitor_sink)
        {
          monitor_controller_->addSink(sink);
        }
      }

      template <typename ScalarT, typename IdxT>
      Overhead<ScalarT, IdxT>::~Overhead()
      {
        delete csr_jac_;
      }

      template <typename ScalarT, typename IdxT>
      int Overhead<ScalarT, IdxT>::allocate()
      {
        size_ = SystemAssembly::totalSize(elements_);

        y_.resize(size_);
        y_.setToZero();
        yp_.resize(size_);
        yp_.setToZero();
        f_.resize(size_);
        f_.setToZero();
        abs_tol_.resize(size_);
        abs_tol_.setToZero();
        tag_.assign(static_cast<size_t>(size_), false);
        gidx_.resize(static_cast<size_t>(size_));
        std::iota(gidx_.begin(), gidx_.end(), IdxT{0});

        SystemAssembly::bindAll(elements_, y_.getData(), yp_.getData(), f_.getData(), gidx_.data());
        for (auto* element : elements_)
        {
          element->setCoordinate(omega_, alpha_);
        }

        initialize();
        tagDifferentiable();
        evaluateResidual();
        evaluateJacobian();
        initializeMonitor();
        startMonitor();

        return 0;
      }

      template <typename ScalarT, typename IdxT>
      int Overhead<ScalarT, IdxT>::initialize()
      {
        for (auto* element : elements_)
        {
          element->initialize();
        }
        return 0;
      }

      template <typename ScalarT, typename IdxT>
      int Overhead<ScalarT, IdxT>::tagDifferentiable()
      {
        for (auto* element : elements_)
        {
          element->tagDifferentiable(tag_);
        }
        return 0;
      }

      template <typename ScalarT, typename IdxT>
      int Overhead<ScalarT, IdxT>::evaluateResidual()
      {
        for (auto* element : elements_)
        {
          element->evaluateResidual();
        }
        return 0;
      }

      template <typename ScalarT, typename IdxT>
      int Overhead<ScalarT, IdxT>::evaluateJacobian()
      {
        for (auto* element : elements_)
        {
          element->evaluateJacobian();
        }

        if (csr_jac_ == nullptr)
        {
          SystemAssembly::buildCsr(elements_, size_, csr_jac_, csr_map_);
          nnz_ = csr_jac_->getNnz();
        }

        SystemAssembly::updateCsrValues(elements_, csr_jac_, csr_map_);

        return 0;
      }

      template <typename ScalarT, typename IdxT>
      void Overhead<ScalarT, IdxT>::updateTime(RealT t, RealT a)
      {
        omega_ = t;
        alpha_ = a;
        for (auto* element : elements_)
        {
          element->setCoordinate(t, a);
        }
      }

      template <typename ScalarT, typename IdxT>
      void Overhead<ScalarT, IdxT>::initializeMonitor()
      {
        if (monitor_initialized_)
        {
          return;
        }

        using Monitor  = Detail::OverheadMonitor<ScalarT, IdxT>;
        using Variable = typename Monitor::Variable;
        auto* monitor  = static_cast<Monitor*>(monitor_.get());
        monitor->setMatrix(Variable::R, series_impedance_.R(), y_.getData());
        monitor->setMatrix(Variable::L, series_impedance_.L(), y_.getData());
        monitor->setMatrix(Variable::G, shunt_admittance_.G(), y_.getData());
        monitor->setMatrix(Variable::C, shunt_admittance_.C(), y_.getData());
        monitor->setComplexMatrix(Variable::Yc, yc_.Gc(), yc_.Bc(), y_.getData());
        monitor->setComplexMatrix(Variable::Zc, zc_.Rc(), zc_.Xc(), y_.getData());

        // Modal variables are observations of gamma, not DAE states:
        // they live in the decomposition's buffer and are refreshed at
        // every emission by printMonitoredVariables().
        if (monitor->selectedAny({Variable::Tv,
                                  Variable::Ti,
                                  Variable::Alpha,
                                  Variable::Beta,
                                  Variable::Tau,
                                  Variable::H}))
        {
          modes_ = std::make_unique<ModalDecomposition<ScalarT, IdxT>>(
              gamma_.alpha().rows, path_.length());
          monitor->setComplexMatrix(Variable::Tv,
                                    modes_->tvReal(),
                                    modes_->tvImag(),
                                    modes_->data());
          monitor->setComplexMatrix(Variable::Ti,
                                    modes_->tiReal(),
                                    modes_->tiImag(),
                                    modes_->data());
          monitor->setVector(Variable::Alpha, modes_->alphaM(), modes_->data());
          monitor->setVector(Variable::Beta, modes_->betaM(), modes_->data());
          monitor->setVector(Variable::Tau, modes_->tau(), modes_->data());
          monitor->setComplexVector(Variable::H,
                                    modes_->hReal(),
                                    modes_->hImag(),
                                    modes_->data());
        }

        monitor_controller_->addMonitor(monitor_.get());
        monitor_initialized_ = true;
      }

      template <typename ScalarT, typename IdxT>
      bool Overhead<ScalarT, IdxT>::monitoring() const
      {
        return monitor_controller_ && !monitor_controller_->empty();
      }

      template <typename ScalarT, typename IdxT>
      void Overhead<ScalarT, IdxT>::startMonitor()
      {
        monitor_controller_->start();
      }

      template <typename ScalarT, typename IdxT>
      void Overhead<ScalarT, IdxT>::stopMonitor()
      {
        monitor_controller_->stop();
      }

      template <typename ScalarT, typename IdxT>
      void Overhead<ScalarT, IdxT>::printMonitoredVariables() const
      {
        if (modes_)
        {
          const int status = modes_->decompose(omega_,
                                               y_.getData()
                                                   + gamma_.alpha().offset,
                                               y_.getData()
                                                   + gamma_.beta().offset);
          if (status != 0)
          {
            throw std::runtime_error(
                "Modal decomposition failed during monitor emission");
          }
        }
        monitor_controller_->print();
      }

      template <typename ScalarT, typename IdxT>
      const Model::VariableMonitorBase* Overhead<ScalarT, IdxT>::getMonitor() const
      {
        return monitor_.get();
      }
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
