#pragma once

#include <chrono>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

#include <GridKit/AutomaticDifferentiation/DependencyTracking/Variable.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentLibrary.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

namespace GridKit
{
  namespace Testing
  {
    using Log = ::GridKit::Utilities::Logger;

    /**
     * @brief Compares and times system Jacobians from Enzyme and DependencyTracking.
     */
    template <class ScalarT, typename IdxT>
    class SystemJacobianTests
    {
    private:
      using RealT      = ScalarT;
      using ClockT     = std::chrono::steady_clock;
      using DurationT  = std::chrono::duration<double, std::micro>;
      using JacobianT  = std::vector<DependencyTracking::Variable::DependencyMap>;
      using ModelDataT = PhasorDynamics::SystemModelData<RealT, IdxT>;

      /// Number of timed evaluations after setup and warm-up
      static constexpr int num_evaluations_ = 100;

      // Timing table, one entry per case (times in microseconds per call, setup once)
      std::vector<std::string> case_names_;
      std::vector<IdxT>        sizes_;
      std::vector<IdxT>        enzyme_nnz_;
      std::vector<double>      enzyme_setup_;
      std::vector<double>      enzyme_residual_;
      std::vector<double>      enzyme_jacobian_;
      std::vector<IdxT>        dependency_tracking_nnz_;
      std::vector<double>      dependency_tracking_setup_;
      std::vector<double>      dependency_tracking_residual_;
      std::vector<double>      dependency_tracking_jacobian_;

    public:
      /**
       * @brief Checks the system Jacobian of the case in `case_file`.
       */
      TestOutcome jacobian(const std::filesystem::path& case_file)
      {
        TestStatus success = true;

        // Looser than the component tests' 10 * eps, for the Math::antiwindup terms
        auto tol = 10000 * std::numeric_limits<RealT>::epsilon();

        const std::string case_name = case_file.stem().stem().string();

        // Timing table entry, filled by EnzymeJacobian() and DependencyTrackingJacobian()
        case_names_.push_back(case_name);
        sizes_.push_back(0);
        enzyme_nnz_.push_back(0);
        enzyme_setup_.push_back(0.0);
        enzyme_residual_.push_back(0.0);
        enzyme_jacobian_.push_back(0.0);
        dependency_tracking_nnz_.push_back(0);
        dependency_tracking_setup_.push_back(0.0);
        dependency_tracking_residual_.push_back(0.0);
        dependency_tracking_jacobian_.push_back(0.0);

        try
        {
          auto data = PhasorDynamics::parseSystemModelData(case_file);

          // No monitor output
          data.monitor_sink.clear();

          // Jacobian via Enzyme
          JacobianT enzyme_jacobian = EnzymeJacobian(data);

          // Jacobian via DependencyTracking
          JacobianT dependency_tracking_jacobian = DependencyTrackingJacobian(data);

          // Compare DependencyTracking dependencies to Enzyme's
          success *= (dependency_tracking_jacobian.size() == enzyme_jacobian.size());
          for (IdxT i = 0; i < dependency_tracking_jacobian.size() && i < enzyme_jacobian.size(); ++i)
          {
            success *= (GridKit::Testing::isEqual(dependency_tracking_jacobian[i], enzyme_jacobian[i], tol));
          }
        }
        catch (const std::exception& error)
        {
          Log::error() << case_name << ": " << error.what() << '\n';
          success = false;
        }

        return success.report(("jacobian " + case_name).c_str());
      }

      /**
       * @brief Prints the timings of all checked cases, sorted by number of non-zeros.
       */
      void printTimings() const
      {
        // Case column fits the longest case name
        IdxT case_width = 4;
        for (const auto& case_name : case_names_)
        {
          if (case_name.size() > case_width)
          {
            case_width = case_name.size();
          }
        }
        case_width += 2;

        auto& output_stream = Log::summary();

        // Restored at the end, the stream is shared
        const auto flags     = output_stream.flags();
        const auto precision = output_stream.precision();

        output_stream << "Times in us per call (setup once); ratio = DependencyTracking / Enzyme.\n";
        output_stream << std::left << std::setw(static_cast<int>(case_width)) << "case"
                      << std::setw(20) << "method" << std::right
                      << std::setw(7) << "size"
                      << std::setw(8) << "nnz"
                      << std::setw(11) << "setup" << std::setw(8) << "ratio"
                      << std::setw(10) << "residual" << std::setw(8) << "ratio"
                      << std::setw(10) << "Jacobian" << std::setw(8) << "ratio" << '\n';

        output_stream << std::fixed << std::setprecision(1);

        std::vector<bool> printed(case_names_.size(), false);
        for (IdxT row = 0; row < case_names_.size(); ++row)
        {
          // Case with the fewest non-zeros among those not printed yet
          IdxT i = case_names_.size();
          for (IdxT j = 0; j < case_names_.size(); ++j)
          {
            if (!printed[j] && (i == case_names_.size() || enzyme_nnz_[j] < enzyme_nnz_[i]))
            {
              i = j;
            }
          }
          printed[i] = true;

          output_stream << std::left << std::setw(static_cast<int>(case_width)) << case_names_[i]
                        << std::setw(20) << "Enzyme" << std::right
                        << std::setw(7) << sizes_[i]
                        << std::setw(8) << enzyme_nnz_[i]
                        << std::setw(11) << enzyme_setup_[i] << std::setw(8) << ""
                        << std::setw(10) << enzyme_residual_[i] << std::setw(8) << ""
                        << std::setw(10) << enzyme_jacobian_[i] << '\n';

          output_stream << std::left << std::setw(static_cast<int>(case_width)) << case_names_[i]
                        << std::setw(20) << "DependencyTracking" << std::right
                        << std::setw(7) << sizes_[i]
                        << std::setw(8) << dependency_tracking_nnz_[i]
                        << std::setw(11) << dependency_tracking_setup_[i]
                        << std::setw(8) << dependency_tracking_setup_[i] / enzyme_setup_[i]
                        << std::setw(10) << dependency_tracking_residual_[i]
                        << std::setw(8) << dependency_tracking_residual_[i] / enzyme_residual_[i]
                        << std::setw(10) << dependency_tracking_jacobian_[i]
                        << std::setw(8) << dependency_tracking_jacobian_[i] / enzyme_jacobian_[i] << '\n';
        }

        output_stream.flags(flags);
        output_stream.precision(precision);
      }

    private:
      JacobianT DependencyTrackingJacobian(const ModelDataT& data)
      {
        PhasorDynamics::SystemModel<DependencyTracking::Variable, IdxT> system(data);
        system.updateTime(0.0, 1.0); // Set alpha to 1.0 to verify d/dy' terms

        // Setup: the residual tracks the dependencies the Jacobian is built from
        auto start = ClockT::now();
        system.allocate();
        system.initialize();
        system.evaluateResidual();
        system.evaluateJacobian();
        auto stop = ClockT::now();

        const DurationT setup = stop - start;

        // Warm-up, untimed
        system.evaluateResidual();
        system.evaluateJacobian();

        start = ClockT::now();
        for (int i = 0; i < num_evaluations_; ++i)
        {
          system.evaluateResidual();
        }
        stop = ClockT::now();

        const DurationT residual = (stop - start) / num_evaluations_;

        start = ClockT::now();
        for (int i = 0; i < num_evaluations_; ++i)
        {
          system.evaluateJacobian();
        }
        stop = ClockT::now();

        const DurationT jacobian = (stop - start) / num_evaluations_;

        auto* system_jacobian = system.getCsrJacobian();

        sizes_.back()                        = system.size();
        dependency_tracking_nnz_.back()      = system_jacobian->getNnz();
        dependency_tracking_setup_.back()    = setup.count();
        dependency_tracking_residual_.back() = residual.count();
        dependency_tracking_jacobian_.back() = jacobian.count();

        auto& output_stream = Log::misc();
        output_stream << "Sparse Csr Matrix: System Jacobian with DependencyTracking\n";
        system_jacobian->print(output_stream);

        return GridKit::Testing::MapFromCsr(system_jacobian);
      }

      JacobianT EnzymeJacobian(const ModelDataT& data)
      {
        PhasorDynamics::SystemModel<ScalarT, IdxT> system(data);
        system.updateTime(0.0, 1.0); // Set alpha to 1.0 to verify d/dy' terms

        // Setup: allocate() also initializes and evaluates the residual and the first Jacobian
        auto start = ClockT::now();
        system.allocate();
        auto stop = ClockT::now();

        const DurationT setup = stop - start;

        // Warm-up, untimed
        system.evaluateResidual();
        system.evaluateJacobian();

        start = ClockT::now();
        for (int i = 0; i < num_evaluations_; ++i)
        {
          system.evaluateResidual();
        }
        stop = ClockT::now();

        const DurationT residual = (stop - start) / num_evaluations_;

        start = ClockT::now();
        for (int i = 0; i < num_evaluations_; ++i)
        {
          system.evaluateJacobian();
        }
        stop = ClockT::now();

        const DurationT jacobian = (stop - start) / num_evaluations_;

        auto* system_jacobian = system.getCsrJacobian();

        sizes_.back()           = system.size();
        enzyme_nnz_.back()      = system_jacobian->getNnz();
        enzyme_setup_.back()    = setup.count();
        enzyme_residual_.back() = residual.count();
        enzyme_jacobian_.back() = jacobian.count();

        auto& output_stream = Log::misc();
        output_stream << "Sparse Csr Matrix: System Jacobian with Enzyme\n";
        system_jacobian->print(output_stream);

        return GridKit::Testing::MapFromCsr(system_jacobian);
      }
    };
  } // namespace Testing
} // namespace GridKit
