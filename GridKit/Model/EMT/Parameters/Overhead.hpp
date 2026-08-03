/**
 * @file Overhead.hpp
 *
 * @brief Aggregate frequency-domain model for an overhead transmission line. It
 * implements Model::Evaluator so IDA can evaluate it over angular frequency,
 * owns the global state, and wires its parameter elements per the Parameters
 * dependency diagram.
 *
 */

#pragma once

#include <cmath>
#include <memory>
#include <vector>

#include <GridKit/Model/EMT/Element.hpp>
#include <GridKit/Model/EMT/Parameters/Effects/Carson/Carson.hpp>
#include <GridKit/Model/EMT/Parameters/Effects/GeometricInductance/GeometricInductance.hpp>
#include <GridKit/Model/EMT/Parameters/Effects/Reduction/ReductionMap.hpp>
#include <GridKit/Model/EMT/Parameters/Effects/Reduction/SeriesReduction.hpp>
#include <GridKit/Model/EMT/Parameters/Effects/Reduction/ShuntReduction.hpp>
#include <GridKit/Model/EMT/Parameters/Effects/SeriesImpedance/SeriesImpedance.hpp>
#include <GridKit/Model/EMT/Parameters/Effects/ShuntAdmittance/ShuntAdmittance.hpp>
#include <GridKit/Model/EMT/Parameters/Effects/ShuntPotential/ShuntPotential.hpp>
#include <GridKit/Model/EMT/Parameters/Effects/SkinEffect/SkinEffect.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Conductor/Conductor.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Path/Path.hpp>
#include <GridKit/Model/EMT/Parameters/Geometry/Tower/Tower.hpp>
#include <GridKit/Model/EMT/Parameters/OverheadData.hpp>
#include <GridKit/Model/EMT/Parameters/Response/Gamma/Gamma.hpp>
#include <GridKit/Model/EMT/Parameters/Response/H/H.hpp>
#include <GridKit/Model/EMT/Parameters/Response/Tau/Tau.hpp>
#include <GridKit/Model/EMT/Parameters/Response/Yc/Yc.hpp>
#include <GridKit/Model/EMT/Parameters/Response/Zc/Zc.hpp>
#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Model/VariableMonitor.hpp>
#include <GridKit/Utilities/Errors.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      template <typename scalar_type, typename index_type>
      class Overhead : public Model::Evaluator<scalar_type, index_type>
      {
      public:
        using ScalarT            = scalar_type;
        using IdxT               = index_type;
        using RealT              = typename Model::Evaluator<ScalarT, IdxT>::RealT;
        using VectorT            = typename Model::Evaluator<ScalarT, IdxT>::VectorT;
        using CsrMatrixT         = typename Model::Evaluator<ScalarT, IdxT>::CsrMatrixT;
        using ElementT           = Element<ScalarT, IdxT>;
        using MonitorControllerT = Model::VariableMonitorController<ScalarT>;

        explicit Overhead(const OverheadData<ScalarT, IdxT>& data);
        ~Overhead() override;

        Overhead(const Overhead&)            = delete;
        Overhead& operator=(const Overhead&) = delete;

        int allocate() override;
        int initialize() override;
        int tagDifferentiable() override;
        int evaluateResidual() override;
        int evaluateJacobian() override;

        void                              updateTime(RealT t, RealT a) override;
        bool                              monitoring() const override;
        void                              startMonitor() override;
        void                              stopMonitor() override;
        void                              printMonitoredVariables() const override;
        const Model::VariableMonitorBase* getMonitor() const override;

        const Tower<ScalarT, IdxT>& tower() const
        {
          return tower_;
        }

        const Conductor<ScalarT, IdxT>& conductor() const
        {
          return conductor_;
        }

        const Path<ScalarT, IdxT>& path() const
        {
          return path_;
        }

        const GeometricInductance<ScalarT, IdxT>& geometricInductance() const
        {
          return geometric_inductance_;
        }

        const SkinEffect<ScalarT, IdxT>& skinEffect() const
        {
          return skin_effect_;
        }

        const Carson<ScalarT, IdxT>& carson() const
        {
          return carson_;
        }

        const SeriesImpedance<ScalarT, IdxT>& seriesImpedance() const
        {
          return series_impedance_;
        }

        const ShuntPotential<ScalarT, IdxT>& shuntPotential() const
        {
          return shunt_potential_;
        }

        const ShuntAdmittance<ScalarT, IdxT>& shuntAdmittance() const
        {
          return shunt_admittance_;
        }

        const ReductionMap<ScalarT, IdxT>& reductionMap() const
        {
          return reduction_map_;
        }

        /// Null when the reduction is the identity and the modal stage
        /// consumes the full-conductor matrices directly.
        const SeriesReduction<ScalarT, IdxT>* seriesReduction() const
        {
          return series_reduction_.get();
        }

        /// Null when the reduction is the identity.
        const ShuntReduction<ScalarT, IdxT>* shuntReduction() const
        {
          return shunt_reduction_.get();
        }

        const Gamma<ScalarT, IdxT>& gamma() const
        {
          return gamma_;
        }

        const Tau<ScalarT, IdxT>& tau() const
        {
          return tau_;
        }

        const H<ScalarT, IdxT>& h() const
        {
          return h_;
        }

        const Yc<ScalarT, IdxT>& yc() const
        {
          return yc_;
        }

        const Zc<ScalarT, IdxT>& zc() const
        {
          return zc_;
        }

        RealT omega() const
        {
          return omega_;
        }

        IdxT size() override
        {
          return size_;
        }

        IdxT nnz() override
        {
          return nnz_;
        }

        bool hasJacobian() override
        {
          return true;
        }

        CsrMatrixT* getCsrJacobian() const override
        {
          return csr_jac_;
        }

        int setAbsoluteTolerance(RealT rel_tol) override
        {
          // Parameter magnitudes span many decades, so a flat absolute
          // tolerance either starves the large variables or demands
          // sub-machine-precision corrections on the small ones. Scale each
          // tolerance to the variable's initialization magnitude instead.
          auto*       tol = abs_tol_.getData();
          const auto* y   = y_.getData();
          for (IdxT i = 0; i < size_; ++i)
          {
            tol[i] = static_cast<ScalarT>(rel_tol) * (1.0 + std::abs(y[i]));
          }
          abs_tol_.setDataUpdated();
          return 0;
        }

        VectorT& absoluteTolerance() override
        {
          return abs_tol_;
        }

        const VectorT& absoluteTolerance() const override
        {
          return abs_tol_;
        }

        VectorT& y() override
        {
          return y_;
        }

        const VectorT& y() const override
        {
          return y_;
        }

        VectorT& yp() override
        {
          return yp_;
        }

        const VectorT& yp() const override
        {
          return yp_;
        }

        std::vector<bool>& tag() override
        {
          return tag_;
        }

        const std::vector<bool>& tag() const override
        {
          return tag_;
        }

        VectorT& getResidual() override
        {
          return f_;
        }

        const VectorT& getResidual() const override
        {
          return f_;
        }

        // The frequency-sweep aggregate does not provide quadrature, adjoint, or
        // parameter sensitivity interfaces.
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

        [[noreturn]] IdxT sizeQuadrature() override
        {
          throw NotImplementedError(__func__);
        }

        [[noreturn]] IdxT sizeParams() override
        {
          throw NotImplementedError(__func__);
        }

        [[noreturn]] VectorT& yB() override
        {
          throw NotImplementedError(__func__);
        }

        [[noreturn]] const VectorT& yB() const override
        {
          throw NotImplementedError(__func__);
        }

        [[noreturn]] VectorT& ypB() override
        {
          throw NotImplementedError(__func__);
        }

        [[noreturn]] const VectorT& ypB() const override
        {
          throw NotImplementedError(__func__);
        }

        [[noreturn]] VectorT& param() override
        {
          throw NotImplementedError(__func__);
        }

        [[noreturn]] const VectorT& param() const override
        {
          throw NotImplementedError(__func__);
        }

        [[noreturn]] VectorT& param_up() override
        {
          throw NotImplementedError(__func__);
        }

        [[noreturn]] const VectorT& param_up() const override
        {
          throw NotImplementedError(__func__);
        }

        [[noreturn]] VectorT& param_lo() override
        {
          throw NotImplementedError(__func__);
        }

        [[noreturn]] const VectorT& param_lo() const override
        {
          throw NotImplementedError(__func__);
        }

        [[noreturn]] VectorT& getIntegrand() override
        {
          throw NotImplementedError(__func__);
        }

        [[noreturn]] const VectorT& getIntegrand() const override
        {
          throw NotImplementedError(__func__);
        }

        [[noreturn]] VectorT& getAdjointResidual() override
        {
          throw NotImplementedError(__func__);
        }

        [[noreturn]] const VectorT& getAdjointResidual() const override
        {
          throw NotImplementedError(__func__);
        }

        [[noreturn]] VectorT& getAdjointIntegrand() override
        {
          throw NotImplementedError(__func__);
        }

        [[noreturn]] const VectorT& getAdjointIntegrand() const override
        {
          throw NotImplementedError(__func__);
        }

      private:
        using NotImplementedError = GridKit::Utilities::NotImplementedError;

        void initializeMonitor();

        Conductor<ScalarT, IdxT>           conductor_;
        Tower<ScalarT, IdxT>               tower_;
        Path<ScalarT, IdxT>                path_;
        GeometricInductance<ScalarT, IdxT> geometric_inductance_;
        SkinEffect<ScalarT, IdxT>          skin_effect_;
        Carson<ScalarT, IdxT>              carson_;
        SeriesImpedance<ScalarT, IdxT>     series_impedance_;
        ShuntPotential<ScalarT, IdxT>      shunt_potential_;
        ShuntAdmittance<ScalarT, IdxT>     shunt_admittance_;

        /// Terminal incidence and the optional reduction stage; absent
        /// when every conductor is its own terminal.
        ReductionMap<ScalarT, IdxT>                     reduction_map_;
        std::unique_ptr<SeriesReduction<ScalarT, IdxT>> series_reduction_;
        std::unique_ptr<ShuntReduction<ScalarT, IdxT>>  shunt_reduction_;

        /// The series pair the modal stage consumes: the reduction when
        /// present, the full-conductor assembly otherwise.
        const SeriesView<ScalarT, IdxT>& modalSeries() const
        {
          if (series_reduction_)
          {
            return *series_reduction_;
          }
          return series_impedance_;
        }

        /// The shunt pair the modal stage consumes.
        const ShuntView<ScalarT, IdxT>& modalShunt() const
        {
          if (shunt_reduction_)
          {
            return *shunt_reduction_;
          }
          return shunt_admittance_;
        }

        Gamma<ScalarT, IdxT>   gamma_;
        Tau<ScalarT, IdxT>     tau_;
        H<ScalarT, IdxT>       h_;
        Yc<ScalarT, IdxT>      yc_;
        Zc<ScalarT, IdxT>      zc_;
        std::vector<ElementT*> elements_;

        VectorT           y_;
        VectorT           yp_;
        VectorT           f_;
        VectorT           abs_tol_;
        std::vector<bool> tag_;
        std::vector<IdxT> gidx_;

        IdxT size_{0};
        IdxT nnz_{0};

        RealT omega_{0.0};
        RealT alpha_{0.0};

        CsrMatrixT*       csr_jac_{nullptr};
        std::vector<IdxT> csr_map_;

        std::unique_ptr<Model::VariableMonitorBase> monitor_;
        std::unique_ptr<MonitorControllerT>         monitor_controller_;
        bool                                        monitor_initialized_{false};
      };
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit

#include <GridKit/Model/EMT/Parameters/OverheadImpl.hpp>
