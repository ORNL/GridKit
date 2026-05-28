/**
 * @file BranchImpl.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Definition of a phasor dynamics branch model.
 *
 * The model uses Cartesian coordinates.
 *
 */

#include <cmath>
#include <complex>
#include <variant>

#include <magic_enum/magic_enum.hpp>

#include <GridKit/Model/PhasorDynamics/Branch/Branch.hpp>
#include <GridKit/Model/PhasorDynamics/Branch/BranchData.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Constructor for a line or off-nominal transformer branch
     *
     * Model size:
     * - Number of equations = 0
     * - Number of internal variables = 0
     */
    template <class ScalarT, typename IdxT>
    Branch<ScalarT, IdxT>::Branch(bus_type* bus1, bus_type* bus2)
      : bus1_(bus1),
        bus2_(bus2),
        R_(0.0),
        X_(0.01),
        G_(0.0),
        B_(0.0),
        tap_(1.0),
        phase_(0.0),
        bus1_id_(0),
        bus2_id_(0)
    {
      size_ = 0;
      setDerivedParams();
    }

    /**
     * @brief Construct a new Branch
     *
     * @tparam ScalarT - scalar type
     * @tparam IdxT    - matrix/vector index type
     * @param bus1 - pointer to bus-1
     * @param bus2 - pointer to bus-2
     * @param R - line series resistance
     * @param X - line series reactance
     * @param G - total shunt conductance
     * @param B - total shunt susceptance
     * @param tap - off-nominal tap magnitude on bus1 side
     * @param phase - phase shift angle in radians
     */
    template <class ScalarT, typename IdxT>
    Branch<ScalarT, IdxT>::Branch(bus_type* bus1,
                                  bus_type* bus2,
                                  RealT     R,
                                  RealT     X,
                                  RealT     G,
                                  RealT     B,
                                  RealT     tap,
                                  RealT     phase)
      : bus1_(bus1),
        bus2_(bus2),
        R_(R),
        X_(X),
        G_(G),
        B_(B),
        tap_(tap),
        phase_(phase),
        bus1_id_(0),
        bus2_id_(0)
    {
      size_ = 0;
      setDerivedParams();
    }

    template <class ScalarT, typename IdxT>
    Branch<ScalarT, IdxT>::Branch(bus_type* bus1, bus_type* bus2, const model_data_type& data)
      : bus1_(bus1),
        bus2_(bus2),
        monitor_(std::make_unique<MonitorT>(data))
    {
      initializeParameters(data);
      initializeMonitor();

      size_ = 0;
      setDerivedParams();
    }

    /**
     * @brief Destroy the Branch
     *
     * @tparam ScalarT
     * @tparam IdxT
     */
    template <class ScalarT, typename IdxT>
    Branch<ScalarT, IdxT>::~Branch()
    {
    }

    /**
     * @brief Set the component ID
     */
    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /*!
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::allocate()
    {
      wb_.resize(2);
      h_.resize(2);

      return 0;
    }

    /**
     * Initialization of the branch model
     *
     */
    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::initialize()
    {
      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::tagDifferentiable()
    {
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::verify() const
    {
      int ret = parameter_error_count_;

      auto check = [&](bool condition, const char* message)
      {
        if (!condition)
        {
          Log::error() << "Branch: " << message << '\n';
          ret += 1;
        }
      };

      check(bus1_ != nullptr, "bus1 pointer is null");
      check(bus2_ != nullptr, "bus2 pointer is null");

      check(std::isfinite(R_), "R must be finite");
      check(std::isfinite(X_), "X must be finite");
      check(std::isfinite(G_), "G must be finite");
      check(std::isfinite(B_), "B must be finite");
      check(std::isfinite(tap_), "tap must be finite");
      check(std::isfinite(phase_), "phase must be finite");
      check(R_ * R_ + X_ * X_ > RealT{0.0}, "R and X cannot both be zero");
      check(tap_ > RealT{0.0}, "tap must be positive");

      return ret;
    }

    template <class ScalarT, typename IdxT>
    __attribute__((always_inline)) inline void Branch<ScalarT, IdxT>::addAdmittanceContribution(
        const AdmittanceBlock& y,
        const ScalarT&         Vr,
        const ScalarT&         Vi,
        ScalarT&               Ir,
        ScalarT&               Ii)
    {
      Ir += y.G * Vr - y.B * Vi;
      Ii += y.B * Vr + y.G * Vi;
    }

    template <class ScalarT, typename IdxT>
    __attribute__((always_inline)) inline void Branch<ScalarT, IdxT>::evaluateAdmittanceBlock(
        const AdmittanceBlock& y,
        const ScalarT*         wb,
        ScalarT*               h)
    {
      const ScalarT Vr = wb[0];
      const ScalarT Vi = wb[1];

      h[0] = y.G * Vr - y.B * Vi;
      h[1] = y.B * Vr + y.G * Vi;
    }

    /**
     * @brief Bus 1 residual contribution from bus 1 variables
     *
     */
    template <class ScalarT, typename IdxT>
    __attribute__((always_inline)) inline int Branch<ScalarT, IdxT>::evaluateBusResidual11(
        [[maybe_unused]] ScalarT* y,
        [[maybe_unused]] ScalarT* yp,
        ScalarT*                  wb,
        ScalarT*                  h)
    {
      evaluateAdmittanceBlock(y11_, wb, h);

      return 0;
    }

    /**
     * @brief Bus 1 residual contribution from bus 2 variables
     *
     */
    template <class ScalarT, typename IdxT>
    __attribute__((always_inline)) inline int Branch<ScalarT, IdxT>::evaluateBusResidual12(
        [[maybe_unused]] ScalarT* y,
        [[maybe_unused]] ScalarT* yp,
        ScalarT*                  wb,
        ScalarT*                  h)
    {
      evaluateAdmittanceBlock(y12_, wb, h);

      return 0;
    }

    /**
     * @brief Bus 2 residual contribution from bus 1 variables
     *
     */
    template <class ScalarT, typename IdxT>
    __attribute__((always_inline)) int Branch<ScalarT, IdxT>::evaluateBusResidual21(
        [[maybe_unused]] ScalarT* y,
        [[maybe_unused]] ScalarT* yp,
        ScalarT*                  wb,
        ScalarT*                  h)
    {
      evaluateAdmittanceBlock(y21_, wb, h);

      return 0;
    }

    /**
     * @brief Bus 2 residual contribution from bus 2 variables
     *
     */
    template <class ScalarT, typename IdxT>
    __attribute__((always_inline)) int Branch<ScalarT, IdxT>::evaluateBusResidual22(
        [[maybe_unused]] ScalarT* y,
        [[maybe_unused]] ScalarT* yp,
        ScalarT*                  wb,
        ScalarT*                  h)
    {
      evaluateAdmittanceBlock(y22_, wb, h);

      return 0;
    }

    /**
     * @brief Residual contribution of the branch is computed and pushed to the terminal buses.
     *
     */
    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::evaluateResidual()
    {
      ScalarT ir1{0.0};
      ScalarT ii1{0.0};
      ScalarT ir2{0.0};
      ScalarT ii2{0.0};

      terminalCurrent1(ir1, ii1);
      terminalCurrent2(ir2, ii2);

      Ir1() += ir1;
      Ii1() += ii1;
      Ir2() += ir2;
      Ii2() += ii2;

      return 0;
    }

    template <class ScalarT, typename IdxT>
    void Branch<ScalarT, IdxT>::terminalCurrent1(ScalarT& Ir, ScalarT& Ii)
    {
      Ir = ScalarT{0.0};
      Ii = ScalarT{0.0};

      addAdmittanceContribution(y11_, Vr1(), Vi1(), Ir, Ii);
      addAdmittanceContribution(y12_, Vr2(), Vi2(), Ir, Ii);
    }

    template <class ScalarT, typename IdxT>
    void Branch<ScalarT, IdxT>::terminalCurrent2(ScalarT& Ir, ScalarT& Ii)
    {
      Ir = ScalarT{0.0};
      Ii = ScalarT{0.0};

      addAdmittanceContribution(y21_, Vr1(), Vi1(), Ir, Ii);
      addAdmittanceContribution(y22_, Vr2(), Vi2(), Ir, Ii);
    }

    template <class ScalarT, typename IdxT>
    void Branch<ScalarT, IdxT>::initializeParameters(const model_data_type& data)
    {
      readRealParameter(data, model_data_type::Parameters::R, R_);
      readRealParameter(data, model_data_type::Parameters::X, X_);
      readRealParameter(data, model_data_type::Parameters::G, G_);
      readRealParameter(data, model_data_type::Parameters::B, B_);
      readRealParameter(data, model_data_type::Parameters::tap, tap_);
      readRealParameter(data, model_data_type::Parameters::phase, phase_);

      if (data.ports.contains(model_data_type::Ports::bus1))
      {
        bus1_id_ = data.ports.at(model_data_type::Ports::bus1);
      }

      if (data.ports.contains(model_data_type::Ports::bus2))
      {
        bus2_id_ = data.ports.at(model_data_type::Ports::bus2);
      }
    }

    template <class ScalarT, typename IdxT>
    bool Branch<ScalarT, IdxT>::readRealParameter(const model_data_type&               data,
                                                  typename model_data_type::Parameters parameter,
                                                  RealT&                               target)
    {
      if (!data.parameters.contains(parameter))
      {
        return false;
      }

      const auto& value = data.parameters.at(parameter);
      if (const auto* real_value = std::get_if<RealT>(&value))
      {
        target = *real_value;
        return true;
      }

      if (const auto* integer_value = std::get_if<IdxT>(&value))
      {
        target = static_cast<RealT>(*integer_value);
        return true;
      }

      Log::error() << "Branch: parameter " << magic_enum::enum_name(parameter)
                   << " must be numeric\n";
      parameter_error_count_ += 1;
      return false;
    }

    template <class ScalarT, typename IdxT>
    const Model::VariableMonitorBase* Branch<ScalarT, IdxT>::getMonitor() const
    {
      return monitor_.get();
    }

    template <class ScalarT, typename IdxT>
    void Branch<ScalarT, IdxT>::initializeMonitor()
    {
      using Variable = typename model_data_type::MonitorableVariables;
      monitor_->set(Variable::ir1, [this]
                    {
                      ScalarT Ir;
                      ScalarT Ii;
                      terminalCurrent1(Ir, Ii);
                      return Ir; });
      monitor_->set(Variable::ii1, [this]
                    {
                      ScalarT Ir;
                      ScalarT Ii;
                      terminalCurrent1(Ir, Ii);
                      return Ii; });
      monitor_->set(Variable::im1, [this]
                    {
                      ScalarT Ir;
                      ScalarT Ii;
                      terminalCurrent1(Ir, Ii);
                      return std::sqrt(Ir * Ir + Ii * Ii); });
      monitor_->set(Variable::p1, [this]
                    {
                      ScalarT Ir;
                      ScalarT Ii;
                      terminalCurrent1(Ir, Ii);
                      return Vr1() * Ir + Vi1() * Ii; });
      monitor_->set(Variable::q1, [this]
                    {
                      ScalarT Ir;
                      ScalarT Ii;
                      terminalCurrent1(Ir, Ii);
                      return Vi1() * Ir - Vr1() * Ii; });
      monitor_->set(Variable::ir2, [this]
                    {
                      ScalarT Ir;
                      ScalarT Ii;
                      terminalCurrent2(Ir, Ii);
                      return Ir; });
      monitor_->set(Variable::ii2, [this]
                    {
                      ScalarT Ir;
                      ScalarT Ii;
                      terminalCurrent2(Ir, Ii);
                      return Ii; });
      monitor_->set(Variable::im2, [this]
                    {
                      ScalarT Ir;
                      ScalarT Ii;
                      terminalCurrent2(Ir, Ii);
                      return std::sqrt(Ir * Ir + Ii * Ii); });
      monitor_->set(Variable::p2, [this]
                    {
                      ScalarT Ir;
                      ScalarT Ii;
                      terminalCurrent2(Ir, Ii);
                      return Vr2() * Ir + Vi2() * Ii; });
      monitor_->set(Variable::q2, [this]
                    {
                      ScalarT Ir;
                      ScalarT Ii;
                      terminalCurrent2(Ir, Ii);
                      return Vi2() * Ir - Vr2() * Ii; });
    }

    /**
     * @brief Derived parameters
     *
     */
    template <class ScalarT, typename IdxT>
    void Branch<ScalarT, IdxT>::setDerivedParams()
    {
      y11_ = {};
      y12_ = {};
      y21_ = {};
      y22_ = {};

      const RealT denom = R_ * R_ + X_ * X_;
      if (denom == RealT{0.0} || tap_ == RealT{0.0})
      {
        return;
      }

      const RealT g       = R_ / denom;
      const RealT b       = -X_ / denom;
      const RealT inv_tap = RealT{1.0} / tap_;

      const std::complex<RealT> ybr{g, b};
      const std::complex<RealT> ysh{G_, B_};
      const std::complex<RealT> rotation{std::cos(phase_), std::sin(phase_)};
      const std::complex<RealT> ydiag = -(ybr + RealT{0.5} * ysh);

      setAdmittanceBlock(y11_, ydiag * inv_tap * inv_tap);
      setAdmittanceBlock(y12_, ybr * rotation * inv_tap);
      setAdmittanceBlock(y21_, ybr * std::conj(rotation) * inv_tap);
      setAdmittanceBlock(y22_, ydiag);
    }

    template <class ScalarT, typename IdxT>
    void Branch<ScalarT, IdxT>::setAdmittanceBlock(AdmittanceBlock&           block,
                                                   const std::complex<RealT>& y)
    {
      block.G = y.real();
      block.B = y.imag();
    }

  } // namespace PhasorDynamics
} // namespace GridKit
