/**
 * @file BranchImpl.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Definition of a phasor dynamics branch model.
 *
 * The model uses Cartesian coordinates.
 *
 */

#include <GridKit/Model/PhasorDynamics/Branch/Branch.hpp>
#include <GridKit/Model/PhasorDynamics/Branch/BranchData.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/VariableMonitorImpl.hpp>
#include <GridKit/Utilities/ConfigurationChecks.hpp>
#include <GridKit/Utilities/ParameterReader.hpp>

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
    template <typename scalar_type, typename index_type>
    Branch<scalar_type, index_type>::Branch(BusT* bus1, BusT* bus2)
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
     * @param bus1 - pointer to bus-1
     * @param bus2 - pointer to bus-2
     * @param R - line series resistance
     * @param X - line series reactance
     * @param G - total line shunt conductance
     * @param B - total line shunt susceptance
     * @param tap - off-nominal tap magnitude on bus1 side
     * @param phase - phase shift angle in radians
     */
    template <typename scalar_type, typename index_type>
    Branch<scalar_type, index_type>::Branch(BusT* bus1,
                                            BusT* bus2,
                                            RealT R,
                                            RealT X,
                                            RealT G,
                                            RealT B,
                                            RealT tap,
                                            RealT phase)
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

    template <typename scalar_type, typename index_type>
    Branch<scalar_type, index_type>::Branch(BusT* bus1, BusT* bus2, const ModelDataT& data)
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
     */
    template <typename scalar_type, typename index_type>
    Branch<scalar_type, index_type>::~Branch()
    {
    }

    /**
     * @brief Set the component ID
     */
    template <typename scalar_type, typename index_type>
    int Branch<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /*!
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <typename scalar_type, typename index_type>
    int Branch<scalar_type, index_type>::allocate()
    {
      if (!allocated_)
      {
        this->allocateVectors(size_);
      }
      auto size = static_cast<std::size_t>(size_);

      variable_indices_.resize(size);
      residual_indices_.resize(size);

      wb_.resize(2);
      h_.resize(2);

      allocated_ = true;
      return 0;
    }

    /**
     * Initialization of the branch model
     *
     */
    template <typename scalar_type, typename index_type>
    int Branch<scalar_type, index_type>::initialize()
    {
      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <typename scalar_type, typename index_type>
    int Branch<scalar_type, index_type>::tagDifferentiable()
    {
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Branch<scalar_type, index_type>::verify() const
    {
      Utilities::ConfigurationChecks checks("Branch");

      checks.check(bus1_ != nullptr, "bus1 pointer is null");
      checks.check(bus2_ != nullptr, "bus2 pointer is null");

      checks.check(std::isfinite(R_), "R must be finite");
      checks.check(std::isfinite(X_), "X must be finite");
      checks.check(std::isfinite(G_), "G must be finite");
      checks.check(std::isfinite(B_), "B must be finite");
      checks.check(std::isfinite(Gmag_), "Gmag must be finite");
      checks.check(std::isfinite(Bmag_), "Bmag must be finite");
      checks.check(std::isfinite(tap_), "tap must be finite");
      checks.check(std::isfinite(phase_), "phase must be finite");
      checks.check(R_ * R_ + X_ * X_ > RealT{0.0}, "R and X cannot both be zero");
      checks.check(tap_ > RealT{0.0}, "tap must be positive");

      return parameter_error_count_ + checks.errorCount();
    }

    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline void Branch<scalar_type, index_type>::addAdmittanceContribution(
        const RealT   G,
        const RealT   B,
        const ScalarT Vr,
        const ScalarT Vi,
        ScalarT&      Ir,
        ScalarT&      Ii)
    {
      Ir += G * Vr - B * Vi;
      Ii += B * Vr + G * Vi;
    }

    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline void Branch<scalar_type, index_type>::evaluateAdmittanceBlock(
        const RealT    G,
        const RealT    B,
        const ScalarT* wb,
        ScalarT*       h)
    {
      const ScalarT Vr = wb[0];
      const ScalarT Vi = wb[1];

      h[0] = G * Vr - B * Vi;
      h[1] = B * Vr + G * Vi;
    }

    /**
     * @brief Compute the absolute tolerance for each variable in the model
     *
     * @param rel_tol The relative tolerance which can be used to pick the
     *        absolute tolerance.
     * @tparam scalar_type Scalar data type
     * @tparam index_type Index data type
     * @return int 0 if successful, non-zero otherwise.
     *
     * This represents a "noise" level close to zero for which pure relative
     * error cannot be used.
     */
    template <typename scalar_type, typename index_type>
    int Branch<scalar_type, index_type>::setAbsoluteTolerance(RealT)
    {
      return 0;
    }

    /**
     * @brief Bus 1 residual contribution from bus 1 variables
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline int Branch<scalar_type, index_type>::evaluateBusResidual11(
        [[maybe_unused]] const ScalarT* y,
        [[maybe_unused]] const ScalarT* yp,
        const ScalarT*                  wb,
        ScalarT*                        h)
    {
      evaluateAdmittanceBlock(g11_, b11_, wb, h);

      return 0;
    }

    /**
     * @brief Bus 1 residual contribution from bus 2 variables
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline int Branch<scalar_type, index_type>::evaluateBusResidual12(
        [[maybe_unused]] const ScalarT* y,
        [[maybe_unused]] const ScalarT* yp,
        const ScalarT*                  wb,
        ScalarT*                        h)
    {
      evaluateAdmittanceBlock(g12_, b12_, wb, h);

      return 0;
    }

    /**
     * @brief Bus 2 residual contribution from bus 1 variables
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int Branch<scalar_type, index_type>::evaluateBusResidual21(
        [[maybe_unused]] const ScalarT* y,
        [[maybe_unused]] const ScalarT* yp,
        const ScalarT*                  wb,
        ScalarT*                        h)
    {
      evaluateAdmittanceBlock(g21_, b21_, wb, h);

      return 0;
    }

    /**
     * @brief Bus 2 residual contribution from bus 2 variables
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int Branch<scalar_type, index_type>::evaluateBusResidual22(
        [[maybe_unused]] const ScalarT* y,
        [[maybe_unused]] const ScalarT* yp,
        const ScalarT*                  wb,
        ScalarT*                        h)
    {
      evaluateAdmittanceBlock(g22_, b22_, wb, h);

      return 0;
    }

    /**
     * @brief Residual contribution of the branch is computed and pushed to the terminal buses.
     *
     */
    template <typename scalar_type, typename index_type>
    int Branch<scalar_type, index_type>::evaluateResidual()
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

      if (bus1_->size() > 0)
      {
        bus1_->getResidual().setDataUpdated();
      }
      if (bus2_->size() > 0)
      {
        bus2_->getResidual().setDataUpdated();
      }

      return 0;
    }

    template <typename scalar_type, typename index_type>
    void Branch<scalar_type, index_type>::terminalCurrent1(ScalarT& Ir, ScalarT& Ii)
    {
      Ir = ScalarT{0.0};
      Ii = ScalarT{0.0};

      addAdmittanceContribution(g11_, b11_, Vr1(), Vi1(), Ir, Ii);
      addAdmittanceContribution(g12_, b12_, Vr2(), Vi2(), Ir, Ii);
    }

    template <typename scalar_type, typename index_type>
    void Branch<scalar_type, index_type>::terminalCurrent2(ScalarT& Ir, ScalarT& Ii)
    {
      Ir = ScalarT{0.0};
      Ii = ScalarT{0.0};

      addAdmittanceContribution(g21_, b21_, Vr1(), Vi1(), Ir, Ii);
      addAdmittanceContribution(g22_, b22_, Vr2(), Vi2(), Ir, Ii);
    }

    template <typename scalar_type, typename index_type>
    void Branch<scalar_type, index_type>::initializeParameters(const ModelDataT& data)
    {
      using Parameter = typename ModelDataT::Parameters;
      using Buses     = typename ModelDataT::Buses;

      Utilities::ConfigurationChecks checks("Branch");
      Utilities::ParameterReader     reader(data, checks);

      reader.loadReal(Parameter::R, R_);
      reader.loadReal(Parameter::X, X_);
      reader.loadReal(Parameter::G, G_);
      reader.loadReal(Parameter::B, B_);
      reader.loadReal(Parameter::Gmag, Gmag_);
      reader.loadReal(Parameter::Bmag, Bmag_);
      reader.loadReal(Parameter::tap, tap_);
      reader.loadReal(Parameter::phase, phase_);

      parameter_error_count_ = checks.errorCount();

      if (data.buses.contains(Buses::bus1))
      {
        bus1_id_ = data.buses.at(Buses::bus1);
      }

      if (data.buses.contains(Buses::bus2))
      {
        bus2_id_ = data.buses.at(Buses::bus2);
      }
    }

    template <typename scalar_type, typename index_type>
    const Model::VariableMonitorBase* Branch<scalar_type, index_type>::getMonitor() const
    {
      return monitor_.get();
    }

    template <typename scalar_type, typename index_type>
    void Branch<scalar_type, index_type>::initializeMonitor()
    {
      using Variable = typename ModelDataT::MonitorableVariables;
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
    template <typename scalar_type, typename index_type>
    void Branch<scalar_type, index_type>::setDerivedParams()
    {
      g11_ = RealT{0.0};
      b11_ = RealT{0.0};
      g12_ = RealT{0.0};
      b12_ = RealT{0.0};
      g21_ = RealT{0.0};
      b21_ = RealT{0.0};
      g22_ = RealT{0.0};
      b22_ = RealT{0.0};

      const RealT denom = R_ * R_ + X_ * X_;
      if (denom == RealT{0.0} || tap_ == RealT{0.0})
      {
        return;
      }

      const RealT g_br    = R_ / denom;
      const RealT b_br    = -X_ / denom;
      const RealT inv_tap = RealT{1.0} / tap_;
      const RealT cos_ph  = std::cos(phase_);
      const RealT sin_ph  = std::sin(phase_);

      const RealT g_diag = -g_br;
      const RealT b_diag = -b_br;

      g11_ = g_diag * inv_tap * inv_tap - RealT{0.5} * G_ - Gmag_;
      b11_ = b_diag * inv_tap * inv_tap - RealT{0.5} * B_ - Bmag_;

      g12_ = (g_br * cos_ph - b_br * sin_ph) * inv_tap;
      b12_ = (b_br * cos_ph + g_br * sin_ph) * inv_tap;

      g21_ = (g_br * cos_ph + b_br * sin_ph) * inv_tap;
      b21_ = (b_br * cos_ph - g_br * sin_ph) * inv_tap;

      g22_ = g_diag - RealT{0.5} * G_;
      b22_ = b_diag - RealT{0.5} * B_;
    }

  } // namespace PhasorDynamics
} // namespace GridKit
