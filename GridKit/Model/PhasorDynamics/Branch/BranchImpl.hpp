/**
 * @file BranchImpl.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Definition of a phasor dynamics branch model.
 *
 * The model uses Cartesian coordinates.
 *
 */

#include <cmath>
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
     * @param G - total shunt conductance
     * @param B - total shunt susceptance
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

      this->allocateExternalVectors(static_cast<IdxT>(BranchExternalVariables::MAXIMUM));
      f_ext_.resize(4);
      residual_indices_ext_.assign(4, INVALID_INDEX<IdxT>);

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
     * @brief External residual contributions to both terminal buses.
     *
     */
    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) inline int Branch<scalar_type, index_type>::evaluateExternalResidual(
        [[maybe_unused]] const ScalarT* y,
        [[maybe_unused]] const ScalarT* yp,
        const ScalarT*                  y_ext,
        ScalarT*                        f_ext)
    {
      const ScalarT Vr1 = y_ext[0];
      const ScalarT Vi1 = y_ext[1];
      const ScalarT Vr2 = y_ext[2];
      const ScalarT Vi2 = y_ext[3];

      f_ext[0] = (g11_ * Vr1 - b11_ * Vi1) + (g12_ * Vr2 - b12_ * Vi2);
      f_ext[1] = (b11_ * Vr1 + g11_ * Vi1) + (b12_ * Vr2 + g12_ * Vi2);
      f_ext[2] = (g21_ * Vr1 - b21_ * Vi1) + (g22_ * Vr2 - b22_ * Vi2);
      f_ext[3] = (b21_ * Vr1 + g21_ * Vi1) + (b22_ * Vr2 + g22_ * Vi2);

      return 0;
    }

    /**
     * @brief The branch owns no internal residual equations.
     *
     */
    template <typename scalar_type, typename index_type>
    int Branch<scalar_type, index_type>::evaluateInternalResidual()
    {
      return 0;
    }

    /**
     * @brief Gather external variables and index maps.
     *
     */
    template <typename scalar_type, typename index_type>
    void Branch<scalar_type, index_type>::gatherExternalVariables()
    {
      y_ext_[0] = Vr1();
      y_ext_[1] = Vi1();
      y_ext_[2] = Vr2();
      y_ext_[3] = Vi2();
      if (bus1_->size() > 0)
      {
        variable_indices_ext_[0] = bus1_->getVariableIndex(0);
        variable_indices_ext_[1] = bus1_->getVariableIndex(1);
        residual_indices_ext_[0] = bus1_->getResidualIndex(0);
        residual_indices_ext_[1] = bus1_->getResidualIndex(1);
      }
      if (bus2_->size() > 0)
      {
        variable_indices_ext_[2] = bus2_->getVariableIndex(0);
        variable_indices_ext_[3] = bus2_->getVariableIndex(1);
        residual_indices_ext_[2] = bus2_->getResidualIndex(0);
        residual_indices_ext_[3] = bus2_->getResidualIndex(1);
      }
    }

    /**
     * @brief External residual contributions to the two terminal buses.
     *
     */
    template <typename scalar_type, typename index_type>
    int Branch<scalar_type, index_type>::evaluateExternalResidual()
    {
      gatherExternalVariables();

      const auto* y  = y_.getData();
      const auto* yp = yp_.getData();
      evaluateExternalResidual(y, yp, y_ext_.data(), f_ext_.data());

      Ir1() += f_ext_[0];
      Ii1() += f_ext_[1];
      Ir2() += f_ext_[2];
      Ii2() += f_ext_[3];

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

    /**
     * @brief Residual contribution of the branch is computed and pushed to the terminal buses.
     *
     */
    template <typename scalar_type, typename index_type>
    int Branch<scalar_type, index_type>::evaluateResidual()
    {
      evaluateInternalResidual();
      return evaluateExternalResidual();
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

      readRealParameter(data, Parameter::R, R_);
      readRealParameter(data, Parameter::X, X_);
      readRealParameter(data, Parameter::G, G_);
      readRealParameter(data, Parameter::B, B_);
      readRealParameter(data, Parameter::tap, tap_);
      readRealParameter(data, Parameter::phase, phase_);

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
    bool Branch<scalar_type, index_type>::readRealParameter(const ModelDataT&               data,
                                                            typename ModelDataT::Parameters parameter,
                                                            RealT&                          target)
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

      const RealT g_diag = -(g_br + RealT{0.5} * G_);
      const RealT b_diag = -(b_br + RealT{0.5} * B_);

      g11_ = g_diag * inv_tap * inv_tap;
      b11_ = b_diag * inv_tap * inv_tap;

      g12_ = (g_br * cos_ph - b_br * sin_ph) * inv_tap;
      b12_ = (b_br * cos_ph + g_br * sin_ph) * inv_tap;

      g21_ = (g_br * cos_ph + b_br * sin_ph) * inv_tap;
      b21_ = (b_br * cos_ph - g_br * sin_ph) * inv_tap;

      g22_ = g_diag;
      b22_ = b_diag;
    }

  } // namespace PhasorDynamics
} // namespace GridKit
