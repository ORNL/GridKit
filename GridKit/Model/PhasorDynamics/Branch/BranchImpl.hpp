/**
 * @file BranchImpl.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Definition of a phasor dynamics branch model.
 *
 * The model uses Cartesian coordinates.
 *
 */

#include <cmath>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/Branch/Branch.hpp>
#include <GridKit/Model/PhasorDynamics/Branch/BranchData.hpp>
#include <GridKit/Model/PhasorDynamics/Bus/Bus.hpp>
#include <GridKit/Model/PhasorDynamics/ConnectedElementImpl.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief Constructor for a pi-model branch
     *
     * Model size:
     * - Number of equations = 0
     * - Number of internal variables = 0
     */
    template <typename ScalarP, typename IdxP>
    Branch<ScalarP, IdxP>::Branch(BusT* bus1, BusT* bus2)
      : bus1_(bus1),
        bus2_(bus2),
        R_(0.0),
        X_(0.01),
        G_(0.0),
        B_(0.0),
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
     * @param G - line shunt conductance
     * @param B - line shunt charging
     */
    template <typename ScalarP, typename IdxP>
    Branch<ScalarP, IdxP>::Branch(BusT* bus1,
                                  BusT* bus2,
                                  RealT R,
                                  RealT X,
                                  RealT G,
                                  RealT B)
      : bus1_(bus1),
        bus2_(bus2),
        R_(R),
        X_(X),
        G_(G),
        B_(B),
        bus1_id_(0),
        bus2_id_(0)
    {
      setDerivedParams();
    }

    template <typename ScalarP, typename IdxP>
    Branch<ScalarP, IdxP>::Branch(BusT* bus1, BusT* bus2, const ModelDataT& data)
      : ConnectedElement<Branch>(data),
        bus1_(bus1),
        bus2_(bus2)
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
    template <typename ScalarP, typename IdxP>
    Branch<ScalarP, IdxP>::~Branch()
    {
      // std::cout << "Destroy Branch..." << std::endl;
    }

    /**
     * @brief Set the component ID
     */
    template <typename ScalarP, typename IdxP>
    int Branch<ScalarP, IdxP>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /*!
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <typename ScalarP, typename IdxP>
    int Branch<ScalarP, IdxP>::allocate()
    {
      // std::cout << "Allocate Branch..." << std::endl;

      wb_.resize(2);
      h_.resize(2);

      return 0;
    }

    /**
     * Initialization of the branch model
     *
     */
    template <typename ScalarP, typename IdxP>
    int Branch<ScalarP, IdxP>::initialize()
    {
      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <typename ScalarP, typename IdxP>
    int Branch<ScalarP, IdxP>::tagDifferentiable()
    {
      return 0;
    }

    /**
     * @brief Bus 1 residual contribution from bus 1 variables
     *
     */
    template <typename ScalarP, typename IdxP>
    __attribute__((always_inline)) inline int Branch<ScalarP, IdxP>::evaluateBusResidual11(
        [[maybe_unused]] ScalarT* y,
        [[maybe_unused]] ScalarT* yp,
        ScalarT*                  wb,
        ScalarT*                  h)
    {
      ScalarT Vr1 = wb[0];
      ScalarT Vi1 = wb[1];
      ScalarT Ir1 = -(g_ + 0.5 * G_) * Vr1 + (b_ + 0.5 * B_) * Vi1;
      ScalarT Ii1 = -(b_ + 0.5 * B_) * Vr1 - (g_ + 0.5 * G_) * Vi1;
      h[0]        = Ir1;
      h[1]        = Ii1;

      return 0;
    }

    /**
     * @brief Bus 1 residual contribution from bus 2 variables
     *
     */
    template <typename ScalarP, typename IdxP>
    __attribute__((always_inline)) inline int Branch<ScalarP, IdxP>::evaluateBusResidual12(
        [[maybe_unused]] ScalarT* y,
        [[maybe_unused]] ScalarT* yp,
        ScalarT*                  wb,
        ScalarT*                  h)
    {
      ScalarT Vr2 = wb[0];
      ScalarT Vi2 = wb[1];
      ScalarT Ir1 = g_ * Vr2 - b_ * Vi2;
      ScalarT Ii1 = b_ * Vr2 + g_ * Vi2;
      h[0]        = Ir1;
      h[1]        = Ii1;

      return 0;
    }

    /**
     * @brief Bus 2 residual contribution from bus 1 variables
     *
     */
    template <typename ScalarP, typename IdxP>
    __attribute__((always_inline)) int Branch<ScalarP, IdxP>::evaluateBusResidual21(
        [[maybe_unused]] ScalarT* y,
        [[maybe_unused]] ScalarT* yp,
        ScalarT*                  wb,
        ScalarT*                  h)
    {
      ScalarT Vr1 = wb[0];
      ScalarT Vi1 = wb[1];
      ScalarT Ir2 = g_ * Vr1 - b_ * Vi1;
      ScalarT Ii2 = b_ * Vr1 + g_ * Vi1;
      h[0]        = Ir2;
      h[1]        = Ii2;

      return 0;
    }

    /**
     * @brief Bus 2 residual contribution from bus 2 variables
     *
     */
    template <typename ScalarP, typename IdxP>
    __attribute__((always_inline)) int Branch<ScalarP, IdxP>::evaluateBusResidual22(
        [[maybe_unused]] ScalarT* y,
        [[maybe_unused]] ScalarT* yp,
        ScalarT*                  wb,
        ScalarT*                  h)
    {
      ScalarT Vr2 = wb[0];
      ScalarT Vi2 = wb[1];
      ScalarT Ir2 = -(g_ + 0.5 * G_) * Vr2 + (b_ + 0.5 * B_) * Vi2;
      ScalarT Ii2 = -(b_ + 0.5 * B_) * Vr2 - (g_ + 0.5 * G_) * Vi2;
      h[0]        = Ir2;
      h[1]        = Ii2;

      return 0;
    }

    /**
     * @brief Residual contribution of the branch is computed and pushed to the terminal buses.
     *
     */
    template <typename ScalarP, typename IdxP>
    int Branch<ScalarP, IdxP>::evaluateResidual()
    {
      wb_[0] = Vr1();
      wb_[1] = Vi1();
      evaluateBusResidual11(y_.data(), yp_.data(), wb_.data(), h_.data());
      Ir1() += h_[0];
      Ii1() += h_[1];
      evaluateBusResidual21(y_.data(), yp_.data(), wb_.data(), h_.data());
      Ir2() += h_[0];
      Ii2() += h_[1];

      wb_[0] = Vr2();
      wb_[1] = Vi2();
      evaluateBusResidual12(y_.data(), yp_.data(), wb_.data(), h_.data());
      Ir1() += h_[0];
      Ii1() += h_[1];
      evaluateBusResidual22(y_.data(), yp_.data(), wb_.data(), h_.data());
      Ir2() += h_[0];
      Ii2() += h_[1];

      return 0;
    }

    template <typename ScalarP, typename IdxP>
    void Branch<ScalarP, IdxP>::initializeParameters(const ModelDataT& data)
    {
      if (data.parameters.contains(ModelDataT::Parameters::R))
      {
        R_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::R));
      }

      if (data.parameters.contains(ModelDataT::Parameters::X))
      {
        X_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::X));
      }

      if (data.parameters.contains(ModelDataT::Parameters::G))
      {
        G_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::G));
      }

      if (data.parameters.contains(ModelDataT::Parameters::B))
      {
        B_ = std::get<RealT>(data.parameters.at(ModelDataT::Parameters::B));
      }

      if (data.ports.contains(ModelDataT::Ports::bus1))
      {
        bus1_id_ = data.ports.at(ModelDataT::Ports::bus1);
      }

      if (data.ports.contains(ModelDataT::Ports::bus2))
      {
        bus2_id_ = data.ports.at(ModelDataT::Ports::bus2);
      }
    }

    template <typename ScalarP, typename IdxP>
    const Model::VariableMonitorBase* Branch<ScalarP, IdxP>::getMonitor() const
    {
      return monitor_.get();
    }

    template <typename ScalarP, typename IdxP>
    void Branch<ScalarP, IdxP>::initializeMonitor()
    {
      using Variable = typename ModelDataT::MonitorableVariables;
      monitor_->set(Variable::ir1, [this]
                    { return Ir1(); });
      monitor_->set(Variable::ii1, [this]
                    { return Ii1(); });
      monitor_->set(Variable::im1, [this]
                    { return std::sqrt(Ir1() * Ir1() + Ii1() * Ii1()); });
      monitor_->set(Variable::p1, [this]
                    { return Vr1() * Ir1() + Vi1() * Ii1(); });
      monitor_->set(Variable::q1, [this]
                    { return Vi1() * Ir1() - Vr1() * Ii1(); });
      monitor_->set(Variable::ir2, [this]
                    { return Ir2(); });
      monitor_->set(Variable::ii2, [this]
                    { return Ii2(); });
      monitor_->set(Variable::im2, [this]
                    { return std::sqrt(Ir2() * Ir2() + Ii2() * Ii2()); });
      monitor_->set(Variable::p2, [this]
                    { return Vr2() * Ir2() + Vi2() * Ii2(); });
      monitor_->set(Variable::q2, [this]
                    { return Vi2() * Ir2() - Vr2() * Ii2(); });
    }

    /**
     * @brief Derived parameters
     *
     */
    template <typename ScalarP, typename IdxP>
    void Branch<ScalarP, IdxP>::setDerivedParams()
    {
      b_ = -X_ / (R_ * R_ + X_ * X_);
      g_ = R_ / (R_ * R_ + X_ * X_);
    }

  } // namespace PhasorDynamics
} // namespace GridKit
