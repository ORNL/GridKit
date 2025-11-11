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
    template <class ScalarT, typename IdxT>
    Branch<ScalarT, IdxT>::Branch(bus_type* bus1, bus_type* bus2)
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
     * @tparam ScalarT - scalar type
     * @tparam IdxT    - matrix/vector index type
     * @param bus1 - pointer to bus-1
     * @param bus2 - pointer to bus-2
     * @param R - line series resistance
     * @param X - line series reactance
     * @param G - line shunt conductance
     * @param B - line shunt charging
     */
    template <class ScalarT, typename IdxT>
    Branch<ScalarT, IdxT>::Branch(bus_type* bus1,
                                  bus_type* bus2,
                                  RealT     R,
                                  RealT     X,
                                  RealT     G,
                                  RealT     B)
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

    template <class ScalarT, typename IdxT>
    Branch<ScalarT, IdxT>::Branch(bus_type* bus1, bus_type* bus2, const model_data_type& data)
      : bus1_(bus1),
        bus2_(bus2)
    {
      if (data.parameters.contains(model_data_type::Parameters::R))
      {
        R_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::R));
      }

      if (data.parameters.contains(model_data_type::Parameters::X))
      {
        X_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::X));
      }

      if (data.parameters.contains(model_data_type::Parameters::G))
      {
        G_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::G));
      }

      if (data.parameters.contains(model_data_type::Parameters::B))
      {
        B_ = std::get<RealT>(data.parameters.at(model_data_type::Parameters::B));
      }

      if (data.ports.contains(model_data_type::Ports::bus1))
      {
        bus1_id_ = data.ports.at(model_data_type::Ports::bus1);
      }

      if (data.ports.contains(model_data_type::Ports::bus2))
      {
        bus2_id_ = data.ports.at(model_data_type::Ports::bus2);
      }

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
      // std::cout << "Destroy Branch..." << std::endl;
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
      // std::cout << "Allocate Branch..." << std::endl;

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
    template <class ScalarT, typename IdxT>
    __attribute__((always_inline)) inline int Branch<ScalarT, IdxT>::evaluateBusResidual12(
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
    template <class ScalarT, typename IdxT>
    __attribute__((always_inline)) int Branch<ScalarT, IdxT>::evaluateBusResidual21(
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
    template <class ScalarT, typename IdxT>
    __attribute__((always_inline)) int Branch<ScalarT, IdxT>::evaluateBusResidual22(
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
    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::evaluateResidual()
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

    /**
     * @brief Derived parameters
     *
     */
    template <class ScalarT, typename IdxT>
    void Branch<ScalarT, IdxT>::setDerivedParams()
    {
      b_ = -X_ / (R_ * R_ + X_ * X_);
      g_ = R_ / (R_ * R_ + X_ * X_);
    }

  } // namespace PhasorDynamics
} // namespace GridKit
