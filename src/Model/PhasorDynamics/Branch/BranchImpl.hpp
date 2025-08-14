/**
 * @file BranchImpl.hpp
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Definition of a phasor dynamics branch model.
 *
 * The model uses Cartesian coordinates.
 *
 */

#include "Branch.hpp"

#include <cmath>
#include <iostream>

#include <Model/PhasorDynamics/Branch/BranchData.hpp>
#include <Model/PhasorDynamics/Bus/Bus.hpp>

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
                                  real_type R,
                                  real_type X,
                                  real_type G,
                                  real_type B)
      : bus1_(bus1),
        bus2_(bus2),
        R_(R),
        X_(X),
        G_(G),
        B_(B),
        bus1_id_(0),
        bus2_id_(0)
    {
    }

    template <class ScalarT, typename IdxT>
    Branch<ScalarT, IdxT>::Branch(bus_type* bus1, bus_type* bus2, const model_data_type& data)
      : bus1_(bus1),
        bus2_(bus2)
    {
      if (data.parameters.contains(model_data_type::Parameters::R))
      {
        R_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::R));
      }

      if (data.parameters.contains(model_data_type::Parameters::X))
      {
        X_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::X));
      }

      if (data.parameters.contains(model_data_type::Parameters::G))
      {
        G_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::G));
      }

      if (data.parameters.contains(model_data_type::Parameters::B))
      {
        B_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::B));
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

    /*!
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::allocate()
    {
      // std::cout << "Allocate Branch..." << std::endl;
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
     * \brief Residual contribution of the branch is pushed to the
     * two terminal buses.
     *
     */
    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::evaluateResidual()
    {
      // std::cout << "Evaluating branch residual ...\n";
      real_type b = -X_ / (R_ * R_ + X_ * X_);
      real_type g = R_ / (R_ * R_ + X_ * X_);

      Ir1() += -(g + 0.5 * G_) * Vr1() + (b + 0.5 * B_) * Vi1() + g * Vr2() - b * Vi2();
      Ii1() += -(b + 0.5 * B_) * Vr1() - (g + 0.5 * G_) * Vi1() + b * Vr2() + g * Vi2();
      Ir2() += g * Vr1() - b * Vi1() - (g + 0.5 * G_) * Vr2() + (b + 0.5 * B_) * Vi2();
      Ii2() += b * Vr1() + g * Vi1() - (b + 0.5 * B_) * Vr2() - (g + 0.5 * G_) * Vi2();

      return 0;
    }

    /**
     * @brief Jacobian evaluation not implemented yet
     *
     * @tparam ScalarT - scalar data type
     * @tparam IdxT    - matrix index data type
     * @return int - error code, 0 = success
     */
    template <class ScalarT, typename IdxT>
    int Branch<ScalarT, IdxT>::evaluateJacobian()
    {
      std::cout << "Evaluate Jacobian for Branch..." << std::endl;
      std::cout << "Jacobian evaluation not implemented!" << std::endl;
      return 0;
    }

  } // namespace PhasorDynamics
} // namespace GridKit
