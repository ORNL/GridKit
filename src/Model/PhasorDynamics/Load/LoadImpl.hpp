#include <cmath>
#include <iostream>

#include "Load.hpp"
#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/Load/LoadData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /*
     * @brief Constructor for a pi-model load
     *
     * System sizes:
     * - Number of equations = 0
     * - Number of independent variables = 0
     */
    template <class ScalarT, typename IdxT>
    Load<ScalarT, IdxT>::Load(bus_type* bus)
      : bus_(bus)
    {
      size_ = 0;
    }

    template <class ScalarT, typename IdxT>
    Load<ScalarT, IdxT>::Load(bus_type* bus,
                              real_type R,
                              real_type X)
      : bus_(bus),
        R_(R),
        X_(X)
    {
    }

    template <class ScalarT, typename IdxT>
    Load<ScalarT, IdxT>::Load(bus_type*              bus,
                              const model_data_type& data)
      : bus_(bus)
    {
      if (data.parameters.contains(model_data_type::Parameters::R))
      {
        R_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::R));
      }

      if (data.parameters.contains(model_data_type::Parameters::X))
      {
        X_ = std::get<real_type>(data.parameters.at(model_data_type::Parameters::X));
      }
    }

    template <class ScalarT, typename IdxT>
    Load<ScalarT, IdxT>::Load(bus_type* bus, IdxT component_id)
      : bus_(bus)
    {
      size_         = 0;
      component_id_ = component_id;
    }

    template <class ScalarT, typename IdxT>
    Load<ScalarT, IdxT>::~Load()
    {
      // std::cout << "Destroy Load..." << std::endl;
    }

    /*!
     * @brief allocate method computes sparsity pattern of the Jacobian.
     */
    template <class ScalarT, typename IdxT>
    int Load<ScalarT, IdxT>::allocate()
    {
      // std::cout << "Allocate Load..." << std::endl;
      return 0;
    }

    /**
     * Initialization of the load model
     *
     */
    template <class ScalarT, typename IdxT>
    int Load<ScalarT, IdxT>::initialize()
    {
      return 0;
    }

    /**
     * \brief Identify differential variables.
     */
    template <class ScalarT, typename IdxT>
    int Load<ScalarT, IdxT>::tagDifferentiable()
    {
      return 0;
    }

    /**
     * @brief Residual contribution computed locally
     *
     */
    template <class ScalarT, typename IdxT>
    int Load<ScalarT, IdxT>::evaluateResidualLocally(ScalarT* y, ScalarT* f)
    {
      real_type b = -X_ / (R_ * R_ + X_ * X_);
      real_type g = R_ / (R_ * R_ + X_ * X_);

      f[0] = -g * y[0] + b * y[1];
      f[1] = -b * y[0] - g * y[1];

      return 0;
    }

    /**
     * @brief Residual contribution of the load is pushed to the bus.
     *
     */
    template <class ScalarT, typename IdxT>
    int Load<ScalarT, IdxT>::evaluateResidual()
    {
      std::vector<ScalarT> y(2);
      std::vector<ScalarT> f(2);
      y[0] = Vr();
      y[1] = Vi();
      evaluateResidualLocally(y.data(), f.data());
      Ir() += f[0];
      Ii() += f[1];

      return 0;
    }
  } // namespace PhasorDynamics
} // namespace GridKit
