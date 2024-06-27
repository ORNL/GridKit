/* Created by Paul Moon 6/26/2024 */
#include <iostream>
#include <cmath>
#include "SlackBus.hpp"

namespace ModelLib {

/*!
 * @brief Constructor for a slack bus
 *
 * Arguments passed to ModelEvaluatorImpl:
 * - Number of equations = 0 (size_)
 * - Number of variables = 0 (size_)
 * - Number of quadratures = 0
 * - Number of optimization parameters = 0
 */
template <class ScalarT, typename IdxT>
SlackBus<ScalarT, IdxT>::SlackBus()
    : PhasorBus<ScalarT, IdxT>(0), Vr_(0.0), Vi_(0.0), P_(0.0), Q_(0.0), PB_(0.0), QB_(0.0)
{
    //std::cout << "Create SlackBus..." << std::endl;
    //std::cout << "Number of equations is " << size_ << std::endl;

    size_ = 0;
}

/*!
 * @brief SlackBus constructor.
 *
 * Arguments passed to ModelEvaluatorImpl:
 * - Number of equations = 0 (size_)
 * - Number of variables = 0 (size_)
 * - Number of quadratures = 0
 * - Number of optimization parameters = 0
 */
template <class ScalarT, typename IdxT>
SlackBus<ScalarT, IdxT>::SlackBus(ScalarT Vr, ScalarT Vi)
    : PhasorBus<ScalarT, IdxT>(0), Vr_(Vr), Vi_(Vi), P_(0.0), Q_(0.0), PB_(0.0), QB_(0.0)
{
    //std::cout << "Create SlackBus..." << std::endl;
    //std::cout << "Number of equations is " << size_ << std::endl;
    P() = 0.0;
    Q() = 0.0;
    size_ = 0;
}

template <class ScalarT, typename IdxT>
SlackBus<ScalarT, IdxT>::SlackBus(BusData& data)
    : PhasorBus<ScalarT, IdxT>(data.bus_i), Vr_(data.Vm*cos(data.Va)), Vi_(data.Vm*sin(data.Va))
{
    //std::cout << "Create SlackBus..." << std::endl;
    //std::cout << "Number of equations is " << size_ << std::endl;
    P() = 0.0;
    Q() = 0.0;
    size_ = 0;
}

template <class ScalarT, typename IdxT>
SlackBus<ScalarT, IdxT>::~SlackBus()
{
}

template <class ScalarT, typename IdxT>
int SlackBus<ScalarT, IdxT>::evaluateResidual()
{
    // std::cout << "Evaluating residual of a slack bus ...\n";
    P() = 0.0;
    Q() = 0.0;
    return 0;
}

template <class ScalarT, typename IdxT>
int SlackBus<ScalarT, IdxT>::evaluateAdjointResidual()
{
    PB() = 0.0;
    QB() = 0.0;
    return 0;
}


// Available template instantiations
template class SlackBus<double, long int>;
template class SlackBus<double, size_t>;


} // namespace ModelLib

