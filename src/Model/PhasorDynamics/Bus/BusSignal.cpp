
#include "BusSignal.hpp"

#include <cmath>
#include <iostream>

#include <Model/PhasorDynamics/Bus/BusData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {


    template <class ScalarT, typename IdxT>
    BusSignal<ScalarT, IdxT>::BusSignal(const DataT& data)
        : BusBase<ScalarT, IdxT>(data.bus_id)
    {
      size_ = 1;
    }


    template <class ScalarT, typename IdxT>
    int BusSignal<ScalarT, IdxT>::allocate()
    {
      // Temporary while we use std::vector in the code
      size_t size = static_cast<size_t>(size_);

      // Resize component model data
      f_.resize(size);
      y_.resize(size);
      yp_.resize(size);
      tag_.resize(size);

      fB_.resize(size);
      yB_.resize(size);
      ypB_.resize(size);

      return 0;
    }

    template <class ScalarT, typename IdxT>
    int BusSignal<ScalarT, IdxT>::tagDifferentiable()
    {
      tag_[0] = false;
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int BusSignal<ScalarT, IdxT>::initialize()
    {
      /// Possible since not being used it causes issue?
      y_[0]  = 0.0;
      yp_[0] = 0.0;

      return 0;
    }

    template <class ScalarT, typename IdxT>
    int BusSignal<ScalarT, IdxT>::evaluateResidual()
    {
      // One side of equality of signal. Reference component
      // must add to this residual after this function has been called.
      f_[0] = - y_[0]; 

      return 0;
    }

    template <class ScalarT, typename IdxT>
    int BusSignal<ScalarT, IdxT>::evaluateJacobian()
    {
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int BusSignal<ScalarT, IdxT>::initializeAdjoint()
    {
      // std::cout << "Initialize BusSignal..." << std::endl;
      yB_[0]  = 0.0;
      ypB_[0] = 0.0;

      return 0;
    }

 
    template <class ScalarT, typename IdxT>
    int BusSignal<ScalarT, IdxT>::evaluateAdjointResidual()
    {
      fB_[0] = 0.0;

      return 0;
    }

    template <class ScalarT, typename IdxT>
    int BusSignal<ScalarT, IdxT>::evaluateIntegrand()
    {
      return 0;
    }

  
    template <class ScalarT, typename IdxT>
    int BusSignal<ScalarT, IdxT>::evaluateAdjointIntegrand()
    {
      return 0;
    }

    // Available template instantiations
    template class BusSignal<double, long int>;
    template class BusSignal<double, size_t>;
    template class BusSignal<DependencyTracking::Variable, long int>;
    template class BusSignal<DependencyTracking::Variable, size_t>;

  } // namespace PhasorDynamics
} // namespace GridKit