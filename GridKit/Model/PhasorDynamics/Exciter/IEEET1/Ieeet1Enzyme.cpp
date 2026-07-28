/**
 * @file Ieeet1Enzyme.cpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "Ieeet1Impl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Exciter
    {
      /**
       * @brief Jacobian evaluation not implemented yet
       *
       * @return int - error code, 0 = success
       */
      template <class scalar_type, typename index_type>
      int Ieeet1<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Ieeet1..." << std::endl;
        Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

        if (J_rows_buffer_ == nullptr)
        {
          // Reserve space for the dense blocks.
          // The size of the buffer is the sum of maximum capacities of the blocks.
          // Enyme will compute the appropriate nnz from sparsification.
          auto size        = static_cast<size_t>(size_);
          auto bus_size    = static_cast<size_t>(bus_->size());
          auto signal_size = static_cast<size_t>(ws_.getSize());
          auto buffer_size = 2 * size * size + size * bus_size + size * signal_size;
          J_rows_buffer_   = new IdxT[buffer_size];
          J_cols_buffer_   = new IdxT[buffer_size];
          J_vals_buffer_   = new RealT[buffer_size];
        }

        nnz_ = 0;

        GridKit::Enzyme::Sparse::DfDy<GridKit::PhasorDynamics::Exciter::Ieeet1<ScalarT, IdxT>,
                                      GridKit::Enzyme::Sparse::MemberFunctions::InternalResidualWithSignal>::eval(this,
                                                                                                                  static_cast<size_t>(f_.getSize()),
                                                                                                                  static_cast<size_t>(y_.getSize()),
                                                                                                                  (this->getResidualIndices()).data(),
                                                                                                                  (this->getVariableIndices()).data(),
                                                                                                                  y_.getData(),
                                                                                                                  yp_.getData(),
                                                                                                                  wb_.getData(),
                                                                                                                  ws_.getData(),
                                                                                                                  J_rows_buffer_,
                                                                                                                  J_cols_buffer_,
                                                                                                                  J_vals_buffer_,
                                                                                                                  nnz_);

        GridKit::Enzyme::Sparse::DfDyp<GridKit::PhasorDynamics::Exciter::Ieeet1<ScalarT, IdxT>,
                                       GridKit::Enzyme::Sparse::MemberFunctions::InternalResidualWithSignal>::eval(this,
                                                                                                                   static_cast<size_t>(f_.getSize()),
                                                                                                                   static_cast<size_t>(y_.getSize()),
                                                                                                                   (this->getResidualIndices()).data(),
                                                                                                                   (this->getVariableIndices()).data(),
                                                                                                                   y_.getData(),
                                                                                                                   yp_.getData(),
                                                                                                                   wb_.getData(),
                                                                                                                   ws_.getData(),
                                                                                                                   alpha(),
                                                                                                                   J_rows_buffer_,
                                                                                                                   J_cols_buffer_,
                                                                                                                   J_vals_buffer_,
                                                                                                                   nnz_);

        GridKit::Enzyme::Sparse::DfDwb<GridKit::PhasorDynamics::Exciter::Ieeet1<ScalarT, IdxT>,
                                       GridKit::Enzyme::Sparse::MemberFunctions::InternalResidualWithSignal>::eval(this,
                                                                                                                   static_cast<size_t>(f_.getSize()),
                                                                                                                   static_cast<size_t>(bus_->size()),
                                                                                                                   (this->getResidualIndices()).data(),
                                                                                                                   (bus_->getVariableIndices()).data(),
                                                                                                                   y_.getData(),
                                                                                                                   yp_.getData(),
                                                                                                                   bus_->y().getData(),
                                                                                                                   ws_.getData(),
                                                                                                                   J_rows_buffer_,
                                                                                                                   J_cols_buffer_,
                                                                                                                   J_vals_buffer_,
                                                                                                                   nnz_);

        GridKit::Enzyme::Sparse::DfDws<GridKit::PhasorDynamics::Exciter::Ieeet1<ScalarT, IdxT>,
                                       GridKit::Enzyme::Sparse::MemberFunctions::InternalResidualWithSignal>::eval(this,
                                                                                                                   static_cast<size_t>(f_.getSize()),
                                                                                                                   static_cast<size_t>(ws_.getSize()),
                                                                                                                   (this->getResidualIndices()).data(),
                                                                                                                   ws_indices_.data(),
                                                                                                                   y_.getData(),
                                                                                                                   yp_.getData(),
                                                                                                                   wb_.getData(),
                                                                                                                   ws_.getData(),
                                                                                                                   J_rows_buffer_,
                                                                                                                   J_cols_buffer_,
                                                                                                                   J_vals_buffer_,
                                                                                                                   nnz_);

        this->constructCoo();

        return 0;
      }

      // Available template instantiations
      template class Ieeet1<double, long int>;
      template class Ieeet1<double, size_t>;

    } // namespace Exciter
  } // namespace PhasorDynamics
} // namespace GridKit
