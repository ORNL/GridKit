/**
 * @file DependentNortonEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Enzyme sparse Jacobian for the Norton source model.
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "DependentNortonImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Source
    {
      /**
       * @brief Assemble voltage and source-current derivatives with Enzyme.
       *
       * @pre evaluateResidual() has run at the current state.
       */
      template <typename scalar_type, typename index_type>
      int DependentNorton<scalar_type, index_type>::evaluateJacobian()
      {
        if (J_rows_buffer_ == nullptr)
        {
          // Reserve space for the dense blocks. Enzyme keeps only structural
          // nonzeros for each differentiated block.
          auto bus_size    = static_cast<size_t>(bus_->size());
          auto signal_size = static_cast<size_t>(ws_.getSize());
          auto buffer_size = bus_size * (bus_size + signal_size);
          J_rows_buffer_   = new IdxT[buffer_size];
          J_cols_buffer_   = new IdxT[buffer_size];
          J_vals_buffer_   = new RealT[buffer_size];
        }

        using DependentNortonT = GridKit::PhasorDynamics::Source::DependentNorton<ScalarT, IdxT>;
        using Fn               = GridKit::Enzyme::Sparse::MemberFunctions;

        nnz_ = 0;

        GridKit::Enzyme::Sparse::DhDwb<DependentNortonT,
                                       Fn::BusResidualWithSignal>::eval(this,
                                                                        static_cast<size_t>(bus_->size()),
                                                                        static_cast<size_t>(bus_->size()),
                                                                        (bus_->getResidualIndices()).data(),
                                                                        (bus_->getVariableIndices()).data(),
                                                                        y_.getData(),
                                                                        yp_.getData(),
                                                                        wb_.getData(),
                                                                        ws_.getData(),
                                                                        J_rows_buffer_,
                                                                        J_cols_buffer_,
                                                                        J_vals_buffer_,
                                                                        nnz_);

        GridKit::Enzyme::Sparse::DhDws<DependentNortonT,
                                       Fn::BusResidualWithSignal>::eval(this,
                                                                        static_cast<size_t>(bus_->size()),
                                                                        static_cast<size_t>(ws_.getSize()),
                                                                        (bus_->getResidualIndices()).data(),
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
      template class DependentNorton<double, long int>;
      template class DependentNorton<double, size_t>;
    } // namespace Source
  } // namespace PhasorDynamics
} // namespace GridKit
