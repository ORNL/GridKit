/**
 * @file RegcaEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Enzyme sparse Jacobian for the REGCA converter model.
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "RegcaImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Converter
    {
      /**
       * @brief Assemble the sparse component Jacobian with Enzyme.
       *
       * @pre evaluateResidual() has run at the current state.
       */
      template <typename scalar_type, typename index_type>
      int Regca<scalar_type, index_type>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Regca..." << std::endl;
        Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

        if (J_rows_buffer_ == nullptr)
        {
          // Reserve space for the dense blocks. Enzyme keeps only structural
          // nonzeros for each differentiated block.
          auto size        = static_cast<size_t>(size_);
          auto bus_size    = static_cast<size_t>(bus_->size());
          auto signal_size = static_cast<size_t>(ws_.getSize());
          auto buffer_size = 2 * size * size + 2 * size * bus_size + size * signal_size;
          J_rows_buffer_   = new IdxT[buffer_size];
          J_cols_buffer_   = new IdxT[buffer_size];
          J_vals_buffer_   = new RealT[buffer_size];
        }

        using RegcaT = GridKit::PhasorDynamics::Converter::Regca<ScalarT, IdxT>;
        using Fn     = GridKit::Enzyme::Sparse::MemberFunctions;

        nnz_ = 0;

        GridKit::Enzyme::Sparse::DfDy<RegcaT,
                                      Fn::InternalResidualWithSignal>::eval(this,
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

        GridKit::Enzyme::Sparse::DfDyp<RegcaT,
                                       Fn::InternalResidualWithSignal>::eval(this,
                                                                             static_cast<size_t>(f_.getSize()),
                                                                             static_cast<size_t>(y_.getSize()),
                                                                             (this->getResidualIndices()).data(),
                                                                             (this->getVariableIndices()).data(),
                                                                             y_.getData(),
                                                                             yp_.getData(),
                                                                             wb_.getData(),
                                                                             ws_.getData(),
                                                                             alpha_,
                                                                             J_rows_buffer_,
                                                                             J_cols_buffer_,
                                                                             J_vals_buffer_,
                                                                             nnz_);

        GridKit::Enzyme::Sparse::DfDwb<RegcaT,
                                       Fn::InternalResidualWithSignal>::eval(this,
                                                                             static_cast<size_t>(f_.getSize()),
                                                                             static_cast<size_t>(bus_->size()),
                                                                             (this->getResidualIndices()).data(),
                                                                             (bus_->getVariableIndices()).data(),
                                                                             y_.getData(),
                                                                             yp_.getData(),
                                                                             wb_.getData(),
                                                                             ws_.getData(),
                                                                             J_rows_buffer_,
                                                                             J_cols_buffer_,
                                                                             J_vals_buffer_,
                                                                             nnz_);

        GridKit::Enzyme::Sparse::DfDws<RegcaT,
                                       Fn::InternalResidualWithSignal>::eval(this,
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

        GridKit::Enzyme::Sparse::DhDy<RegcaT,
                                      Fn::BusResidual>::eval(this,
                                                             static_cast<size_t>(bus_->size()),
                                                             static_cast<size_t>(y_.getSize()),
                                                             (bus_->getResidualIndices()).data(),
                                                             (this->getVariableIndices()).data(),
                                                             y_.getData(),
                                                             yp_.getData(),
                                                             wb_.getData(),
                                                             J_rows_buffer_,
                                                             J_cols_buffer_,
                                                             J_vals_buffer_,
                                                             nnz_);
        this->constructCoo();

        return 0;
      }

      template class Regca<double, long int>;
      template class Regca<double, size_t>;
    } // namespace Converter
  } // namespace PhasorDynamics
} // namespace GridKit
