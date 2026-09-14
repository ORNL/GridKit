/**
 * @file IeeestEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Enzyme sparse Jacobian for the IEEEST stabilizer model.
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "IeeestImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Stabilizer
    {
      /**
       * @brief Evaluate the sparse residual Jacobian.
       *
       * @tparam scalar_type Scalar data type.
       * @tparam index_type Index data type.
       * @tparam order Notch-denominator degree.
       * @return int - error code, 0 = success
       */
      template <typename scalar_type, typename index_type, size_t order>
      int Ieeest<scalar_type, index_type, order>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Ieeest..." << std::endl;
        Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

        if (J_rows_buffer_ == nullptr)
        {
          // Reserve space for the dense blocks.
          // The size of the buffer is the sum of maximum capacities of the blocks.
          // Enzyme will compute the appropriate nnz from sparsification.
          auto size        = static_cast<size_t>(size_);
          auto signal_size = static_cast<size_t>(ws_.getSize());
          auto buffer_size = 2 * size * size + size * signal_size;
          J_rows_buffer_   = new IdxT[buffer_size];
          J_cols_buffer_   = new IdxT[buffer_size];
          J_vals_buffer_   = new RealT[buffer_size];
        }

        using ModelT = Ieeest<ScalarT, IdxT, order>;
        using Fn     = GridKit::Enzyme::Sparse::MemberFunctions;

        nnz_ = 0;

        GridKit::Enzyme::Sparse::DfDy<ModelT, Fn::InternalResidualWithSignal>::eval(this,
                                                                                    static_cast<size_t>(f_.getSize()),
                                                                                    static_cast<size_t>(y_.getSize()),
                                                                                    (this->getResidualIndices()).data(),
                                                                                    (this->getVariableIndices()).data(),
                                                                                    y_.getData(),
                                                                                    yp_.getData(),
                                                                                    nullptr,
                                                                                    ws_.getData(),
                                                                                    J_rows_buffer_,
                                                                                    J_cols_buffer_,
                                                                                    J_vals_buffer_,
                                                                                    nnz_);

        GridKit::Enzyme::Sparse::DfDyp<ModelT, Fn::InternalResidualWithSignal>::eval(this,
                                                                                     static_cast<size_t>(f_.getSize()),
                                                                                     static_cast<size_t>(y_.getSize()),
                                                                                     (this->getResidualIndices()).data(),
                                                                                     (this->getVariableIndices()).data(),
                                                                                     y_.getData(),
                                                                                     yp_.getData(),
                                                                                     nullptr,
                                                                                     ws_.getData(),
                                                                                     alpha_,
                                                                                     J_rows_buffer_,
                                                                                     J_cols_buffer_,
                                                                                     J_vals_buffer_,
                                                                                     nnz_);

        GridKit::Enzyme::Sparse::DfDws<ModelT, Fn::InternalResidualWithSignal>::eval(this,
                                                                                     static_cast<size_t>(f_.getSize()),
                                                                                     static_cast<size_t>(ws_.getSize()),
                                                                                     (this->getResidualIndices()).data(),
                                                                                     ws_indices_.data(),
                                                                                     y_.getData(),
                                                                                     yp_.getData(),
                                                                                     nullptr,
                                                                                     ws_.getData(),
                                                                                     J_rows_buffer_,
                                                                                     J_cols_buffer_,
                                                                                     J_vals_buffer_,
                                                                                     nnz_);

        this->constructCoo();

        return 0;
      }

      // Available template instantiations
      template class Ieeest<double, long int, 0>;
      template class Ieeest<double, long int, 1>;
      template class Ieeest<double, long int, 2>;
      template class Ieeest<double, long int, 3>;
      template class Ieeest<double, long int, 4>;
      template class Ieeest<double, size_t, 0>;
      template class Ieeest<double, size_t, 1>;
      template class Ieeest<double, size_t, 2>;
      template class Ieeest<double, size_t, 3>;
      template class Ieeest<double, size_t, 4>;

    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
