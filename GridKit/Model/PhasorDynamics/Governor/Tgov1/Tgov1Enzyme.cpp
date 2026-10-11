/**
 * @file Tgov1Enzyme.cpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/Jacobian.hpp>

#include "Tgov1Impl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      /**
       * @brief Jacobian evaluation for the Enzyme build.
       *
       * @return int - error code, 0 = success
       */
      template <typename scalar_type, typename index_type>
      int Tgov1<scalar_type, index_type>::evaluateJacobian()
      {
        gatherExternalVariables();

        const auto n     = static_cast<size_t>(size_);
        const auto n_ext = static_cast<size_t>(y_ext_.getSize());

        if (J_rows_buffer_ == nullptr)
        {
          // Each (row, column) of [f | f_ext] x [y | y_ext] is stored at most once.
          const auto buffer_size = (n + n_ext) * (n + n_ext);
          J_rows_buffer_         = new IdxT[buffer_size];
          J_cols_buffer_         = new IdxT[buffer_size];
          J_vals_buffer_         = new RealT[buffer_size];
        }

        nnz_ = 0;

        GridKit::Enzyme::Sparse::Jacobian<Tgov1>::eval(this,
                                                       n,
                                                       n_ext,
                                                       y_.getData(),
                                                       yp_.getData(),
                                                       y_ext_.getData(),
                                                       alpha_,
                                                       residual_indices_.data(),
                                                       variable_indices_.data(),
                                                       residual_indices_ext_.data(),
                                                       variable_indices_ext_.data(),
                                                       J_rows_buffer_,
                                                       J_cols_buffer_,
                                                       J_vals_buffer_,
                                                       nnz_);

        this->constructCoo();

        return 0;
      }

      // Available template instantiations
      template class Tgov1<double, long int>;
      template class Tgov1<double, size_t>;

    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
