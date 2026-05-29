/**
 * @file HygovEnzyme.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Enzyme sparse Jacobian for the HYGOV governor model.
 */

#include <GridKit/AutomaticDifferentiation/Enzyme/SparseJacobians.hpp>

#include "HygovImpl.hpp"

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Governor
    {
      namespace
      {
        template <typename RealT>
        RealT deadband1Derivative(RealT x, RealT lower, RealT upper)
        {
          static constexpr RealT MU = static_cast<RealT>(240.0);

          const RealT lower_sig = Math::sigmoid(lower - x);
          const RealT upper_sig = Math::sigmoid(x - upper);
          return lower_sig + upper_sig
                 + x * MU
                       * (upper_sig * (ONE<RealT> - upper_sig)
                          - lower_sig * (ONE<RealT> - lower_sig));
        }
      } // namespace

      template <class ScalarT, typename IdxT>
      int Hygov<ScalarT, IdxT>::evaluateJacobian()
      {
        Log::misc() << "Evaluate Jacobian for Hygov..." << std::endl;
        Log::misc() << "Jacobian evaluation is experimental!" << std::endl;

        J_.zeroMatrix();
        if (J_rows_buffer_ == nullptr)
        {
          J_rows_buffer_ = new IdxT[static_cast<size_t>(size_) * static_cast<size_t>(size_)];
          J_cols_buffer_ = new IdxT[static_cast<size_t>(size_) * static_cast<size_t>(size_)];
          J_vals_buffer_ = new RealT[static_cast<size_t>(size_) * static_cast<size_t>(size_)];
        }

        using ModelT = GridKit::PhasorDynamics::Governor::Hygov<ScalarT, IdxT>;
        using Fn     = GridKit::Enzyme::Sparse::MemberFunctions;

        GridKit::Enzyme::Sparse::DfDy<ModelT,
                                      Fn::InternalResidualWithSignal,
                                      ScalarT,
                                      IdxT>::eval(this,
                                                  f_.size(),
                                                  y_.size(),
                                                  (this->getResidualIndices()).data(),
                                                  (this->getVariableIndices()).data(),
                                                  y_.data(),
                                                  yp_.data(),
                                                  wb_.data(),
                                                  ws_.data(),
                                                  alpha_,
                                                  J_rows_buffer_,
                                                  J_cols_buffer_,
                                                  J_vals_buffer_,
                                                  J_);

        const auto OMEGADB = static_cast<size_t>(HygovInternalVariables::OMEGADB);
        const auto EF      = static_cast<size_t>(HygovInternalVariables::EF);
        const auto PMECH   = static_cast<size_t>(HygovInternalVariables::PMECH);
        const auto G       = static_cast<size_t>(HygovInternalVariables::G);

        const auto OMEGA = static_cast<size_t>(HygovExternalVariables::OMEGA);
        const auto PREF  = static_cast<size_t>(HygovExternalVariables::PREF);
        const auto PAUX  = static_cast<size_t>(HygovExternalVariables::PAUX);

        IdxT nnz                   = 0;
        auto storeSignalDerivative = [&](size_t row, size_t col, RealT value)
        {
          if (ws_indices_[col] == INVALID_INDEX<IdxT>)
          {
            return;
          }

          J_rows_buffer_[nnz] = residual_indices_.at(row);
          J_cols_buffer_[nnz] = ws_indices_[col];
          J_vals_buffer_[nnz] = value;
          ++nnz;
        };

        storeSignalDerivative(OMEGADB,
                              OMEGA,
                              deadband1Derivative(static_cast<RealT>(ws_[OMEGA]),
                                                  -db1_,
                                                  db1_));
        storeSignalDerivative(EF, PREF, ONE<RealT>);
        storeSignalDerivative(EF, PAUX, ONE<RealT>);
        storeSignalDerivative(PMECH,
                              OMEGA,
                              -Dturb_ * static_cast<RealT>(y_[G]) * Trate_ / pmech_mva_base_);

        J_.setValues(1.0, J_rows_buffer_, J_cols_buffer_, J_vals_buffer_, nnz);
        return 0;
      }

      template class Hygov<double, long int>;
      template class Hygov<double, size_t>;
    } // namespace Governor
  } // namespace PhasorDynamics
} // namespace GridKit
