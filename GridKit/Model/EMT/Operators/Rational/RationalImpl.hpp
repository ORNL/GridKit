/**
 * @file RationalImpl.hpp
 * @brief Implementation of the EMT Rational model.
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <GridKit/Model/EMT/Operators/Rational/Rational.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    typename Rational<scalar_type, index_type>::IdxT Rational<scalar_type, index_type>::rows() const
    {
      return static_cast<IdxT>(rows_);
    }

    template <typename scalar_type, typename index_type>
    typename Rational<scalar_type, index_type>::IdxT Rational<scalar_type, index_type>::cols() const
    {
      return static_cast<IdxT>(cols_);
    }

    template <typename scalar_type, typename index_type>
    void Rational<scalar_type, index_type>::attachInput(const std::vector<SignalT*>& input)
    {
      if (input.size() != cols_)
      {
        throw std::invalid_argument("Rational: input dimension mismatch");
      }
      input_ = input;
    }

    template <typename scalar_type, typename index_type>
    void Rational<scalar_type, index_type>::attachOutput(const std::vector<SignalT*>& output)
    {
      if (output.size() != rows_)
      {
        throw std::invalid_argument("Rational: output dimension mismatch");
      }
      output_ = output;
    }

    template <typename scalar_type, typename index_type>
    void Rational<scalar_type, index_type>::attachInput(SignalT* a, SignalT* b, SignalT* c)
    {
      attachInput(std::vector<SignalT*>{a, b, c});
    }

    template <typename scalar_type, typename index_type>
    void Rational<scalar_type, index_type>::attachOutput(SignalT* a, SignalT* b, SignalT* c)
    {
      attachOutput(std::vector<SignalT*>{a, b, c});
    }

    template <typename scalar_type, typename index_type>
    bool Rational<scalar_type, index_type>::hasInputDerivative(IdxT k) const
    {
      for (size_t n = 0; n < rows_; ++n)
      {
        if (E_[n][static_cast<size_t>(k)] != RealT{0})
        {
          return true;
        }
      }
      return false;
    }

    template <typename scalar_type, typename index_type>
    bool Rational<scalar_type, index_type>::hasFeedthroughDerivative() const
    {
      for (size_t k = 0; k < cols_; ++k)
      {
        if (hasInputDerivative(static_cast<IdxT>(k)))
        {
          return true;
        }
      }
      return false;
    }

    template <typename scalar_type, typename index_type>
    int Rational<scalar_type, index_type>::setGridKitComponentID(IdxT id)
    {
      gridkit_component_id_ = id;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Rational<scalar_type, index_type>::verify() const
    {
      int errors = errors_;
      for (size_t k = 0; k < cols_; ++k)
      {
        if (input_[k] == nullptr)
        {
          ++errors;
          continue;
        }
        if (hasInputDerivative(static_cast<IdxT>(k)) && input_[k]->computed())
        {
          ++errors;
        }
        if (coupling_allocated_ && (!input_[k]->linked() || (hasInputDerivative(static_cast<IdxT>(k)) && !input_[k]->derivativeLinked())))
        {
          ++errors;
        }
      }
      for (auto* output : output_)
      {
        if (output && coupling_allocated_ && !output->residualLinked())
        {
          ++errors;
        }
      }
      return errors;
    }

    template <typename scalar_type, typename index_type>
    int Rational<scalar_type, index_type>::allocate()
    {
      if (verify() != 0)
      {
        return 1;
      }
      if (!allocated_)
      {
        this->allocateVectors(size_);
      }
      const auto size = static_cast<size_t>(size_);
      tag_.resize(size);
      variable_indices_.resize(size);
      residual_indices_.resize(size);
      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, j);
        this->setResidualIndex(j, j);
      }
      this->allocateExternalVectors(cols(), rows());
      for (IdxT k = 0; k < cols(); ++k)
      {
        this->setExternalVariableSignal(k, input_[static_cast<size_t>(k)]);
      }
      for (IdxT n = 0; n < rows(); ++n)
      {
        this->setExternalResidualSignal(n, output_[static_cast<size_t>(n)]);
      }
      allocated_          = true;
      coupling_allocated_ = true;
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Rational<scalar_type, index_type>::initializeState(const std::map<std::string, RealT>& values)
    {
      if (!values.empty())
        throw std::invalid_argument("Rational states are initialized from their inputs");
      return initialize();
    }

    template <typename scalar_type, typename index_type>
    typename Component<scalar_type, index_type>::InitializationPortsT Rational<scalar_type, index_type>::initializationPorts()
    {
      return {input_, {}, {}};
    }

    template <typename scalar_type, typename index_type>
    int Rational<scalar_type, index_type>::initializeSteadyState(RealT omega)
    {
      std::vector<RealT> u(cols_), udot(cols_);
      for (size_t k = 0; k < cols_; ++k)
      {
        u[k]    = static_cast<RealT>(input_[k]->read());
        udot[k] = static_cast<RealT>(input_[k]->readDerivative());
      }
      return initializeSteadyState(omega, u, udot);
    }

    template <typename scalar_type, typename index_type>
    int Rational<scalar_type, index_type>::initialize()
    {
      if (!coupling_allocated_ || verify() != 0)
      {
        return 1;
      }
      if (size_ == 0)
      {
        return 0;
      }
      auto* y  = y_.getData();
      auto* yp = yp_.getData();
      for (IdxT j = 0; j < size_; ++j)
      {
        y[j]  = RealT{0};
        yp[j] = RealT{0};
      }
      // Zero memory is a valid transient start; its initial slope is B u.
      evaluateInternalResidual();
      for (IdxT j = 0; j < size_; ++j)
      {
        yp[j] = f_.getData()[j];
      }
      y_.setDataUpdated();
      yp_.setDataUpdated();
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Rational<scalar_type, index_type>::initializeSteadyState(RealT omega, std::span<const RealT> u, std::span<const RealT> udot)
    {
      if (!coupling_allocated_ || errors_ || !std::isfinite(omega) || u.size() != cols_ || udot.size() != cols_)
      {
        return 1;
      }
      for (size_t k = 0; k < cols_; ++k)
      {
        if (!std::isfinite(u[k]) || !std::isfinite(udot[k]) || (omega == RealT{0} && udot[k] != RealT{0}))
        {
          return 1;
        }
      }
      // Validate the forcing frequency before changing any state.
      for (const auto& section : sections_)
      {
        if (ComplexT{0, omega} == ComplexT{section.a, section.w}
            || (section.pair && ComplexT{0, omega} == ComplexT{section.a, -section.w}))
        {
          return 1;
        }
      }
      std::vector<RealT> y(static_cast<size_t>(size_)), yp(y.size());
      size_t             offset = 0;
      for (const auto& section : sections_)
      {
        for (size_t j = 0; j < section.order; ++j)
        {
          ComplexT forcing_plus{}, forcing_minus{};
          for (size_t k = 0; k < cols_; ++k)
          {
            const ComplexT U{u[k], omega == RealT{0} ? RealT{0} : -udot[k] / omega};
            const ComplexT B{section.Br[j][k], section.Bi[j][k]};
            forcing_plus  += B * U;
            forcing_minus += std::conj(B) * U;
          }
          const ComplexT plus = forcing_plus / ComplexT{-section.a, omega - section.w};
          ComplexT       W    = plus, V{};
          if (section.pair)
          {
            const ComplexT minus = forcing_minus / ComplexT{-section.a, omega + section.w};
            W                    = (plus + minus) / RealT{2};
            V                    = (plus - minus) / ComplexT{0, 2};
          }
          y[offset + j]  = W.real();
          yp[offset + j] = -omega * W.imag();
          if (section.pair)
          {
            y[offset + section.order + j]  = V.real();
            yp[offset + section.order + j] = -omega * V.imag();
          }
          if (!std::isfinite(W.real()) || !std::isfinite(W.imag()) || !std::isfinite(V.real()) || !std::isfinite(V.imag()))
          {
            return 1;
          }
        }
        offset += section.order * (section.pair ? 2 : 1);
      }
      for (size_t j = 0; j < y.size(); ++j)
      {
        y_.getData()[j]  = y[j];
        yp_.getData()[j] = yp[j];
      }
      y_.setDataUpdated();
      yp_.setDataUpdated();
      return 0;
    }

    template <typename scalar_type, typename index_type>
    void Rational<scalar_type, index_type>::transfer(RealT omega, MatrixT& re, MatrixT& im) const
    {
      if (!std::isfinite(omega))
      {
        throw std::invalid_argument("Rational: frequency must be finite");
      }
      re = MatrixT(rows_, cols_);
      im = MatrixT(rows_, cols_);
      for (size_t n = 0; n < rows_; ++n)
      {
        for (size_t k = 0; k < cols_; ++k)
        {
          ComplexT H{D_[n][k], omega * E_[n][k]};
          for (const auto& section : sections_)
          {
            const ComplexT denominator{-section.a, omega - section.w};
            const ComplexT partner{-section.a, omega + section.w};
            if (denominator == ComplexT{} || (section.pair && partner == ComplexT{}))
            {
              throw std::domain_error("Rational: transfer requested at a pole");
            }
            for (size_t j = 0; j < section.order; ++j)
            {
              const ComplexT residue = ComplexT{section.Cr[n][j], section.Ci[n][j]}
                                       * ComplexT{section.Br[j][k], section.Bi[j][k]};
              if (section.pair)
              {
                H += (residue / denominator + std::conj(residue) / partner) / RealT{2};
              }
              else
              {
                H += residue / denominator;
              }
            }
          }
          re[n][k] = H.real();
          im[n][k] = H.imag();
        }
      }
    }

    template <typename scalar_type, typename index_type>
    void Rational<scalar_type, index_type>::transfer(RealT omega, ABCMatrix<RealT>& re, ABCMatrix<RealT>& im) const
    {
      if (rows_ != 3 || cols_ != 3)
      {
        throw std::invalid_argument("Rational: three-phase transfer dimension mismatch");
      }
      MatrixT dynamic_re, dynamic_im;
      transfer(omega, dynamic_re, dynamic_im);
      for (size_t n = 0; n < 3; ++n)
      {
        for (size_t k = 0; k < 3; ++k)
        {
          re[n][k] = dynamic_re[n][k];
          im[n][k] = dynamic_im[n][k];
        }
      }
    }

    template <typename scalar_type, typename index_type>
    int Rational<scalar_type, index_type>::setAbsoluteTolerance(RealT tolerance)
    {
      size_t offset = 0;
      for (const auto& section : sections_)
      {
        const RealT floor = tolerance / std::max(RealT{1}, std::hypot(section.a, section.w));
        for (size_t j = 0; j < section.order * (section.pair ? 2 : 1); ++j)
        {
          abs_tol_.getData()[offset++] = floor;
        }
      }
      abs_tol_.setDataUpdated();
      return 0;
    }

    template <typename scalar_type, typename index_type>
    typename Rational<scalar_type, index_type>::ScalarT Rational<scalar_type, index_type>::output(IdxT n) const
    {
      auto* model = const_cast<Rational*>(this);
      model->gatherExternalVariables();
      return evaluateOutput(n, size_ == 0 ? nullptr : y_.getData(), y_ext_.data(), yp_ext_.data());
    }

    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) typename Rational<scalar_type, index_type>::ScalarT Rational<scalar_type, index_type>::evaluateOutput(IdxT n, const ScalarT* y, const ScalarT* y_ext, const ScalarT* yp_ext) const
    {
      ScalarT      value{0};
      const size_t row = static_cast<size_t>(n);
      for (size_t k = 0; k < cols_; ++k)
      {
        if (D_[row][k] != RealT{0})
        {
          value += D_[row][k] * y_ext[k];
        }
        if (E_[row][k] != RealT{0})
        {
          value += E_[row][k] * yp_ext[k];
        }
      }
      size_t offset = 0;
      for (const auto& section : sections_)
      {
        for (size_t j = 0; j < section.order; ++j)
        {
          if (section.Cr[row][j] != RealT{0})
          {
            value += section.Cr[row][j] * y[offset + j];
          }
          if (section.pair && section.Ci[row][j] != RealT{0})
          {
            value -= section.Ci[row][j] * y[offset + section.order + j];
          }
        }
        offset += section.order * (section.pair ? 2 : 1);
      }
      return value;
    }

    template <typename scalar_type, typename index_type>
    typename Rational<scalar_type, index_type>::ScalarT Rational<scalar_type, index_type>::outputDerivative(IdxT n) const
    {
      if (hasFeedthroughDerivative())
        throw std::logic_error("Rational: output derivative requires a proper transfer");
      auto* model = const_cast<Rational*>(this);
      for (size_t k = 0; k < cols_; ++k)
        if (D_[static_cast<size_t>(n)][k] != ZERO<RealT>)
          model->yp_ext_[k] = input_[k]->readDerivative();
      return evaluateOutput(n, size_ == 0 ? nullptr : yp_.getData(), yp_ext_.data(), nullptr);
    }

    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int Rational<scalar_type, index_type>::evaluateInternalResidual(const ScalarT* y,
                                                                                                   const ScalarT* yp,
                                                                                                   const ScalarT* y_ext,
                                                                                                   const ScalarT*,
                                                                                                   ScalarT* f)
    {
      size_t offset = 0;
      for (const auto& section : sections_)
      {
        for (size_t j = 0; j < section.order; ++j)
        {
          const size_t w = offset + j, v = w + section.order;
          ScalarT      fw = -yp[w] + section.a * y[w];
          ScalarT      fv{0};
          if (section.pair)
          {
            fw -= section.w * y[v];
            fv  = -yp[v] + section.w * y[w] + section.a * y[v];
          }
          for (size_t k = 0; k < cols_; ++k)
          {
            if (section.Br[j][k] != RealT{0})
            {
              fw += section.Br[j][k] * y_ext[k];
            }
            if (section.pair && section.Bi[j][k] != RealT{0})
            {
              fv += section.Bi[j][k] * y_ext[k];
            }
          }
          f[w] = fw;
          if (section.pair)
            f[v] = fv;
        }
        offset += section.order * (section.pair ? 2 : 1);
      }
      return 0;
    }

    template <typename scalar_type, typename index_type>
    __attribute__((always_inline)) int Rational<scalar_type, index_type>::evaluateExternalResidual(const ScalarT* y,
                                                                                                   const ScalarT*,
                                                                                                   const ScalarT* y_ext,
                                                                                                   const ScalarT* yp_ext,
                                                                                                   ScalarT*       f_ext)
    {
      for (IdxT n = 0; n < rows(); ++n)
        f_ext[n] = evaluateOutput(n, y, y_ext, yp_ext);
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int Rational<scalar_type, index_type>::evaluateInternalResidual()
    {
      if (!coupling_allocated_ || errors_)
        return 1;
      if (size_ == 0)
        return 0;
      this->gatherExternalVariables();
      const int status = evaluateInternalResidual(y_.getData(), yp_.getData(), y_ext_.data(), yp_ext_.data(), f_.getData());
      f_.setDataUpdated();
      return status;
    }

    template <typename scalar_type, typename index_type>
    int Rational<scalar_type, index_type>::evaluateExternalResidual()
    {
      if (!coupling_allocated_ || errors_)
        return 1;
      this->gatherExternalVariables();
      const int status = evaluateExternalResidual(size_ == 0 ? nullptr : y_.getData(), nullptr, y_ext_.data(), yp_ext_.data(), f_ext_.data());
      this->scatterExternalResidual();
      return status;
    }

    template <typename scalar_type, typename index_type>
    int Rational<scalar_type, index_type>::evaluateResidual()
    {
      const int status = evaluateInternalResidual();
      return status == 0 ? evaluateExternalResidual() : status;
    }

    template <typename scalar_type, typename index_type>
    typename Rational<scalar_type, index_type>::IdxT Rational<scalar_type, index_type>::jacobianCapacity() const
    {
      size_t capacity = 2 * rows_ * cols_;
      for (const auto& section : sections_)
      {
        capacity += section.order * (section.pair ? 2 : 1) * (3 + cols_ + rows_);
      }
      return static_cast<IdxT>(capacity);
    }

    template <typename scalar_type, typename index_type>
    Rational<scalar_type, index_type>::Rational(size_t rows, size_t cols, const MatrixT& D, const MatrixT& E, RealT scale, int errors)
      : rows_(rows), cols_(cols), D_(rows, cols), E_(rows, cols), input_(cols), output_(rows), errors_(errors)
    {
      if (!std::isfinite(scale))
      {
        ++errors_;
      }
      if (errors_ != 0)
      {
        return;
      }
      for (size_t n = 0; n < rows_; ++n)
      {
        for (size_t k = 0; k < cols_; ++k)
        {
          D_[n][k] = scale * D[n][k];
          E_[n][k] = scale * E[n][k];
        }
      }
    }

    template <typename scalar_type, typename index_type>
    void Rational<scalar_type, index_type>::finishSections()
    {
      for (const auto& section : sections_)
      {
        size_ += static_cast<IdxT>(section.order * (section.pair ? 2 : 1));
        for (const auto* matrix : {&section.Br, &section.Bi, &section.Cr, &section.Ci})
        {
          for (const auto& row : *matrix)
          {
            for (const auto value : row)
            {
              if (!std::isfinite(value))
              {
                ++errors_;
              }
            }
          }
        }
      }
      for (const auto* matrix : {&D_, &E_})
      {
        for (const auto& row : *matrix)
        {
          for (const auto value : row)
          {
            if (!std::isfinite(value))
            {
              ++errors_;
            }
          }
        }
      }
    }
  } // namespace EMT
} // namespace GridKit
