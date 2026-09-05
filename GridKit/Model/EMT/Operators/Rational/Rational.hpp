#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <span>
#include <stdexcept>
#include <vector>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Operators/Rational/RationalMatrix.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Real realization shared by the matrix-residue and factorized operators.
     *
     * Each section realizes x' = p x + B u with real storage for real poles
     * and real/imaginary storage for a conjugate pair. C maps its states into
     * the output rows. Section factors already include the consumer scale and
     * the pair factor of two. VectorFit uses K states per pole; StateSpace uses one.
     */
    template <typename scalar_type, typename index_type>
    class Rational : public Component<scalar_type, index_type>
    {
    public:
      using ScalarT  = scalar_type;
      using IdxT     = index_type;
      using Base     = Component<ScalarT, IdxT>;
      using RealT    = typename Base::RealT;
      using SignalT  = Signal<ScalarT, IdxT>;
      using Port3T   = Port3<ScalarT, IdxT>;
      using MatrixT  = RationalMatrix<RealT>;
      using ComplexT = std::complex<RealT>;

      IdxT rows() const
      {
        return static_cast<IdxT>(rows_);
      }

      IdxT cols() const
      {
        return static_cast<IdxT>(cols_);
      }

      void attachInput(const std::vector<SignalT*>& input)
      {
        if (input.size() != cols_)
        {
          throw std::invalid_argument("Rational: input dimension mismatch");
        }
        input_ = input;
      }

      void attachOutput(const std::vector<SignalT*>& output)
      {
        if (output.size() != rows_)
        {
          throw std::invalid_argument("Rational: output dimension mismatch");
        }
        output_ = output;
      }

      void attachInput(SignalT* a, SignalT* b, SignalT* c)
      {
        attachInput(std::vector<SignalT*>{a, b, c});
      }

      void attachOutput(SignalT* a, SignalT* b, SignalT* c)
      {
        attachOutput(std::vector<SignalT*>{a, b, c});
      }

      void attachInput(Port3T* port)
      {
        attachInput(port->a(), port->b(), port->c());
      }

      void attachOutput(Port3T* port)
      {
        attachOutput(port->a(), port->b(), port->c());
      }

      bool hasInputDerivative(IdxT k) const
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

      bool hasFeedthroughDerivative() const
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

      int setGridKitComponentID(IdxT id) override
      {
        this->gridkit_component_id_ = id;
        return 0;
      }

      int verify() const override
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
          if (output == nullptr || (coupling_allocated_ && !output->residualLinked()))
          {
            ++errors;
          }
        }
        return errors;
      }

      int allocate() override
      {
        if (verify() != 0)
        {
          return 1;
        }
        if (!this->allocated_)
        {
          this->allocateVectors(this->size_);
        }
        const auto size = static_cast<size_t>(this->size_);
        this->tag_.resize(size);
        this->variable_indices_.resize(size);
        this->residual_indices_.resize(size);
        for (IdxT j = 0; j < this->size_; ++j)
        {
          this->setVariableIndex(j, j);
          this->setResidualIndex(j, j);
        }
        this->allocateExternalVectors(cols(), rows());
        for (IdxT k = 0; k < cols(); ++k)
        {
          this->setExternalVariableSignal(k, input_[static_cast<size_t>(k)]);
          if (hasInputDerivative(k))
          {
            input_[static_cast<size_t>(k)]->markDerivativeCoupling();
          }
        }
        for (IdxT n = 0; n < rows(); ++n)
        {
          this->setExternalResidualSignal(n, output_[static_cast<size_t>(n)]);
        }
        this->allocated_    = true;
        coupling_allocated_ = true;
        return 0;
      }

      int initialize() override
      {
        if (!coupling_allocated_ || verify() != 0)
        {
          return 1;
        }
        if (this->size_ == 0)
        {
          return 0;
        }
        auto* y  = this->y_.getData();
        auto* yp = this->yp_.getData();
        for (IdxT j = 0; j < this->size_; ++j)
        {
          y[j]  = RealT{0};
          yp[j] = RealT{0};
        }
        // Zero memory is a valid transient start; its initial slope is B u.
        evaluateInternalResidual();
        for (IdxT j = 0; j < this->size_; ++j)
        {
          yp[j] = this->f_.getData()[j];
        }
        this->y_.setDataUpdated();
        this->yp_.setDataUpdated();
        return 0;
      }

      /** Initialize sinusoidal steady state, or constant equilibrium at omega=0. */
      int initializeSteadyState(RealT omega, std::span<const RealT> u, std::span<const RealT> udot)
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
        std::vector<RealT> y(static_cast<size_t>(this->size_)), yp(y.size());
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
          this->y_.getData()[j]  = y[j];
          this->yp_.getData()[j] = yp[j];
        }
        this->y_.setDataUpdated();
        this->yp_.setDataUpdated();
        return 0;
      }

      void transfer(RealT omega, MatrixT& re, MatrixT& im) const
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

      void transfer(RealT omega, ABCMatrix<RealT>& re, ABCMatrix<RealT>& im) const
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

      int tagDifferentiable() override
      {
        std::fill(this->tag_.begin(), this->tag_.end(), true);
        return 0;
      }

      int setAbsoluteTolerance(RealT tolerance) override
      {
        size_t offset = 0;
        for (const auto& section : sections_)
        {
          const RealT floor = tolerance / std::max(RealT{1}, std::hypot(section.a, section.w));
          for (size_t j = 0; j < section.order * (section.pair ? 2 : 1); ++j)
          {
            this->abs_tol_.getData()[offset++] = floor;
          }
        }
        this->abs_tol_.setDataUpdated();
        return 0;
      }

      ScalarT output(IdxT n) const
      {
        ScalarT      value{0};
        const size_t row = static_cast<size_t>(n);
        for (size_t k = 0; k < cols_; ++k)
        {
          if (D_[row][k] != RealT{0})
          {
            value += D_[row][k] * input_[k]->read();
          }
          if (E_[row][k] != RealT{0})
          {
            value += E_[row][k] * input_[k]->readDerivative();
          }
        }
        size_t      offset = 0;
        const auto* y      = this->size_ == 0 ? nullptr : this->y_.getData();
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

      int evaluateInternalResidual() override
      {
        if (!coupling_allocated_ || errors_)
        {
          return 1;
        }
        this->gatherExternalVariables();
        if (this->size_ == 0)
        {
          return 0;
        }
        const auto* y      = this->y_.getData();
        const auto* yp     = this->yp_.getData();
        auto*       f      = this->f_.getData();
        size_t      offset = 0;
        for (const auto& section : sections_)
        {
          for (size_t j = 0; j < section.order; ++j)
          {
            const size_t w = offset + j, v = w + section.order;
            f[w] = -yp[w] + section.a * y[w];
            if (section.pair)
            {
              f[w] -= section.w * y[v];
              f[v]  = -yp[v] + section.w * y[w] + section.a * y[v];
            }
            for (size_t k = 0; k < cols_; ++k)
            {
              if (section.Br[j][k] != RealT{0})
              {
                f[w] += section.Br[j][k] * this->y_ext_[k];
              }
              if (section.pair && section.Bi[j][k] != RealT{0})
              {
                f[v] += section.Bi[j][k] * this->y_ext_[k];
              }
            }
          }
          offset += section.order * (section.pair ? 2 : 1);
        }
        this->f_.setDataUpdated();
        return 0;
      }

      int evaluateExternalResidual() override
      {
        if (!coupling_allocated_ || errors_)
        {
          return 1;
        }
        for (IdxT n = 0; n < rows(); ++n)
        {
          this->f_ext_[static_cast<size_t>(n)] = output(n);
        }
        this->scatterExternalResidual();
        return 0;
      }

      int evaluateResidual() override
      {
        const int status = evaluateInternalResidual();
        return status == 0 ? evaluateExternalResidual() : status;
      }

      IdxT jacobianCapacity() const
      {
        size_t capacity = 2 * rows_ * cols_;
        for (const auto& section : sections_)
        {
          capacity += section.order * (section.pair ? 2 : 1) * (3 + cols_ + rows_);
        }
        return static_cast<IdxT>(capacity);
      }

      /** Exact J = dF/dy + alpha dF/dyp, including computed-input chain rules. */
      int evaluateJacobian() override
      {
        if (!coupling_allocated_ || errors_)
        {
          return 1;
        }
        this->gatherExternalVariables();
        std::vector<typename SignalT::GradientT> gradients(cols_);
        for (size_t k = 0; k < cols_; ++k)
        {
          input_[k]->appendGradient(gradients[k]);
        }
        const size_t capacity = std::max(size_t{1}, static_cast<size_t>(jacobianCapacity()) * this->externalJacobianExpansion());
        // Gradient structure may change with a computed input; resize and rebuild safely.
        if (this->hasComputedInputs() || capacity > capacity_)
        {
          this->resetJacobianStructure();
        }
        if (capacity > capacity_)
        {
          delete[] this->J_rows_buffer_;
          delete[] this->J_cols_buffer_;
          delete[] this->J_vals_buffer_;
          this->J_rows_buffer_ = new IdxT[capacity];
          this->J_cols_buffer_ = new IdxT[capacity];
          this->J_vals_buffer_ = new RealT[capacity];
          capacity_            = capacity;
        }
        this->nnz_  = 0;
        auto append = [&](IdxT row, IdxT col, RealT value)
        {
          if (row == INVALID_INDEX<IdxT> || col == INVALID_INDEX<IdxT>)
          {
            return;
          }
          const auto j            = this->nnz_++;
          this->J_rows_buffer_[j] = row;
          this->J_cols_buffer_[j] = col;
          this->J_vals_buffer_[j] = value;
        };
        auto input = [&](IdxT row, size_t k, RealT value)
        {
          for (const auto& [col, coefficient] : gradients[k])
          {
            append(row, col, value * coefficient);
          }
        };
        for (size_t n = 0; n < rows_; ++n)
        {
          for (size_t k = 0; k < cols_; ++k)
          {
            const auto row = this->residual_indices_ext_[n];
            if (D_[n][k] != RealT{0})
            {
              input(row, k, D_[n][k]);
            }
            if (E_[n][k] != RealT{0})
            {
              append(row, this->variable_indices_ext_[k], this->alpha_ * E_[n][k]);
            }
          }
        }
        size_t offset = 0;
        for (const auto& section : sections_)
        {
          for (size_t j = 0; j < section.order; ++j)
          {
            const size_t w = offset + j, v = w + section.order;
            const auto   rw = this->residual_indices_[w], cw = this->variable_indices_[w];
            append(rw, cw, section.a - this->alpha_);
            if (section.pair)
            {
              append(rw, this->variable_indices_[v], -section.w);
              append(this->residual_indices_[v], cw, section.w);
              append(this->residual_indices_[v], this->variable_indices_[v], section.a - this->alpha_);
            }
            for (size_t k = 0; k < cols_; ++k)
            {
              if (section.Br[j][k] != RealT{0})
              {
                input(rw, k, section.Br[j][k]);
              }
              if (section.pair && section.Bi[j][k] != RealT{0})
              {
                input(this->residual_indices_[v], k, section.Bi[j][k]);
              }
            }
            for (size_t n = 0; n < rows_; ++n)
            {
              if (section.Cr[n][j] != RealT{0})
              {
                append(this->residual_indices_ext_[n], cw, section.Cr[n][j]);
              }
              if (section.pair && section.Ci[n][j] != RealT{0})
              {
                append(this->residual_indices_ext_[n], this->variable_indices_[v], -section.Ci[n][j]);
              }
            }
          }
          offset += section.order * (section.pair ? 2 : 1);
        }
        this->constructCoo();
        return 0;
      }

    protected:
      struct Section
      {
        RealT   a{}, w{};
        bool    pair{false};
        size_t  order{};
        MatrixT Br, Bi, Cr, Ci;

        Section(size_t states, size_t rows, size_t cols)
          : order(states), Br(states, cols), Bi(states, cols), Cr(rows, states), Ci(rows, states)
        {
        }
      };

      Rational(size_t rows, size_t cols, const MatrixT& D, const MatrixT& E, RealT scale, int errors)
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

      void finishSections()
      {
        for (const auto& section : sections_)
        {
          this->size_ += static_cast<IdxT>(section.order * (section.pair ? 2 : 1));
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

      size_t                rows_, cols_;
      MatrixT               D_, E_;
      std::vector<SignalT*> input_, output_;
      std::vector<Section>  sections_;
      int                   errors_{};
      size_t                capacity_{};
      bool                  coupling_allocated_{false};
    };
  } // namespace EMT
} // namespace GridKit
