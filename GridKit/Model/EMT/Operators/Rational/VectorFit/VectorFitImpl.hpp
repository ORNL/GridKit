#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>

#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFit.hpp>
#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFitData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief Constructor for a three-phase rational operator
     *
     * System sizes:
     * - Number of equations = 3 Qr + 6 Qc
     * - Number of independent variables = 3 Qr + 6 Qc
     *
     * The poles are partitioned into real sections and conjugate-pair
     * sections by the specification's adjacency rule, and the consumer scale
     * factor, along with the pair factor of two, is folded into the stored
     * coefficient matrices. Complex values do not survive construction.
     */
    template <typename scalar_type, typename index_type>
    VectorFit<scalar_type, index_type>::VectorFit(const ModelDataT& data, RealT scale)
    {
      pole_error_count_ = static_cast<IdxT>(data.validate());

      for (size_t n = 0; n < 3; ++n)
      {
        for (size_t k = 0; k < 3; ++k)
        {
          feedthrough_.D[n][k] = scale * data.D[n][k];
          feedthrough_.E[n][k] = scale * data.E[n][k];
        }
      }

      size_t q = 0;
      while (q < data.poles.size() && q < data.residues.size())
      {
        if (data.poles[q].imag() == 0.0)
        {
          auto& section = real_sections_.emplace_back();
          section.a     = data.poles[q].real();
          for (size_t n = 0; n < 3; ++n)
          {
            for (size_t k = 0; k < 3; ++k)
            {
              section.A[n][k] = scale * data.residues[q][n][k].real();
            }
          }
          ++q;
          continue;
        }

        if (q + 1 >= data.poles.size())
        {
          break;
        }

        auto& section = complex_sections_.emplace_back();
        section.a     = data.poles[q].real();
        section.w     = data.poles[q].imag();
        for (size_t n = 0; n < 3; ++n)
        {
          for (size_t k = 0; k < 3; ++k)
          {
            section.A[n][k] = TWO<RealT> * scale * data.residues[q][n][k].real();
            section.B[n][k] = TWO<RealT> * scale * data.residues[q][n][k].imag();
          }
        }
        q += 2;
      }

      size_ = static_cast<IdxT>(3 * real_sections_.size() + 6 * complex_sections_.size());
    }

    template <typename scalar_type, typename index_type>
    VectorFit<scalar_type, index_type>::~VectorFit()
    {
    }

    /**
     * @brief Set the component ID
     */
    template <typename scalar_type, typename index_type>
    int VectorFit<scalar_type, index_type>::setGridKitComponentID(IdxT component_id)
    {
      gridkit_component_id_ = component_id;
      return 0;
    }

    /**
     * @brief Whether the operator reads the input derivative.
     */
    template <typename scalar_type, typename index_type>
    bool VectorFit<scalar_type, index_type>::hasFeedthroughDerivative() const
    {
      bool nonzero = false;
      for (size_t n = 0; n < 3; ++n)
      {
        for (size_t k = 0; k < 3; ++k)
        {
          if (feedthrough_.E[n][k] != 0.0)
          {
            nonzero = true;
          }
        }
      }
      return nonzero;
    }

    /**
     * @brief Jacobian triplet capacity of this operator.
     */
    template <typename scalar_type, typename index_type>
    index_type VectorFit<scalar_type, index_type>::jacobianCapacity() const
    {
      // Dense block capacities per section: a real section carries 9 entries
      // in each of its four blocks; a pair section carries 36 + 18 + 18 + 36;
      // the feedthrough carries 9 + 9.
      return static_cast<IdxT>(36 * real_sections_.size()
                               + 108 * complex_sections_.size() + 18);
    }

    /*!
     * @brief allocate method resizes local storage and registers coupling signals.
     */
    template <typename scalar_type, typename index_type>
    int VectorFit<scalar_type, index_type>::allocate()
    {
      if (!allocated_)
      {
        this->allocateVectors(size_);
      }

      auto size = static_cast<size_t>(size_); // avoid compiler warnings

      tag_.resize(size);

      variable_indices_.resize(size);
      residual_indices_.resize(size);

      // Default variable and residual index mapping to local index
      for (IdxT j = 0; j < size_; ++j)
      {
        this->setVariableIndex(j, j);
        this->setResidualIndex(j, j);
      }

      // Resize coupling data
      this->allocateExternalVectors(3, 3);
      for (IdxT n = 0; n < 3; ++n)
      {
        this->setExternalVariableSignal(n, input_[static_cast<size_t>(n)]);
        this->setExternalResidualSignal(n, output_[static_cast<size_t>(n)]);
      }

      // The feedthrough reads the input derivative, so the input variables
      // become differential.
      if (hasFeedthroughDerivative())
      {
        for (IdxT n = 0; n < 3; ++n)
        {
          input_[static_cast<size_t>(n)]->markDerivativeCoupling();
        }
      }

      allocated_ = true;
      return 0;
    }

    /**
     * @brief Check model correctness
     *
     * @return Number of model configuration errors found
     */
    template <typename scalar_type, typename index_type>
    int VectorFit<scalar_type, index_type>::verify() const
    {
      int error_count = static_cast<int>(pole_error_count_);

      for (size_t n = 0; n < 3; ++n)
      {
        if (input_[n] == nullptr || output_[n] == nullptr)
        {
          Log::error() << "VectorFit: the input and output signals must be attached\n";
          ++error_count;
          break;
        }
      }

      return error_count;
    }

    /**
     * Initialization of the rational operator model
     *
     * The memory states start de-energized; consumers with a sinusoidal
     * steady state call initializeSteadyState afterwards.
     */
    template <typename scalar_type, typename index_type>
    int VectorFit<scalar_type, index_type>::initialize()
    {
      auto* y  = y_.getData();
      auto* yp = yp_.getData();

      for (IdxT j = 0; j < size_; ++j)
      {
        y[j]  = 0.0;
        yp[j] = 0.0;
      }

      y_.setDataUpdated();
      yp_.setDataUpdated();

      return 0;
    }

    /**
     * @brief Set the memory states and their derivatives to the sinusoidal
     * steady state.
     *
     * The instantaneous input pair at the initialization instant determines
     * the rotating input phasor per phase; each section's state phasors
     * follow in closed form, and the states and derivatives are its samples.
     */
    template <typename scalar_type, typename index_type>
    int VectorFit<scalar_type, index_type>::initializeSteadyState(RealT                   omega0,
                                                                  const ABCVector<RealT>& u,
                                                                  const ABCVector<RealT>& u_dot)
    {
      auto* y  = y_.getData();
      auto* yp = yp_.getData();

      IdxT offset = 0;
      for (const auto& section : real_sections_)
      {
        // W = U / (j omega0 - a), evaluated in real pairs
        const RealT d_re  = -section.a;
        const RealT d_im  = omega0;
        const RealT d_mag = d_re * d_re + d_im * d_im;
        for (IdxT n = 0; n < 3; ++n)
        {
          const RealT u_re = u[static_cast<size_t>(n)];
          const RealT u_im = -u_dot[static_cast<size_t>(n)] / omega0;
          const RealT w_re = (u_re * d_re + u_im * d_im) / d_mag;
          const RealT w_im = (u_im * d_re - u_re * d_im) / d_mag;

          y[offset + n]  = w_re;
          yp[offset + n] = -omega0 * w_im;
        }
        offset += 3;
      }

      for (const auto& section : complex_sections_)
      {
        // With g = j omega0 - a and det = g^2 + omega_q^2:
        //   W = g U / det and V = omega_q U / det
        const RealT g_re    = -section.a;
        const RealT g_im    = omega0;
        const RealT det_re  = g_re * g_re - g_im * g_im + section.w * section.w;
        const RealT det_im  = TWO<RealT> * g_re * g_im;
        const RealT det_mag = det_re * det_re + det_im * det_im;
        for (IdxT n = 0; n < 3; ++n)
        {
          const RealT u_re = u[static_cast<size_t>(n)];
          const RealT u_im = -u_dot[static_cast<size_t>(n)] / omega0;

          const RealT gu_re = g_re * u_re - g_im * u_im;
          const RealT gu_im = g_re * u_im + g_im * u_re;

          const RealT w_re = (gu_re * det_re + gu_im * det_im) / det_mag;
          const RealT w_im = (gu_im * det_re - gu_re * det_im) / det_mag;
          const RealT v_re = section.w * (u_re * det_re + u_im * det_im) / det_mag;
          const RealT v_im = section.w * (u_im * det_re - u_re * det_im) / det_mag;

          y[offset + n]      = w_re;
          yp[offset + n]     = -omega0 * w_im;
          y[offset + 3 + n]  = v_re;
          yp[offset + 3 + n] = -omega0 * v_im;
        }
        offset += 6;
      }

      y_.setDataUpdated();
      yp_.setDataUpdated();

      return 0;
    }

    /**
     * @brief The frequency response at one angular frequency, in real pairs.
     *
     * The stored coefficient matrices already carry the consumer scale and
     * the pair factor of two, so the response matches the assembled residual
     * contributions exactly.
     */
    template <typename scalar_type, typename index_type>
    void VectorFit<scalar_type, index_type>::transfer(RealT             omega0,
                                                      ABCMatrix<RealT>& H_re,
                                                      ABCMatrix<RealT>& H_im) const
    {
      for (size_t n = 0; n < 3; ++n)
      {
        for (size_t k = 0; k < 3; ++k)
        {
          H_re[n][k] = feedthrough_.D[n][k];
          H_im[n][k] = omega0 * feedthrough_.E[n][k];
        }
      }

      for (const auto& section : real_sections_)
      {
        // c = 1 / (j omega0 - a)
        const RealT d_re  = -section.a;
        const RealT d_im  = omega0;
        const RealT d_mag = d_re * d_re + d_im * d_im;
        const RealT c_re  = d_re / d_mag;
        const RealT c_im  = -d_im / d_mag;
        for (size_t n = 0; n < 3; ++n)
        {
          for (size_t k = 0; k < 3; ++k)
          {
            H_re[n][k] += section.A[n][k] * c_re;
            H_im[n][k] += section.A[n][k] * c_im;
          }
        }
      }

      for (const auto& section : complex_sections_)
      {
        // The pair sum R/(s - p) + R*/(s - p*) at s = j omega0 reduces to
        // (A g - B omega_q) / det with g = j omega0 - a and
        // det = g^2 + omega_q^2, matching the state realization exactly
        const RealT g_re    = -section.a;
        const RealT g_im    = omega0;
        const RealT det_re  = g_re * g_re - g_im * g_im + section.w * section.w;
        const RealT det_im  = TWO<RealT> * g_re * g_im;
        const RealT det_mag = det_re * det_re + det_im * det_im;
        for (size_t n = 0; n < 3; ++n)
        {
          for (size_t k = 0; k < 3; ++k)
          {
            const RealT num_re  = section.A[n][k] * g_re - section.B[n][k] * section.w;
            const RealT num_im  = section.A[n][k] * g_im;
            H_re[n][k]         += (num_re * det_re + num_im * det_im) / det_mag;
            H_im[n][k]         += (num_im * det_re - num_re * det_im) / det_mag;
          }
        }
      }
    }

    /**
     * \brief Identify differential variables.
     */
    template <typename scalar_type, typename index_type>
    int VectorFit<scalar_type, index_type>::tagDifferentiable()
    {
      for (size_t j = 0; j < static_cast<size_t>(size_); ++j)
      {
        tag_[j] = true;
      }

      return 0;
    }

    /**
     * @brief Compute the absolute tolerance for each variable in the model
     *
     * @param rel_tol The relative tolerance which can be used to pick the
     *        absolute tolerance.
     * @tparam scalar_type Scalar data type
     * @tparam index_type Index data type
     * @return int 0 if successful, non-zero otherwise.
     *
     * A memory state scales like the input over the pole distance, so fast
     * poles carry proportionally smaller states and receive proportionally
     * smaller absolute floors.
     */
    template <typename scalar_type, typename index_type>
    int VectorFit<scalar_type, index_type>::setAbsoluteTolerance(RealT rel_tol)
    {
      if (size_ == 0)
      {
        return 0;
      }

      auto* abs_tol = abs_tol_.getData();

      IdxT offset = 0;
      for (const auto& section : real_sections_)
      {
        const RealT floor = rel_tol / std::max(ONE<RealT>, std::abs(section.a));
        for (IdxT n = 0; n < 3; ++n)
        {
          abs_tol[offset + n] = static_cast<ScalarT>(floor);
        }
        offset += 3;
      }

      for (const auto& section : complex_sections_)
      {
        const RealT magnitude = std::sqrt(section.a * section.a + section.w * section.w);
        const RealT floor     = rel_tol / std::max(ONE<RealT>, magnitude);
        for (IdxT n = 0; n < 6; ++n)
        {
          abs_tol[offset + n] = static_cast<ScalarT>(floor);
        }
        offset += 6;
      }

      abs_tol_.setDataUpdated();

      return 0;
    }

    /**
     * @brief Evaluate the operator output from the current states.
     */
    template <typename scalar_type, typename index_type>
    scalar_type VectorFit<scalar_type, index_type>::output(IdxT n) const
    {
      const auto* y   = y_.getData();
      const auto* ye  = y_ext_.data();
      const auto* ype = yp_ext_.data();

      ScalarT value{0.0};

      const auto row = static_cast<size_t>(n);
      for (size_t k = 0; k < 3; ++k)
      {
        value += feedthrough_.D[row][k] * ye[k] + feedthrough_.E[row][k] * ype[k];
      }

      IdxT offset = 0;
      for (const auto& section : real_sections_)
      {
        for (size_t k = 0; k < 3; ++k)
        {
          value += section.A[row][k] * y[offset + static_cast<IdxT>(k)];
        }
        offset += 3;
      }

      for (const auto& section : complex_sections_)
      {
        for (size_t k = 0; k < 3; ++k)
        {
          value += section.A[row][k] * y[offset + static_cast<IdxT>(k)]
                   - section.B[row][k] * y[offset + 3 + static_cast<IdxT>(k)];
        }
        offset += 6;
      }

      return value;
    }

    template <typename scalar_type, typename index_type>
    int VectorFit<scalar_type, index_type>::evaluateInternalResidual()
    {
      this->gatherExternalVariables();

      auto*       y   = y_.getData();
      auto*       yp  = yp_.getData();
      auto*       f   = f_.getData();
      const auto* ye  = y_ext_.data();
      const auto* ype = yp_ext_.data();

      IdxT offset = 0;
      for (auto& section : real_sections_)
      {
        section.evaluateInternalResidual(y + offset, yp + offset, ye, ype, f + offset);
        offset += 3;
      }
      for (auto& section : complex_sections_)
      {
        section.evaluateInternalResidual(y + offset, yp + offset, ye, ype, f + offset);
        offset += 6;
      }

      f_.setDataUpdated();

      return 0;
    }

    /**
     * @brief External residual contributions to the output rows.
     *
     * Section contributions accumulate locally and scatter once.
     */
    template <typename scalar_type, typename index_type>
    int VectorFit<scalar_type, index_type>::evaluateExternalResidual()
    {
      auto*       y   = y_.getData();
      auto*       yp  = yp_.getData();
      const auto* ye  = y_ext_.data();
      const auto* ype = yp_ext_.data();

      ScalarT contribution[3]{};

      feedthrough_.evaluateExternalResidual(y, yp, ye, ype, f_ext_.data());

      IdxT offset = 0;
      for (auto& section : real_sections_)
      {
        section.evaluateExternalResidual(y + offset, yp + offset, ye, ype, contribution);
        f_ext_[0] += contribution[0];
        f_ext_[1] += contribution[1];
        f_ext_[2] += contribution[2];
        offset    += 3;
      }
      for (auto& section : complex_sections_)
      {
        section.evaluateExternalResidual(y + offset, yp + offset, ye, ype, contribution);
        f_ext_[0] += contribution[0];
        f_ext_[1] += contribution[1];
        f_ext_[2] += contribution[2];
        offset    += 6;
      }

      this->scatterExternalResidual();

      return 0;
    }

    /**
     * @brief Residual contribution of the operator is pushed to the output rows.
     *
     */
    template <typename scalar_type, typename index_type>
    int VectorFit<scalar_type, index_type>::evaluateResidual()
    {
      evaluateInternalResidual();
      return evaluateExternalResidual();
    }

    /**
     * @brief Invalidate the Jacobian structure and the linear cache.
     */
    template <typename scalar_type, typename index_type>
    void VectorFit<scalar_type, index_type>::resetJacobianStructure()
    {
      jacobian_cached_ = false;
      Component<scalar_type, index_type>::resetJacobianStructure();
    }

  } // namespace EMT
} // namespace GridKit
