/**
 * @file VectorFit.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of the EMT vector-fitted rational operator model.
 *
 */

#pragma once

#include <array>
#include <vector>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFitData.hpp>

// Forward declarations.
namespace GridKit
{
  namespace EMT
  {
    template <typename real_type, typename index_type>
    struct VectorFitData;
  } // namespace EMT
} // namespace GridKit

namespace GridKit
{
  namespace EMT
  {
    /**
     * @brief First-order section of a rational operator, one real pole.
     *
     * The section owns the memory-state triple `y[0..2]` and reads the input
     * triple through `y_ext[0..2]`. It is a plain differentiation model for
     * `SparseJacobian`, not a component; the containing `VectorFit` supplies
     * storage, indices, and orchestration.
     */
    template <typename scalar_type, typename index_type>
    struct VectorFitRealSection
    {
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename GridKit::ScalarTraits<ScalarT>::RealT;

      RealT            a{0.0}; ///< \f$a_q\f$
      ABCMatrix<RealT> A{};    ///< Scaled \f$\mathbf{A}_q\f$

      __attribute__((always_inline)) inline int evaluateInternalResidual(
          const ScalarT*                  y,
          const ScalarT*                  yp,
          const ScalarT*                  y_ext,
          [[maybe_unused]] const ScalarT* yp_ext,
          ScalarT*                        f)
      {
        f[0] = -yp[0] + a * y[0] + y_ext[0];
        f[1] = -yp[1] + a * y[1] + y_ext[1];
        f[2] = -yp[2] + a * y[2] + y_ext[2];

        return 0;
      }

      __attribute__((always_inline)) inline int evaluateExternalResidual(
          const ScalarT*                  y,
          [[maybe_unused]] const ScalarT* yp,
          [[maybe_unused]] const ScalarT* y_ext,
          [[maybe_unused]] const ScalarT* yp_ext,
          ScalarT*                        f_ext)
      {
        f_ext[0] = A[0][0] * y[0] + A[0][1] * y[1] + A[0][2] * y[2];
        f_ext[1] = A[1][0] * y[0] + A[1][1] * y[1] + A[1][2] * y[2];
        f_ext[2] = A[2][0] * y[0] + A[2][1] * y[1] + A[2][2] * y[2];

        return 0;
      }
    };

    /**
     * @brief Second-order section of a rational operator, one conjugate pair.
     *
     * The section owns the real memory states `y[0..2]` and the imaginary
     * memory states `y[3..5]`. The pair factor of two is folded into the
     * stored residue matrices at construction.
     */
    template <typename scalar_type, typename index_type>
    struct VectorFitComplexSection
    {
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename GridKit::ScalarTraits<ScalarT>::RealT;

      RealT            a{0.0}; ///< \f$a_q = \mathrm{Re}\,p_q\f$
      RealT            w{0.0}; ///< \f$\omega_q = \mathrm{Im}\,p_q\f$
      ABCMatrix<RealT> A{};    ///< Scaled and doubled \f$\mathbf{A}_q\f$
      ABCMatrix<RealT> B{};    ///< Scaled and doubled \f$\mathbf{B}_q\f$

      __attribute__((always_inline)) inline int evaluateInternalResidual(
          const ScalarT*                  y,
          const ScalarT*                  yp,
          const ScalarT*                  y_ext,
          [[maybe_unused]] const ScalarT* yp_ext,
          ScalarT*                        f)
      {
        f[0] = -yp[0] + a * y[0] - w * y[3] + y_ext[0];
        f[1] = -yp[1] + a * y[1] - w * y[4] + y_ext[1];
        f[2] = -yp[2] + a * y[2] - w * y[5] + y_ext[2];
        f[3] = -yp[3] + w * y[0] + a * y[3];
        f[4] = -yp[4] + w * y[1] + a * y[4];
        f[5] = -yp[5] + w * y[2] + a * y[5];

        return 0;
      }

      __attribute__((always_inline)) inline int evaluateExternalResidual(
          const ScalarT*                  y,
          [[maybe_unused]] const ScalarT* yp,
          [[maybe_unused]] const ScalarT* y_ext,
          [[maybe_unused]] const ScalarT* yp_ext,
          ScalarT*                        f_ext)
      {
        f_ext[0] = A[0][0] * y[0] + A[0][1] * y[1] + A[0][2] * y[2]
                   - B[0][0] * y[3] - B[0][1] * y[4] - B[0][2] * y[5];
        f_ext[1] = A[1][0] * y[0] + A[1][1] * y[1] + A[1][2] * y[2]
                   - B[1][0] * y[3] - B[1][1] * y[4] - B[1][2] * y[5];
        f_ext[2] = A[2][0] * y[0] + A[2][1] * y[1] + A[2][2] * y[2]
                   - B[2][0] * y[3] - B[2][1] * y[4] - B[2][2] * y[5];

        return 0;
      }
    };

    /**
     * @brief Feedthrough section of a rational operator.
     *
     * The section owns no states and contributes the constant and linear
     * coefficient terms directly from the input triple and its derivative.
     */
    template <typename scalar_type, typename index_type>
    struct VectorFitFeedthroughSection
    {
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename GridKit::ScalarTraits<ScalarT>::RealT;

      ABCMatrix<RealT> D{}; ///< Scaled \f$\mathbf{D}\f$
      ABCMatrix<RealT> E{}; ///< Scaled \f$\mathbf{E}\f$

      __attribute__((always_inline)) inline int evaluateExternalResidual(
          [[maybe_unused]] const ScalarT* y,
          [[maybe_unused]] const ScalarT* yp,
          const ScalarT*                  y_ext,
          const ScalarT*                  yp_ext,
          ScalarT*                        f_ext)
      {
        f_ext[0] = D[0][0] * y_ext[0] + D[0][1] * y_ext[1] + D[0][2] * y_ext[2]
                   + E[0][0] * yp_ext[0] + E[0][1] * yp_ext[1] + E[0][2] * yp_ext[2];
        f_ext[1] = D[1][0] * y_ext[0] + D[1][1] * y_ext[1] + D[1][2] * y_ext[2]
                   + E[1][0] * yp_ext[0] + E[1][1] * yp_ext[1] + E[1][2] * yp_ext[2];
        f_ext[2] = D[2][0] * y_ext[0] + D[2][1] * y_ext[1] + D[2][2] * y_ext[2]
                   + E[2][0] * yp_ext[0] + E[2][1] * yp_ext[1] + E[2][2] * yp_ext[2];

        return 0;
      }
    };

    /*!
     * @brief Implementation of a three-phase vector-fitted rational operator.
     *
     * The operator is a parallel bank of fixed-size sections: one first-order
     * section per real pole, one second-order section per conjugate pair, and
     * one feedthrough. Every differentiated kernel is a straight-line section
     * kernel; the runtime pole count lives only in orchestration loops. The
     * output terms accumulate into the residual rows of the attached output
     * port, and their summation over sections happens through the residual
     * scatter and the system Jacobian deduplication.
     *
     * The model is consumed as a registered submodel of a component that
     * attaches the input and output ports, but presents the full `Component`
     * surface and can be assembled and tested standalone.
     */
    template <typename scalar_type, typename index_type>
    class VectorFit : public Component<scalar_type, index_type>
    {
      using Component<scalar_type, index_type>::gridkit_component_id_;
      using Component<scalar_type, index_type>::size_;
      using Component<scalar_type, index_type>::nnz_;
      using Component<scalar_type, index_type>::time_;
      using Component<scalar_type, index_type>::alpha_;
      using Component<scalar_type, index_type>::y_;
      using Component<scalar_type, index_type>::yp_;
      using Component<scalar_type, index_type>::abs_tol_;
      using Component<scalar_type, index_type>::tag_;
      using Component<scalar_type, index_type>::y_ext_;
      using Component<scalar_type, index_type>::yp_ext_;
      using Component<scalar_type, index_type>::variable_indices_ext_;
      using Component<scalar_type, index_type>::residual_indices_ext_;
      using Component<scalar_type, index_type>::f_ext_;
      using Component<scalar_type, index_type>::f_;
      using Component<scalar_type, index_type>::J_rows_buffer_;
      using Component<scalar_type, index_type>::J_cols_buffer_;
      using Component<scalar_type, index_type>::J_vals_buffer_;
      using Component<scalar_type, index_type>::variable_indices_;
      using Component<scalar_type, index_type>::residual_indices_;
      using Component<scalar_type, index_type>::allocated_;

    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename Component<ScalarT, IdxT>::RealT;
      using ModelDataT = VectorFitData<RealT, IdxT>;
      using SignalT    = Signal<ScalarT, IdxT>;
      using Port3T     = Port3<ScalarT, IdxT>;

      VectorFit(const ModelDataT& data, RealT scale);
      virtual ~VectorFit();

      virtual int setGridKitComponentID(IdxT) override final;
      virtual int allocate() override final;
      virtual int verify() const override final;
      virtual int initialize() override final;
      virtual int tagDifferentiable() override final;
      virtual int setAbsoluteTolerance(RealT) override final;
      virtual int evaluateInternalResidual() override final;
      virtual int evaluateResidual() override final;
      virtual int evaluateExternalResidual() override final;
      virtual int evaluateJacobian() override final;

      virtual void resetJacobianStructure() override final;

      /**
       * @brief Attach the input signals exposing the input triple.
       *
       * The signals supply the input values, their derivatives, and the global
       * variable indices.
       */
      void attachInput(SignalT* a, SignalT* b, SignalT* c)
      {
        input_ = {a, b, c};
      }

      void attachInput(Port3T* input)
      {
        attachInput(input->a(), input->b(), input->c());
      }

      /**
       * @brief Attach the output signals whose residual rows receive the
       * operator output terms.
       */
      void attachOutput(SignalT* a, SignalT* b, SignalT* c)
      {
        output_ = {a, b, c};
      }

      void attachOutput(Port3T* output)
      {
        attachOutput(output->a(), output->b(), output->c());
      }

      /**
       * @brief Whether the operator reads the input derivative.
       *
       * A nonzero linear coefficient makes the input differential and marks
       * derivative coupling on the input signals.
       */
      bool hasFeedthroughDerivative() const;

      /**
       * @brief Jacobian triplet capacity of this operator.
       *
       * A consuming component adds this into its own buffer sizing before
       * appending the operator Jacobian.
       */
      IdxT jacobianCapacity() const;

      /**
       * @brief Evaluate the operator output from the current states.
       *
       * Orchestration only, for monitors and diagnostics.
       */
      ScalarT output(IdxT n) const;

      /**
       * @brief The frequency response at one angular frequency, in real pairs.
       */
      void transfer(RealT omega0, ABCMatrix<RealT>& H_re, ABCMatrix<RealT>& H_im) const;

      /**
       * @brief Set the memory states and their derivatives to the sinusoidal
       * steady state driven by the instantaneous input pair at the
       * initialization instant.
       */
      int initializeSteadyState(RealT                   omega0,
                                const ABCVector<RealT>& u,
                                const ABCVector<RealT>& u_dot);

    private:
      /* Input parameters */
      IdxT pole_error_count_{0};

      /* Derivied parameters */
      std::vector<VectorFitRealSection<ScalarT, IdxT>>    real_sections_;
      std::vector<VectorFitComplexSection<ScalarT, IdxT>> complex_sections_;
      VectorFitFeedthroughSection<ScalarT, IdxT>          feedthrough_;

      std::array<SignalT*, 3> input_{};
      std::array<SignalT*, 3> output_{};

      /* Linear time-invariant Jacobian cache */
      bool               jacobian_cached_{false};
      IdxT               alpha_range_begin_{0};
      IdxT               alpha_range_end_{0};
      std::vector<RealT> unit_vals_;
    };

  } // namespace EMT
} // namespace GridKit
