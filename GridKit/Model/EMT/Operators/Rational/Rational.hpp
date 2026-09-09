/**
 * @file Rational.hpp
 * @brief Declaration of the EMT Rational model.
 */
#pragma once

#include <complex>
#include <span>
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
    protected:
      using Component<scalar_type, index_type>::J_cols_buffer_;
      using Component<scalar_type, index_type>::J_rows_buffer_;
      using Component<scalar_type, index_type>::J_vals_buffer_;
      using Component<scalar_type, index_type>::abs_tol_;
      using Component<scalar_type, index_type>::allocated_;
      using Component<scalar_type, index_type>::f_;
      using Component<scalar_type, index_type>::f_ext_;
      using Component<scalar_type, index_type>::gridkit_component_id_;
      using Component<scalar_type, index_type>::nnz_;
      using Component<scalar_type, index_type>::residual_indices_;
      using Component<scalar_type, index_type>::residual_indices_ext_;
      using Component<scalar_type, index_type>::size_;
      using Component<scalar_type, index_type>::tag_;
      using Component<scalar_type, index_type>::variable_indices_;
      using Component<scalar_type, index_type>::variable_indices_ext_;
      using Component<scalar_type, index_type>::y_;
      using Component<scalar_type, index_type>::y_ext_;
      using Component<scalar_type, index_type>::yp_;
      using Component<scalar_type, index_type>::yp_ext_;

    public:
      using ScalarT  = scalar_type;
      using IdxT     = index_type;
      using RealT    = typename Component<ScalarT, IdxT>::RealT;
      using SignalT  = Signal<ScalarT, IdxT>;
      using MatrixT  = RationalMatrix<RealT>;
      using ComplexT = std::complex<RealT>;

      IdxT rows() const;
      IdxT cols() const;

      void attachInput(const std::vector<SignalT*>& input);
      void attachOutput(const std::vector<SignalT*>& output);
      void attachInput(SignalT* a, SignalT* b, SignalT* c);
      void attachOutput(SignalT* a, SignalT* b, SignalT* c);

      bool hasInputDerivative(IdxT k) const;
      bool hasFeedthroughDerivative() const;

      int setGridKitComponentID(IdxT id) override;
      int allocate() override;
      int verify() const override;

      int initialize();
      int initializeState(const std::map<std::string, RealT>& values) override;
      int initializeSteadyState(RealT omega);
      int initializeSteadyState(RealT omega, std::span<const RealT> u, std::span<const RealT> udot);

      typename Component<ScalarT, IdxT>::InitializationPortsT initializationPorts() override;

      void transfer(RealT omega, MatrixT& re, MatrixT& im) const;
      void transfer(RealT omega, ABCMatrix<RealT>& re, ABCMatrix<RealT>& im) const;

      ScalarT output(IdxT n) const;
      /// Time derivative of a proper output at an accepted solver state.
      ScalarT outputDerivative(IdxT n) const;
      void    appendOutputGradient(IdxT n, typename SignalT::GradientT& gradient, RealT scale) const;

      int setAbsoluteTolerance(RealT tolerance) override;
      int evaluateInternalResidual() override;
      int evaluateExternalResidual() override;
      int evaluateResidual() override;
      int assembleJacobian(RealT y_scale, RealT yp_scale) override;

      IdxT jacobianCapacity() const;

      __attribute__((always_inline)) inline ScalarT evaluateOutput(
          IdxT n, const ScalarT* y, const ScalarT* y_ext, const ScalarT* yp_ext) const;

      __attribute__((always_inline)) inline int evaluateInternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

      __attribute__((always_inline)) inline int evaluateExternalResidual(
          const ScalarT*, const ScalarT*, const ScalarT*, const ScalarT*, ScalarT*);

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

      Rational(size_t rows, size_t cols, const MatrixT& D, const MatrixT& E, RealT scale, int errors);
      void finishSections();

      size_t                rows_, cols_;
      MatrixT               D_, E_;
      std::vector<SignalT*> input_, output_;
      std::vector<Section>  sections_;
      int                   errors_{};
      size_t                capacity_{};
      bool                  coupling_allocated_{false};

    private:
      /// Local output partials depend only on the fixed realization coefficients.
      mutable std::vector<std::vector<RealT>> output_partials_;
      mutable std::vector<ScalarT>            output_direction_;
    };
  } // namespace EMT
} // namespace GridKit
