/**
 * @file VectorFit.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Declaration of an EMT VectorFit operator.
 */

#pragma once

#include <cstddef>
#include <string_view>
#include <vector>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFitData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Internal variables of a `VectorFit`
    enum class VectorFitInternalVariables : size_t
    {
      w, ///< \f$\mathbf{w}_q\f$ Real memory states [s [u]]
      v, ///< \f$\mathbf{v}_q\f$ Imaginary memory states [s [u]]
    };

    /**
     * @brief Vector-fitted rational operator
     *
     * The output dimension N and the input dimension K are those of the
     * first of `D`, `E` and the residues that the data gives. `w` has K
     * entries per real pole, then per conjugate pair, and `v` has K entries
     * per conjugate pair.
     */
    template <typename scalar_type, typename index_type>
    class VectorFit final : public Component<scalar_type, index_type>
    {
      using ComponentT = Component<scalar_type, index_type>;
      using IdxT       = index_type;
      using RealT      = typename ComponentT::RealT;

    public:
      using ModelDataT         = VectorFitData<RealT, IdxT>;
      using InternalVariablesT = VectorFitInternalVariables;
      using TerminalsT         = Utilities::EmptyEnum;
      using InputsT            = VectorFitInputs;
      using OutputsT           = VectorFitOutputs;

      static constexpr std::string_view CLASS_NAME = "VectorFit";

      explicit VectorFit(const ModelDataT& data);

      /// Output dimension N
      IdxT outputs() const
      {
        return N_;
      }

      /// Input dimension K
      IdxT inputs() const
      {
        return K_;
      }

      Model::ConfigurationChecks verify() const override;
      int                        setAbsoluteTolerance(RealT rel_tol) override;
      int                        tagDifferentiable() override;

      int evaluateOutputs() override;

      template <typename T>
      void internalResidual(const T* y, const T* yp, const T* u, const T* up, T* f) const;

      template <typename T>
      void output(OutputsT output, const T* y, const T* yp, const T* u, const T* up, T* f) const;

    private:
      std::vector<RealT> coefficients(const Matrix<RealT>& matrix) const;

      void setDerivedParameters(const ModelDataT& data);

      IdxT N_{0};
      IdxT K_{0};

      std::vector<RealT> D_; ///< N by K
      std::vector<RealT> E_; ///< N by K

      std::vector<RealT> ar_; ///< Real poles
      std::vector<RealT> Ar_; ///< Their residues, Qr by N by K
      std::vector<RealT> ac_; ///< Real part of the first pole of each conjugate pair
      std::vector<RealT> wc_; ///< Its imaginary part
      std::vector<RealT> Ac_; ///< Real part of its residues, Qc by N by K
      std::vector<RealT> Bc_; ///< Imaginary part of its residues, Qc by N by K
      bool               residues_shaped_{true};
      bool               real_residues_{true};
      bool               conjugate_pairs_{true};
    };
  } // namespace EMT
} // namespace GridKit
