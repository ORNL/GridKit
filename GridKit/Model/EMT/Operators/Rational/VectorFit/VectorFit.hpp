#pragma once

#include <GridKit/Model/EMT/Operators/Rational/Rational.hpp>
#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFitData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /** Matrix-residue realization with K real states per pole. */
    template <typename scalar_type, typename index_type>
    class VectorFit : public Rational<scalar_type, index_type>
    {
    public:
      using Base       = Rational<scalar_type, index_type>;
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename Base::RealT;
      using ModelDataT = VectorFitData<RealT, IdxT>;

      explicit VectorFit(const ModelDataT& data, RealT scale = RealT{1})
        : Base(data.rows > 0 ? static_cast<size_t>(data.rows) : 0,
               data.cols > 0 ? static_cast<size_t>(data.cols) : 0,
               data.D,
               data.E,
               scale,
               data.validate())
      {
        if (this->errors_)
        {
          return;
        }
        // Keep real sections first to preserve the established state ordering.
        for (bool pair : {false, true})
        {
          for (size_t q = 0; q < data.poles.size(); ++q)
          {
            const bool complex = data.poles[q].imag() != RealT{0};
            if (complex == pair)
            {
              auto& section = this->sections_.emplace_back(this->cols_, this->rows_, this->cols_);
              section.a     = data.poles[q].real();
              section.w     = data.poles[q].imag();
              section.pair  = pair;
              for (size_t k = 0; k < this->cols_; ++k)
              {
                section.Br[k][k] = RealT{1};
              }
              for (size_t n = 0; n < this->rows_; ++n)
              {
                for (size_t k = 0; k < this->cols_; ++k)
                {
                  const auto residue = scale * (pair ? RealT{2} : RealT{1}) * data.residues[q][n][k];
                  section.Cr[n][k]   = residue.real();
                  section.Ci[n][k]   = residue.imag();
                }
              }
            }
            if (complex)
            {
              ++q;
            }
          }
        }
        this->finishSections();
      }
    };
  } // namespace EMT
} // namespace GridKit
