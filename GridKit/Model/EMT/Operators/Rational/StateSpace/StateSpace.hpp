#pragma once

#include <GridKit/Model/EMT/Operators/Rational/Rational.hpp>
#include <GridKit/Model/EMT/Operators/Rational/StateSpace/StateSpaceData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /** Factorized rational realization with one real state per pole. */
    template <typename scalar_type, typename index_type>
    class StateSpace : public Rational<scalar_type, index_type>
    {
    public:
      using Base       = Rational<scalar_type, index_type>;
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename Base::RealT;
      using ModelDataT = StateSpaceData<RealT, IdxT>;

      explicit StateSpace(const ModelDataT& data, RealT scale = RealT{1})
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
        for (bool pair : {false, true})
        {
          for (size_t q = 0; q < data.poles.size(); ++q)
          {
            const bool complex = data.poles[q].imag() != RealT{0};
            if (complex == pair)
            {
              auto& section = this->sections_.emplace_back(1, this->rows_, this->cols_);
              section.a     = data.poles[q].real();
              section.w     = data.poles[q].imag();
              section.pair  = pair;
              for (size_t k = 0; k < this->cols_; ++k)
              {
                section.Br[0][k] = data.B[q][k].real();
                section.Bi[0][k] = data.B[q][k].imag();
              }
              for (size_t n = 0; n < this->rows_; ++n)
              {
                const auto factor = scale * (pair ? RealT{2} : RealT{1}) * data.C[n][q];
                section.Cr[n][0]  = factor.real();
                section.Ci[n][0]  = factor.imag();
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
