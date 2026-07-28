#pragma once

#include <algorithm>
#include <vector>

#include <GridKit/Model/PhasorDynamics/Component.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /**
     * @brief The constant complex bus admittance matrix of the network.
     *
     * The linear part of the phasor network is time invariant: branches and
     * constant-impedance loads contribute a fixed admittance to the bus current
     * balance. Assembling that once as a sparse matrix lets the residual
     * evaluate it as a single sparse matrix-vector product instead of walking
     * one object per network element.
     *
     * Rows are the buses that own a residual slice, in the order they were bound
     * to system storage. Columns are stored as global system indices of the real
     * voltage component, so evaluation needs no index translation.
     *
     * Entries are complex, held as interleaved (g, b) pairs. Each expands to the
     * real 2x2 block [[g, -b], [b, g]] acting on the interleaved (Vr, Vi) bus
     * slice, which is why one complex entry covers four real matrix entries.
     */
    template <typename scalar_type, typename index_type>
    class NetworkAdmittance
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename Model::Evaluator<ScalarT, IdxT>::RealT;
      using StampT  = AdmittanceStamp<RealT, IdxT>;

      /**
       * @brief Assemble CSR storage from unordered stamps.
       *
       * Stamps that land on the same (row, column) are summed here, once, rather
       * than at every residual evaluation. Rows carrying no stamp are kept as
       * empty rows so that evaluation still writes their zero.
       *
       * @param[in,out] stamps     - Stamps to assemble; reordered in place.
       * @param[in]     row_offset - Global residual index of each row bus's Ir
       *                             entry, ascending.
       *
       * @pre Every stamp row matches some entry of `row_offset`.
       */
      void assemble(std::vector<StampT>& stamps, std::vector<IdxT> row_offset)
      {
        row_offset_ = std::move(row_offset);

        const std::size_t rows = row_offset_.size();

        row_ptr_.assign(rows + 1, IdxT{0});
        col_.clear();
        gb_.clear();

        if (stamps.empty())
        {
          return;
        }

        // Translate each stamp's global row index to its dense row position.
        // row_offset_ is ascending, so a binary search suffices.
        for (auto& stamp : stamps)
        {
          const auto it = std::lower_bound(row_offset_.begin(), row_offset_.end(), stamp.row);
          stamp.row     = static_cast<IdxT>(it - row_offset_.begin());
        }

        std::sort(stamps.begin(),
                  stamps.end(),
                  [](const StampT& a, const StampT& b)
                  {
                    return (a.row != b.row) ? (a.row < b.row) : (a.col < b.col);
                  });

        col_.reserve(stamps.size());
        gb_.reserve(2 * stamps.size());

        for (std::size_t i = 0; i < stamps.size();)
        {
          std::size_t j = i + 1;
          RealT       g = stamps[i].g;
          RealT       b = stamps[i].b;
          while (j < stamps.size() && stamps[j].row == stamps[i].row && stamps[j].col == stamps[i].col)
          {
            g += stamps[j].g;
            b += stamps[j].b;
            ++j;
          }

          col_.push_back(stamps[i].col);
          gb_.push_back(g);
          gb_.push_back(b);
          ++row_ptr_[static_cast<std::size_t>(stamps[i].row) + 1];

          i = j;
        }

        for (std::size_t r = 0; r < rows; ++r)
        {
          row_ptr_[r + 1] += row_ptr_[r];
        }
      }

      /**
       * @brief Write the network current injection of every row bus.
       *
       * Overwrites, so it also establishes the zero that bus rows without any
       * network incidence need, and it leaves every bus row ready for the
       * nonlinear components to accumulate onto.
       *
       * @param[in]  y - System state vector data.
       * @param[out] f - System residual vector data.
       */
      void multiply(const ScalarT* y, ScalarT* f) const
      {
        const std::size_t rows = row_offset_.size();

        for (std::size_t r = 0; r < rows; ++r)
        {
          ScalarT ir{0.0};
          ScalarT ii{0.0};

          const auto end = static_cast<std::size_t>(row_ptr_[r + 1]);
          for (auto k = static_cast<std::size_t>(row_ptr_[r]); k < end; ++k)
          {
            const auto    c  = static_cast<std::size_t>(col_[k]);
            const ScalarT vr = y[c];
            const ScalarT vi = y[c + 1];
            const RealT   g  = gb_[2 * k];
            const RealT   b  = gb_[2 * k + 1];

            ir += g * vr - b * vi;
            ii += b * vr + g * vi;
          }

          const auto row = static_cast<std::size_t>(row_offset_[r]);
          f[row]         = ir;
          f[row + 1]     = ii;
        }
      }

      bool empty() const
      {
        return row_offset_.empty();
      }

      IdxT rowCount() const
      {
        return static_cast<IdxT>(row_offset_.size());
      }

      IdxT nnz() const
      {
        return static_cast<IdxT>(col_.size());
      }

    private:
      /// Global residual index of each row bus's Ir entry
      std::vector<IdxT>  row_offset_;
      std::vector<IdxT>  row_ptr_;
      /// Global variable index of each entry's column bus Vr entry
      std::vector<IdxT>  col_;
      /// Interleaved (g, b) per entry
      std::vector<RealT> gb_;
    };

  } // namespace PhasorDynamics
} // namespace GridKit
