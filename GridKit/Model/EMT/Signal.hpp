/**
 * @file Signal.hpp
 *
 * @brief Shared EMT wiring contract. A Signal is durable metadata naming a
 * shaped block of the global state or state derivative. A Slice is a non-owning
 * shaped pointer view into bound state, residual, or signal arrays.
 *
 */

#pragma once

namespace GridKit
{
  namespace EMT
  {
    template <typename value_type, typename index_type>
    struct Slice
    {
      using ValueT = value_type;
      using IdxT   = index_type;

      ValueT*     value{nullptr};
      const IdxT* index{nullptr};
      IdxT        rows{0};
      IdxT        cols{1};
      IdxT        base_offset{0};

      ValueT& operator()(IdxT r, IdxT c = 0)
      {
        return value[offset(r, c)];
      }

      const ValueT& operator()(IdxT r, IdxT c = 0) const
      {
        return value[offset(r, c)];
      }

      ValueT& operator[](IdxT i)
      {
        return value[base_offset + i];
      }

      const ValueT& operator[](IdxT i) const
      {
        return value[base_offset + i];
      }

      IdxT global(IdxT r, IdxT c = 0) const
      {
        return index[offset(r, c)];
      }

      IdxT column(IdxT r, IdxT c = 0) const
      {
        return global(r, c);
      }

      IdxT size() const
      {
        return rows * cols;
      }

    private:
      IdxT offset(IdxT r, IdxT c = 0) const
      {
        return base_offset + r * cols + c;
      }
    };

    template <typename index_type>
    struct Signal
    {
      using IdxT = index_type;

      enum class Storage
      {
        State,
        Derivative
      };

      IdxT    offset{0};
      IdxT    rows{0};
      IdxT    cols{1};
      Storage storage{Storage::State};

      IdxT size() const
      {
        return rows * cols;
      }
    };

    template <typename scalar_type, typename index_type>
    using SignalView = Slice<scalar_type, index_type>;
  } // namespace EMT
} // namespace GridKit
