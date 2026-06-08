#pragma once

#include <cassert>
#include <cstddef>

namespace GridKit
{
  namespace LinearAlgebra
  {
    /**
     * @brief Non-owning view of contiguous vector data.
     *
     * VectorView does not allocate, free, or synchronize memory. It only
     * provides indexed access to storage owned elsewhere.
     */
    template <typename ScalarT, typename IdxT>
    class VectorView
    {
    public:
      VectorView() = default;

      VectorView(ScalarT* data, std::size_t size)
        : data_(data),
          size_(size)
      {
        assert(data_ != nullptr || size_ == 0);
      }

      void bind(ScalarT* data, std::size_t size)
      {
        assert(data != nullptr || size == 0);
        data_ = data;
        size_ = size;
      }

      ScalarT* data()
      {
        return data_;
      }

      const ScalarT* data() const
      {
        return data_;
      }

      std::size_t size() const
      {
        return size_;
      }

      IdxT getSize() const
      {
        return static_cast<IdxT>(size_);
      }

      bool empty() const
      {
        return size_ == 0;
      }

      ScalarT& operator[](std::size_t i)
      {
        assert(i < size_);
        return data_[i];
      }

      const ScalarT& operator[](std::size_t i) const
      {
        assert(i < size_);
        return data_[i];
      }

    private:
      ScalarT*    data_{nullptr};
      std::size_t size_{0};
    };
  } // namespace LinearAlgebra
} // namespace GridKit
