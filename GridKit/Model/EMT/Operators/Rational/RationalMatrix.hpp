#pragma once

#include <vector>

#include <GridKit/Model/EMT/ComponentData.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Rectangular coefficient storage; three-phase data retain their usual syntax.
    template <typename T>
    class RationalMatrix
    {
    public:
      RationalMatrix()
        : RationalMatrix(3, 3)
      {
      }

      RationalMatrix(size_t rows, size_t cols)
        : values_(rows, std::vector<T>(cols))
      {
      }

      RationalMatrix(const ABCMatrix<T>& values)
      {
        *this = values;
      }

      RationalMatrix& operator=(const ABCMatrix<T>& values)
      {
        values_.resize(3);
        for (size_t n = 0; n < 3; ++n)
        {
          values_[n].assign(values[n].begin(), values[n].end());
        }
        return *this;
      }

      auto& operator[](size_t n)
      {
        return values_[n];
      }

      const auto& operator[](size_t n) const
      {
        return values_[n];
      }

      size_t size() const
      {
        return values_.size();
      }

      auto begin()
      {
        return values_.begin();
      }

      auto end()
      {
        return values_.end();
      }

      auto begin() const
      {
        return values_.begin();
      }

      auto end() const
      {
        return values_.end();
      }

      bool hasShape(size_t rows, size_t cols) const
      {
        if (values_.size() != rows)
        {
          return false;
        }
        for (const auto& row : values_)
        {
          if (row.size() != cols)
          {
            return false;
          }
        }
        return true;
      }

    private:
      std::vector<std::vector<T>> values_;
    };
  } // namespace EMT
} // namespace GridKit
