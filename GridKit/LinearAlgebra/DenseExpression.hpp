/**
 * @file DenseExpression.hpp
 *
 * @brief Small lazy dense-expression helpers for readable EMT residuals.
 *
 */

#pragma once

#include <cassert>
#include <concepts>

namespace GridKit
{
  namespace LinearAlgebra
  {
    template <typename T>
    concept DenseLike = requires(T a, typename T::IdxT i) {
      typename T::ValueT;
      typename T::IdxT;
      a.rows;
      a.cols;
      a(i, i);
    };

    template <typename T>
    concept SliceLike = DenseLike<T> && requires(T a) {
      a.value;
      a.base_offset;
    };

    template <SliceLike SliceT>
    __attribute__((always_inline)) typename SliceT::ValueT eval(SliceT                s,
                                                                typename SliceT::IdxT i,
                                                                typename SliceT::IdxT j = 0)
    {
      return s.value[s.base_offset + i * s.cols + j];
    }

    template <DenseLike ExprT>
      requires(!SliceLike<ExprT>)
    __attribute__((always_inline)) typename ExprT::ValueT eval(ExprT                expr,
                                                               typename ExprT::IdxT i,
                                                               typename ExprT::IdxT j = 0)
    {
      return expr(i, j);
    }

    template <SliceLike SliceT>
    __attribute__((always_inline)) void store(SliceT                  out,
                                              typename SliceT::IdxT   i,
                                              typename SliceT::IdxT   j,
                                              typename SliceT::ValueT value)
    {
      out.value[out.base_offset + i * out.cols + j] = value;
    }

    template <DenseLike A, DenseLike B>
    struct Add
    {
      using ValueT = typename A::ValueT;
      using IdxT   = typename A::IdxT;

      A    a;
      B    b;
      IdxT rows;
      IdxT cols;

      __attribute__((always_inline)) ValueT operator()(IdxT i, IdxT j = 0) const
      {
        return eval(a, i, j) + eval(b, i, j);
      }
    };

    template <DenseLike A, DenseLike B>
    struct Subtract
    {
      using ValueT = typename A::ValueT;
      using IdxT   = typename A::IdxT;

      A    a;
      B    b;
      IdxT rows;
      IdxT cols;

      __attribute__((always_inline)) ValueT operator()(IdxT i, IdxT j = 0) const
      {
        return eval(a, i, j) - eval(b, i, j);
      }
    };

    template <DenseLike A>
    struct Negate
    {
      using ValueT = typename A::ValueT;
      using IdxT   = typename A::IdxT;

      A    a;
      IdxT rows;
      IdxT cols;

      __attribute__((always_inline)) ValueT operator()(IdxT i, IdxT j = 0) const
      {
        return -eval(a, i, j);
      }
    };

    template <DenseLike A, DenseLike B>
    struct Product
    {
      using ValueT = typename A::ValueT;
      using IdxT   = typename A::IdxT;

      A    a;
      B    b;
      IdxT rows;
      IdxT cols;
      IdxT inner;

      __attribute__((always_inline)) ValueT operator()(IdxT i, IdxT j = 0) const
      {
        ValueT value = 0.0;
        for (IdxT k = 0; k < inner; ++k)
        {
          value += eval(a, i, k) * eval(b, k, j);
        }
        return value;
      }
    };

    template <typename S, DenseLike A>
    struct Scale
    {
      using ValueT = typename A::ValueT;
      using IdxT   = typename A::IdxT;

      S    s;
      A    a;
      IdxT rows;
      IdxT cols;

      __attribute__((always_inline)) ValueT operator()(IdxT i, IdxT j = 0) const
      {
        return s * eval(a, i, j);
      }
    };

    template <DenseLike A, typename S>
    struct Divide
    {
      using ValueT = typename A::ValueT;
      using IdxT   = typename A::IdxT;

      A    a;
      S    s;
      IdxT rows;
      IdxT cols;

      __attribute__((always_inline)) ValueT operator()(IdxT i, IdxT j = 0) const
      {
        return eval(a, i, j) / s;
      }
    };

    template <typename value_type, typename index_type, typename Fn>
    struct Matrix
    {
      using ValueT = value_type;
      using IdxT   = index_type;

      IdxT rows;
      IdxT cols;
      Fn   fn;

      __attribute__((always_inline)) ValueT operator()(IdxT i, IdxT j = 0) const
      {
        return fn(i, j);
      }
    };

    template <typename value_type, typename index_type>
    struct Identity
    {
      using ValueT = value_type;
      using IdxT   = index_type;

      IdxT rows;
      IdxT cols;

      __attribute__((always_inline)) ValueT operator()(IdxT i, IdxT j = 0) const
      {
        return static_cast<ValueT>(i == j);
      }
    };

    template <typename value_type, typename index_type>
    struct Constant
    {
      using ValueT = value_type;
      using IdxT   = index_type;

      IdxT   rows;
      IdxT   cols;
      ValueT value;

      __attribute__((always_inline)) ValueT operator()(IdxT, IdxT = 0) const
      {
        return value;
      }
    };

    template <DenseLike V>
    struct Diagonal
    {
      using ValueT = typename V::ValueT;
      using IdxT   = typename V::IdxT;

      V    v;
      IdxT rows;
      IdxT cols;

      __attribute__((always_inline)) ValueT operator()(IdxT i, IdxT j = 0) const
      {
        return (i == j) ? eval(v, i) : 0.0;
      }
    };

    template <DenseLike A, DenseLike B>
    struct Hadamard
    {
      using ValueT = typename A::ValueT;
      using IdxT   = typename A::IdxT;

      A    a;
      B    b;
      IdxT rows;
      IdxT cols;

      __attribute__((always_inline)) ValueT operator()(IdxT i, IdxT j = 0) const
      {
        return eval(a, i, j) * eval(b, i, j);
      }
    };

    template <DenseLike A, DenseLike B>
    __attribute__((always_inline)) auto operator+(A a, B b)
    {
      assert(a.rows == b.rows);
      assert(a.cols == b.cols);
      return Add<A, B>{a, b, a.rows, a.cols};
    }

    template <DenseLike A, DenseLike B>
    __attribute__((always_inline)) auto operator-(A a, B b)
    {
      assert(a.rows == b.rows);
      assert(a.cols == b.cols);
      return Subtract<A, B>{a, b, a.rows, a.cols};
    }

    template <DenseLike A>
    __attribute__((always_inline)) auto operator-(A a)
    {
      return Negate<A>{a, a.rows, a.cols};
    }

    template <DenseLike A, DenseLike B>
    __attribute__((always_inline)) auto operator*(A a, B b)
    {
      assert(a.cols == b.rows);
      return Product<A, B>{a, b, a.rows, b.cols, a.cols};
    }

    template <typename S, DenseLike A>
      requires(!DenseLike<S>)
    __attribute__((always_inline)) auto operator*(S s, A a)
    {
      return Scale<S, A>{s, a, a.rows, a.cols};
    }

    template <DenseLike A, typename S>
      requires(!DenseLike<S>)
    __attribute__((always_inline)) auto operator*(A a, S s)
    {
      return Scale<S, A>{s, a, a.rows, a.cols};
    }

    template <DenseLike A, typename S>
      requires(!DenseLike<S>)
    __attribute__((always_inline)) auto operator/(A a, S s)
    {
      return Divide<A, S>{a, s, a.rows, a.cols};
    }

    template <typename ValueT, typename IdxT, typename Fn>
    __attribute__((always_inline)) auto matrix(IdxT rows, IdxT cols, Fn fn)
    {
      return Matrix<ValueT, IdxT, Fn>{rows, cols, fn};
    }

    template <typename IdxT, typename Fn>
    __attribute__((always_inline)) auto matrix(IdxT rows, IdxT cols, Fn fn)
    {
      using ValueT = decltype(fn(IdxT{}, IdxT{}));
      return Matrix<ValueT, IdxT, Fn>{rows, cols, fn};
    }

    template <typename ValueT, typename IdxT, typename Fn>
    __attribute__((always_inline)) auto vector(IdxT rows, Fn fn)
    {
      auto matrix_fn = [fn](IdxT i, IdxT)
      {
        return fn(i);
      };
      return Matrix<ValueT, IdxT, decltype(matrix_fn)>{rows, 1, matrix_fn};
    }

    template <typename IdxT, typename Fn>
    __attribute__((always_inline)) auto vector(IdxT rows, Fn fn)
    {
      auto matrix_fn = [fn](IdxT i, IdxT)
      {
        return fn(i);
      };
      using ValueT = decltype(fn(IdxT{}));
      return Matrix<ValueT, IdxT, decltype(matrix_fn)>{rows, 1, matrix_fn};
    }

    template <typename ValueT, typename IdxT>
    __attribute__((always_inline)) auto identity(IdxT n)
    {
      return Identity<ValueT, IdxT>{n, n};
    }

    template <DenseLike A>
    __attribute__((always_inline)) auto identity(A a)
    {
      assert(a.rows == a.cols);
      return Identity<typename A::ValueT, typename A::IdxT>{a.rows, a.cols};
    }

    template <DenseLike A>
    __attribute__((always_inline)) auto constant(A a, typename A::ValueT value)
    {
      return Constant<typename A::ValueT, typename A::IdxT>{a.rows, a.cols, value};
    }

    template <DenseLike V>
    __attribute__((always_inline)) auto diag(V v)
    {
      return Diagonal<V>{v, v.rows, v.rows};
    }

    template <DenseLike A, DenseLike B>
    __attribute__((always_inline)) auto hadamard(A a, B b)
    {
      assert(a.rows == b.rows);
      assert(a.cols == b.cols);
      return Hadamard<A, B>{a, b, a.rows, a.cols};
    }

    template <SliceLike Out, DenseLike ExprT>
    __attribute__((always_inline)) void lower(Out out, ExprT expr)
    {
      assert(out.rows == expr.rows);
      assert(out.cols == expr.cols);
      for (typename Out::IdxT i = 0; i < out.rows; ++i)
      {
        for (typename Out::IdxT j = 0; j < out.cols; ++j)
        {
          store(out, i, j, eval(expr, i, j));
        }
      }
    }

    template <SliceLike Out>
    struct EquationProxy
    {
      Out out;

      template <DenseLike ExprT>
      __attribute__((always_inline)) void operator=(ExprT expr)
      {
        lower(out, expr);
      }
    };

    template <SliceLike Out>
    __attribute__((always_inline)) auto equation(Out out)
    {
      return EquationProxy<Out>{out};
    }

    template <SliceLike SliceT>
    __attribute__((always_inline)) void zero(SliceT out)
    {
      for (typename SliceT::IdxT i = 0; i < out.rows; ++i)
      {
        for (typename SliceT::IdxT j = 0; j < out.cols; ++j)
        {
          store(out, i, j, 0.0);
        }
      }
    }
  } // namespace LinearAlgebra
} // namespace GridKit
