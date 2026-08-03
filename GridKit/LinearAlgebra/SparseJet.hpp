/**
 * @file SparseJet.hpp
 *
 * @brief Small sparse forward-mode scalar used by EMT equation residuals.
 *
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>
#include <utility>
#include <vector>

namespace GridKit
{
  namespace LinearAlgebra
  {
    template <typename RealT, typename IdxT>
    struct SparseJet
    {
      using ValueT = RealT;
      using IdxT_  = IdxT;

      RealT                               value{0.0};
      std::vector<std::pair<IdxT, RealT>> tangent;

      SparseJet() = default;

      SparseJet(RealT v)
        : value(v)
      {
      }

      static SparseJet variable(RealT value, IdxT col, RealT scale = RealT{1.0})
      {
        SparseJet jet(value);
        jet.tangent.emplace_back(col, scale);
        return jet;
      }
    };

    template <typename T>
    struct IsSparseJet
    {
      static constexpr bool value = false;
    };

    template <typename RealT, typename IdxT>
    struct IsSparseJet<SparseJet<RealT, IdxT>>
    {
      static constexpr bool value = true;
    };

    template <typename T>
    inline constexpr bool IsSparseJetV = IsSparseJet<T>::value;

    template <typename RealT, typename IdxT>
    void addTangent(SparseJet<RealT, IdxT>& out, IdxT col, RealT value)
    {
      auto found = std::ranges::find_if(out.tangent,
                                        [col](const auto& term)
                                        {
                                          return term.first == col;
                                        });
      if (found != out.tangent.end())
      {
        found->second += value;
      }
      else
      {
        out.tangent.emplace_back(col, value);
      }
    }

    template <typename RealT, typename IdxT>
    void addScaledTangent(SparseJet<RealT, IdxT>&       out,
                          const SparseJet<RealT, IdxT>& in,
                          RealT                         scale)
    {
      for (const auto& [col, value] : in.tangent)
      {
        addTangent(out, col, scale * value);
      }
    }

    template <typename RealT, typename IdxT>
    SparseJet<RealT, IdxT> operator+(const SparseJet<RealT, IdxT>& a,
                                     const SparseJet<RealT, IdxT>& b)
    {
      SparseJet<RealT, IdxT> out(a.value + b.value);
      addScaledTangent(out, a, RealT{1.0});
      addScaledTangent(out, b, RealT{1.0});
      return out;
    }

    template <typename RealT, typename IdxT, typename S>
      requires(std::is_arithmetic_v<S>)
    SparseJet<RealT, IdxT> operator+(const SparseJet<RealT, IdxT>& a, S b)
    {
      return a + SparseJet<RealT, IdxT>(static_cast<RealT>(b));
    }

    template <typename S, typename RealT, typename IdxT>
      requires(std::is_arithmetic_v<S>)
    SparseJet<RealT, IdxT> operator+(S a, const SparseJet<RealT, IdxT>& b)
    {
      return SparseJet<RealT, IdxT>(static_cast<RealT>(a)) + b;
    }

    template <typename RealT, typename IdxT>
    SparseJet<RealT, IdxT> operator-(const SparseJet<RealT, IdxT>& a,
                                     const SparseJet<RealT, IdxT>& b)
    {
      SparseJet<RealT, IdxT> out(a.value - b.value);
      addScaledTangent(out, a, RealT{1.0});
      addScaledTangent(out, b, RealT{-1.0});
      return out;
    }

    template <typename RealT, typename IdxT, typename S>
      requires(std::is_arithmetic_v<S>)
    SparseJet<RealT, IdxT> operator-(const SparseJet<RealT, IdxT>& a, S b)
    {
      return a - SparseJet<RealT, IdxT>(static_cast<RealT>(b));
    }

    template <typename S, typename RealT, typename IdxT>
      requires(std::is_arithmetic_v<S>)
    SparseJet<RealT, IdxT> operator-(S a, const SparseJet<RealT, IdxT>& b)
    {
      return SparseJet<RealT, IdxT>(static_cast<RealT>(a)) - b;
    }

    template <typename RealT, typename IdxT>
    SparseJet<RealT, IdxT> operator-(const SparseJet<RealT, IdxT>& a)
    {
      SparseJet<RealT, IdxT> out(-a.value);
      addScaledTangent(out, a, RealT{-1.0});
      return out;
    }

    template <typename RealT, typename IdxT, typename S>
      requires(std::is_arithmetic_v<S>)
    SparseJet<RealT, IdxT> operator*(const SparseJet<RealT, IdxT>& a, S b)
    {
      return a * SparseJet<RealT, IdxT>(static_cast<RealT>(b));
    }

    template <typename S, typename RealT, typename IdxT>
      requires(std::is_arithmetic_v<S>)
    SparseJet<RealT, IdxT> operator*(S a, const SparseJet<RealT, IdxT>& b)
    {
      return SparseJet<RealT, IdxT>(static_cast<RealT>(a)) * b;
    }

    template <typename RealT, typename IdxT>
    SparseJet<RealT, IdxT> operator*(const SparseJet<RealT, IdxT>& a,
                                     const SparseJet<RealT, IdxT>& b)
    {
      SparseJet<RealT, IdxT> out(a.value * b.value);
      addScaledTangent(out, a, b.value);
      addScaledTangent(out, b, a.value);
      return out;
    }

    template <typename RealT, typename IdxT, typename S>
      requires(std::is_arithmetic_v<S>)
    SparseJet<RealT, IdxT> operator/(const SparseJet<RealT, IdxT>& a, S b)
    {
      return a / SparseJet<RealT, IdxT>(static_cast<RealT>(b));
    }

    template <typename S, typename RealT, typename IdxT>
      requires(std::is_arithmetic_v<S>)
    SparseJet<RealT, IdxT> operator/(S a, const SparseJet<RealT, IdxT>& b)
    {
      return SparseJet<RealT, IdxT>(static_cast<RealT>(a)) / b;
    }

    template <typename RealT, typename IdxT>
    SparseJet<RealT, IdxT>& operator+=(SparseJet<RealT, IdxT>&       a,
                                       const SparseJet<RealT, IdxT>& b)
    {
      a = a + b;
      return a;
    }

    template <typename RealT, typename IdxT>
    SparseJet<RealT, IdxT>& operator-=(SparseJet<RealT, IdxT>&       a,
                                       const SparseJet<RealT, IdxT>& b)
    {
      a = a - b;
      return a;
    }

    template <typename RealT, typename IdxT>
    SparseJet<RealT, IdxT>& operator*=(SparseJet<RealT, IdxT>&       a,
                                       const SparseJet<RealT, IdxT>& b)
    {
      a = a * b;
      return a;
    }

    template <typename RealT, typename IdxT>
    SparseJet<RealT, IdxT>& operator/=(SparseJet<RealT, IdxT>&       a,
                                       const SparseJet<RealT, IdxT>& b)
    {
      a = a / b;
      return a;
    }

    template <typename RealT, typename IdxT>
    SparseJet<RealT, IdxT> operator/(const SparseJet<RealT, IdxT>& a,
                                     const SparseJet<RealT, IdxT>& b)
    {
      SparseJet<RealT, IdxT> out(a.value / b.value);
      addScaledTangent(out, a, RealT{1.0} / b.value);
      addScaledTangent(out, b, -a.value / (b.value * b.value));
      return out;
    }

    template <typename RealT, typename IdxT>
    SparseJet<RealT, IdxT> log(const SparseJet<RealT, IdxT>& a)
    {
      SparseJet<RealT, IdxT> out(std::log(a.value));
      addScaledTangent(out, a, RealT{1.0} / a.value);
      return out;
    }

    template <typename RealT, typename IdxT>
    SparseJet<RealT, IdxT> sqrt(const SparseJet<RealT, IdxT>& a)
    {
      const RealT            root = std::sqrt(a.value);
      SparseJet<RealT, IdxT> out(root);
      addScaledTangent(out, a, RealT{0.5} / root);
      return out;
    }

    template <typename RealT, typename IdxT>
    SparseJet<RealT, IdxT> exp(const SparseJet<RealT, IdxT>& a)
    {
      const RealT            value = std::exp(a.value);
      SparseJet<RealT, IdxT> out(value);
      addScaledTangent(out, a, value);
      return out;
    }

    template <typename RealT, typename IdxT>
    SparseJet<RealT, IdxT> sin(const SparseJet<RealT, IdxT>& a)
    {
      SparseJet<RealT, IdxT> out(std::sin(a.value));
      addScaledTangent(out, a, std::cos(a.value));
      return out;
    }

    template <typename RealT, typename IdxT>
    SparseJet<RealT, IdxT> cos(const SparseJet<RealT, IdxT>& a)
    {
      SparseJet<RealT, IdxT> out(std::cos(a.value));
      addScaledTangent(out, a, -std::sin(a.value));
      return out;
    }

    template <typename RealT, typename IdxT>
    SparseJet<RealT, IdxT> atan2(const SparseJet<RealT, IdxT>& y,
                                 const SparseJet<RealT, IdxT>& x)
    {
      const RealT            denom = x.value * x.value + y.value * y.value;
      SparseJet<RealT, IdxT> out(std::atan2(y.value, x.value));
      addScaledTangent(out, y, x.value / denom);
      addScaledTangent(out, x, -y.value / denom);
      return out;
    }
  } // namespace LinearAlgebra
} // namespace GridKit
