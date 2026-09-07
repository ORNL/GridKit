#pragma once

#include <span>

#include <GridKit/Model/EMT/Component.hpp>
#include <GridKit/Model/EMT/Operators/Rational/RationalMatrix.hpp>

namespace GridKit
{
  namespace EMT
  {
    /// Common port and frequency-response contract for transport operators.
    template <typename scalar_type, typename index_type>
    class Shift : public Component<scalar_type, index_type>
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using Base    = Component<ScalarT, IdxT>;
      using RealT   = typename Base::RealT;
      using SignalT = typename Base::SignalT;
      using MatrixT = RationalMatrix<RealT>;

      virtual void     attachInput(const std::vector<SignalT*>& input)                                                                  = 0;
      virtual SignalT& outputSignal(IdxT channel)                                                                                       = 0;
      virtual void     transfer(RealT omega, MatrixT& re, MatrixT& im) const                                                            = 0;
      virtual int      initializeSteadyState(RealT omega, std::span<const RealT> u, std::span<const RealT> udot, RealT time = RealT{0}) = 0;

      int setGridKitComponentID(IdxT id) override
      {
        this->gridkit_component_id_ = id;
        return 0;
      }
    };
  } // namespace EMT
} // namespace GridKit
