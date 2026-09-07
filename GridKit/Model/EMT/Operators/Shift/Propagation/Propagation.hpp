/**
 * @file Propagation.hpp
 * @brief Sum of independently fitted and delayed matrix modes.
 */
#pragma once

#include <memory>

#include <GridKit/Model/EMT/Operators/Rational/VectorFit/VectorFit.hpp>
#include <GridKit/Model/EMT/Operators/Shift/Delay/Delay.hpp>
#include <GridKit/Model/EMT/Operators/Shift/Propagation/PropagationData.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    class Propagation : public Shift<scalar_type, index_type>
    {
    public:
      using Base       = Shift<scalar_type, index_type>;
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename Base::RealT;
      using SignalT    = typename Base::SignalT;
      using MatrixT    = typename Base::MatrixT;
      using ModelDataT = PropagationData<RealT, IdxT>;
      using VectorFitT = VectorFit<ScalarT, IdxT>;
      using DelayT     = Delay<ScalarT, IdxT>;

      explicit Propagation(const ModelDataT& data);
      void     attachInput(const std::vector<SignalT*>& input) override final;
      SignalT& outputSignal(IdxT channel) override final;
      int      allocate() override final;
      int      verify() const override final;
      int      initializeSteadyState(RealT omega, std::span<const RealT> u, std::span<const RealT> udot, RealT time = ZERO<RealT>) override final;
      void     transfer(RealT omega, MatrixT& re, MatrixT& im) const override final;
      int      setAbsoluteTolerance(RealT tolerance) override final;
      int      evaluateInternalResidual() override final;
      int      evaluateResidual() override final;
      int      assembleJacobian(RealT y_scale, RealT yp_scale) override final;

    private:
      static DelayData<RealT, IdxT> delayData(const ModelDataT& data);

      ModelDataT                               data_;
      std::vector<std::unique_ptr<VectorFitT>> fits_;
      std::vector<SignalT>                     filtered_, output_;
      DelayT                                   delay_;
      size_t                                   capacity_{0};
    };
  } // namespace EMT
} // namespace GridKit
