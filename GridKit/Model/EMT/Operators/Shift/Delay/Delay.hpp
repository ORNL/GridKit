/**
 * @file Delay.hpp
 * @brief Constant transport delays with accepted-step history.
 */
#pragma once

#include <deque>
#include <functional>
#include <set>

#include <GridKit/Model/EMT/Operators/Shift/Delay/DelayData.hpp>
#include <GridKit/Model/EMT/Operators/Shift/Shift.hpp>

namespace GridKit
{
  namespace EMT
  {
    template <typename scalar_type, typename index_type>
    class Delay : public Shift<scalar_type, index_type>
    {
    public:
      using Base       = Shift<scalar_type, index_type>;
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename Base::RealT;
      using SignalT    = typename Base::SignalT;
      using MatrixT    = typename Base::MatrixT;
      using ModelDataT = DelayData<RealT, IdxT>;
      using Prehistory = std::function<std::pair<RealT, RealT>(size_t, RealT)>;

      explicit Delay(const ModelDataT& data);
      void     attachInput(const std::vector<SignalT*>& input) override final;
      void     setInputDerivative(std::function<RealT(size_t)> derivative);
      SignalT& outputSignal(IdxT channel) override final;
      void     setPrehistory(RealT time, Prehistory prehistory);
      int      initializeSteadyState(RealT omega, std::span<const RealT> u, std::span<const RealT> udot, RealT time = ZERO<RealT>) override final;
      int      initialize();
      int      allocate() override final;
      int      verify() const override final;
      int      setAbsoluteTolerance(RealT tolerance) override final;
      void     transfer(RealT omega, MatrixT& re, MatrixT& im) const override final;
      void     resetHistory() override final;
      void     acceptStep(RealT time) override final;
      RealT    maximumStepSize() const override final;
      RealT    nextDiscontinuityTime(RealT after) const override final;
      void     beginDiscontinuity(RealT time) override final;
      void     updateTime(RealT time, RealT alpha) override final;
      int      evaluateInternalResidual() override final;
      int      evaluateResidual() override final;
      int      assembleJacobian(RealT y_scale, RealT yp_scale) override final;

    private:
      struct Knot
      {
        RealT              time;
        std::vector<RealT> value, slope;
        bool               smooth;
      };

      RealT history(size_t k, RealT time) const;
      void  reserveJacobian(size_t capacity);
      void  appendJacobian(IdxT row, IdxT column, RealT value);

      size_t jacobian_capacity_{0};

      ModelDataT                               data_;
      std::vector<SignalT*>                    input_;
      std::vector<SignalT>                     output_;
      std::vector<RealT>                       constant_, coefficient_;
      std::vector<typename SignalT::GradientT> gradients_;
      std::deque<Knot>                         knots_;
      Prehistory                               prehistory_;
      std::function<RealT(size_t)>             derivative_;
      RealT                                    origin_{0.0};
      RealT                                    tau_max_{0.0};
      RealT                                    right_limit_time_{0.0};
      std::set<RealT>                          discontinuities_;
    };
  } // namespace EMT
} // namespace GridKit
