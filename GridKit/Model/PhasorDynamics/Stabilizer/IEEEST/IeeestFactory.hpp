/**
 * @file IeeestFactory.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Construct the IEEEST specialization for the notch order of its data.
 */

#pragma once

#include <GridKit/Model/ParameterReader.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeSet.hpp>
#include <GridKit/Model/PhasorDynamics/Stabilizer/IEEEST/Ieeest.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace Stabilizer
    {
      template <typename scalar_type, typename index_type>
      class IeeestFactory
      {
      public:
        using ScalarT        = scalar_type;
        using IdxT           = index_type;
        using ComponentT     = Component<ScalarT, IdxT>;
        using RealT          = typename ComponentT::RealT;
        using ModelDataT     = IeeestData<RealT, IdxT>;
        using SignalNodeSetT = SignalNodeSet<ScalarT, IdxT>;

        IeeestFactory() = delete;

        /// Construct and connect the specialization for the notch order of
        /// the data; the caller owns the result.
        static ComponentT* create(const ModelDataT& data, SignalNodeSetT& signal_nodes)
        {
          switch (readNotchOrder(data))
          {
          case 0:
            return construct<0>(data, signal_nodes);
          case 1:
            return construct<1>(data, signal_nodes);
          case 2:
            return construct<2>(data, signal_nodes);
          case 3:
            return construct<3>(data, signal_nodes);
          default:
            return construct<4>(data, signal_nodes);
          }
        }

        /// Read the notch order from the denominator coefficients.
        static size_t readNotchOrder(const ModelDataT& data)
        {
          using Params = typename ModelDataT::Parameters;

          Model::ParameterReader reader(data, "Ieeest");

          const RealT A1 = reader.readReal(Params::A1).value_or(0);
          const RealT A2 = reader.readReal(Params::A2).value_or(0);
          const RealT A3 = reader.readReal(Params::A3).value_or(0);
          const RealT A4 = reader.readReal(Params::A4).value_or(0);
          return notchOrder(A1, A2, A3, A4);
        }

      private:
        template <size_t order>
        static ComponentT* construct(const ModelDataT& data, SignalNodeSetT& signal_nodes)
        {
          auto* stabilizer = new Ieeest<ScalarT, IdxT, order>(data);
          stabilizer->getPorts().connect(data, signal_nodes);
          return stabilizer;
        }
      };
    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
