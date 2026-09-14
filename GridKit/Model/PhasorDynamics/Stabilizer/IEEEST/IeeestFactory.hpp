/**
 * @file IeeestFactory.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Construct the IEEEST specialization implied by its notch coefficients.
 */

#pragma once

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <variant>

#include <magic_enum/magic_enum.hpp>

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

        /// Construct and connect the derived order; the caller owns the result.
        static ComponentT* create(const ModelDataT& data, SignalNodeSetT& signal_nodes)
        {
          using Params   = typename ModelDataT::Parameters;
          const RealT A1 = readCoefficient(data, Params::A1);
          const RealT A2 = readCoefficient(data, Params::A2);
          const RealT A3 = readCoefficient(data, Params::A3);
          const RealT A4 = readCoefficient(data, Params::A4);

          switch (ieeestNotchOrder(A1, A2, A3, A4))
          {
          case 0:
            return createIeeest<0>(data, signal_nodes);
          case 1:
            return createIeeest<1>(data, signal_nodes);
          case 2:
            return createIeeest<2>(data, signal_nodes);
          case 3:
            return createIeeest<3>(data, signal_nodes);
          default:
            // The sum of two quadratic degrees cannot exceed four.
            return createIeeest<4>(data, signal_nodes);
          }
        }

      private:
        static RealT readCoefficient(const ModelDataT& data, IeeestParameters parameter)
        {
          const auto entry = data.parameters.find(parameter);
          if (entry == data.parameters.end())
          {
            return RealT{0};
          }

          RealT value{};
          if (const auto* real = std::get_if<RealT>(&entry->second))
          {
            value = *real;
          }
          else if (const auto* integer = std::get_if<IdxT>(&entry->second))
          {
            value = static_cast<RealT>(*integer);
          }
          else
          {
            throw std::invalid_argument("Ieeest: parameter '" + std::string(magic_enum::enum_name(parameter))
                                        + "' must be numeric");
          }
          if (!std::isfinite(value))
          {
            throw std::invalid_argument("Ieeest: parameter '" + std::string(magic_enum::enum_name(parameter))
                                        + "' must be finite");
          }
          return value;
        }

        template <size_t order>
        static ComponentT* createIeeest(const ModelDataT& data, SignalNodeSetT& signal_nodes)
        {
          auto stabilizer = std::make_unique<Ieeest<ScalarT, IdxT, order>>(data);
          stabilizer->getPorts().connect(data, signal_nodes);
          return stabilizer.release();
        }
      };
    } // namespace Stabilizer
  } // namespace PhasorDynamics
} // namespace GridKit
