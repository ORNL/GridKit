/**
 * @file Views.hpp
 *
 * @brief Signal views over the per-unit-length parameter matrices.
 *
 * Consumers of the series and shunt pairs depend on these views rather
 * than on the element that assembles them, so a reduction stage can
 * substitute for the full-conductor assembly without the consumers
 * changing. Signals carry global offsets that exist only after element
 * binding, so the accessors are late calls rather than stored values.
 *
 */

#pragma once

#include <GridKit/Model/EMT/Signal.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace Parameters
    {
      /// Series pair Z(omega) = R + j omega L as matrix signals.
      template <typename scalar_type, typename index_type>
      class SeriesView
      {
      public:
        using SignalT = Signal<index_type>;

        virtual ~SeriesView() = default;

        virtual SignalT R() const = 0;
        virtual SignalT L() const = 0;
      };

      /// Shunt pair Y(omega) = G + j omega C as matrix signals.
      template <typename scalar_type, typename index_type>
      class ShuntView
      {
      public:
        using SignalT = Signal<index_type>;

        virtual ~ShuntView() = default;

        virtual SignalT G() const = 0;
        virtual SignalT C() const = 0;
      };
    } // namespace Parameters
  } // namespace EMT
} // namespace GridKit
