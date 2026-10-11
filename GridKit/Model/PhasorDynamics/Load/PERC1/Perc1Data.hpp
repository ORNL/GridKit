#pragma once

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// PERC1 parameters. Powers use system base; current limits use component base.
    enum class Perc1Parameters : size_t
    {
      Pnom,    ///< Initial active-power consumption on system base
      Qnom,    ///< Initial reactive-power consumption on system base
      Lfm,     ///< Loading factor defining the component power base
      QPratio, ///< Controlled reactive-to-active power ratio
      Dbfl,    ///< Lower frequency deadband threshold
      Dbfh,    ///< Upper frequency deadband threshold
      Kdroop,  ///< Frequency-response gain
      Kvp,     ///< Active-power voltage-washout coefficient
      Tvp,     ///< Active-power washout time constant
      Kvq,     ///< Reactive-power voltage-washout coefficient
      Tvq,     ///< Reactive-power washout time constant
      Tap,     ///< Active-power lead time constant
      Tbp,     ///< Active-power lag time constant
      Taq,     ///< Reactive-power lead time constant
      Tbq,     ///< Reactive-power lag time constant
      nP,      ///< Active-power voltage exponent
      nQ,      ///< Reactive-power voltage exponent
      Ipmax,   ///< Upper active-current limit on component base
      Ipmin,   ///< Lower active-current limit on component base
      Iqmax,   ///< Upper controlled reactive-current limit on component base
      Iqmin,   ///< Lower controlled reactive-current limit on component base
      Fcease,  ///< Fraction subject to cessation
      Vcease,  ///< Undervoltage detection threshold
      Tcease,  ///< Undervoltage qualification duration
      Tdelay,  ///< Committed delay before cessation
      Vrecon,  ///< Voltage recovery threshold
      Trecon,  ///< Recovery qualification duration
      Tramp,   ///< Reconnection ramp duration
      Frecon,  ///< Reconnected amount relative to the ceased fraction
      Tt,      ///< Output-current lag time constant
      Tv,      ///< Voltage measurement time constant
      Tf,      ///< Frequency measurement time constant
    };

    enum class Perc1Buses : size_t
    {
      bus, ///< Terminal bus
    };

    enum class Perc1SignalInputs : size_t
    {
      freq, ///< Optional absolute per-unit frequency, known at initialization
    };

    enum class Perc1SignalOutputs : size_t
    {
    };

    enum class Perc1MonitorableVariables : size_t
    {
      ir,
      ii,
      im,
      p,
      q,
      ip,
      iq,
      fracon,
    };

    /// Data for the power-electronic reconnecting and ceasing load model.
    template <typename real_type, typename index_type>
    using Perc1Data = ComponentData<real_type,
                                    index_type,
                                    Perc1Parameters,
                                    Perc1Buses,
                                    Perc1SignalInputs,
                                    Perc1SignalOutputs,
                                    Perc1MonitorableVariables>;
  } // namespace PhasorDynamics
} // namespace GridKit
