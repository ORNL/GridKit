# Phasor Dynamics Events

## Overview

This document describes the runtime event system for PhasorDynamics
components. A runtime event mutates the state of a single component during simulation. Events are scheduled in a solver file
(see [PDSim](../../../application/PhasorDynamics/README.md)) and delivered to components by `SystemModel` at their scheduled times. Buses are addressed for events by their existing `number`; devices by their
existing `id`.

## Event

An event is a triple of fields:

Name     | Description
---------|------------------------------------------------------
`time`   | Simulation time at which the event fires, in seconds
`target` | Routing key for the component receiving the event
`action` | The mutation to perform, drawn from the action vocabulary

For buses, the `target` is the bus's `number` rendered as a string
(e.g. `"1001"`). For devices, it is the device's `id` field from the case
file (e.g. `"BR_1001_1064_1"`). Buses and devices share a single routing
namespace; collisions between a bus number and a device id are an error at registration.


```cpp
struct Event {
    double      time;
    std::string target;
    Action      action;
};
```

## Action vocabulary

`Action` is a `std::variant` of action structs in
`GridKit::PhasorDynamics::Actions`. Each struct has a canonical name string
exposed as `static constexpr std::string_view name`.

C++ type         | JSON name | Params              | Applies to
-----------------|-----------|---------------------|----------------
`Actions::Open`  | `open`    | none                | `Branch`
`Actions::Close` | `close`   | none                | `Branch`
`Actions::Fault` | `fault`   | `R`, `X`, `percent` | `Bus`, `Branch`
`Actions::Clear` | `clear`   | none                | `Bus`, `Branch`

The semantics of each action are documented in the receiving component's
README. For `fault`, the `percent` parameter (position along the line as a
percent in `[0, 100]`) is used by `Branch` and ignored by `Bus`.

## Dispatch

`Component::apply(const Action&)` is virtual on the base class. The default
throws. Components that handle events override it with a `std::visit` over
the variant, dispatching to per-action handlers and falling through to a
generic catch-all for unhandled action types:

```cpp
void Bus::apply(const Action& a) override {
  std::visit(overloaded{
    [this](const Actions::Fault& f) {},
    [this](const Actions::Clear&)   {},
    [this](const auto& other) {
      using A = std::decay_t<decltype(other)>;
      throw std::runtime_error("Bus does not handle action '" + std::string(A::name) + "'");
    }
  }, a);
}
```

`SystemModel::apply(const Event&)` looks up the target by string id and
forwards the action to the component.

## Schedule

The `schedule` field of the solver file is an array of events. The parser
stable-sorts the schedule by `time` on read. If the input was not already
sorted, an info line is logged and parsing continues.

Multiple events at the same time are applied sequentially in listing order in a single re-init of the integrator (`IDACalcIC`). When two same-time events conflict on a single target, the last-listed event wins.

The same-time grouping uses exact `double` equality on `time`. Both sides of the comparison are read from the parsed event records, so events with the same JSON literal compare equal.

## Example schedule

A schedule that faults bus `2` at `t=1.0`, clears the fault at `t=1.1`, and
simultaneously opens line `L23` at `t=1.1`:

```json
"schedule": [
  { "time": 1.0, "target": "2",   "action": "fault", "params": { "R": 0.0, "X": 1.0e-5 } },
  { "time": 1.1, "target": "2",   "action": "clear" },
  { "time": 1.1, "target": "L23", "action": "open"  }
]
```

## Errors

Failure                        | Source              | Message
-------------------------------|---------------------|--------
Unknown target id              | `SystemModel`       | `No event target with id '<id>'`
Unhandled action for component | `Component::apply`  | `Component '<class>' does not handle action '<name>'`
Unknown action string          | parser              | `unknown action '<str>' (valid: open, close, fault, clear)`
Missing required `params`      | parser              | `action '<name>' requires params field with <fields>`
Schedule entry past `tmax`     | parser              | `schedule entry at t=<t> exceeds tmax=<tmax>`

## Adding a new action

1. Add a struct to `GridKit::PhasorDynamics::Actions` with a
   `static constexpr std::string_view name` and any payload fields.
2. Add the struct to the `Action` variant alias.
3. Add a parser case mapping the JSON name to the struct.
4. Override `apply` in the receiving component(s) to handle the new action.
5. Add a row to the action vocabulary table above and an entry in each
   receiving component's README Actions section.
