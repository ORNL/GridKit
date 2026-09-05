# EMT state format specification

## Overview

This document describes the JSON data format for EMT operating points. A
case file carries parameters and topology only; a state file carries the
operating point. The state format is the shared
[Model::StateData](../StateData.hpp) container with the instantaneous
phase-coordinate fields, so partial states are legal: any missing record or
field keeps its default.

The state carries instantaneous SI values at the initialization instant.
Synthesizing instantaneous values from an RMS or phasor description is an
upstream tool concern and is outside the EMT model specification.

## Format

The root object may contain `header`, `buses`, and `devices`.

### Header

Contained in the optional `header` key is an object with the following items:

   Name          | Value
 ----------------|-------------------------------------------------------
  `version`      | Optional integer identifying the state format version
  `time`         | Optional floating-point model time of the state
  `created`      | Optional string with the wall-clock creation time
  `description`  | Optional string describing the state

### Buses

Contained in the `buses` key is an object mapping Bus component paths from the
case file to bus states. A root Bus uses its local ID; a nested Bus uses its
dot-qualified path, for example `left.network.bus`:

   Name  | Value
 --------|-------------------------------------------------------
  `va`   | Optional instantaneous phase a voltage in volts
  `vb`   | Optional instantaneous phase b voltage in volts
  `vc`   | Optional instantaneous phase c voltage in volts

### Devices

Contained in the `devices` key is an object mapping component paths from the
case file to device states. Nested paths are qualified by each containing
Container, for example `plant.machine`:

   Name   | Value
 ---------|-------------------------------------------------------
  `p`     | Optional active power injection in watts
  `q`     | Optional reactive power injection in vars
  `open`  | Optional Boolean switch command, true is open

## Application

[StateDataAdapter](StateDataAdapter.hpp) applies a parsed state to parsed
system model data before the system model is constructed, so model
constructors remain the single ingestion path: bus voltages land in the bus
initial conditions, machine dispatch lands in the machine `p0` and `q0`
parameters, and switch commands land in the switch `open` parameter.
