# Co-simulation example

The client and server exchange voltage and current over ZeroMQ on TCP port
8800. The interface buses connect directly to signal nodes in the case files:

- The client uses `BusSignalVoltageOut` (JSON class `SignalVoltageOut`). It
  sends voltage on signals 1 and 2 (`vr`, `vi`) and receives current on
  signals 3 and 4 (`ir`, `ii`).
- The server uses `BusSignalVoltageIn` (JSON class `SignalVoltageIn`). It
  receives voltage on signals 3 and 4 (`vr`, `vi`) and sends the accumulated
  current injections on signals 1 and 2 (`ir`, `ii`).

Note that this example code is temporary and meant simply as an initial
step towards implementing multi-instance co-simulation with GridKit. It
will be generalized into an application that can be used for other
co-simulation cases.
