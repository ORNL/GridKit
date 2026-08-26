# Three-bus co-simulation example case with constant signal source

This example duplicates the behavior of the ThreeBusConstantSource example, but
with the `ConstantSignalSource` component removed and its behavior reproduced in
the app (external to the system model) with currents being received over zmq
from another app.

Note that this example code is temporary and meant simply as an initial step
towards implementing multi-instance co-simulation with GridKit. It will be
generalized into an application that can be used for other co-simulation cases.
