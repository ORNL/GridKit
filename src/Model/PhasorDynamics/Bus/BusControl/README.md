# Control Bus

A control bus is an auxillary channel representing a control 
signal or the exchange of shared information between models.

## Types of Control Buses

- Signal Bus (See [BusSignal](./BusSignal/))

# Signal Bus

A signal bus is allows multiple models to read/access information.


## Application

A signal bus allows a model to interface with any other model simply by publishing
its accessor method to the bus signal. Only one model is allowed to set the bus signal.

## Example

An example use case is demonstrated below. The accessor function of a generator
is passed as a callback function to the set_source function of the signal
which allows the generator speed to be accessed by any object that has
permissions to poll the signal object.

signal.set_source(gen.speed)

To read/poll the signal

signal.poll()

A full implementation may look like this

gen = ...
gov = ...

signal.set_source(gen.speed)
gov.set_speed_signal(signal)

Then inside the governor residual function

...

ScalarT speed = speed_signal->poll()
...
