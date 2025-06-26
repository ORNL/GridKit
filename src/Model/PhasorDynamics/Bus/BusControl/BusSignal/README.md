# Signal Bus

A signal bus is allows multiple models to read/access information.


## Application

A signal bus allows a model to interface with any other model simply by publishing
its accessor method to the bus signal. Only one model is allowed to set the bus signal.
