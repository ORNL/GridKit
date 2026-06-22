# Bus-to-Signal Adapter

This component enables signals to send and receive bus variables. It has five
ports:

## Bus Port
- `bus` for the bus whose variables are managed by the adapter

## Input Ports
- `ir` ($I_r$)
- `ii` ($I_i$)

External current injections are read from input signal nodes added to currents
on the bus.

## Output Ports
- `vr` ($V_r$)
- `vi` ($V_i$)

Voltages read from the bus are made available to signal nodes.
