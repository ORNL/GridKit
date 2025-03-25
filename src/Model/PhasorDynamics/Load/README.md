# Load Model

Load modeling is one of the more complex aspects of power system dynamics.
The simplest model, which is used for this challenge problem, is to model
the load as a complex shunt impedance with the impedance given by:
``` math
Z = R + jX
```
where $`R`$ is the load resistance, $`X`$ is the load reactance. The current
drawn by the load is then obtained as
```math
I_{\mathrm{load}} = \frac{V_{\mathrm{bus}}}{Z},
```
where $`V_{bus}`$ is the voltage on the bus to which the load is connected.

After some algebra, one obtains expressions for real and imaginary components
for the currents entering the bus:
```math
I_{r} = -g V_{r} - b V_{i} 
```

```math
I_{i} = b V_{r} - g V_{i}
```
where
```math
g = \frac{R}{R^2+X^2} ~~~\mathrm{and}~~~ b = \frac{-X}{R^2+X^2}.
```
