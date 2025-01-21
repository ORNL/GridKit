# Load Model


Load models represent the relationship between the power and the voltage on the load bus. Depending on the type of studies different models are used. Generally, models can be divided into two categories: static and dynamic load models. The two most common models used for the static analysis are constant power (1) and constant impedance (2). 

```math
S_{L}=P_{L}+jQ_{L} 
```

```math
S_{L}=\frac{\vert V \vert ^2}{Z_{L}} 
```

Constant power loads are included in PF as a negative bus injection. The constant impedance load can be modeled as shunt element connected to a bus. 
The dispatchable load can be modeled as a negative generator.
For more advanced studies, frequency and voltage dependence should be considered too.  


# Shunt Model


Besides the main network elements listed above, the power grid also includes other devices such as shunt capacitors, reactors, and power electronic reactive control devices. These devices are used to control reactive injection into a given bus (thus control the voltage). There are passive (capacitors and reactors) as well as active shunts (advance power electronic devices that can vary the reactive output of the devices independent on the bus voltage). 
Passive elements are included in the network models as a fixed impedance to ground at a bus. 

## Passive 
```math
Y_{SH}=G_{SH}+jB_{SH}
```

## Active 

to be added