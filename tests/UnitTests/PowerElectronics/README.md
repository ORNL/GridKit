## HIRES Partitioning Test Problem

The HIRES test problem is a simple ODE system with eight variables. To
demonstrate the partitioning machinery in GridKit, the system is divided into
three components. Equations 1--3 belong to **Component 1**, equations 4--5
belong to **Component 2**, which is modeled as a bus, and equations 6--8 belong
to **Component 3**.

> **Note:** HIRES is not a circuit problem. In this example, it is modeled to
> resemble GridKit circuit components so that the existing partitioning
> machinery can be used. The example also provides a simple test problem for
> verifying the order of co-simulation methods.

The full HIRES system is

$$
\begin{aligned}
f_1 &= \frac{dy_1}{dt} +1.71y_1 -0.43y_2 -8.32y_3 -0.0007, \\
f_2 &= \frac{dy_2}{dt} -1.71y_1 +8.75y_2, \\
f_3 &= \frac{dy_3}{dt} +10.03y_3 -0.43y_4 -0.035y_5, \\
f_4 &= \frac{dy_4}{dt} -8.32y_2 -1.71y_3 +1.12y_4, \\
f_5 &= \frac{dy_5}{dt} +1.745y_5 -0.43y_6 -0.43y_7, \\
f_6 &= \frac{dy_6}{dt} +280y_6y_8 -0.69y_4 -1.71y_5
       +0.43y_6 -0.69y_7, \\
f_7 &= \frac{dy_7}{dt} -280y_6y_8 +1.81y_7, \\
f_8 &= \frac{dy_8}{dt} +280y_6y_8 -1.81y_7.
\end{aligned}
$$

For the component representation, equations $f_4$ and $f_5$ are decomposed
to expose the contributions from each component:

$$
\begin{aligned}
f_4
&= \frac{dy_4}{dt} -8.32y_2 -1.71y_3 +1.12y_4 \qquad
&\longrightarrow
\left(\frac{dy_4}{dt}+y_4\right)
+\left(0.1y_4-8.32y_2-1.71y_3\right)
+\left[0.02y_4\right],
\\ \\
f_5
&= \frac{dy_5}{dt}+1.745y_5-0.43y_6-0.43y_7 \qquad
&\longrightarrow
\left(\frac{dy_5}{dt}+y_5\right)
+\left(0.7y_5\right)
+\left[0.045y_5-0.43y_6-0.43y_7\right].
\end{aligned}
$$

The decomposition does not change the HIRES system. It separates the terms in
the bus equations (conviniently choosen to be equation 4 and 5) according to the component responsible for each contribution.

### Component 1

Component 1 has three internal equations:

$$
\begin{aligned}
f_1 &= \frac{dy_1}{dt} +1.71y_1 -0.43y_2 -8.32y_3 -0.0007, \\
f_2 &= \frac{dy_2}{dt} -1.71y_1 +8.75y_2, \\
f_3 &= \frac{dy_3}{dt} +10.03y_3 -0.43y_4 -0.035y_5.
\end{aligned}
$$

It also contributes the following terms to the bus equations as its external contribution:

$$
\begin{aligned}
f_4^{(1)} &= 0.1y_4 -8.32y_2 -1.71y_3, \\
f_5^{(1)} &= 0.7y_5.
\end{aligned}
$$

### Component 2 (HiresBus)

Component 2 represents the bus and owns the following contributions to
equations 4 and 5:

$$
\begin{aligned}
f_4^{(2)} &= \frac{dy_4}{dt} + y_4, \\
f_5^{(2)} &= \frac{dy_5}{dt} + y_5.
\end{aligned}
$$

The remaining terms in these equations are supplied by the components
connected to the bus.

### Component 3

Component 3 has three internal equations:

$$
\begin{aligned}
f_6 &= \frac{dy_6}{dt} -280y_6y_8 +0.69y_4 +1.71y_5
       -0.43y_6 +0.69y_7, \\
f_7 &= \frac{dy_7}{dt} +280y_6y_8 -1.81y_7, \\
f_8 &= \frac{dy_8}{dt} -280y_6y_8 +1.81y_7.
\end{aligned}
$$

It also contributes the following terms to the bus equations as its external contribution:

$$
\begin{aligned}
f_4^{(3)} &= 0.02y_4, \\
f_5^{(3)} &= 0.045y_5 -0.43y_6 -0.43y_7.
\end{aligned}
$$

### Bus Residual Assembly

The complete bus residuals are obtained by adding the contributions from
Components 1, 2, and 3:

$$
f_4 = f_4^{(1)} + f_4^{(2)} + f_4^{(3)},
$$

$$
f_5 = f_5^{(1)} + f_5^{(2)} + f_5^{(3)}.
$$

This reconstruction gives exactly the corresponding equations in the original
HIRES system.

### Partitioning

The full HIRES system in component form looks like this:

```text
Component 1  --------  Component 2 (HiresBus)  --------  Component 3
```

The system is divided between **Component 2 (HiresBus)** and **Component 3**, and
a `BusPartitionInterface` is added to the first partition. The resulting
partition residuals are then evaluated independently and are then compared with the full-system residual to verify
that the partitioned evaluation reproduces the original system.
