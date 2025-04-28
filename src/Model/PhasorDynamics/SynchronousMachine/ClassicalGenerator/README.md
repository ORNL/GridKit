Differential equations:

```math 
\dot{\delta} = \omega \cdot \omega _0 
```
```math
\dot{\omega} = \frac{1}{2H}\bigg( \frac{P_{mech} - D\omega _0}{1 + \omega}   - T_{elec}\bigg)
```

Algebraic Equations: 

```math
    T_{elec} = \frac{1}{1+\omega}\bigg( g E_p^2 - E_p \bigg((gV_r - bV_i)cos\,\delta + (bV_r + gV_i)sin\,\delta \bigg)\bigg)
```

Network Interface Equations:

```math
I_r = -gV_r + bV_i + E_p(g \cos \delta - b \sin \delta)
```
```math
I_i = -gV_r - bV_i + E_p(b \cos \delta + g \sin \delta)
```

Intialization notes: <br>
To initialize the model, given $V_r$, $V_i$, $P$ and $Q$, we use following equations:
<br><br>
```math
I_r = \frac{PV_r + QV_i}{V_r^2 + V_i^2} 
```
```math
I_i = \frac{PV_i - QV_r}{V_r^2 + V_i^2} 
```
```math
E_r = \frac{ g(I_r + gV_r - bV_i) + b (I_i + bV_r + gV_i)   }{g^2 + b^2}
```
```math
E_i = \frac{ -b(I_r + gV_r - bV_i) + g (I_i + bV_r + gV_i)   }{g^2 + b^2}
```
```math
E_p = \sqrt{E_r^2 + E_i^2}
```
```math
\delta = atan2(E_i, E_r) 
```
```math
\omega = 0 
```
```math
T_{elec} = gE_p^2 - E_p \bigg( \bigg(gV_r -  bV_i \bigg) \cos \delta + \bigg(bV_r + gV_i \bigg)\sin \delta\bigg)
```
```math
P_{mech} = T_{elec}
```