# Reference Frame Operators

Reference frame operators track the electrical reference angle and frequency and
transform real, instantaneous three-component signals. Clarke and Park use
power-invariant normalization. The transformations support the inverse
direction through the `inverse` parameter.

## Models

- [PLL](PLL/README.md): bus voltage to reference angle and frequency.
- [Clarke](Clarke/README.md): $abc$ to $\alpha\beta0$; specification only.
- [Park](Park/README.md): $abc$ to $dq0$.
- [Rotation](Rotation/README.md): $\alpha\beta0$ to $dq0$; specification only.
