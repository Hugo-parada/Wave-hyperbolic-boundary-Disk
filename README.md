# Spectrum of the wave equation with hyperbolic boundary condition on the disk

Python script for the spectral computation in:

> H. Parada and N. Vanspranghe, **Optimal stabilization rate for the wave equation with hyperbolic boundary condition**.  
> HAL: https://hal.science/hal-05426803

`spectrum_P_disk.py` computes high-frequency eigenvalues for the radial mode and compares them with the asymptotic formula from the paper.

## Spectral equation

On the unit disk, separation of variables gives radial eigenfunctions

```text
u(r,theta)=J_n(lambda*r)e^{i n theta}.
```

The boundary condition gives the characteristic equation

```text
(n^2+i*lambda-lambda^2) J_n(lambda)
    +lambda J_n'(lambda)=0.
```

For each root `lambda`,

```text
z=i*lambda
```

is the corresponding eigenvalue of the quadratic pencil `P`.

For the radial mode `n=0`, the equation becomes

```text
(i*lambda-lambda^2)J_0(lambda)+lambda J_0'(lambda)=0.
```

## What the code does

The solver starts near the positive zeros `alpha_k` of `J_0`.

For `k_start <= k <= k_end`, it solves

```text
Re F(lambda)=0,
Im F(lambda)=0,
```

where

```text
F(lambda)=(i*lambda-lambda^2)J_0(lambda)+lambda J_0'(lambda).
```

The initial guess comes from the asymptotic formula in the paper. The script sets

```text
z_k=i*lambda_k
```

and plots the resulting eigenvalues in the complex plane.

## Relation with the asymptotic formula in the paper

For

```text
phi_k=k*pi-pi/4,
```

the paper proves

```text
lambda_k = phi_k
           +1/(8*phi_k)
           +phi_k/(1+phi_k^2)
           +i/(1+phi_k^2)
           +O(k^-3).
```

Therefore

```text
Re z_k=-Im lambda_k=-1/(1+phi_k^2)+O(k^-3),
Im z_k= Re lambda_k ~ k*pi.
```

In particular,

```text
|Re z_k| ~ (Im z_k)^(-2).
```

The paper uses this behaviour to show that the resolvent bound and the `1/t` energy decay rate are sharp.

The numerical comparison has two parts:

1. A fit of `|Re z_k| = C*(Im z_k)^p` on the high-frequency tail, with expected values `p=-2` and `C=1`.
2. A comparison with the asymptotic formula above, using the scaled error `k^3 |lambda_k-lambda_k^asym|`.

Typical output with the default range `k=8,...,60`:

```text
p = -1.999664,
C = 0.998137,
maximum root residual = 2.1e-11,
max k^3 |lambda_k-lambda_k^asym| on the fitted tail = 6.70e-2.
```

## Requirements

```text
Python 3
NumPy
SciPy
Matplotlib
```

Run with

```bash
python spectrum_P_disk.py
```
