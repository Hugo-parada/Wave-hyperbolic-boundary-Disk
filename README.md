# Spectrum of the wave equation with hyperbolic boundary condition on the disk

Python code related to the spectral construction used in:

> H. Parada and N. Vanspranghe, **Optimal stabilization rate for the wave equation with hyperbolic boundary condition**.  
> HAL: https://hal.science/hal-05426803

The script `spectrum_P_disk.py` is a numerical check of the high-frequency eigenvalues used in the optimality part of the paper.

## Spectral equation

On the unit disk, separation of variables gives radial eigenfunctions

```text
u(r,theta)=J_n(lambda*r)e^{i n theta}.
```

For a complex number `lambda`, the boundary condition gives the characteristic equation

```text
(n^2+i*lambda-lambda^2) J_n(lambda)
    +lambda J_n'(lambda)=0.
```

If `lambda` is a root, then

```text
z=i*lambda
```

is an eigenvalue of the quadratic pencil `P` considered in the paper.

The optimality argument only needs the radial mode `n=0`, for which the equation becomes

```text
(i*lambda-lambda^2)J_0(lambda)+lambda J_0'(lambda)=0.
```

## What the code does

The code starts from the positive zeros `alpha_k` of `J_0`, because the complex roots of the dissipative problem are small perturbations of them at high frequency.

For every `k` between `k_start` and `k_end`, it solves the two real equations

```text
Re F(lambda)=0,
Im F(lambda)=0,
```

where

```text
F(lambda)=(i*lambda-lambda^2)J_0(lambda)+lambda J_0'(lambda).
```

The initial guess uses the asymptotic displacement proved in the paper. The script then sets

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

This is the spectral behaviour used to show that the resolvent estimate, and consequently the `1/t` energy decay rate, cannot be improved.

The script checks this numerically in two ways:

1. it fits `|Re z_k| = C*(Im z_k)^p` on the high-frequency tail and should obtain `p` close to `-2` and `C` close to `1`;
2. it compares the computed roots with the complete asymptotic expression above and reports the scaled error `k^3 |lambda_k-lambda_k^asym|`.

With the default values `k=8,...,60`, the code gives approximately

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
