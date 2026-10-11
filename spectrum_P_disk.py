import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize, special

# Numerical check of the spectrum in the disk.
# See the section "Focusing modes in the disk and optimality of the decay rate".

n = 0
k_start = 8
k_end = 60
fit_tail = 20
tol = 1e-11
save_figure = True


def dispersion(lam, mode_n):
    """Equation (n^2+i*lambda-lambda^2)J_n(lambda)+lambda J_n'(lambda)=0."""
    return ((mode_n**2 + 1j*lam - lam**2)*special.jv(mode_n, lam)
            + lam*special.jvp(mode_n, lam, 1))


def solve_root(alpha, mode_n, root_tol):
    """Root close to a positive zero alpha of J_n."""

    # For n=0 the paper gives, up to O(k^-3),
    # lambda-alpha = alpha/(1+alpha^2) + i/(1+alpha^2).
    # This is also a good initial guess for the numerical solver.
    eps0 = alpha/(1+alpha**2) + 1j/(1+alpha**2)

    def root_system(xy):
        lam = xy[0] + 1j*xy[1]
        val = dispersion(lam, mode_n)
        return [val.real, val.imag]

    guesses = [
        alpha + eps0,
        alpha + 1j/(1+alpha**2),
        alpha + 0.5*eps0,
        alpha + 2*eps0,
    ]

    best_lam = guesses[0]
    best_res = abs(dispersion(best_lam, mode_n))

    for guess in guesses:
        out = optimize.root(
            root_system,
            x0=[guess.real, guess.imag],
            method="hybr",
            tol=root_tol,
        )

        lam = out.x[0] + 1j*out.x[1]
        res = abs(dispersion(lam, mode_n))

        if res < best_res:
            best_lam = lam
            best_res = res

        if out.success and res < 1e-9:
            return lam

    if best_res > 1e-7:
        raise RuntimeError(
            f"Root solver failed near alpha={alpha:.8f}; residual={best_res:.3e}"
        )

    return best_lam


print(f"Computing roots for n={n}, k={k_start},...,{k_end}")

k_values = np.arange(k_start, k_end+1)
alphas = special.jn_zeros(n, k_end)[k_start-1:k_end]
lambda_vals = np.array([solve_root(alpha, n, tol) for alpha in alphas])
z_vals = 1j*lambda_vals

residuals = np.abs([dispersion(lam, n) for lam in lambda_vals])
print(f"Maximum root residual: {residuals.max():.3e}")

# Asymptotic formula in the paper for n=0.
phi = k_values*np.pi - np.pi/4
lambda_asym = (phi + 1/(8*phi) + phi/(1+phi**2)
               + 1j/(1+phi**2))
z_asym = 1j*lambda_asym
error = np.abs(lambda_vals-lambda_asym)

omega = z_vals.imag
sigma = np.abs(z_vals.real)

tail = min(fit_tail, len(omega))
slope, b = np.polyfit(np.log(omega[-tail:]), np.log(sigma[-tail:]), 1)
cte = np.exp(b)

print("\nHigh-frequency behaviour")
print("  |Re z_k| ~ C*(Im z_k)^p")
print(f"  p = {slope:.6f}   (theoretical exponent: -2)")
print(f"  C = {cte:.6e}      (asymptotic constant: 1)")
print(f"  max k^3*|lambda_k-lambda_k^asym| on the fitted tail: "
      f"{np.max(k_values[-tail:]**3*error[-tail:]):.6e}")

fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))

ax[0].axvline(0.0, color="k", linestyle="--", linewidth=1.0)
ax[0].scatter(z_vals.real, z_vals.imag, s=26, color="#b22222")
ax[0].set_xlabel(r"$\Re z$")
ax[0].set_ylabel(r"$\Im z$")
ax[0].set_title(f"Eigenvalues of the pencil $P$ (mode n={n})")
ax[0].grid(True, alpha=0.3)

fit_curve = cte*omega**slope
ax[1].loglog(omega, sigma, "o", label=r"$|\Re z_k|$")
ax[1].loglog(omega, fit_curve, "--",
             label=rf"fit: ${cte:.2e}\,\omega^{{{slope:.2f}}}$")
ax[1].loglog(z_asym.imag, np.abs(z_asym.real), "-",
             label=r"asymptotic formula in the paper")
ax[1].set_xlabel(r"$\omega_k=\Im z_k$")
ax[1].set_ylabel(r"$|\Re z_k|$")
ax[1].set_title("Approach to the imaginary axis")
ax[1].legend()
ax[1].grid(True, which="both", alpha=0.3)

fig.tight_layout()
if save_figure:
    fig.savefig("spectrum_P_disk.pdf", bbox_inches="tight")
plt.show()
