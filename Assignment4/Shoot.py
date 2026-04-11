import os
import numpy as np
import matplotlib.pyplot as plt

def rk2_solve(x_l, x_r, dx, V, E, u_l=0.0, v_l=1e-5):
    """
    Solves the 2nd order ODE: u'' = 2(V(x) - E)u from x_l to x_r using RK2.
    """
    # Determine the number of steps and ensure dx lands exactly on boundaries
    N = int(np.abs((x_r - x_l) / dx)) + 1
    x_vals = np.linspace(x_l, x_r, N)
    actual_dx = x_vals[1] - x_vals[0]
    u = np.zeros(N)
    v = np.zeros(N)
    u[0] = u_l
    v[0] = v_l
    for i in range(N - 1):
        xi = x_vals[i]
        ui = u[i]
        vi = v[i]
        
        # Euler prediction
        k1_u = vi
        k1_v = 2.0 * (V(xi) - E) * ui
        x_mid = xi + actual_dx / 2.0
        u_mid = ui + actual_dx * k1_u / 2.0
        v_mid = vi + actual_dx * k1_v / 2.0
        
        # Midpoint correction
        k2_u = v_mid
        k2_v = 2.0 * (V(x_mid) - E) * u_mid
        
        # Update states
        u[i+1] = ui + actual_dx * k2_u
        v[i+1] = vi + actual_dx * k2_v
        
    return x_vals, u, v

def numerov_solve(x_l, x_r, dx, V, E, u_l=0.0, u_l_next=1e-5):
    N = int(np.abs((x_r - x_l) / dx)) + 1
    x_vals = np.linspace(x_l, x_r, N)
    actual_dx = x_vals[1] - x_vals[0]
    u = np.zeros(N)
    u[0] = u_l
    u[1] = u_l_next
    k = 2.0 * (V(x_vals) - E)
    factor = (actual_dx**2) / 12.0
    for i in range(1, N - 1):
        c0 = 1.0 - factor * k[i+1]
        c1 = 2.0 + 10.0 * factor * k[i]
        c2 = 1.0 - factor * k[i-1]
        u[i+1] = (c1 * u[i] - c2 * u[i-1]) / c0
    return x_vals, u

def shooting_method(x_l, x_r, dx, V, E_vec, v_l=1e-5, v_r=-1e-5, x_m=None, method='numerov'):
    E_vec = np.asarray(E_vec)
    num_E = len(E_vec)
    if x_m is None: x_m = 0.5 * (x_l + x_r)
    
    N_left = int(np.abs(x_m - x_l) / dx) + 1
    N_right = int(np.abs(x_r - x_m) / dx) + 1
    x_left = np.linspace(x_l, x_m, N_left)
    x_right = np.linspace(x_r, x_m, N_right)
    dx_L = x_left[1] - x_left[0]
    dx_R = x_right[1] - x_right[0]
    factor_L = (dx_L**2) / 12.0
    factor_R = (dx_R**2) / 12.0

    u_up = np.zeros(num_E)
    if method == 'rk2':
        v_up = np.full(num_E, v_l)
    elif method == 'numerov':
        u_up_prev = np.zeros(num_E)
        u_up = np.full(num_E, v_l * dx_L) 
        v_up = np.zeros(num_E) 

    for i in range(N_left - 1):
        xi = x_left[i]
        if method == 'rk2':
            k1_u = v_up
            k1_v = 2.0 * (V(xi) - E_vec) * u_up
            x_mid, u_mid, v_mid = xi + dx_L / 2.0, u_up + dx_L * k1_u / 2.0, v_up + dx_L * k1_v / 2.0
            k2_u = v_mid
            k2_v = 2.0 * (V(x_mid) - E_vec) * u_mid
            u_up, v_up = u_up + dx_L * k2_u, v_up + dx_L * k2_v
        elif method == 'numerov':
            xi_prev = x_left[i-1] if i > 0 else xi - dx_L
            xi_next = x_left[i+1]
            k_prev, k_curr, k_next = 2.0*(V(xi_prev)-E_vec), 2.0*(V(xi)-E_vec), 2.0*(V(xi_next)-E_vec)
            c0, c1, c2 = 1.0 - factor_L*k_next, 2.0 + 10.0*factor_L*k_curr, 1.0 - factor_L*k_prev
            u_next = (c1 * u_up - c2 * u_up_prev) / c0         
            u_up_prev, u_up = u_up.copy(), u_next
            v_up = (u_up - u_up_prev) / dx_L 

    u_down = np.zeros(num_E)
    if method == 'rk2':
        v_down = np.full(num_E, v_r)
    elif method == 'numerov':
        u_down_prev = np.zeros(num_E)
        u_down = np.full(num_E, v_r * dx_R) 
        v_down = np.zeros(num_E)

    for i in range(N_right - 1):
        xi = x_right[i]
        if method == 'rk2':
            k1_u = v_down
            k1_v = 2.0 * (V(xi) - E_vec) * u_down
            x_mid, u_mid, v_mid = xi + dx_R / 2.0, u_down + dx_R * k1_u / 2.0, v_down + dx_R * k1_v / 2.0
            k2_u = v_mid
            k2_v = 2.0 * (V(x_mid) - E_vec) * u_mid
            u_down, v_down = u_down + dx_R * k2_u, v_down + dx_R * k2_v
        elif method == 'numerov':
            xi_prev = x_right[i-1] if i > 0 else xi - dx_R
            xi_next = x_right[i+1]
            k_prev, k_curr, k_next = 2.0*(V(xi_prev)-E_vec), 2.0*(V(xi)-E_vec), 2.0*(V(xi_next)-E_vec)
            c0, c1, c2 = 1.0 - factor_R*k_next, 2.0 + 10.0*factor_R*k_curr, 1.0 - factor_R*k_prev
            u_next = (c1 * u_down - c2 * u_down_prev) / c0
            u_down_prev, u_down = u_down.copy(), u_next
            v_down = (u_down - u_down_prev) / dx_R 

    f_E = u_up * v_down - u_down * v_up
    sign_changes = np.where(np.diff(np.sign(f_E)))[0]
    roots = [E_vec[i] - f_E[i] * (E_vec[i+1] - E_vec[i]) / (f_E[i+1] - f_E[i]) for i in sign_changes]
    return f_E, roots

def electrostatic_potential_eff(r:(np.ndarray|float), l=0, mu=0, lamb=0, Z=1):
    return (lamb*r*r) + (l*(l+1)/(2*r*r)) - (Z*np.exp(-mu*r)/r)

def solve_radial_atom(l=0, mu=0.0, lamb=0.0, Z=1.0, r_max=100.0, dr=0.001, num_E=1000, E_min=-1.5, E_max=-0.01, return_f_E=True):
    """
    Solves for the eigen-energies of the general effective radial potential.
    """    
    def V_eff(r):
        return (lamb * r**2) + (l * (l + 1) / (2 * r**2)) - (Z * np.exp(-mu * r) / r)
    r_min = dr 
    E_vec = np.linspace(E_min, E_max, num_E)
    f_E, roots = shooting_method(r_min, r_max, dr, V_eff, E_vec, x_m = (0.9 * r_min + 0.1 * r_max))
    if return_f_E:
        return roots, E_vec, f_E
    return roots

def visualise_f_E(roots, E_vec, f_E, condition="",
    fig=None, ax=None, final=True):
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(E_vec, f_E, label=f"$f(E)$ {condition}")
    if roots:ax.scatter(roots, np.zeros_like(roots), zorder=5, label=f"roots {condition}")        
    if final:
        ax.axhline(0, color='black', linestyle='--', alpha=0.7)
        y_limit = np.median(np.abs(f_E)) * 2
        if roots:
            y_limit = np.median(np.abs(f_E)[E_vec + 0.01 > roots[0]]) * 2
            ax.set_xlim(roots[0]-0.01, roots[-1]+0.01)
        ax.set_ylim(-y_limit, y_limit)
        ax.set_yscale('symlog', linthresh=1e-3)
        ax.set_xlabel('$E$')
        ax.set_ylabel('Mismatch $f(E)$')
        ax.set_title(f'Radial Mismatch Function')
        ax.grid(True, alpha=0.3)
        ax.legend()
    return fig, ax

def E2n(E):return int(round((-2*E)**-0.5))

def visualise_eigenstates(roots, V_func, x_l, x_r, dx=0.001, condition="", 
                          fig=None, ax=None, final=True, method='numerov'):
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    condition = condition.strip()
    if condition != "" and not condition.startswith("(") and not condition.endswith(")"):
        condition = "(" + condition + ")"
        
    for E in roots:
        if method=='rk2':
            x_vals, u, v = rk2_solve(x_l, x_r, dx, V_func, E)
        elif method=='numerov':
            x_vals, u = numerov_solve(x_l, x_r, dx, V_func, E)
        else:raise ValueError(f"unknown method {method}")        
        u_norm = u / np.max(np.abs(u))
        n = E2n(E)
        ax.plot(x_vals, u_norm, label=f"$E_{n}={E:.4f}$ {condition}")
        
    if final:
        ax.axhline(0, color='black', linestyle='--', alpha=0.7)
        ax.set_xlabel('Coordinate ($r$)')
        ax.set_ylabel('Normalised Wavefunction $u(r)$')
        if not ax.get_title():ax.set_title('Eigenstates $u(r)$')
        ax.legend()
        
    return fig, ax

def visualise_eigenstates_stitched(roots, V_func, x_l, x_r, dx=0.001, condition="", 
                                   fig=None, ax=None, final=True,method='rk2'):
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
        
    N_total = int(np.abs(x_r - x_l) / dx) + 1
    x_vals = np.linspace(x_l, x_r, N_total)    
    V_arr = V_func(x_vals)

    condition = condition.strip()
    if condition != "" and not condition.startswith("(") and not condition.endswith(")"):
        condition = "(" + condition + ")"
    
    for E in roots:
        n = E2n(E)
        allowed_indices = np.where(E > V_arr)[0]
        if len(allowed_indices) > 0: match_idx = allowed_indices[-1] 
        else: match_idx = N_total // 2
        x_m = x_vals[match_idx]
        if method=='rk2':
            x_fwd, u_fwd, _ = rk2_solve(x_l, x_m, dx, V_func, E)
            x_bwd, u_bwd, _ = rk2_solve(x_r, x_m, -dx, V_func, E)
        elif method=='numerov':
            x_fwd, u_fwd = numerov_solve(x_l, x_m, dx, V_func, E)
            x_bwd, u_bwd = numerov_solve(x_r, x_m, -dx, V_func, E)
        else:raise ValueError(f"unknown method {method}")
        
        x_bwd = x_bwd[::-1]
        u_bwd = u_bwd[::-1]
        
        scale_factor = u_fwd[-1] / u_bwd[0]
        u_bwd_scaled = u_bwd * scale_factor
        
        x_stitched = np.concatenate((x_fwd, x_bwd[1:]))
        u_stitched = np.concatenate((u_fwd, u_bwd_scaled[1:]))
        u_norm = u_stitched / np.max(np.abs(u_stitched))
    
        ax.plot(x_stitched, u_norm, label=f"$E_{n}={E:.4f}$ {condition}")
        
    if final:
        ax.axhline(0, color='black', linestyle='--', alpha=0.7)
        ax.set_xlabel('$r$')
        ax.set_ylabel('Normalised $u(r)$')
        if not ax.get_title():ax.set_title('Eigenstates')
        ax.legend()
        
    return fig, ax

def compute_expectation(g, r, u):
    """
    Computes the expectation value <g(r)> for a given radial wavefunction u(r).
    As per the derived notes, this is given by:
    <g> = integral(g(r) * u^2 * dr) / integral(u^2 * dr)
    """
    # Calculate the integrand for the numerator and denominator
    numerator_integrand = g(r) * (u ** 2)
    denominator_integrand = u ** 2
    
    # Integrate using the trapezoidal rule
    numerator = np.trapz(numerator_integrand, r)
    denominator = np.trapz(denominator_integrand, r)
    
    return numerator / denominator

def MiniProject2():
    """
    Text results are saved to 'MP2.md', and plots are saved to the 'images/' folder.
    """
    if not os.path.exists('images'):
        os.makedirs('images')
        
    results_file = open("MP2.md", "w")
    
    def log(msg):
        print(msg)
        results_file.write(msg + "\n")
        
    # Common parameters
    r_max = 100.0  
    dr = 0.001
    E_min, E_max = -0.6, -0.01
    num_E = 2000

    log("## Part (a) - Plotting Effective Potentials")
    r_vals = np.linspace(0.1, 20.0, 500)
    plt.figure(figsize=(8, 5))
    for l in [0, 1, 2]:
        V_eff = electrostatic_potential_eff(r_vals, l=l, Z=1.0)
        plt.plot(r_vals, V_eff, label=f"$l={l}$")
    plt.ylim(-1.5, 0.5)
    plt.xlim(0, 20)
    plt.axhline(0, color='black', linestyle='--', alpha=0.5)
    plt.xlabel('$r$ (atomic units)')
    plt.ylabel(r'$V_{eff}(r)$')
    plt.title('Effective Potential for Hydrogen Atom')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig('images/P2a_V_eff.png')
    plt.close()
    log("* Saved images/P2a_V_eff.png\n")

    log("## Computing Normal Hydrogen states (l=0, 1, 2)")
    roots_l0, E_vec, f_E_l0 = solve_radial_atom(l=0, r_max=r_max, dr=dr, num_E=num_E, E_min=E_min, E_max=E_max)
    roots_l1, _, f_E_l1 = solve_radial_atom(l=1, r_max=r_max, dr=dr, num_E=num_E, E_min=E_min, E_max=E_max)
    roots_l2, _, f_E_l2 = solve_radial_atom(l=2, r_max=r_max, dr=dr, num_E=num_E, E_min=E_min, E_max=E_max)
    
    V_0 = lambda r: electrostatic_potential_eff(r, l=0, Z=1.0)
    V_1 = lambda r: electrostatic_potential_eff(r, l=1, Z=1.0)
    V_2 = lambda r: electrostatic_potential_eff(r, l=2, Z=1.0)

    fig1, axes1 = plt.subplots(3, 2, figsize=(12, 12), sharex='col')
    
    # l=0
    if len(roots_l0) > 0:
        visualise_eigenstates_stitched(roots_l0[:3], V_0, x_l=0.001, x_r=r_max, dx=dr, method='numerov', fig=fig1, ax=axes1[0,0], final=True)
    visualise_f_E(roots_l0, E_vec, f_E_l0, fig=fig1, ax=axes1[0,1], final=True)
    axes1[0,0].set_title('Eigenstates ($l=0$)')
    axes1[0,1].set_title('$f(E)$ Mismatch Function ($l=0$)')
    axes1[0,0].set_xlim(0, 30)

    # l=1
    if len(roots_l1) > 0:
        visualise_eigenstates_stitched(roots_l1[:2], V_1, x_l=0.001, x_r=r_max, dx=dr, method='numerov', fig=fig1, ax=axes1[1,0], final=True)
    visualise_f_E(roots_l1, E_vec, f_E_l1, fig=fig1, ax=axes1[1,1], final=True)
    axes1[1,0].set_title('Eigenstates ($l=1$)')
    axes1[1,1].set_title('$f(E)$ Mismatch Function ($l=1$)')
    axes1[1,0].set_xlim(0, 30)

    # l=2
    if len(roots_l2) > 0:
        visualise_eigenstates_stitched(roots_l2[:1], V_2, x_l=0.001, x_r=r_max, dx=dr, method='numerov',  fig=fig1, ax=axes1[2,0], final=True)
    visualise_f_E(roots_l2, E_vec, f_E_l2, fig=fig1, ax=axes1[2,1], final=True)
    axes1[2,0].set_title('Eigenstates ($l=2$)')
    axes1[2,1].set_title('$f(E)$ Mismatch Function ($l=2$)')
    axes1[2,0].set_xlim(0, 30)

    plt.tight_layout()
    fig1.savefig('images/P2_normal_states_grid.png')
    plt.close(fig1)
    log("* Saved image to images/P2_normal_states_grid.png")

    log("## Part (f) - Breaking Degeneracy")
    
    lam_val = 1e-5
    roots_lam0, _, f_E_lam0 = solve_radial_atom(l=0, lamb=lam_val, r_max=r_max, dr=dr, num_E=num_E, E_min=E_min, E_max=E_max)
    roots_lam1, _, f_E_lam1 = solve_radial_atom(l=1, lamb=lam_val, r_max=r_max, dr=dr, num_E=num_E, E_min=E_min, E_max=E_max)
    V_lam_0 = lambda r: electrostatic_potential_eff(r, l=0, lamb=lam_val, Z=1.0)
    V_lam_1 = lambda r: electrostatic_potential_eff(r, l=1, lamb=lam_val, Z=1.0)

    mu_val = 1e-2
    roots_mu0, _, f_E_mu0 = solve_radial_atom(l=0, mu=mu_val, r_max=r_max, dr=dr, num_E=num_E, E_min=E_min, E_max=E_max)
    roots_mu1, _, f_E_mu1 = solve_radial_atom(l=1, mu=mu_val, r_max=r_max, dr=dr, num_E=num_E, E_min=E_min, E_max=E_max)
    V_mu_0 = lambda r: electrostatic_potential_eff(r, l=0, mu=mu_val, Z=1.0)
    V_mu_1 = lambda r: electrostatic_potential_eff(r, l=1, mu=mu_val, Z=1.0)

    fig2, axes2 = plt.subplots(3, 2, figsize=(12, 12), sharex='col')

    # Row 0: Normal Potential
    visualise_eigenstates_stitched(roots_l0[1:3], V_0, x_l=0.001, x_r=r_max, dx=dr, condition="($s$)", method='numerov', fig = fig2, ax=axes2[0,0], final=False)
    visualise_eigenstates_stitched(roots_l1[0:2], V_1, x_l=0.001, x_r=r_max, dx=dr, condition="($p$)", method='numerov', fig = fig2, ax=axes2[0,0], final=True)
    visualise_f_E(roots_l0, E_vec, f_E_l0, condition="($l=0$)", fig = fig2, ax=axes2[0,1], final=False)
    visualise_f_E(roots_l1, E_vec, f_E_l1, condition="($l=1$)", fig = fig2, ax=axes2[0,1], final=True)
    axes2[0,0].set_title('Normal Potential (Degenerate)')
    axes2[0,0].set_xlim(0, 30)

    # Row 1: Non-zero lambda
    if len(roots_lam0) > 1: visualise_eigenstates_stitched(roots_lam0[1:3], V_lam_0, x_l=0.001, x_r=r_max, dx=dr, condition="($s_{\lambda})$", method='numerov', fig = fig2, ax=axes2[1,0], final=False)
    if len(roots_lam1) > 0: visualise_eigenstates_stitched(roots_lam1[0:2], V_lam_1, x_l=0.001, x_r=r_max, dx=dr, condition="($p_{\lambda})$", method='numerov', fig = fig2, ax=axes2[1,0], final=True)
    visualise_f_E(roots_lam0, E_vec, f_E_lam0, condition="($l=0$)", fig = fig2, ax=axes2[1,1], final=False)
    visualise_f_E(roots_lam1, E_vec, f_E_lam1, condition="($l=1$)", fig = fig2, ax=axes2[1,1], final=True)
    axes2[1,0].set_title(f'Harmonic Perturbation ($\lambda={lam_val}$)')
    axes2[1,0].set_xlim(0, 30)

    # Row 2: Non-zero mu
    if len(roots_mu0) > 1: visualise_eigenstates_stitched(roots_mu0[1:3], V_mu_0, x_l=0.001, x_r=r_max, dx=dr, condition="($s_{\mu}$)", method='numerov', fig = fig2, ax=axes2[2,0], final=False)
    if len(roots_mu1) > 0: visualise_eigenstates_stitched(roots_mu1[0:2], V_mu_1, x_l=0.001, x_r=r_max, dx=dr, condition="($p_{\mu}$)", method='numerov', fig = fig2, ax=axes2[2,0], final=True)
    visualise_f_E(roots_mu0, E_vec, f_E_mu0, condition="($l=0$)", fig = fig2, ax=axes2[2,1], final=False)
    visualise_f_E(roots_mu1, E_vec, f_E_mu1, condition="($l=1$)", fig = fig2, ax=axes2[2,1], final=True)
    axes2[2,0].set_title(f'Yukawa Screening ($\mu={mu_val}$)')
    axes2[2,0].set_xlim(0, 30)

    plt.tight_layout()
    fig2.savefig('images/P2_degeneracy_breaking_grid.png')
    plt.close(fig2)
    log("* Saved image to images/P2_degeneracy_breaking_grid.png")
    
    # Log results of the degeneracy breaking
    E_2s = roots_l0[1] if len(roots_l0) > 1 else None
    E_2p = roots_l1[0] if len(roots_l1) > 0 else None
    E_2s_lam = roots_lam0[1] if len(roots_lam0) > 1 else None
    E_2p_lam = roots_lam1[0] if len(roots_lam1) > 0 else None
    E_2s_mu = roots_mu0[1] if len(roots_mu0) > 1 else None
    E_2p_mu = roots_mu1[0] if len(roots_mu1) > 0 else None
    
    if E_2s and E_2p:           log(f"* Normal :\t\t|E_2s - E_2p| = {abs(E_2s - E_2p):.2e}")
    if E_2s_lam and E_2p_lam:   log(f"* $\lambda$ = {lam_val} :\t\t|E_2s - E_2p| = {abs(E_2s_lam - E_2p_lam):.2e}")
    if E_2s_mu and E_2p_mu:     log(f"* $\mu$ = {mu_val}$ :\t\t|E_2s - E_2p| = {abs(E_2s_mu - E_2p_mu):.2e}")

    log("\n## Expectation Values for " + r"$\braket{r}$ and $\braket{1/r}$")
    log(r"State | $\braket{r}$ | $\braket{1/r}$")
    log("-" * 45)
    
    def get_stitched_u(E, V_func):
        x_l, x_r = 0.001, r_max
        N_total = int(np.abs(x_r - x_l) / dr) + 1
        x_vals = np.linspace(x_l, x_r, N_total)
        V_arr = V_func(x_vals)
        allowed_indices = np.where(E > V_arr)[0]
        match_idx = allowed_indices[-1] if len(allowed_indices) > 0 else N_total // 2
        x_m = x_vals[match_idx]
        
        x_fwd, u_fwd = numerov_solve(x_l, x_m, dr, V_func, E)
        x_bwd, u_bwd = numerov_solve(x_r, x_m, -dr, V_func, E)
        
        x_bwd, u_bwd = x_bwd[::-1], u_bwd[::-1]
        scale_factor = u_fwd[-1] / u_bwd[0]
        u_bwd_scaled = u_bwd * scale_factor
        
        r_stitched = np.concatenate((x_fwd, x_bwd[1:]))
        u_stitched = np.concatenate((u_fwd, u_bwd_scaled[1:]))
        return r_stitched, u_stitched

    g_r = lambda r: r
    g_1_over_r = lambda r: 1.0 / r

    states_to_compute = [
        ("1s", roots_l0, 0, V_0),
        ("2s", roots_l0, 1, V_0),
        ("2p", roots_l1, 0, V_1),
        ("3s", roots_l0, 2, V_0),
        ("3p", roots_l1, 1, V_1),
        ("3d", roots_l2, 0, V_2)
    ]

    for name, roots, idx, V_func in states_to_compute:
        if len(roots) > idx:
            E_val = roots[idx]
            r_arr, u_arr = get_stitched_u(E_val, V_func)
            exp_r = compute_expectation(g_r, r_arr, u_arr)
            exp_1_over_r = compute_expectation(g_1_over_r, r_arr, u_arr)
            log(f"{name:<10} | {exp_r:<15.5f} | {exp_1_over_r:<15.5f}")
    results_file.close()

if __name__ == "__main__":
    MiniProject2()