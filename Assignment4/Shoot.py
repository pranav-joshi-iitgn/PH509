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
        
        # RK2 Step 1 (Euler predictor)
        k1_u = vi
        k1_v = 2.0 * (V(xi) - E) * ui
        
        x_mid = xi + actual_dx / 2.0
        u_mid = ui + actual_dx * k1_u / 2.0
        v_mid = vi + actual_dx * k1_v / 2.0
        
        # RK2 Step 2 (Midpoint corrector)
        k2_u = v_mid
        k2_v = 2.0 * (V(x_mid) - E) * u_mid
        
        # Update states
        u[i+1] = ui + actual_dx * k2_u
        v[i+1] = vi + actual_dx * k2_v
        
    return x_vals, u, v

def numerov_solve(x_l, x_r, dx, V, E, u_l=0.0, u_l_next=1e-5):
    """
    Solves u'' = k(x)u using Numerov's Method for O(dx^4) precision.
    """
    N = int(np.abs((x_r - x_l) / dx)) + 1
    x_vals = np.linspace(x_l, x_r, N)
    actual_dx = x_vals[1] - x_vals[0]
    
    u = np.zeros(N)
    u[0] = u_l
    u[1] = u_l_next
    
    k = 2.0 * (V(x_vals) - E)
    factor = (actual_dx**2) / 12.0
    
    for i in range(1, N - 1):
        # THE FIX: Corrected signs for the Numerov coefficients
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
    dx_R = x_right[1] - x_right[0] # Negative for stepping backwards

    factor_L = (dx_L**2) / 12.0
    factor_R = (dx_R**2) / 12.0

    # ---- Forward Integration (x_l to x_m) ----
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
            
            # THE FIX: Corrected signs
            c0, c1, c2 = 1.0 - factor_L*k_next, 2.0 + 10.0*factor_L*k_curr, 1.0 - factor_L*k_prev
            u_next = (c1 * u_up - c2 * u_up_prev) / c0
            
            u_up_prev, u_up = u_up.copy(), u_next
            v_up = (u_up - u_up_prev) / dx_L 

    # ---- Backward Integration (x_r to x_m) ----
    u_down = np.zeros(num_E)
    if method == 'rk2':
        v_down = np.full(num_E, v_r)
    elif method == 'numerov':
        u_down_prev = np.zeros(num_E)
        # THE FIX: Removed the incorrect minus sign here (dx_R is already negative)
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
            
            # THE FIX: Corrected signs
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
    print("Eigen-energies:",roots)
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

def visualise_eigenstates(roots, V_func, x_l, x_r, dx=0.001, condition="", 
                          fig=None, ax=None, final=True, method='numerov'):
    """
    Runs the RK2 solver for each eigen-energy (root) and plots the normalised wavefunctions.
    """
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
        
    for E in roots:
        # Run the RK2 solver from x_l to x_r
        if method=='rk2':
            x_vals, u, v = rk2_solve(x_l, x_r, dx, V_func, E)
        elif method=='numerov':
            x_vals, u = numerov_solve(x_l, x_r, dx, V_func, E)
        else:raise ValueError(f"unknown method {method}")
        
        # Normalise the wave-function (L-infinity norm) so they stack nicely on the same plot
        u_norm = u / np.max(np.abs(u))
        
        ax.plot(x_vals, u_norm, label=f"$E={E:.4f}$ {condition}")
        
    if final:
        ax.axhline(0, color='black', linestyle='--', alpha=0.7)
        ax.set_xlabel('Coordinate ($r$)')
        ax.set_ylabel('Normalised Wavefunction $u(r)$')
        
        if not ax.get_title():
            ax.set_title('Eigenstates $u(r)$')
            
        ax.grid(True, alpha=0.3)
        ax.legend()
        
    return fig, ax

def visualise_eigenstates_stitched(roots, V_func, x_l, x_r, dx=0.001, condition="", 
                                   fig=None, ax=None, final=True,method='rk2'):
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
        
    N_total = int(np.abs(x_r - x_l) / dx) + 1
    x_vals = np.linspace(x_l, x_r, N_total)
    
    # Pre-calculate the potential across the grid
    V_arr = V_func(x_vals)
    
    for E in roots:
        # --- DYNAMIC MATCHING POINT ---
        # Find the classical turning point (where E > V). 
        # We take the last index where this is true (the outer turning point).
        allowed_indices = np.where(E > V_arr)[0]
        if len(allowed_indices) > 0:
            match_idx = allowed_indices[-1] 
        else:
            match_idx = N_total // 2 # Fallback if no classical region is found
            
        x_m = x_vals[match_idx]
        
        if method=='rk2':
            # 1. Forward integration (x_l to x_m)
            x_fwd, u_fwd, _ = rk2_solve(x_l, x_m, dx, V_func, E)
            # 2. Backward integration (x_r to x_m)
            x_bwd, u_bwd, _ = rk2_solve(x_r, x_m, -dx, V_func, E)
        elif method=='numerov':
            # 1. Forward integration (x_l to x_m)
            x_fwd, u_fwd = numerov_solve(x_l, x_m, dx, V_func, E)
            # 2. Backward integration (x_r to x_m)
            x_bwd, u_bwd = numerov_solve(x_r, x_m, -dx, V_func, E)
        else:raise ValueError(f"unknown method {method}")
        
        x_bwd = x_bwd[::-1]
        u_bwd = u_bwd[::-1]
        
        # 3. Stitch them together
        scale_factor = u_fwd[-1] / u_bwd[0]
        u_bwd_scaled = u_bwd * scale_factor
        
        x_stitched = np.concatenate((x_fwd, x_bwd[1:]))
        u_stitched = np.concatenate((u_fwd, u_bwd_scaled[1:]))
        u_norm = u_stitched / np.max(np.abs(u_stitched))
    
        ax.plot(x_stitched, u_norm, label=f"$E={E:.4f}$ {condition}")
        
    if final:
        ax.axhline(0, color='black', linestyle='--', alpha=0.7)
        ax.set_xlabel('Coordinate ($r$)')
        ax.set_ylabel('Normalised Wavefunction $u(r)$')
        
        if not ax.get_title():
            ax.set_title('Eigenstates $u(r)$ (Stitched at Turning Point)')
            
        ax.grid(True, alpha=0.3)
        ax.legend()
        
    return fig, ax

if __name__ == "__main__":
    
    # 1. Define the parameters
    Z, mu, lamb = 1.0, 0.0, 0.0
    r_min, r_max, dr = 0.001, 100.0, 0.001
    
    # 2. Get the roots for l=0
    l_0 = 0
    roots_0, E_0, f_E_0 = solve_radial_atom(l=l_0, mu=mu, lamb=lamb, Z=Z, r_max=r_max, dr=dr, num_E=1000)
    
    # Define the potential function for l=0
    V_0 = lambda r: electrostatic_potential_eff(r, l=l_0, mu=mu, lamb=lamb, Z=Z)
    
    # 3. Get the roots for l=1
    l_1 = 1
    roots_1, E_1, f_E_1 = solve_radial_atom(l=l_1, mu=mu, lamb=lamb, Z=Z, r_max=r_max, dr=dr, num_E=1000)
    
    # Define the potential function for l=1
    V_1 = lambda r: electrostatic_potential_eff(r, l=l_1, mu=mu, lamb=lamb, Z=Z)
    
    # 4. Visualise the wavefunctions
    # We only plot the first 2 roots of each to avoid overcrowding the plot
    fig, ax = visualise_eigenstates(roots_0[:2], V_0, x_l=r_max, x_r=r_min, dx=-dr, condition="($l=0$)", final=False)
    fig, ax = visualise_eigenstates(roots_1[:2], V_1, x_l=r_max, x_r=r_min, dx=-dr, condition="($l=1$)", fig=fig, ax=ax, final=True)
    
    # (Optional: zooming in to see the wave structure since the tails are long)
    # ax.set_xlim(0, 30) 
    
    plt.show()

    # # 1. Define the parameters for pure Hydrogen
    # Z, mu, lamb = 1.0, 0.0, 0.0
    # r_min, r_max, dr = 0.001, 100.0, 0.001
    
    # # 2. Solve for l=0 
    # # (Assuming solve_radial_atom returns: roots, E_vec, f_E)
    # roots_0, E_0, f_E_0 = solve_radial_atom(l=0, mu=mu, lamb=lamb, Z=Z, r_max=r_max, dr=dr, num_E=1000)
    
    # # 3. Visualise the newly stabilised f(E)
    # fig, ax = plt.subplots(figsize=(10, 6))
    # fig, ax = visualise_f_E(roots_0, E_0, f_E_0, condition="($l=0$, Dynamic Turning Point)", fig=fig, ax=ax, final=True)
    
    # # Optional: Plot l=1 on the same graph to see how the mismatch functions compare
    # roots_1, E_1, f_E_1 = solve_radial_atom(l=1, mu=mu, lamb=lamb, Z=Z, r_max=r_max, dr=dr, num_E=1000)
    # fig, ax = visualise_f_E(roots_1, E_1, f_E_1, condition="($l=1$, Dynamic Turning Point)", fig=fig, ax=ax, final=True)
    
    # plt.title("Stabilised Mismatch Function using Dynamic Turning Points")
    # plt.show()