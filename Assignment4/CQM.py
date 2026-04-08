import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from scipy.linalg import solve_banded
from IPython.display import HTML
import matplotlib as mpl
mpl.rcParams['animation.embed_limit'] = 50.0 # Set limit to 50 MB

def init_segment(dx, x_left, x_right, psi_func):
    """
    Initializes the spatial grid and the wave function.
    Using np.arange ensures consistent step sizes matching dx.
    """
    x = np.arange(x_left, x_right + dx, dx)
    psi = psi_func(x)
    return x, psi.astype(complex)

def gaussian_potential(x, amplitude=1.0, center=0.0, width=1.0):
    """Returns a Gaussian potential barrier/well."""
    return amplitude * np.exp(-((x - center) ** 2) / (2 * width ** 2))

def step_potential(x, step_pos=0.0, v_left=0.0, v_right=1.0):
    """Returns a step potential using numpy's vectorized where()."""
    return np.where(x < step_pos, v_left, v_right)

def harmonic_potential(x, omega=1.0):
    return 0.5 * omega**2 * x**2

def harmonic_force(x, omega=1.0):
    return -omega**2 * x

def linear_potential(x, F=1.0):
    return F * x

def linear_force(x, F=1.0):
    return -F * np.ones_like(x)

def compute_a(dx, dt, V_x):
    """Computes the main diagonal array for the A matrix (implicit step)."""
    return ( 4j * (dx**2 / dt) ) - ( 2 * (dx**2) * V_x ) - 2

def compute_b(dx, dt, V_x):
    """Computes the main diagonal array for the B matrix (explicit step)."""
    return ( 4j * (dx**2 / dt) ) + ( 2 * (dx**2) * V_x ) + 2

def compute_matrix_A_compact(a_x):
    """
    Constructs the A matrix in compact form for scipy.linalg.solve_banded.
    Row 0: Upper diagonal (padded with 0 at index 0)
    Row 1: Main diagonal a(x)
    Row 2: Lower diagonal (padded with 0 at index N-1)
    """
    N = len(a_x)
    A_banded = np.zeros((3, N), dtype=complex)
    A_banded[0, 1:] = 1.0   # Upper diagonal
    A_banded[1, :] = a_x    # Main diagonal
    A_banded[2, :-1] = 1.0  # Lower diagonal
    return A_banded

def compute_matrix_B_compact(b_x):
    """
    Constructs the B matrix in compact form.
    Note: When computing B * psi_0, you usually just do standard matrix 
    multiplication, but keeping it banded can help if you write a custom sparse dot product.
    """
    N = len(b_x)
    B_banded = np.zeros((3, N), dtype=complex)
    B_banded[0, 1:] = -1.0   # Upper diagonal
    B_banded[1, :] = b_x     # Main diagonal
    B_banded[2, :-1] = -1.0  # Lower diagonal
    return B_banded

def multiply_B_psi(B_banded, psi):
    """
    Multiplies a compact tridiagonal matrix B (3xN) with a vector psi (N,).
    Uses fast numpy array slicing to avoid explicit loops.
    """
    
    # 1. Main diagonal contribution
    result = B_banded[1, :] * psi
    
    # 2. Upper diagonal contribution
    # B_{i, i+1} * psi_{i+1} -> affects result[0] to result[N-2]
    # Note: B_banded[0, 0] is padding, the upper diagonal starts at index 1
    result[:-1] += B_banded[0, 1:] * psi[1:]
    
    # 3. Lower diagonal contribution
    # B_{i, i-1} * psi_{i-1} -> affects result[1] to result[N-1]
    # Note: B_banded[2, N-1] is padding, the lower diagonal ends at index N-2
    result[1:] += B_banded[2, :-1] * psi[:-1]
    
    return result

def test_1():
    # Parameters
    dx = 0.1
    dt = 0.05
    x_left = -10.0
    x_right = 10.0

    # Define an initial Gaussian wave packet (real-valued for simplicity here)
    def initial_psi(x):
        return np.exp(-(x + 3)**2 / 2) + 0j

    # Init grid
    x, psi = init_segment(dx, x_left, x_right, initial_psi)

    # Compute potentials
    V_gauss = gaussian_potential(x, amplitude=5.0, center=2.0, width=1.5)
    V_step = step_potential(x, step_pos=0.0, v_left=0.0, v_right=3.0)

    # Compute coefficients
    a_gauss = compute_a(dx, dt, V_gauss)
    b_gauss = compute_b(dx, dt, V_gauss)
    a_step = compute_a(dx, dt, V_step)
    b_step = compute_b(dx, dt, V_step)

    # Compute matrices
    A_gauss = compute_matrix_A_compact(a_gauss)
    B_gauss = compute_matrix_B_compact(b_gauss)
    A_step = compute_matrix_A_compact(a_step)
    B_step = compute_matrix_B_compact(b_step)

    # Output sanity checks
    print("=== Testing Gaussian Potential Matrices ===")
    print("A_compact (first 4 columns):\n", np.round(A_gauss[:, :4], 2))
    print("\nB_compact (first 4 columns):\n", np.round(B_gauss[:, :4], 2))

    print("\n=== Testing Step Potential Matrices ===")
    print("A_compact (first 4 columns):\n", np.round(A_step[:, :4], 2))
    print("\nB_compact (first 4 columns):\n", np.round(B_step[:, :4], 2))

    # --- NEW: Test multiply_B_psi ---
    print("\n=== Testing multiply_B_psi function ===")
    
    # Construct a standard dense NxN matrix from the compact B_gauss matrix using np.diag
    B_dense = (np.diag(B_gauss[1, :], k=0) +       # Main diagonal
               np.diag(B_gauss[0, 1:], k=1) +      # Upper diagonal
               np.diag(B_gauss[2, :-1], k=-1))     # Lower diagonal
    
    # Standard matrix-vector multiplication
    expected_result = B_dense @ psi
    
    # Optimized multiplication
    optimized_result = multiply_B_psi(B_gauss, psi)
    
    # Check if they are numerically equivalent
    is_correct = np.allclose(expected_result, optimized_result)
    print(f"Optimized multiply_B_psi matches dense dot product: {is_correct}")
    if not is_correct:
        print("Max absolute difference:", np.max(np.abs(expected_result - optimized_result)))

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(x, np.abs(psi)**2, label="Initial $|\psi(x)|^2$", color="blue")
    plt.plot(x, V_gauss, label="Gaussian $V(x)$", linestyle="--", color="red")
    plt.plot(x, V_step, label="Step $V(x)$", linestyle=":", color="green")
    
    plt.title("Test 1: Initial Wave Packet and Potentials")
    plt.xlabel("Position (x)")
    plt.ylabel("Amplitude / Potential Energy")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

def psi_to_phi(x, psi):
    """
    Converts position wave function psi(x) to wave-number wave function phi(k).
    Applies the 1/sqrt(2*pi) normalization and phase shift for discrete grids.
    """
    dx = x[1] - x[0]
    N = len(x)
    
    # Generate frequencies and shift the zero-frequency component to the center
    k = np.fft.fftshift(np.fft.fftfreq(N, d=dx)) * 2 * np.pi
    
    # Compute the discrete Fourier transform
    phi_unshifted = np.fft.fft(psi) * dx / np.sqrt(2 * np.pi)
    
    # Shift array and apply phase correction for the grid's starting offset
    phi = np.fft.fftshift(phi_unshifted) * np.exp(-1j * k * x[0])
    
    return k, phi

def phi_to_psi(k, phi, x):
    """
    Converts wave-number wave function phi(k) back to position wave function psi(x).
    """
    dx = x[1] - x[0]
    
    # Remove phase correction
    phi_shifted = phi * np.exp(1j * k * x[0])
    
    # Undo the shift and compute the Inverse FFT
    phi_unshifted = np.fft.ifftshift(phi_shifted)
    psi = np.fft.ifft(phi_unshifted) / dx * np.sqrt(2 * np.pi)
    
    return psi

def update_psi(psi, A_banded, B_banded, method="split-operator", dt=0.01):
    """
    Advances the wave function by a single time step.
    Currently supports the implicit split-operator method via Crank-Nicolson.
    """
    if method == "split-operator":
        # 1. Compute RHS = B * psi_0
        RHS = multiply_B_psi(B_banded, psi)
        
        # 2. Solve A * psi_t = RHS
        # solve_banded takes (number of lower diagonals, number of upper diagonals)
        psi_next = solve_banded((1, 1), A_banded, RHS)
        
        return psi_next
    elif method == "SDLF":
        return step_sdlf(psi, A_banded, dt)
    else:
        raise ValueError(f"Method '{method}' is not implemented.")

def simulate_1D(
    V_func, psi_0_func, dx, dt, x_left, x_right, frames=100, methods=["split-operator"],
    F_func=None, x0_classical=None, p0_classical=None, show_animation=True, output_format='widget',
    animation_frame_scaling=1
    ):
    """
    Simulates the TDSE and returns a dictionary with history metrics and figures.
    """
    # 1. Initialize grid and potentials
    x, psi_initial = init_segment(dx, x_left, x_right, psi_0_func)
    V_x = V_func(x)
    
    # 2. Precompute matrices for the methods
    matrices = {}
    for method in methods:
        if method == "split-operator":
            a_x = compute_a(dx, dt, V_x)
            b_x = compute_b(dx, dt, V_x)
            matrices[method] = {
                'A': compute_matrix_A_compact(a_x),
                'B': compute_matrix_B_compact(b_x),
                'psi': psi_initial.copy()
            }
        elif method == "SDLF":
            matrices[method] = {
                "A":compute_F_matrix_compact(dx, V_x),
                "B":None,
                'psi':psi_initial.copy()
            }

    # --- UPDATED: Metrics Helper Function ---
    def get_metrics(psi, phi, k_arr, x_arr):
        dk = k_arr[1] - k_arr[0]
        norm_x = np.sum(np.abs(psi)**2) * dx
        norm_k = np.sum(np.abs(phi)**2) * dk
        exp_x = np.sum(x_arr * np.abs(psi)**2) * dx
        exp_p = np.sum(k_arr * np.abs(phi)**2) * dk
        exp_H = np.sum(0.5 * k_arr**2 * np.abs(phi)**2) * dk + np.sum(V_x * np.abs(psi)**2) * dx
        
        # Most probable values
        x_p = x_arr[np.argmax(np.abs(psi)**2)]
        k_p = k_arr[np.argmax(np.abs(phi)**2)]
        return norm_x, norm_k, exp_x, exp_p, exp_H, x_p, k_p

    primary_method = methods[0]
    k_init, phi_init = psi_to_phi(x, matrices[primary_method]['psi'])
    norm_x_0, norm_k_0, exp_x_0, exp_p_0, exp_H_0, x_p_0, k_p_0 = get_metrics(matrices[primary_method]['psi'], phi_init, k_init, x)

    steps = frames * animation_frame_scaling

    # --- NEW: History and Classical Init ---
    history = {
        't': np.zeros(steps), 'norm_x': np.zeros(steps), 'norm_k': np.zeros(steps),
        'exp_x': np.zeros(steps), 'exp_p': np.zeros(steps), 'exp_H': np.zeros(steps),
        'x_p': np.zeros(steps), 'k_p': np.zeros(steps), 'x_c': np.zeros(steps)
    }
    
    classical_state = {
        'x': x0_classical if x0_classical is not None else exp_x_0,
        'p': p0_classical if p0_classical is not None else exp_p_0
    }

    # 3. Setup Figure and Axes for Animation
    fig = plt.figure(figsize=(8, 8))
    row1, row2 = fig.subplots(2, 2)
    ax1, ax3 = row1
    ax2, ax4 = row2
    
    # Dictionary to store line objects for each method
    lines_psi = {}
    lines_phi = {}
    
    for method in methods:
        line_psi, = ax1.plot(x, np.abs(matrices[method]['psi'])**2, label=f"|$\psi$|² ({method})")
        lines_psi[method] = line_psi
        
        k, phi = psi_to_phi(x, matrices[method]['psi'])
        line_phi, = ax3.plot(k, np.abs(phi)**2, label=f"|$\phi$|² ({method})")
        lines_phi[method] = line_phi

    if show_animation:
        # Setup ax1 (Position Space Wavefunction)
        ax1.set_title("Position Wave Function $|\psi(x,t)|^2$")
        ax1.set_ylabel("Probability Density")
        ax1.set_ylim(0, np.max(np.abs(psi_initial)**2) * 1.5)
        vline_x = ax1.axvline(exp_x_0, color='purple', linestyle='--', label=r'$\langle x \rangle$')
        vline_xp = ax1.axvline(x_p_0, color='red', linestyle=':', label=r'$x_p$') # <-- ADDED
        ax1.legend(loc="upper right")
        ax1.grid(True, alpha=0.3)

        # Setup ax2 (Potential Energy)
        ax2.plot(x, V_x, color="black", label="V(x)")
        ax2.set_title("Potential Energy $V(x)$")
        ax2.set_ylabel("Energy")
        ax2.legend(loc="upper right")
        ax2.grid(True, alpha=0.3)
        
        # Setup ax3 (Momentum Space Wavefunction)
        ax3.set_title("Wave-Number Wave Function $|\phi(k,t)|^2$")
        ax3.set_xlabel("Wave number (k)")
        # ax3.set_ylabel("Probability Density")
        ax3.set_xlim(-15, 15)
        ax3.set_ylim(0, np.max(np.abs(phi_init)**2) * 1.5)
        vline_p = ax3.axvline(exp_p_0, color='purple', linestyle='--', label=r'$\langle p \rangle$')
        vline_kp = ax3.axvline(k_p_0, color='red', linestyle=':', label=r'$k_p$') # <-- ADDED
        ax3.legend(loc="upper right")
        ax3.grid(True, alpha=0.3)

        # Setup ax4 (Bar Chart for Invariants)
        # ax4.set_title("Percentage Change of Invariants")
        bar_labels = [r'$\int |\psi|^2 dx$', r'$\int |\phi|^2 dk$', r'$\langle H \rangle$']
        bars = ax4.bar(bar_labels, [0, 0, 0], color=['blue', 'orange', 'green'])
        ax4.set_ylim(-1, 1) # Small initial bounds, will scale dynamically
        ax4.set_ylabel("% Change")
        ax4.yaxis.set_label_position('right')
        ax4.yaxis.tick_right()
        ax4.grid(axis='y', alpha=0.3)

    # 4. Computation and Update Step
    def compute_step(frame):
        # 1. Advance Quantum Time Step
        if frame > 0:
            for method in methods:
                matrices[method]['psi'] = update_psi(
                    matrices[method]['psi'], matrices[method]['A'], matrices[method]['B'], 
                    dt=dt, method=method)
        
        psi_curr = matrices[primary_method]['psi']
        k_curr, phi_curr = psi_to_phi(x, psi_curr)
        
        # 2. Extract Metrics
        n_x, n_k, e_x, e_p, e_H, x_p, k_p = get_metrics(psi_curr, phi_curr, k_curr, x)
        
        # 3. Classical Position Verlet Step (Atomic Units: m=1)
        x_c = classical_state['x']
        p_c = classical_state['p']
        if frame > 0 and F_func is not None:
            x_mid = x_c + 0.5 * p_c * dt
            a_mid = F_func(x_mid)
            p_next = p_c + a_mid * dt
            x_next = x_mid + 0.5 * p_next * dt
            
            classical_state['x'] = x_next
            classical_state['p'] = p_next
            x_c, p_c = x_next, p_next

        # 4. Log to History
        history['t'][frame] = frame * dt
        history['norm_x'][frame], history['norm_k'][frame] = n_x, n_k
        history['exp_x'][frame], history['exp_p'][frame], history['exp_H'][frame] = e_x, e_p, e_H
        history['x_p'][frame], history['k_p'][frame] = x_p, k_p
        history['x_c'][frame] = x_c
        
        return psi_curr, phi_curr, n_x, n_k, e_x, e_p, e_H, x_p, k_p

    def animate(frame):
        print(f"\rProcessing frame {frame + 1}/{frames}...", end="", flush=True)
        for step in range(animation_frame_scaling-1):compute_step(frame*animation_frame_scaling + step)
        psi_curr, phi_curr, n_x, n_k, e_x, e_p, e_H, x_p, k_p = compute_step((frame+1)*animation_frame_scaling -1)
        
        # Update plotting data
        for method in methods:
            lines_psi[method].set_ydata(np.abs(matrices[method]['psi'])**2)
            k_arr, current_phi = psi_to_phi(x, matrices[method]['psi'])
            lines_phi[method].set_ydata(np.abs(current_phi)**2)
            
        vline_x.set_xdata([e_x, e_x])
        vline_xp.set_xdata([x_p, x_p])
        vline_p.set_xdata([e_p, e_p])
        vline_kp.set_xdata([k_p, k_p])
        
        # Update bar chart
        pct_n_x = (n_x - norm_x_0) / norm_x_0 * 100
        pct_n_k = (n_k - norm_k_0) / norm_k_0 * 100
        pct_H = (e_H - exp_H_0) / exp_H_0 * 100 if exp_H_0 != 0 else 0
        pct_changes = [pct_n_x, pct_n_k, pct_H]
        for bar, h in zip(bars, pct_changes):
            bar.set_height(h)
            
        max_abs = max(np.abs(pct_changes))
        if max_abs > ax4.get_ylim()[1] or max_abs == 0:
            scale = max_abs * 1.5 if max_abs > 0 else 1.0
            ax4.set_ylim(-scale, scale)
            
        fig.suptitle(f"Frame: {frame} | $\\langle x \\rangle$: {e_x:.2f} | $\\langle p \\rangle$: {e_p:.2f}", fontsize=12)
        return list(lines_psi.values()) + list(lines_phi.values()) + [vline_x, vline_xp, vline_p, vline_kp] + list(bars)

    # 5. Execute Simulation
    anim_widget = None
    ani = None
    if show_animation:
        print("\nSimulation started (with animation)!")
        ani = FuncAnimation(fig, animate, frames=frames, interval=50, blit=True)
        if output_format == 'widget':
            plt.close(fig)
            anim_widget = HTML(ani.to_jshtml())
        elif output_format == 'video':
            plt.close(fig)
            anim_widget = HTML(ani.to_html5_video())
    else:
        print("\nRunning fast simulation (no animation)...")
        for frame in range(frames):
            for step in range(animation_frame_scaling):
                compute_step(frame*animation_frame_scaling + step)
        print("Done!")

    # Return Output Dictionary
    return {
        'history': history,
        'animation': anim_widget,
        'raw_ani': ani if show_animation else None,
        'final_psi': matrices[primary_method]['psi'],
        'x_grid': x
    }

def plot_static_results(history, F_func=None):
    """Generates the static plots from a populated history dictionary."""
    fig_static, (ax_s1, ax_s2) = plt.subplots(2, 1, figsize=(8, 8))
    
    # Trajectories
    ax_s1.plot(history['t'], history['exp_x'], label=r'$\langle x \rangle$ (Quantum)', color='blue')
    ax_s1.plot(history['t'], history['x_p'], label=r'$x_p$', color='green', linestyle='--')
    if F_func is not None:
        ax_s1.plot(history['t'], history['x_c'], label=r'$x_c$ (Classical)', color='red', linestyle='-.')
    ax_s1.set_title("Position Trajectories over Time")
    ax_s1.set_xlabel("Time (t)")
    ax_s1.set_ylabel("Position (x)")
    ax_s1.legend()
    ax_s1.grid(True, alpha=0.3)

    # Normalization Error
    ax_s2.plot(history['t'], history['norm_x'], label=r'$N(t)$')
    ax_s2.set_title("Wavefunction Normalization over Time")
    ax_s2.set_xlabel("Time (t)")
    ax_s2.set_ylabel("Norm")
    ax_s2.grid(True, alpha=0.3)
    plt.tight_layout()
    
    return fig_static

def test_1D():
    # Setup for Harmonic Oscillator (Assignment parameters)
    def initial_packet(x):
        sigma = 0.7
        x0 = -3.0
        # No initial momentum, purely real Gaussian
        return (1 / (sigma * np.sqrt(np.pi)))**0.5 * np.exp(-(x - x0)**2 / (2 * sigma**2)) + 0j

    # Run Simulation
    res = simulate_1D(
        V_func=harmonic_potential, 
        F_func=harmonic_force,     # <--- Passes the Force to enable classical tracking
        psi_0_func=initial_packet, 
        dx=0.1, 
        dt=0.001, 
        x_left=-10.0, 
        x_right=10.0, 
        frames=100,
        animation_frame_scaling=200,
        # methods=["split-operator"],
        methods=["SDLF","split-operator"],
        show_animation=True,
        output_format='save'       # Avoids compiling the heavy JS widget in terminal
    )

    if res['raw_ani'] is not None:
        writer = PillowWriter(fps=10)
        res['raw_ani'].save("harmonic_simulation.gif", writer=writer)
        print("\nGIF saved successfully!")

    # 3. Generate and Save the Static Figure NOW
    fig_static = plot_static_results(res['history'], F_func=harmonic_force)
    fig_static.savefig("static_results.png", dpi=300, bbox_inches='tight')
    print("Saved static plot as 'static_results.png'")

    plt.close('all')

def compute_F_matrix_compact(dx, V_x):
    """
    Constructs the spatial finite-difference matrix F in a 3xN compact banded format.
    F = -H, where H is the Hamiltonian in Atomic Units.
    """
    N = len(V_x)
    F = np.zeros((3, N), dtype=float)
    
    # Middle row (index 1): Main diagonal = -2 * (1 + V_x * dx^2)
    F[1, :] = -2.0 * (1.0 + V_x * dx**2)
    
    # Top row (index 0): Upper diagonal (shifted for scipy/banded convention)
    F[0, 1:] = 1.0
    
    # Bottom row (index 2): Lower diagonal (shifted for scipy/banded convention)
    F[2, :-1] = 1.0
    
    # Scale entire matrix by 1 / (2 * dx^2)
    F *= 1.0 / (2.0 * dx**2)
    
    return F

def step_sdlf(psi, F_banded, dt):
    """
    Advances the wave function by one time step 'dt' using the 
    Space Discretised Leap Frog (Position-Verlet) method.
    Assumes `multiply_B_psi` is available in the global scope.
    """
    R = np.real(psi)
    I = np.imag(psi)
    
    # Step 1: Advance Imaginary part by half a step
    # dI/dt = F * R
    F_R = multiply_B_psi(F_banded, R)
    I_half = I + (dt / 2.0) * F_R
    
    # Step 2: Advance Real part by a full step using the half-step Imaginary part
    # dR/dt = -F * I
    F_I_half = multiply_B_psi(F_banded, I_half)
    R_next = R - dt * F_I_half
    
    # Step 3: Advance Imaginary part by the remaining half step using the new Real part
    F_R_next = multiply_B_psi(F_banded, R_next)
    I_next = I_half + (dt / 2.0) * F_R_next
    
    # Combine back into complex wave function
    return R_next + 1j * I_next


#####################################################################################
################## 2D case ##########################################################
#####################################################################################

def init_grid_2d_old(dx, x_left, x_right, y_down, y_up):
    """
    Initializes the 2D spatial grid.
    Uses a single dx for both dimensions, assuming delta y = delta x.
    """
    x = np.arange(x_left, x_right + dx, dx)
    y = np.arange(y_down, y_up + dx, dx)
    
    # Create 2D meshgrid matrices
    X, Y = np.meshgrid(x, y)
    
    return x, y, X, Y

def init_grid_2d(dx, x_left, x_right, y_down, y_up):
    """
    Initializes the 2D spatial grid.
    Uses np.linspace to avoid floating-point accumulation errors and ensure symmetry.
    """
    # Calculate the exact number of points
    Nx = int(round((x_right - x_left) / dx)) + 1
    Ny = int(round((y_up - y_down) / dx)) + 1
    
    # Generate exactly symmetric arrays
    x = np.linspace(x_left, x_right, Nx)
    y = np.linspace(y_down, y_up, Ny)
    
    # Create 2D meshgrid matrices
    X, Y = np.meshgrid(x, y)
    
    return x, y, X, Y

def init_moving_gaussian_2d(X, Y, x0=-8.0, y0=0.0, sigma=1.0, kx=2.0, ky=0.0):
    """
    Initializes a 2D Gaussian wave packet with initial momentum.
    
    Parameters:
    - X, Y: Spatial meshgrid matrices
    - x0, y0: Initial center of the wave packet
    - sigma: Spatial spread (standard deviation)
    - kx, ky: Initial wave numbers (momentum in atomic units)
    """
    # 2D Normalization constant so that integral of |psi|^2 dx dy = 1
    norm = 1.0 / (sigma * np.sqrt(np.pi))
    
    # Real Gaussian envelope
    envelope = np.exp(-((X - x0)**2 + (Y - y0)**2) / (2 * sigma**2))
    
    # Complex phase for momentum p = hbar * k (in atomic units, hbar = 1)
    phase = np.exp(1j * (kx * X + ky * Y))
    
    # Final wave function
    psi_0 = norm * envelope * phase
    
    return psi_0.astype(complex)

def gaussian_potential_2d(X, Y, amplitude=1.0, center_x=0.0, center_y=0.0, width_x=1.0, width_y=1.0):
    """Returns a 2D Gaussian potential barrier/well."""
    return amplitude * np.exp(-((X - center_x)**2 / (2 * width_x**2) + (Y - center_y)**2 / (2 * width_y**2)))

def step_potential_2d(X, Y, step_pos_x=0.0, v_left=0.0, v_right=1.0):
    """
    Returns a 2D step potential along the x-axis. 
    (Can be easily modified to step along y or an arbitrary line).
    """
    return np.where(X < step_pos_x, v_left, v_right)

def disk_potential_2d(X, Y, R=1.5, center_x=0.0, center_y=0.0, v_out=0.0, v_in=1.0):
    """
    Returns a 2D disk potential with radius R 
    """
    epsilon = 1e-10
    return np.where(((X-center_x)**2 + (Y-center_y)**2) < R**2 + epsilon, v_out, v_in)


def harmonic_potential_2d(X, Y, omega_x=1.0, omega_y=1.0):
    """Returns a 2D harmonic oscillator potential."""
    return 0.5 * (omega_x**2 * X**2 + omega_y**2 * Y**2)

def linear_potential_2d(X, Y, F_x=1.0, F_y=1.0):
    """Returns a 2D linear potential (constant force)."""
    return F_x * X + F_y * Y

def compute_a_x(dx, dt, V_2d):
    """
    Computes the main diagonal array for the A_x matrices (implicit x-step).
    Uses dt/2 for the time step and V(x,y)/2 for the potential split.
    """
    return (8j * (dx**2 / dt)) - (dx**2 * V_2d) - 2

def compute_b_x(dx, dt, V_2d):
    """
    Computes the main diagonal array for the B_x matrices (explicit x-step).
    """
    return (8j * (dx**2 / dt)) + (dx**2 * V_2d) + 2

def compute_a_y(dy, dt, V_2d):
    return (8j * (dy**2 / dt)) - (dy**2 * V_2d) - 2

def compute_b_y(dy, dt, V_2d):
    return (8j * (dy**2 / dt)) + (dy**2 * V_2d) + 2

def compute_matrices_A_x(a_x_2d):
    """
    Constructs the A_x matrices for all rows (y-coordinates).
    Returns a 3D array of shape (Ny, 3, Nx) where A_x[y] is the compact 
    banded matrix for the y-th row.
    """
    Ny, Nx = a_x_2d.shape
    A_x = np.zeros((Ny, 3, Nx), dtype=complex)
    A_x[:, 0, 1:] = 1.0     # Upper diagonal (padded at index 0)
    A_x[:, 1, :] = a_x_2d   # Main diagonal
    A_x[:, 2, :-1] = 1.0    # Lower diagonal (padded at index Nx-1)
    return A_x

def compute_matrices_B_x(b_x_2d):
    """
    Constructs the B_x matrices for all rows (y-coordinates).
    Returns a 3D array of shape (Ny, 3, Nx) where B_x[y] is the compact
    banded matrix for the y-th row.
    """
    Ny, Nx = b_x_2d.shape
    B_x = np.zeros((Ny, 3, Nx), dtype=complex)
    B_x[:, 0, 1:] = -1.0    
    B_x[:, 1, :] = b_x_2d   
    B_x[:, 2, :-1] = -1.0   
    return B_x

def compute_matrices_A_y(a_y_2d):
    """
    Constructs the A_y matrices for all columns (x-coordinates).
    Returns a 3D array of shape (Nx, 3, Ny) where A_y[x] is the compact 
    banded matrix for the x-th column.
    Note: We transpose a_y_2d so the y-elements align as rows for solve_banded.
    """
    Ny, Nx = a_y_2d.shape
    A_y = np.zeros((Nx, 3, Ny), dtype=complex)
    A_y[:, 0, 1:] = 1.0     
    A_y[:, 1, :] = a_y_2d.T # Transpose to iterate over columns
    A_y[:, 2, :-1] = 1.0    
    return A_y

def compute_matrices_B_y(b_y_2d):
    """
    Constructs the B_y matrices for all columns (x-coordinates).
    Returns a 3D array of shape (Nx, 3, Ny) where B_y[x] is the compact
    banded matrix for the x-th column.
    """
    Ny, Nx = b_y_2d.shape
    B_y = np.zeros((Nx, 3, Ny), dtype=complex)
    B_y[:, 0, 1:] = -1.0    
    B_y[:, 1, :] = b_y_2d.T 
    B_y[:, 2, :-1] = -1.0   
    return B_y

def multiply_banded_matrix_vector(B_banded, vec):
    """
    Multiplies a compact tridiagonal matrix (shape 3 x N) by a vector of length N.
    B_banded[0] is the upper diagonal.
    B_banded[1] is the main diagonal.
    B_banded[2] is the lower diagonal.
    """
    res = np.zeros_like(vec, dtype=complex)
    
    # Main diagonal
    res += B_banded[1] * vec 
    
    # Upper diagonal (shifts vector left)
    res[:-1] += B_banded[0, 1:] * vec[1:] 
    
    # Lower diagonal (shifts vector right)
    res[1:] += B_banded[2, :-1] * vec[:-1]
    
    return res

def step_adi_2d(psi_0, A_x, B_x, A_y, B_y):
    """
    Performs a single time-step evolution of the 2D TDSE using the ADI method.
    
    Parameters:
    - psi_0: 2D array of the wave function at time t, shape (Ny, Nx), mapped as Psi[y,x].
    - A_x, B_x: 3D arrays of compact banded matrices for the x-sweep, shape (Ny, 3, Nx).
    - A_y, B_y: 3D arrays of compact banded matrices for the y-sweep, shape (Nx, 3, Ny).
    
    Returns:
    - psi_t: 2D array of the wave function at time t + dt.
    """
    Ny, Nx = psi_0.shape
    
    # Temporary array to hold the half-step wave function Psi_{t/2}
    psi_half = np.zeros_like(psi_0, dtype=complex)
    
    # Final array to hold the full-step wave function Psi_t
    psi_t = np.zeros_like(psi_0, dtype=complex)

    # ---------------------------------------------------------
    # Step 1: X-Sweep
    # Solving A_x(y) * Psi_{t/2}[y] = B_x(y) * Psi_0[y] for all y
    # ---------------------------------------------------------
    for y in range(Ny):
        # Extract the y-th row (Psi_0[y])
        row_vec_0 = psi_0[y, :]
        
        # Calculate RHS: B_x(y) * Psi_0[y]
        rhs_x = multiply_banded_matrix_vector(B_x[y], row_vec_0)
        
        # Solve for the half-step row (Psi_{t/2}[y])
        # (1, 1) indicates 1 lower diagonal and 1 upper diagonal
        psi_half[y, :] = solve_banded((1, 1), A_x[y], rhs_x)

    # ---------------------------------------------------------
    # Step 2: Y-Sweep
    # Solving A_y(x) * Psi_t^T[x] = B_y(x) * Psi_{t/2}^T[x] for all x
    # ---------------------------------------------------------
    for x in range(Nx):
        # Extract the x-th column. 
        # Slicing [:, x] from the matrix naturally accesses Psi^T[x]
        col_vec_half = psi_half[:, x]
        
        # Calculate RHS: B_y(x) * Psi_{t/2}^T[x]
        rhs_y = multiply_banded_matrix_vector(B_y[x], col_vec_half)
        
        # Solve for the full-step column (Psi_t^T[x])
        psi_t[:, x] = solve_banded((1, 1), A_y[x], rhs_y)

    return psi_t

def test_2():
    """
    Tests the 2D grid generation, wave packet initialization, and 2D potentials.
    Displays the X matrix, Y matrix, initial probability density, and the 
    Harmonic potential in a 2x2 subplot figure.
    """
    # 1. Grid Parameters
    dx = 0.5
    x_left, x_right = -10.0, 10.0
    y_down, y_up = -10.0, 10.0

    # Initialize Grid
    x, y, X, Y = init_grid_2d(dx, x_left, x_right, y_down, y_up)

    # 2. Wave Packet Parameters
    x0, y0 = -4.0, -4.0
    sigma = 1.5
    kx, ky = 2.0, 1.0

    # Initialize Particle
    psi = init_moving_gaussian_2d(X, Y, x0, y0, sigma, kx, ky)
    probability_density = np.abs(psi)**2

    # 3. Initialize 2D Potential (Harmonic Oscillator)
    # Using an anisotropic oscillator (different omegas) for visual interest
    V_harmonic = harmonic_potential_2d(X, Y, omega_x=0.5, omega_y=0.8)

    # 4. Plotting
    # Create a 2x2 subplot layout
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Flatten axes array for easy indexing
    ax_flat = axes.flatten()
    
    # Set the extent for imshow
    extent = [x_left, x_right, y_down, y_up]

    # Plot X Mesh Matrix
    im0 = ax_flat[0].imshow(X, extent=extent, origin='lower', cmap='coolwarm')
    ax_flat[0].set_title("X Mesh Matrix")
    ax_flat[0].set_xlabel("x")
    ax_flat[0].set_ylabel("y")
    fig.colorbar(im0, ax=ax_flat[0], fraction=0.046, pad=0.04)

    # Plot Y Mesh Matrix
    im1 = ax_flat[1].imshow(Y, extent=extent, origin='lower', cmap='coolwarm')
    ax_flat[1].set_title("Y Mesh Matrix")
    ax_flat[1].set_xlabel("x")
    ax_flat[1].set_ylabel("y")
    fig.colorbar(im1, ax=ax_flat[1], fraction=0.046, pad=0.04)

    # Plot Probability Density
    im2 = ax_flat[2].imshow(probability_density, extent=extent, origin='lower', cmap='inferno')
    ax_flat[2].set_title(r"Initial Probability Density $|\psi(x,y,0)|^2$")
    ax_flat[2].set_xlabel("x")
    ax_flat[2].set_ylabel("y")
    fig.colorbar(im2, ax=ax_flat[2], fraction=0.046, pad=0.04)

    # Plot 2D Harmonic Potential
    im3 = ax_flat[3].imshow(V_harmonic, extent=extent, origin='lower', cmap='viridis')
    ax_flat[3].set_title("2D Harmonic Potential $V(x,y)$")
    ax_flat[3].set_xlabel("x")
    ax_flat[3].set_ylabel("y")
    # Add contour lines to make the potential shape clearer
    ax_flat[3].contour(X, Y, V_harmonic, colors='white', alpha=0.3, levels=10)
    fig.colorbar(im3, ax=ax_flat[3], fraction=0.046, pad=0.04)

    plt.tight_layout()
    plt.show()

def psi_to_phi_2d(x, y, psi):
    """
    Converts 2D position wave function psi(x,y) to momentum wave function phi(k_x,k_y).
    """
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    Ny, Nx = psi.shape
    
    # Generate 1D frequency arrays and shift them
    kx = np.fft.fftshift(np.fft.fftfreq(Nx, d=dx)) * 2 * np.pi
    ky = np.fft.fftshift(np.fft.fftfreq(Ny, d=dy)) * 2 * np.pi
    
    # Create 2D momentum grids
    KX, KY = np.meshgrid(kx, ky)
    
    # Compute the 2D discrete Fourier transform
    phi_unshifted = np.fft.fft2(psi) * dx * dy / (2 * np.pi)
    
    # Shift array and apply phase correction for the grid's starting offset
    phi = np.fft.fftshift(phi_unshifted) * np.exp(-1j * (KX * x[0] + KY * y[0]))
    
    return kx, ky, KX, KY, phi

def compute_P_theta(psi, X, Y, dx, r_m, delta_r, theta_array, delta_theta):
    """
    Evaluates the angular probability distribution P(theta) within a specified radial band.
    
    Parameters:
    - psi: 2D array of the wave function.
    - X, Y: 2D meshgrid arrays of spatial coordinates.
    - dx: Spatial step size (assumes dx = dy).
    - r_m: Mean radius of the detection band.
    - delta_r: Half-width of the radial band.
    - theta_array: 1D array of angles (in radians) at which to evaluate P(theta).
    - delta_theta: Half-width of the angular sector for integration.
    
    Returns:
    - P_theta: 1D array of probabilities matching the shape of theta_array.
    """
    # 1. Convert Cartesian mesh to Polar mesh
    R = np.sqrt(X**2 + Y**2)
    Theta = np.arctan2(Y, X) # Angles in the range [-pi, pi]
    
    # 2. Probability density multiplied by Cartesian area element (dx * dy)
    # This is equivalent to integrating with r * dr * dtheta in polar coordinates
    prob_density_area = np.abs(psi)**2 * (dx**2)
    
    # 3. Pre-compute the radial mask (this is constant for all angles)
    r_mask = (R >= r_m - delta_r) & (R <= r_m + delta_r)
    
    P_theta = np.zeros_like(theta_array, dtype=float)
    
    # 4. Sweep over target angles and integrate
    for i, target_theta in enumerate(theta_array):
        
        # Calculate shortest angular distance to handle -pi/pi wrap-around
        dTheta = Theta - target_theta
        dTheta = (dTheta + np.pi) % (2 * np.pi) - np.pi
        
        # Create angular mask for the current sector
        theta_mask = np.abs(dTheta) <= delta_theta
        
        # Combine masks to find grid points in the target region
        combined_mask = r_mask & theta_mask
        
        # Sum the probability area elements to approximate the integral
        P_theta[i] = np.sum(prob_density_area[combined_mask])
        
    return P_theta

def simulate_2D(
    V_func, psi_0_func, dx, dt, x_left, x_right, y_down, y_up, 
    frames=100, method="ADI", show_animation=True, output_format='widget',
    animation_frame_scaling=1,
    r_m = 13.0,                    # Mean radius for detection (adjust as needed)
    delta_r = 7.0,                 # Radial band half-width
    delta_theta = np.pi / 36,      # 10 degrees half-width for angular sector
    ):
    """
    Simulates the 2D TDSE and returns an animation/history.
    """
    # 1. Initialize grid and potentials
    x, y, X, Y = init_grid_2d(dx, x_left, x_right, y_down, y_up)
    psi_initial = psi_0_func(X, Y)
    V_2d = V_func(X, Y)
    
    # 2. Precompute matrices for ADI
    if method == "ADI":
        A_x = compute_matrices_A_x(compute_a_x(dx, dt, V_2d))
        B_x = compute_matrices_B_x(compute_b_x(dx, dt, V_2d))
        A_y = compute_matrices_A_y(compute_a_y(dx, dt, V_2d))
        B_y = compute_matrices_B_y(compute_b_y(dx, dt, V_2d))
    else:
        raise ValueError(f"Method '{method}' is not implemented for 2D.")

    psi_curr = psi_initial.copy()

    # 3. Metrics Helper Function (Updated to track most probable values)
    def get_metrics_2d(psi, phi, KX, KY, dkx, dky):
        # Norms and expectations
        norm_x = np.sum(np.abs(psi)**2) * dx * dx 
        norm_k = np.sum(np.abs(phi)**2) * dkx * dky
        exp_x = np.sum(X * np.abs(psi)**2) * dx * dx
        exp_y = np.sum(Y * np.abs(psi)**2) * dx * dx
        exp_px = np.sum(KX * np.abs(phi)**2) * dkx * dky
        exp_py = np.sum(KY * np.abs(phi)**2) * dkx * dky
        exp_H = np.sum(0.5 * (KX**2 + KY**2) * np.abs(phi)**2) * dkx * dky + \
                np.sum(V_2d * np.abs(psi)**2) * dx * dx
        
        # Most probable position coordinates
        idx_y, idx_x = np.unravel_index(np.argmax(np.abs(psi)**2), psi.shape)
        x_p, y_p = X[idx_y, idx_x], Y[idx_y, idx_x]
        
        # Most probable momentum coordinates
        idx_ky, idx_kx = np.unravel_index(np.argmax(np.abs(phi)**2), phi.shape)
        kx_p, ky_p = KX[idx_ky, idx_kx], KY[idx_ky, idx_kx]
        
        return norm_x, norm_k, exp_x, exp_y, exp_px, exp_py, exp_H, x_p, y_p, kx_p, ky_p

    # Calculate initial momentum space and metrics
    kx, ky, KX, KY, phi_init = psi_to_phi_2d(x, y, psi_curr)
    dkx, dky = kx[1] - kx[0], ky[1] - ky[0]
    
    n_x_0, n_k_0, e_x_0, e_y_0, e_px_0, e_py_0, e_H_0, x_p_0, y_p_0, kx_p_0, ky_p_0 = get_metrics_2d(psi_curr, phi_init, KX, KY, dkx, dky)

    # History Dictionary
    steps = frames*animation_frame_scaling
    history = {
        't': np.zeros(steps),
        'norm_x': np.zeros(steps), 'norm_k': np.zeros(steps),
        'exp_x': np.zeros(steps), 'exp_y': np.zeros(steps),
        'exp_px': np.zeros(steps), 'exp_py': np.zeros(steps),
        'exp_H': np.zeros(steps),
        'x_p': np.zeros(steps), 'y_p': np.zeros(steps),
        'kx_p': np.zeros(steps), 'ky_p': np.zeros(steps)
    }

    # # 4. Setup Figure and Axes for Animation
    # fig = plt.figure(figsize=(10, 10))
    # ((ax1, ax2), (ax3, ax4)) = fig.subplots(2, 2)

    # 4. Setup Figure and Axes for Animation
    fig = plt.figure(figsize=(12, 10)) # Slightly wider to accommodate the polar plot
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224, projection='polar') # Set to polar projection

    #-----------------------------
    
    extent_r = [x_left, x_right, y_down, y_up]
    extent_k = [kx[0], kx[-1], ky[0], ky[-1]]

    if show_animation:
        # Position Space
        im_psi = ax1.imshow(np.abs(psi_curr)**2, extent=extent_r, origin='lower', cmap='inferno', vmin=0)
        ax1.set_title(r"Position Space $|\psi(x,y)|^2$")
        ax1.set_xlabel("x")
        ax1.set_ylabel("y")
        fig.colorbar(im_psi, ax=ax1, fraction=0.046, pad=0.04)
        
        # Position Scatter Overlays
        scatter_exp_r = ax1.scatter([e_x_0], [e_y_0], color='cyan', marker='x', s=100, label=r'$\langle \mathbf{r} \rangle$')
        scatter_prob_r = ax1.scatter([x_p_0], [y_p_0], color='purple', marker='+', s=100, label=r'Most Probable $\mathbf{r}$')
        ax1.legend(loc='upper right')

        # Measurement Band Circles 
        circle_inner = plt.Circle((0, 0), r_m - delta_r, color='white', fill=False, linewidth=0.8, linestyle='--')
        circle_outer = plt.Circle((0, 0), r_m + delta_r, color='white', fill=False, linewidth=0.8, linestyle='--')
        ax1.add_patch(circle_inner)
        ax1.add_patch(circle_outer)

        # Potential Energy
        im_V = ax2.imshow(V_2d, extent=extent_r, origin='lower', cmap='viridis')
        ax2.contour(X, Y, V_2d, colors='white', alpha=0.3, levels=10)
        ax2.set_title(r"Potential Energy $V(x,y)$")
        ax2.set_xlabel("x")
        ax2.set_ylabel("y")
        fig.colorbar(im_V, ax=ax2, fraction=0.046, pad=0.04)

        # Momentum Space
        im_phi = ax3.imshow(np.abs(phi_init)**2, extent=extent_k, origin='lower', cmap='magma', vmin=0)
        ax3.set_title(r"Momentum Space $|\phi(k_x, k_y)|^2$")
        ax3.set_xlabel("$k_x$")
        ax3.set_ylabel("$k_y$")
        
        # Momentum Scatter Overlays
        scatter_exp_k = ax3.scatter([e_px_0], [e_py_0], color='cyan', marker='x', s=100, label=r'$\langle \mathbf{p} \rangle$')
        scatter_prob_k = ax3.scatter([kx_p_0], [ky_p_0], color='purple', marker='+', s=100, label=r'Most Probable $\mathbf{p}$')
        ax3.legend(loc='upper right')
        
        k_bound = max(abs(e_px_0), abs(e_py_0)) + 4.0 
        ax3.set_xlim(-k_bound, k_bound)
        ax3.set_ylim(-k_bound, k_bound)
        fig.colorbar(im_phi, ax=ax3, fraction=0.046, pad=0.04)

        theta_array = np.arange(-np.pi, np.pi + 1e-15, delta_theta)
        P_theta_init = compute_P_theta(psi_curr, X, Y, dx, r_m, delta_r, theta_array, delta_theta)
        r_plot_init = np.sqrt(P_theta_init)
        line_P_theta, = ax4.plot(theta_array, r_plot_init, '-', color='green')
        ax4.set_title(fr"$\sqrt{{P(\theta)}}$ at $r={r_m} \pm {delta_r}$", y=-0.15)
        ax4.set_rmin(0)
        ax4.set_rmax(1)
        ax4.set_xticks((theta_array % (2*np.pi)))
        ax4.set_xticklabels((
            (f"${theta/np.pi:.2f}\pi$" if (i%5 == 0) else "") 
            for i,theta in enumerate(theta_array % (2*np.pi))
            ))

        # plt.show()

        # if input() == "": exit()

        # ----------------------------------------

    # 5. Computation and Update Step
    def animate(frame):
        nonlocal psi_curr
        
        print('doing frame',frame+1,'/',frames,"(computation)",end='\r',flush=True)
        
        for step in range(animation_frame_scaling):

            global_step = frame*animation_frame_scaling + step
            # Advance Step
            if global_step > 0:
                if method == "ADI":
                    psi_curr = step_adi_2d(psi_curr, A_x, B_x, A_y, B_y)
            

            # Get momentum space & metrics
            _, _, _, _, phi_curr = psi_to_phi_2d(x, y, psi_curr)
            n_x, n_k, e_x, e_y, e_px, e_py, e_H, x_p, y_p, kx_p, ky_p = get_metrics_2d(psi_curr, phi_curr, KX, KY, dkx, dky)
            
            # Record into history
            history['t'][global_step] = global_step * dt
            history['norm_x'][global_step], history['norm_k'][global_step] = n_x, n_k
            history['exp_x'][global_step], history['exp_y'][global_step] = e_x, e_y
            history['exp_px'][global_step], history['exp_py'][global_step] = e_px, e_py
            history['exp_H'][global_step] = e_H
            history['x_p'][global_step], history['y_p'][global_step] = x_p, y_p
            history['kx_p'][global_step], history['ky_p'][global_step] = kx_p, ky_p

        if show_animation:
            print('doing frame',frame+1,'/',frames,"(plotting)   ",end='\r',flush=True)
            prob_density_x = np.abs(psi_curr)**2
            prob_density_k = np.abs(phi_curr)**2
            im_psi.set_data(prob_density_x)
            im_phi.set_data(prob_density_k)
            
            im_psi.set_clim(vmax=np.max(prob_density_x) or 1.0)
            im_phi.set_clim(vmax=np.max(prob_density_k) or 1.0)
            
            # Update scatter dot positions
            scatter_exp_r.set_offsets([[e_x, e_y]])
            scatter_prob_r.set_offsets([[x_p, y_p]])
            scatter_exp_k.set_offsets([[e_px, e_py]])
            scatter_prob_k.set_offsets([[kx_p, ky_p]])

            P_theta_curr = compute_P_theta(psi_curr, X, Y, dx, r_m, delta_r, theta_array, delta_theta)
            r_plot_curr = np.sqrt(P_theta_curr)
            line_P_theta.set_ydata(r_plot_curr)
                
            fig.suptitle(f"Frame: {frame} | $\\langle x \\rangle$: {e_x:.2f} | $\\langle y \\rangle$: {e_y:.2f}", fontsize=12)
            
            # Return updated artists
            return [im_psi, im_phi, scatter_exp_r, scatter_prob_r, scatter_exp_k, scatter_prob_k, line_P_theta]

            #--------------------------------------------------------

    # 6. Execute Simulation
    anim_widget = None
    ani = None
    if show_animation:
        print(f"\nSimulation started (with animation) for {frames} frames!")
        ani = FuncAnimation(fig, animate, frames=frames, interval=50, blit=True)
        if output_format == 'widget':
            plt.close(fig)
            anim_widget = HTML(ani.to_jshtml())
        elif output_format == 'video':
            plt.close(fig)
            anim_widget = HTML(ani.to_html5_video())
    else:
        print(f"\nRunning fast simulation (no animation) for {frames} frames...")
        for frame in range(frames):
            animate(frame)
        print("\nDone!")

    return {
        'history': history,          # NEW: Added history to output
        'animation': anim_widget,
        'raw_ani': ani if show_animation else None,
        'final_psi': psi_curr,
        'X_grid': X,
        'Y_grid': Y
    }

def plot_static_results_2d(history, cond='', axes=None, fig_static=None, final=True, show_most_prob=True):
    """Generates the static plots from a populated 2D history dictionary."""
    if fig_static is None:
        fig_static, axes = plt.subplots(4, 2, figsize=(8, 12), sharex=True)
    fts = 5

    # 1. Position Trajectories
    axes[0][0].plot(history['t'], history['exp_x'], label=r'$\langle x \rangle$' + f" ({cond})")
    if show_most_prob : axes[0][0].plot(history['t'], history['x_p'], label=r'$x_p$' + f" ({cond})", linestyle='--')
    axes[0][1].plot(history['t'], history['exp_y'], label=r'$\langle y \rangle$' + f" ({cond})")
    if show_most_prob : axes[0][1].plot(history['t'], history['y_p'], label=r'$y_p$' + f" ({cond})", linestyle='--')
    
    if final:
        axes[0][0].set_title("x-coordinate over time")
        axes[0][0].set_ylabel("$x$")
        axes[0][0].legend(loc='upper right', fontsize=fts)
        axes[0][0].grid(True, alpha=0.3)
        axes[0][1].set_title("y-coordinate over time")
        axes[0][1].set_ylabel("$y$")
        axes[0][1].legend(loc='upper right', fontsize=fts)
        axes[0][1].grid(True, alpha=0.3)

        axes[0][1].yaxis.set_label_position('right')
        axes[0][1].yaxis.tick_right()

    # 2. Momentum Trajectories
    axes[1][0].plot(history['t'], history['exp_px'], label=r'$\langle p_x \rangle$' + f" ({cond})")
    if show_most_prob : axes[1][0].plot(history['t'], history['kx_p'], label=r'$p_{xp}$'+ f" ({cond})", linestyle='--')
    axes[1][1].plot(history['t'], history['exp_py'], label=r'$\langle p_y \rangle$' + f" ({cond})")
    if show_most_prob : axes[1][1].plot(history['t'], history['ky_p'], label=r'$p_{yp}$'+ f" ({cond})", linestyle='--')
    if final:
        axes[1][0].set_title("$p_x$ over time")
        axes[1][0].set_ylabel("$p_x$")
        axes[1][0].legend(loc='upper right', fontsize=fts)
        axes[1][0].grid(True, alpha=0.3)
        axes[1][1].set_title("$p_y$ over time")
        axes[1][1].set_ylabel("$p_y$")
        axes[1][1].legend(loc='upper right', fontsize=fts)
        axes[1][1].grid(True, alpha=0.3)

        axes[1][1].yaxis.set_label_position('right')
        axes[1][1].yaxis.tick_right()

    # 3. Normalization and Energy
    axes[2][0].plot(history['t'], history['norm_x'], label=r'$\int |\psi|^2 dx dy$'+ f" ({cond})")
    axes[2][1].plot(history['t'], history['norm_k'], label=r'$\int |\phi|^2 dk_x dk_y$'+ f" ({cond})")
    if final:
        axes[2][0].set_title(r"$\langle\psi|\psi\rangle$ over time")
        axes[2][0].set_ylabel(r"$\langle\psi|\psi\rangle$")
        axes[2][0].legend(loc='upper right', fontsize=fts)
        axes[2][0].grid(True, alpha=0.3)
        axes[2][1].set_title(r"$\langle\phi|\phi\rangle$ over time")
        axes[2][1].set_ylabel(r"$\langle\phi|\phi\rangle$")
        axes[2][1].legend(loc='upper right', fontsize=fts)
        axes[2][1].grid(True, alpha=0.3)

        axes[2][1].yaxis.set_label_position('right')
        axes[2][1].yaxis.tick_right()
    
    axes[3][0].plot(history['t'], history['exp_H'], label=r'$\langle H \rangle$'+ f" ({cond})")
    if final:
        axes[3][0].set_ylabel("Energy")
        axes[3][0].legend(loc='upper right', fontsize=fts)
        axes[3][0].set_xlabel("Time (t)")
        # plt.tight_layout()
    return fig_static

def MiniProject3():
    """
    Generates an animation and saves comprehensive static line plots of history.
    """

    param_list = [
        {"R":3, "k_x":2,"sigma":1, "V_0":1}, # primary setup
        {"R":1, "k_x":2,"sigma":1, "V_0":1},  # same size for disk and packet
        {"R":0.3, "k_x":2,"sigma":1, "V_0":1},  # larger packet than disk
        {"R":3, "k_x":2,"sigma":1, "V_0":-1},  # attractive potential
        {"R":3, "k_x":2,"sigma":1, "V_0":0},  # free particle
        {"R":3, "k_x":5,"sigma":1, "V_0":1},  # faster packet
    ]

    for params in param_list[::-1]:

        print("parameters :",params)
        param_str = "__".join([f"{key}_{val}" for key,val in params.items()])

        def initial_packet(X, Y):
            sigma = params['sigma']
            x0 = -8.0
            y0 = 0.0
            kx = params['k_x']
            ky = 0.0
            return init_moving_gaussian_2d(X, Y, x0, y0, sigma, kx, ky)

        def disk(X,Y):
            return disk_potential_2d(X,Y, params['R'], v_in=params['V_0'])

        print("Initializing 2D Simulation...")
        res = simulate_2D(
            V_func=disk,
            psi_0_func=initial_packet,
            dx=0.1,
            dt=0.003,
            x_left=-20.0,
            x_right=20.0,
            y_down=-20.0,
            y_up=20.0,
            frames=50,  
            animation_frame_scaling=100,             
            method="ADI",
            show_animation=True,
            output_format='save'      
        )

        if res['raw_ani'] is not None:
            anim_file = "images/scatterer_2d__" + param_str + ".gif"
            writer = PillowWriter(fps=10)
            res['raw_ani'].save(anim_file, writer=writer)
            print(f"\nGIF saved successfully as '{anim_file}'!")

        # Generate and Save the Static History Figure
        stat_file = 'images/scatterer_2d_summary__' + param_str + '.png'
        print("Generating static history figures...")
        fig_static = plot_static_results_2d(res['history'])
        fig_static.savefig(stat_file, dpi=300, bbox_inches='tight')
        print(f"Saved static plot as '{stat_file}'")


    plt.close('all')
    fig_static, axes = plt.subplots(4, 2, figsize=(8, 12), sharex=True)
    plt.delaxes(axes[3][1])

    dx_dt_list = [
        # dx  dt
        (0.5, 0.1),
        (0.2, 0.1),
        (0.1, 0.1),
        (0.1, 0.02),
        (0.1, 0.003),
    ]
    for i,dx_dt in enumerate(dx_dt_list[::-1]):
            dx,dt = dx_dt
            print('\ndt=',dt,'and dx=',dx)
            res = simulate_2D(
                V_func=disk_potential_2d,
                psi_0_func=init_moving_gaussian_2d,
                dx=dx,dt=dt,
                x_left=-20.0,
                x_right=20.0,
                y_down=-20.0,
                y_up=20.0,
                frames=int(round(15/dt)),  
                method="ADI",
                show_animation=False,
                output_format='save'      
            )
            history = res['history']
            plot_static_results_2d(
                history,
                f"$\Delta t={dt}, \, \Delta x ={dx}$",
                axes, fig_static, 
                (len(dx_dt_list)==i+1),
                False
                )
    stat_file = "images/scatterer_2d_convergence.png"
    fig_static.savefig(stat_file, dpi=300, bbox_inches='tight')
    print(f"Saved static plot as '{stat_file}'")
    plt.close()


if __name__ == "__main__": 
    print("-"*15 + "   starting simulations   " + "-"*15)
    # MiniProject3()
    print("-"*15 + " all simulations finished " + "-"*15)