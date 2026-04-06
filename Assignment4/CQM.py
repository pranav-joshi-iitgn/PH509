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

def update_psi(psi, A_banded, B_banded, method="split-operator"):
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
    else:
        raise ValueError(f"Method '{method}' is not implemented.")

def simulate_1D(
    V_func, psi_0_func, dx, dt, x_left, x_right, frames=100, methods=["split-operator"],
    F_func=None, x0_classical=None, p0_classical=None, show_animation=True, output_format='widget'
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

    # --- NEW: History and Classical Init ---
    history = {
        't': np.zeros(frames), 'norm_x': np.zeros(frames), 'norm_k': np.zeros(frames),
        'exp_x': np.zeros(frames), 'exp_p': np.zeros(frames), 'exp_H': np.zeros(frames),
        'x_p': np.zeros(frames), 'k_p': np.zeros(frames), 'x_c': np.zeros(frames)
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
        vline_xp = ax1.axvline(x_p_0, color='red', linestyle=':', label=r'$x_p$ (Most Probable)') # <-- ADDED
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
        vline_kp = ax3.axvline(k_p_0, color='red', linestyle=':', label=r'$k_p$ (Most Probable)') # <-- ADDED
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
                matrices[method]['psi'] = update_psi(matrices[method]['psi'], matrices[method]['A'], matrices[method]['B'], method=method)
        
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
        psi_curr, phi_curr, n_x, n_k, e_x, e_p, e_H, x_p, k_p = compute_step(frame)
        
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
            compute_step(frame)
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
    ax_s1.plot(history['t'], history['x_p'], label=r'$x_p$ (Most Probable)', color='green', linestyle='--')
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

if __name__ == "__main__":
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
        dx=0.05, 
        dt=0.01, 
        x_left=-10.0, 
        x_right=10.0, 
        frames=300, 
        methods=["split-operator"],
        show_animation=True,
        output_format='save'       # Avoids compiling the heavy JS widget in terminal
    )

    if res['raw_ani'] is not None:
        writer = PillowWriter(fps=20)
        res['raw_ani'].save("harmonic_simulation.gif", writer=writer)
        print("\nGIF saved successfully!")

    # 3. Generate and Save the Static Figure NOW
    fig_static = plot_static_results(res['history'], F_func=harmonic_force)
    fig_static.savefig("static_results.png", dpi=300, bbox_inches='tight')
    print("Saved static plot as 'static_results.png'")

    plt.close('all')