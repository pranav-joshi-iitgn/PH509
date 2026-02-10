import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from IPython.display import HTML


last_method = None
method_params = None
def step(i,t,x,dydt,method,dt):
    global last_method, method_params
    if method=="Euler":
        k1 = dydt(t[i],x[i])*dt
        last_method = method
        return x[i] + k1
    elif method=="RK2":
        k1 = dydt(t[i],x[i])*dt
        k2 = dydt(t[i]+dt,x[i]+k1)*dt
        last_method = method
        return x[i] + (k1+k2)/2
    elif method.startswith("RK2("):
        if method == last_method:
            a1,a2,p,q = method_params
        else:
            params = method[4:-1].strip().split(",")
            ARGS = []
            KWARGS = {"a1":None,"a2":None,"p":None,"q":None}
            for x in params:
                if "=" in x :
                    P,V = x.split("=")
                    KWARGS[P.strip().lower()] = float(V)
                else:
                    ARGS.append(float(x))
            ARGS = ARGS[::-1]
            for x in KWARGS:
                if KWARGS[x] is None:
                    KWARGS[x] = ARGS.pop()
            if KWARGS["a2"] is None: a2 = 1 - a1
            else: a2 = KWARGS["a2"]
            if KWARGS["a1"] is None: a1 = 1 - a2
            else: a1 = KWARGS["a1"]
            if KWARGS["p"] is None: p = 1/(2*a2)
            else: p = KWARGS["p"]
            if KWARGS["q"] is None: q = 1/(2*a2)
            else: q = KWARGS["q"]
            last_method = method
            method_params = (a1,a2,p,q)

        k1 = dydt(t[i],x[i])*dt
        k2 = dydt(t[i]+(p*dt), x[i]+(k1*q*dt))*dt
        return x[i] + (a1*k1 + a2*k2)

    elif method == "EC":
        k1_theta, k1_omega = dydt(t[i], x[i])*dt
        omega_new = x[i][1] + k1_omega
        theta_new = x[i][0] + dt*omega_new
        return (theta_new,omega_new)

    elif method == "LF": # Leap Frog
        theta, omega = x[i]
        theta_intermediate = theta + dt*omega/2
        k1_theta,k1_omega = dydt(t[i],np.array(theta_intermediate,omega))*dt
        omega_new = omega + k1_omega
        theta_new = theta_intermediate + omega_new*dt/2
        return (theta_new,omega_new)


def Numerical_Solver(
    x_0, # initial value. Can be a single value or numpy array
    dydt, # a function mapping x,t to the derivative
    dt=None,  # Interval \Delat t
    t_ini=0,t_max=None, # real values  
    tau = None,
    precision_level = 100,
    methods=["Euler"],
    x_analytical=None, # analytical solution (function). If avaliable, compute errors too
    visualise=False, # if True, does some plotting
    variable_name=None,
    compute_error=True,
    compute_one_step_error=True,
) -> (tuple[np.ndarray,dict]|tuple[np.ndarray,dict,HTML]):
    if tau is not None:
        if dt is None: dt = tau/precision_level
        if t_max is None: t_max = 10*tau
    t = np.arange(t_ini, t_max, dt)
    N = t.shape[0] -1 # excluding the start value
    if isinstance(x_0,np.ndarray):
        d = x_0.shape[0] # dimension
    elif isinstance(x_0,list) or isinstance(x_0,tuple):
        d = len(x_0)
        x_0 = np.array(x_0)
    else:d =1
    X = {}
    if x_analytical is not None:
        try: x_analytical = x_analytical(t)
        except:x_analytical = np.array([x_analytical(t_i) for t_i in t])
        X["analytical"] = x_analytical
    for method in methods:
        if d > 1:x = np.zeros((N+1,d))
        else:x = np.zeros(N+1)
        x[0] = x_0
        for i in range(0,N):x[i+1] = step(i,t,x,dydt,method,dt)
        X[method] = x
        if x_analytical is not None:
            if compute_error:
                err = x_analytical -x
                X[method + " error"] = err
            if compute_one_step_error:
                x_one_step = x_analytical.copy() 
                for i in range(0,N):
                    x_one_step[i+1] = step(i,t,x_analytical,dydt,method,dt)
                err_one_step = x_analytical - x_one_step
                X[method + " one step error"] = err_one_step
    if visualise and d == 1:
        plt.figure()
        for method in methods:
            plt.plot(t, X[method], label=method)
        if x_analytical is not None:
            plt.plot(t, x_analytical, label="analytical")
        plt.xlabel("$t$")
        if variable_name is not None:plt.ylabel(variable_name)
        plt.title("Numerical Solution (" + ",".join(methods) + ")")
        plt.legend()
        if compute_error and (x_analytical is not None):
            plt.figure()
            for method in methods:
                lbl = method + " error"
                plt.plot(t, X[lbl],label=lbl)
            plt.xlabel("$t$")
            plt.ylabel("error")
            plt.legend()
        if compute_one_step_error and (x_analytical is not None):
            plt.figure()
            for method in methods:
                lbl = method + " one step error"
                plt.plot(t, X[lbl],label=lbl)
            plt.xlabel("$t$")
            plt.ylabel("error")
            plt.legend()
        plt.show()
    elif visualise=="PS" and d == 2:
        fig = plt.figure()
        ax = plt.subplot()
        all_x_min = min([X[m][:,0].min() for m in methods])
        all_y_min = min([X[m][:,1].min() for m in methods])
        all_x_max = max([X[m][:,0].max() for m in methods])
        all_y_max = max([X[m][:,1].max() for m in methods])
        Dx = all_x_max - all_x_min
        Dy = all_y_max - all_y_min
        ax.set_xlim(all_x_min - Dx/10, all_x_max + Dx/10)
        ax.set_ylim(all_y_min - Dy/10, all_y_max + Dy/10)

        LINES = {method : ax.plot([],[],'-',label=method)[0] for method in methods}
        POINTS = {method : ax.plot([],[],'bo')[0] for method in methods}
        ax.set_xlabel(variable_name[0])
        ax.set_ylabel(variable_name[1])
        ax.legend()
        def update(i):
            for method in methods:
                x,y = X[method][:i+1,0],X[method][:i+1,1]
                LINES[method].set_data(x,y)
                POINTS[method].set_data([x[i]],[y[i]])
            return [x for x in LINES.values()] + [x for x in POINTS.values()]
        ani = FuncAnimation(fig, update, frames=len(t), interval=50, blit=False)
        plt.close(fig)
        ani = HTML(ani.to_html5_video())
        return t,X,ani
    elif visualise and d >= 2:
        fig = plt.figure()
        AX:list[plt.Axes] = fig.subplots(d,1,sharex=True)
        for i in range(d) : 
            for method in methods:
                AX[i].plot(t, X[method][:,i], label=method)
            if x_analytical is not None:
                AX[i].plot(t, x_analytical[:,i], label="analytical")
            if variable_name is not None:AX[i].set_ylabel(variable_name[i])
            AX[i].legend()
        AX[-1].set_xlabel("$t$")
        fig.suptitle("Numerical Solution (" + ",".join(methods) + ")")
        if compute_error and (x_analytical is not None):
            fig = plt.figure()
            AX:list[plt.Axes] = fig.subplots(d,1,sharex=True)
            for i in range(d):
                for method in methods:
                    lbl = method + " error"
                    AX[i].plot(t, X[lbl][:,i],label=lbl)
                if variable_name is not None: vn = variable_name[i]
                else: vn = f"var{i}"
                AX[i].set_ylabel("error in " + vn)
                AX[i].legend()
            AX[-1].set_xlabel("$t$")
        if compute_one_step_error and (x_analytical is not None):
            fig = plt.figure()
            AX:list[plt.Axes] = fig.subplots(d,1,sharex=True)
            for i in range(d):
                for method in methods:
                    lbl = method + " one step error"
                    AX[i].plot(t, X[lbl][:,i],label=lbl)
                if variable_name is not None: vn = variable_name[i]
                else: vn = f"var{i}"
                AX[i].set_ylabel("error in " + vn)
                AX[i].legend()
            AX[-1].set_xlabel("$t$")
        plt.show()
    return t,X


def RadioactiveDecay(N_0=1, tau=1, t_max=5, **kwargs):
    def dydt(t,N):return - N / tau
    def N(t):return N_0 * np.exp(-t/tau)
    return Numerical_Solver(
        N_0,dydt,tau=tau,
        x_analytical=N,
        variable_name="$N(t)$",
        **kwargs)

def Aharmonic(y_0=[np.pi/4,0], g=10, l=10, **kwargs):
    def dydt(t,y):
        theta,omega = y[0],y[1]
        return np.array([omega, -(g/l)* np.sin(theta)])
    tau = 2 * np.pi * (l/g)**0.5
    return Numerical_Solver(
        y_0,dydt,tau=tau,
        variable_name=[r"$\theta(t)$",r"$\omega(t)$"],
        **kwargs)

def CoupledRadioactiveDecay(tau_A=1, tau_B=2, N_A_0=1, N_B_0=0,**kwargs):
    def dydt(t,N):
        N_A,N_B = N[0],N[1]
        return np.array([
            -N_A/tau_A,
            N_A/tau_A - N_B/tau_B
            ])
    def N(t):
        assert (not isinstance(t,np.ndarray)) or (t.shape == (1,))
        return np.array([
            np.exp(-t/tau_A),
            (np.exp(-t/tau_A) -np.exp(-t/tau_B))/((tau_A/tau_B)-1)
        ])
    return Numerical_Solver(
        np.array([N_A_0,N_B_0]),dydt,tau=min(tau_A,tau_B),
        x_analytical=N, # analytical solution 
        variable_name=["$N_A(t)$","$N_B(t)$"],
        **kwargs)

def LF_step(q,p,a,dt,i,**kwargs):
    q_mid = q[i] + p[i]*dt/2
    p[i+1] = p[i] + a(q_mid,**kwargs)*dt
    q[i+1] = q_mid + p[i+1]*dt/2

def LFTT_step(W,t,q,p,a,dt,i,**kwargs):
    q_mid = q[i] + p[i]*(dt/2)/W[i]
    t_mid = t[i] + (dt/2)/W[i]
    q_mid_mag = np.linalg.norm(q_mid)
    p[i+1] = p[i] + a(q_mid,**kwargs)*dt*q_mid_mag
    W[i+1] = W[i] - np.dot(q_mid, ((p[i+1] + p[i])/2))*dt*(q_mid_mag**-2)
    q[i+1] = q_mid + p[i+1]*(dt/2)/W[i+1]
    t[i+1] = t_mid + (dt/2)/W[i+1]


GM_s = 4*(np.pi**2) # the GM values for sun in AU^3 year^-2
plt.rcParams['ytick.labelsize'] = 8


def gravity(r,GM=GM_s,delta=0,lamb=0):
    if lamb == 0: return - GM*r/(np.linalg.norm(r)**(3+delta))
    else: 
        r_mag = np.linalg.norm(r)
        return -GM*r*(r_mag**-(3+delta) + lamb*r_mag**-(5+delta))

def Orbit(
    GM=GM_s,
    r_init=np.array([1.017, 0.0]),
    v_init=np.array([0.0,6.179]),
    methods = {"LF":0.001},
    deltas =[0,0.01],
    lambs = [0],
    T_max = 50, # years
    a = gravity,
    vis=False
    ):
    if vis:
        fig = plt.figure(figsize=(10,20))
        AX = fig.subplots(13,1,sharex=True)
        fig2 = plt.figure()
        AX2 = fig2.subplots(1,1)
        colors = ["blue","red","green","black","purple","yellow","brown"][::-1]

    for delta in deltas:
        for method, dt in methods.items():
            for lamb in lambs:
                c = colors.pop()
                t = np.arange(0,T_max+dt/2,dt)
                r = np.zeros((t.shape[0],2))
                v = np.zeros((t.shape[0],2))
                r[0] = r_init
                v[0] = v_init
                if method=="LF":
                    for i in range(t.shape[0]-1):
                        LF_step(r,v,gravity,dt,i,lamb=lamb,delta=delta)
                elif method=="LFTT":
                    W = np.zeros_like(t)
                    W[0] = 1/np.linalg.norm(r_init)
                    for i in range(t.shape[0]-1):
                        LFTT_step(W,t,r,v,gravity,dt,i,lamb=lamb,delta=delta)
                else: raise ValueError(f"unknown method: {method}")
                x = r[:,0]
                y = r[:,1]
                v_x = v[:,0]
                v_y = v[:,1]
                theta = np.arctan2(y,x)
                theta = theta + 2*np.pi*(theta <  0)
                r = (x**2 + y**2)**0.5
                v = (v_x**2 + v_y**2)**0.5
                E = 0.5*v**2 - GM/(r*(1+delta))
                L = x*v_y - y*v_x
                A_x,A_y = [L*v_y -GM*x/r,-L*v_x -GM*y/r]
                A = (A_x**2 + A_y**2)**0.5
                phi = np.arctan2(A_y,A_x)
                phi = 2*np.pi*(phi < 0) + phi
                v_r = (x*v_x + y*v_y)/r
                v_theta = (-y*v_x + x*v_y)/r
                omega = v_theta/r
                a = ((2/r) - ((v**2)/GM))**-1
                T = 2*np.pi*((a**3)/GM)**0.5 
                r_prev = np.concatenate((r[:1],r[:1],r[:-2]))
                r_next = np.concatenate((r[2:],r[-1:],r[-1:]))
                i_r = np.flatnonzero((r < r_prev) & (r <= r_next))
                theta_prev = np.concatenate((theta[:1],theta[:-1]))
                theta_next = np.concatenate((theta[1:],theta[-1:]))
                i_theta = np.flatnonzero((theta < theta_prev) & (theta <= theta_next))
                n_r = np.zeros_like(t,dtype=int)
                for i in i_r :
                    n_r = n_r + (np.arange(t.shape[0]) >= i)
                n_theta = np.zeros_like(t,dtype=int)
                for i in i_theta :
                    n_theta = n_theta + (np.arange(t.shape[0]) >= i)
                if not vis:continue
                cond = (
                    f"({method}, " +
                    r"$\delta={" + str(delta)+ 
                    r"}$, $\Delta t=" + str(dt) + 
                    r"$, $\lambda=" + str(lamb)+"$)")
                AX[0].plot(t,v_x,label=cond,color=c) 
                AX[1].plot(t,v_y,label=cond,color=c)
                AX[2].plot(t,theta,label=cond,color=c)
                AX[2].plot(t[i_theta],theta[i_theta],'o',label="minima " + cond,color=c)
                AX[3].plot(t,r,label=cond,color=c)
                AX[3].plot(t[i_r],r[i_r],'o',label="minima " + cond,color=c)
                AX[4].plot(t,v,label=cond,color=c)
                AX[5].plot(t,v_r,label=cond,color=c)
                AX[6].plot(t,omega,label=cond,color=c)
                AX[7].plot(t,n_r,label=cond,color=c)
                AX[8].plot(t,n_theta,label=cond,color=c)
                AX[9].plot(t,E,label=cond,color=c)
                AX[10].plot(t,L,label=cond,color=c)
                AX[11].plot(t,A,label=cond,color=c)
                AX[12].plot(t,phi,label=cond,color=c)

                AX2.plot(x,y,label=cond,color=c)
                AX2.plot(x[0],y[0],'o',label="initial position" + cond,color=c)

    if vis:
        AX[2].set_ylim(0,2*np.pi)
        AX[2].set_yticks([0,2*np.pi],[0,r"$2\pi$"])
        AX[2].legend()
        # AX[3].legend()
        AX[7].legend()
        AX[8].legend()
        AX[12].set_ylim(0,2*np.pi)
        AX[12].set_yticks([0,2*np.pi],[0,r"$2\pi$"])
        ylabels = ["v_x","v_y",r"\theta","r","v","v_r",r"\omega",r"n_r",r"n_\theta",r"E/m",r"L/m",r"A/m^2",r"\phi"]
        for ax,lbl in zip(AX,ylabels): 
            ax.set_ylabel("$" + lbl + "$")
            # ax.legend()

        AX[12].set_xlabel("$t$ (years)")
        AX[12].legend()

        AX2.set_xlabel("$x$ (AU)")
        AX2.set_ylabel("$y$ (AU)")
        AX2.plot([0],[0],'o',label="sun",color="yellow",markersize=10)
        AX2.axhline(0)
        AX2.axvline(0)
        AX2.set_aspect('equal')
        AX2.set_xlim(-0.1-np.linalg.norm(r_init),0.1+2*np.linalg.norm(r_init))
        AX2.set_ylim(-0.1-np.linalg.norm(r_init),0.1+2*np.linalg.norm(r_init))
        AX2.legend()

Orbit(vis=True,T_max=10,deltas=[0],methods={"LF":0.01,"LFTT":0.01})
plt.show()