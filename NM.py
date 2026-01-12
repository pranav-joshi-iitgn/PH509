import numpy as np
import matplotlib.pyplot as plt

def step(i,t,x,dydt,method,dt):
    if method=="Euler":
        k1 = dydt(t[i],x[i])*dt
        return x[i] + k1
    elif method=="RK2":
        k1 = dydt(t[i],x[i])*dt
        k2 = dydt(t[i]+dt,x[i]+k1)*dt
        return x[i] + (k1+k2)/2

def Numerical_Solver(
    x_0, # initial value. Can be a single value or numpy array
    dt,  # Interval \Delat t
    t_ini,t_max, # real values  
    dydt, # a function mapping x,t to the derivative
    methods=["Euler"],
    x_analytical=None, # analytical solution (function). If avaliable, compute errors too
    visualise=False, # if True, does some plotting
    variable_name=None,
    compute_error=True,
    compute_one_step_error=True,
):
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
                # x_one_step = x_analytical.copy() 
                # for i in range(0,N):
                #     x_one_step[i+1] = step(i,t,x_analytical,dydt,method)
                # err_one_step = x_analytical - x_one_step
                err_one_step = - np.array([0] + [step(i,t,x_analytical,dydt,method,dt) for i in range(0,N)])
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
    return t,X

def RadioactiveDecay(
    N_0=1, tau=1, t_max=5, dt = None,
    precision_level = 100,
    visualise = False,
    methods = ["Euler"],
    compute_error = False,
    compute_one_step_error = False
    ):
    def dydt(t,N):return - N / tau
    def N(t):return N_0 * np.exp(-t/tau)
    if dt is None:dt = tau/precision_level
    t,results = Numerical_Solver(
        N_0,dt,0,t_max,dydt,
        methods=methods,
        x_analytical=N, # analytical solution (function). If avaliable, compute errors too
        visualise=visualise, # if True, does some plotting
        variable_name="$N(t)$",
        compute_error=compute_error,
        compute_one_step_error=compute_one_step_error
    )
    return t,results
