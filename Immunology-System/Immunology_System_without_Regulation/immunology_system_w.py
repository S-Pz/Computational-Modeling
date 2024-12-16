import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def ode_system(t, u, constants):

    """
    system of first order differential equations
    t: discrete time step value
    y: state vector (vetor das variaveis/populacoes do modelo)
    """
    alpha = constants[0]
    betha = constants[1]
    betha_2 =constants[2]
    gama = constants[3]
    m = constants[4]
    
    N = u[0]
    TD = u[1]
    CH = u[2]

    dNdt = (-alpha*N)+(betha*CH)+(betha_2*TD)
    dTDdt = alpha*N
    dCHdt = gama*N - m *CH

    return np.array([dNdt,dTDdt,dCHdt])

def rk4(f, tk, _uk, _dt=0.01, **kwargs):
    """
    single-step fourth-order numerical integration (RK4) method
    f: system of first order ODEs
    tk: current time step
    _uk: current state vector [y1, y2, y3, ...]
    _dt: discrete time step size
    **kwargs: additional parameters for ODE system
    returns: y evaluated at time k+1
    """
    # evaluate derivative at several stages within time interval
    k1 = f(tk, _uk, **kwargs)
    k2 = f(tk + _dt / 2, _uk + (k1 * (_dt / 2)), **kwargs)
    k3 = f(tk + _dt / 2, _uk + (k2 * (_dt / 2)), **kwargs)
    k4 = f(tk + _dt, _uk + (k3 * _dt), **kwargs)

    # return an average of the derivative over tk, tk + dt
    return _uk + (_dt / 6) * (k1 + (2 * k2) + (2 * k3) + k4)

def euler(f, tk, _uk, _dt=0.001, **kwargs):
    """
    single-step explicit method
    func: system of first order ODEs
    tk: current time step
    _uk: current state vector [y1, y2, y3, ...]
    _dt: discrete time step size
    **kwargs: additional parameters for ODE system
    returns: y evaluated at time k+1
    """
    return _uk + _dt*f(tk, _uk, **kwargs)

def solveSystem(time, dt, y0, method):
    t = 0
    yk = y0
    state = []
    
    #Adicionar a ordem dos parâmetros
    parameters = [0.1,0.1,0.1,0.1,0.1]

    if method == "euler":
        for t in time:
            state.append(yk)
            yk = euler(ode_system, t, yk, dt, constants=parameters)

    elif method == "rk4":
        for t in time:
            state.append(yk)
            yk = rk4(ode_system, t, yk, dt, constants=parameters)

    return np.array(state)

def save(time, state, names):
    df = pd.DataFrame(state, columns = names)
    df.insert(0, 'time', time)
    df.to_csv('results.csv', float_format='%.5f', sep=',')

def plot(time, state, names):
    fig, ax = plt.subplots()
    fig.set_size_inches(8, 6)
    for i in range(len(names)):
        ax.plot(time, state[:, i], label=names[i], linewidth='2')

    #time2 = np.arange(0, tfinal + 0.001, 0.001)
    #sol = np.exp(r*time2)
    #ax.plot(time2, sol, label='solucao analitica', linewidth='2')
    ax.set(xlabel='dias', ylabel='populacao')
    plt.legend(loc='best')
    fig.savefig('exponencial.png', format='png')
    plt.show()

if __name__ == "__main__":
    names = ['N', 'TD', 'CH']
    dt = 0.01
    tfinal = 50
    time = np.arange(0, tfinal + dt, dt)
    initial_condition = np.array([0,10,0])
    result = solveSystem(time, dt, initial_condition, "rk4")
    save(time, result, names)
    plot(time, result, names)