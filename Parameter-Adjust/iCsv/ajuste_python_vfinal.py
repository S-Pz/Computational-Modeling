import numpy as np
import pandas as pd
from scipy.integrate import  solve_ivp
from scipy.optimize import differential_evolution
import math
import matplotlib.pyplot as plt 

data_path = './I.csv'

dt = 0.01
tfinal = 50
times = np.arange(0,tfinal+dt,dt)

S0 = 90082.0
I0 = 203.0 #203
R0 = 0.0

def odeSystem(t, u, beta, alpha):

    S, I, R = u
    dS_dt = - beta*S*I 
    dI_dt = beta*S*I - alpha*I
    dR_dt = alpha*I

    return [dS_dt, dI_dt, dR_dt]  

def isReferenceTime(times, ct):
    for t in times: 
        if (abs(ct - t) <= 10**(-5)):
            return True 
    return False

def solve(x):
    global data, reference_times

    u = [S0, I0, R0]

    beta = x[0]
    alpha = x[1]
    params = (beta, alpha)
    
    def solveOde(t, y):
        return odeSystem(t, y, *params)

    results = solve_ivp(solveOde,(0, tfinal), u, t_eval=times, method='RK45')    

    u = results.y[1,:]
    error = 0
    sumobs = 0
    i = 0
    j = 0    
    for t in times:
        if isReferenceTime(reference_times,t):
            p_data = data["I"][i]
            error += (u[j] - p_data)*(u[j] - p_data) 
            sumobs += p_data*p_data
            i = i + 1
        j = j + 1

    error = math.sqrt(error/sumobs) #Erro norma 2                  
    return error 

if __name__ == "__main__":
    
    global data, reference_times 
    data = pd.read_csv(data_path, delimiter=',') if data_path.endswith('.csv') else pd.read_excel(data_path)
    
    reference_times = data["Semana"]
    dados_I = data["I"]
    #plota os dados experimentais 
    fig = plt.figure()
    fig.set_size_inches(8, 6)
    plt.scatter(reference_times, dados_I, marker='o', color='black', label='dados')
    
    bounds = [
        (0.00001, 0.01), (0.01, 0.9)
    ]

    #chama evolução diferencial, result contém o melhor individuo
    solucao = differential_evolution(solve, bounds, strategy='best2bin', maxiter=50, popsize=40,atol=10**(-3), tol=10**(-3), mutation=0.8, recombination=0.5, disp=True, workers=4)
    
    print(solucao.x)
    #saving the best offspring...
    np.savetxt('solucao_ajuste.txt',solucao.x, fmt='%.6f')        
    best = solucao.x
    error = solve(best)
    #print("ERROR ", error)
    print(solucao.population)
    print(solucao.population_energies)
    
    u = [S0, I0, R0]
    result_best = solve_ivp(odeSystem,(0, tfinal+dt), u, t_eval=times, args=best, method='RK45')
    plt.plot(result_best.t, result_best.y[1,:], color='red', label='I')
    plt.legend(loc='best')    
    fig.savefig('I.png', format='png')
    plt.show()    

    fig = plt.figure()
    fig.set_size_inches(8, 6)
    plt.plot(result_best.t, result_best.y[0,:], color='blue', label='S')
    plt.legend(loc='best')
    fig.savefig('S.png', format='png')
    plt.show()
    
    fig = plt.figure()
    fig.set_size_inches(8, 6)
    plt.plot(result_best.t, result_best.y[2,:], color='green', label='R')
    plt.legend(loc='best')
    fig.savefig('R.png', format='png')
    plt.show()
