import math
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Saída
pdf = PdfPages("resposta_inflamatoria_estocastica.pdf")

# Tempo final 
tf = 500

# Parâmetros 
kpg = 0.3
sm = 0.005
kpn = 1.8
micro_n = 0.05
sc = 0.0125
micro_c = 0.1
knp = 0.1
pinf = 20
micro_m = 0.002
snr = 0.08
kdn = 0.35
kcn = 0.04
cinf = 0.28
knd = 0.02
kpm = 0.6
kmp = 0.01
micro_nr = 0.12
micro_d = 0.02
kcnd = 48
xdn = 0.06
knn = 0.01
muj = 0.12
mud = 0.02
mu_c = 0.1
sigma = 0.02

def f(x, CA=0.125, cinf=0.28):
    return x / (1 + (CA/cinf)**2)

def fs(x, xdn=0.06):
    return x**6 / (xdn**6 + x**6)

def R(N_star, P, D, knn=0.01, kdn=0.35, knp=0.1):
    return knn * N_star + knp * P + kdn * D


if __name__ == "__main__":

    num_simulations = 10

    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    fig4, ax4 = plt.subplots()

    ax1.set_title('População de Patógenos')
    ax2.set_title('População de Fagócitos Ativados')
    ax3.set_title('Dano Tecidual')
    ax4.set_title('Mediadores Anti-inflamatórios')

    for k in range(num_simulations):
        
        random.seed(k)
        
        times = []
        pathogen_result = []
        phagocyte_result = []
        damage_result = []
        mediator_result = []
        
        # Condições iniciais
        P = 1 
        D = 0
        N_star = 0.125
        CA = 0
        t = 0
        
        while t < tf:
            if P <= 0 and D <= 0:
                break
            
            # Ruído
            tau = random.uniform(0.1, 0.5)
            dW = np.random.normal(0, sigma)
            
            # EDOs
            dP = (kpg * P * (1 - P / pinf) - 
                    (kpm * sm * P) / (micro_m + kmp * P) - 
                    kpn * f(N_star) * P + dW)
            
            dN_star = ((snr * R(N_star, P, D)) / (muj + R(N_star, P, D)) - 
                        muj * N_star + dW)
            
            dD = kdn * fs(f(N_star)) - mud * D + dW
            
            dCA = (sc + kcn * f(N_star + kcnd * D) / 
                    (1 + f(N_star + kcnd * D)) - mu_c * CA + dW)
            
            # Atualização das variáveis
            P = max(0, P + dP)
            N_star = max(0, N_star + dN_star)
            D = max(0, D + dD)
            CA = max(0, CA + dCA)
            
            
            # Atualização do tempo
            t += tau
            
            times.append(t)
            pathogen_result.append(P)
            phagocyte_result.append(N_star)
            damage_result.append(D)
            mediator_result.append(CA)
        
        ax1.plot(times, pathogen_result, label=f"Run {k}")
        ax2.plot(times, phagocyte_result, label=f"Run {k}")
        ax3.plot(times, damage_result, label=f"Run {k}")
        ax4.plot(times, mediator_result, label=f"Run {k}")

    pdf.savefig(fig1)
    pdf.savefig(fig2)
    pdf.savefig(fig3)
    pdf.savefig(fig4)

    pdf.close()
    plt.show()

