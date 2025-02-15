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
kpg = 0.3  # Taxa de crescimento da população de patógenos
sm = 0.005  # Parâmetro relacionado à resistência do organismo (possivelmente relacionado à imunidade ou à resposta do sistema)
kpn = 1.8  # Taxa de proliferação dos patógenos em função da população de fagócitos
micro_n = 0.05  # Taxa de produção de fagócitos
sc = 0.0125  # Taxa de aumento dos mediadores anti-inflamatórios
micro_c = 0.1  # Taxa de produção de mediadores anti-inflamatórios
knp = 0.1  # Taxa de proliferação dos patógenos devido à interação com os fagócitos
pinf = 20  # População máxima de patógenos
micro_m = 0.002  # Taxa de crescimento de mediadores
snr = 0.08  # Taxa de interação entre fagócitos e mediadores
kdn = 0.35  # Taxa de produção de dano tecidual devido à interação entre fagócitos e mediadores
kcn = 0.04  # Taxa de modulação de mediadores
cinf = 0.28  # Nível máximo de mediadores anti-inflamatórios
knd = 0.02  # Taxa de modulação do dano tecidual
kpm = 0.6  # Taxa de modulação do patógeno pelos mediadores
kmp = 0.01  # Taxa de modulação do mediador pelos patógenos
micro_nr = 0.12  # Taxa de produção de células fagocíticas
micro_d = 0.02  # Taxa de produção de dano tecidual
kcnd = 48  # Parâmetro de modulação entre dano tecidual e mediadores
xdn = 0.06  # Parâmetro relacionado ao dano tecidual
knn = 0.01  # Taxa de interação entre fagócitos e mediadores
muj = 0.12  # Taxa de degradação de fagócitos
mud = 0.02  # Taxa de degradação do dano tecidual
mu_c = 0.1  # Taxa de degradação dos mediadores
sigma = 0.02  # Desvio padrão do ruído estocástico aplicado no modelo

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

