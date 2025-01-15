import math
import random

import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages

pdf = PdfPages("pred_pray.pdf")

tf = 25 #Tempo final

r = 0.2 #Taxa reprodução presa
a = 0.008 #Taxa de predação
m = 0.01 #Taxa de mortalidade dos predadores

fig1 = plt.figure(figsize=(6, 4))
plt.title('Presas')
plt.xlabel('Tempo')
plt.ylabel('População')
ax1 = fig1.gca()

fig2 = plt.figure(figsize=(6, 4))
plt.title('Predadores')
plt.xlabel('Tempo')
plt.ylabel('População')
ax2 = fig2.gca()

k = 0
colorindex = 0

while k < 10:
    random.seed(k)
    times = []
    h_result = []
    p_result = []

    H = 50 # Presas
    P = 5 # Predadores
    t = 0

    h_result.append(H)
    p_result.append(P)
    times.append(0)

    while t < tf:
      
        if H == 0.0 and P == 0.0:
            break

        p1 = r*H
        p2 = a*H*P
        p3 = m*P
        S = p1 + p2 + p3
        
        ran_tau = random.uniform(0, 1)
        tau = -math.log(ran_tau) / S
        t = t + tau

        # Escolha do evento
        p = random.uniform(0, 1) 

        if p <= p1/S:
            H = H + 1
        elif p <= (p1 + p2)/S:
            H = H - 1
            P = P + 1
        else:
            P = P - 1

        h_result.append(H)
        p_result.append(P)
        times.append(t)
    
    ax1.plot(times, h_result, color='C'+str(colorindex), label="H")
    ax2.plot(times, p_result, color='C'+str(colorindex) , label="P")
    
    colorindex += 1
    k += 1

pdf.savefig(fig1)
pdf.savefig(fig2)
pdf.close()
plt.show()