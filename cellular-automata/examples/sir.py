try:
    from ca import *

except ModuleNotFoundError:
    from sys import path
    path.insert(0, '..')
    from ca import *
    
from time import sleep
from sys import argv
from random import random 

#0: Infectado 
#1: Suscetivel
#2: Recuperado 
# neighbors8(self, x, y, old=True, pos=False)
INFECTADO = 0
SUSCETIVEL = 1
IMUNE = 2

class SIR(CA):

    def rule(self, x, y):
        s = self[x, y]
        n = neighbors8(self, x, y)

        infectados = [i for i in neighbors8(self, x, y) 
                          if i == INFECTADO]

        p_infeccao = 0.5 
        p_recuperacao = 0.3
        p = random()

        if s == SUSCETIVEL:
            if len(infectados) >= 1 and p < p_infeccao: 
                return INFECTADO 
            else: 
                return SUSCETIVEL
        elif s == INFECTADO and p < p_recuperacao: 
            return IMUNE
        else:
            return s 

c = SIR(30, values=[SUSCETIVEL]*50 + [INFECTADO]*2 + [IMUNE])

plot(c, N=50, out='sir.pdf', graphic=True, vmax=2,
    colors=['red', 'green', 'orange'],
    names=['infected', 'susceptible', 'immune']
)