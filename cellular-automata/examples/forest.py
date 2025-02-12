try:
    from ca import *

except ModuleNotFoundError:
    from sys import path
    path.insert(0, '..')
    from ca import *
    
from time import sleep
from sys import argv
from random import random  

from random import randint

#0: Infectado 
#1: Suscetivel
#2: Recuperado 
# neighbors8(self, x, y, old=True, pos=False)
VAZIO = 0
FLORESTA = 1
QUEIMANDO = 2
QUEIMADO = 3

pf = 0.8 #fire probability
#p= 0 
#p= 1 fogo se espalha (com certeza)

class FOREST(CA):

    def rule(self, x, y):
        s = self[x, y]
        n = neighbors8(self, x, y)

        queimando = [i for i in neighbors8(self, x, y) 
                          if i == QUEIMANDO]
        
        p = random()

        if s == FLORESTA :
            if len(queimando) >= 1 and p < pf: 
                return QUEIMANDO 
            else: 
                return FLORESTA
            
        elif s == QUEIMANDO: 
            return QUEIMADO
        
        else:
            return s

tam_dim = 50
arvores = 0.5 * tam_dim * tam_dim    # 80% das celulas com arvore
#arvores = 0.8 * tam_dim * tam_dim 

c = FOREST(50, values=4, random_values=False)

for i in range(int(arvores)):
    c.add(value=FLORESTA, points=[(randint(0, tam_dim), randint(0, tam_dim))])

c.add(value=QUEIMANDO, points=[(10, 10)], size=(1, 1))
c.add(value=QUEIMANDO, points=[(40, 30)], size=(1, 1))

plot(c, N=50, out='forest.pdf', graphic=True, vmax=3,
    colors=['white','green', 'red', 'black'],
    names=['VAZIO','FLORESTA','QUEIMANDO', 'QUEIMADO']
)