try:
    from ca import *

except ModuleNotFoundError:
    from sys import path
    path.insert(0, '..')
    from ca import *


class MeuCA(CA):
    def rule(self, x, y):
    ...

c = MeuCA(30,values=1,random_values=True)

plot(c, N=50, out="meuca.pdf")
