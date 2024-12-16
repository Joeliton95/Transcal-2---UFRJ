from matplotlib.colors import Normalize
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from scipy.integrate import solve_bvp
import time

H = 12 #milimetros
#npoints = 90
#ne = npoints-1
#dx = (H/ne)/1000 #metros

# parametros

rhob = 1000 #kg/m3
cb = 4200 # J/(kg·K)
wb = 0.0005
Ql = 3380 #W/m3
h = 100 #W/(m2·K)
K = 0.5
Te = 20.0   # temperatura da esquerda
Ta = 37.0 # temperatura da direita

msize = [0.1, 0.7, 0.8, 2.0,8.0] #espessura das camadas
rho = [1200, 1200, 1200, 1000, 1085, 1030, 1060]  # densidade do material ρ, kg/m3
c = [3589, 3300, 3300, 2674, 3800, 3852, 3700]  # c, J/kg·K
k = [0.235, 0.445, 0.445, 0.185, 0.51, 0.558]
wbas = [0.0, 0.0002, 0.0013, 0.0001, 0.0027, 0.0063]
Qbas = [0.0, 368.1, 368.1, 368.3, 684.2, 3700]
q = 2

nIter = 800
dt = 1 # dt < dx*dx/alpha


X1 = np.linspace(0.0, 0.1, 5)
X2 = np.linspace(0.1, 0.7, 5)
X3 = np.linspace(0.7, 0.8, 5)
X4 = np.linspace(0.8, 2, 5)
X5 = np.linspace(2, H, 15)

X = np.concatenate((X1,X2,X3,X4,X5))
X = X/1000
npoints = len(X)
ne = npoints-1
dx = (H/ne)/1000

print(X)
m = [] # representa as camadas da pele

#associa cada ponto do dominio a uma camada m
for aux in range(len(X)):
 if 0<X[aux]<=0.1:
  m.append(1)
 elif 0.1<X[aux]<=0.7:
  m.append(2)
 elif 0.7<X[aux]<=0.8:
  m.append(3)
 elif 0.8<X[aux]<=2.0:
  m.append(4)
 elif 2.0<X[aux]<=H:
  m.append(5)


print(m)
print(dx)

print(X[1]-X[0])

#temperatura inicial
T = np.ones( (npoints),dtype='float')*37
# condicao de contorno

plt.ion()
A = np.zeros( (npoints,npoints),dtype='float')
b = np.zeros( (npoints),dtype='float')

for n in range(0,nIter):

 # matrix A
 #A = np.zeros( (npoints,npoints),dtype='float')
 #b = np.zeros( (npoints),dtype='float')
 # condicoes de contorno em A e b
 A[0,0] = 1
 A[-1,-1] = 1
 b[0] = (2*h*dx*Te - K*T[2] + 4*K*T[1])/(3*K + 2*h*dx)
 b[-1] = Ta

 for i in range(1,npoints-1): # loop no interior da malha

  W = wbas[m[i]-1]*q**((T[i]-Ta)/10)
  Q = Qbas[m[i]-1]*q**((T[i]-Ta)/10)


  A[i,i-1] = (-k[m[i]-1]*dt)/dx**2
  A[i,i]   = rho[m[i]-1]*c[m[i]-1] + (2*k[m[i]-1]*dt)/dx**2
  A[i,i+1] = (-k[m[i]-1]*dt)/dx**2
  b[i] = rho[m[i]-1]*c[m[i]-1]*T[i] + dt*rhob*cb*wb*(Ta-T[i]) + Ql

 # solucao do sistema linear
 T = np.linalg.solve(A,b)


 plt.plot(X,T,'ko-')
 plt.pause(0.00000005)
 plt.clf()


plt.ioff()
plt.plot(X,T,'ko-')
plt.show()