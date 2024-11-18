import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

#constantes da equação
Tout = 100 # temperatura fora da esfera
Rf = 200 # raio de f
Rc = Rf+100 # raio de c
kf = 0.478  # condutividade térmica f
kc = 0.56165 # condutividade térmica c
Sn0 = 11.95*(10**6) #taxa volumétrica de geração de calor
b = 1 #constante positiva sem dimensão
rpoints = 13 #pontos no raio

#pontos em cada coordenada esférica
r = np.linspace(0,Rc,rpoints)
phi = np.linspace(0,np.pi, 3*rpoints)
theta = np.linspace(0,(3/2)*np.pi, 3*rpoints)

xf = []
yf = []
zf = []
xc = []
yc = []
zc = []
Tf = []
Tc = []
Traio = []

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# converte os pontos de coordenada esferica para cilindrica e calcula a temperatura no regime permanente
for idx, i in enumerate(r):
    for j in phi:
        for k in theta:
            if i<=Rf:    
                xf.append(i*np.sin(j)*np.cos(k))
                yf.append(i*np.sin(j)*np.sin(k))
                zf.append(i*np.cos(j))
                
                T1 = (((Sn0*Rf**2)/(6*kf))*((1-(i/Rf)**2) + (3/10)*b*(1-(i/Rf)**4)) + ((Sn0*Rf**2)/(3*kc))*(1+(3/5)*b)*(1 - Rf/Rc) + Tout)
                Tf.append(T1)
            else:
                xc.append(i*np.sin(j)*np.cos(k))
                yc.append(i*np.sin(j)*np.sin(k))
                zc.append(i*np.cos(j))

                T2 = ((Sn0*Rf**2)/(3*kc))*(1+(3/5)*b)*(Rf/i - Rf/Rc) + Tout
                Tc.append(T2)


    if i<=Rf:    
        Traio.append(((Sn0*Rf**2)/(6*kf))*((1-(i/Rf)**2) + (3/10)*b*(1-(i/Rf)**4)) + ((Sn0*Rf**2)/(3*kc))*(1+(3/5)*b)*(1 - Rf/Rc) + Tout)
    else:
        Traio.append(((Sn0*Rf**2)/(3*kc))*(1+(3/5)*b)*(Rf/i - Rf/Rc) + Tout)
    
            
Tt = np.concatenate([Tf, Tc])
xt = np.concatenate([xf, xc])
yt = np.concatenate([yf, yc])
zt = np.concatenate([zf, zc])

ax.set_box_aspect([1, 1, 1])  # Manter proporções iguais para todos os eixos
ax.set_xticks([])             # Remover números do eixo X
ax.set_yticks([])             # Remover números do eixo Y
ax.set_zticks([])             # Remover números do eixo Z
ax.grid(False)  

# plota o gráfico 3D, com os pontos da esfera e o grafiente de temperatura

scat = ax.scatter(xt, yt, zt, c = Tt, s=1, cmap="coolwarm")
ax.set_xlim([-Rc, Rc])
ax.set_ylim([-Rc, Rc])
ax.set_zlim([-Rc, Rc])
ax.set_box_aspect([1, 1, 1])
plt.colorbar(scat, ax=ax)

#rotaciona a esfera e cria a animação
elev = 10
def calculate_temperatures(frame):

   ax.view_init(elev, -frame)

ani = FuncAnimation(fig, calculate_temperatures, frames=np.arange(0, 90), interval=100)

# plota grafico 2D da temperatura pelo raio
plt.figure() 
plt.plot(r, Traio, color="black")  # Gráfico da função seno
plt.xlabel('Raio')             # Rótulo do eixo X
plt.ylabel('Temperatura') 

plt.show()