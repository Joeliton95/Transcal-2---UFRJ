from matplotlib.colors import Normalize
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from scipy.integrate import solve_bvp
import time

#constantes
#T0= 100 
Tout = 100
Rf = 200
Rc = Rf+100
kf = 0.478
kc = 0.56165
Sn0 = 11.95*(10**6) #taxa volumétrica de geração de calor
b = 1 #constante positiva sem dimensão
rpoints = 9
rfpoints = 0

#pontos em cada coordenada esférica
r = np.linspace(0,Rc,rpoints)
phi = np.linspace(0,np.pi, 3*rpoints)
theta = np.linspace(0,2*np.pi, 3*rpoints)

xf = []
yf = []
zf = []
xc = []
yc = []
zc = []
Tf = []
Tc = []

# Cria o gráfico 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_box_aspect([1, 1, 1])  # Manter proporções iguais para todos os eixos
ax.set_xticks([])             # Remover números do eixo X
ax.set_yticks([])             # Remover números do eixo Y
ax.set_zticks([])             # Remover números do eixo Z
ax.grid(False)   
aux = 0
aux2 = []

# converte os pontos de coordenada esferica para cilindrica e calcula a temperatura no regime permanente
for i in r:
    for j in phi:
        for k in theta:
            if i<=Rf:    
                xf.append(i*np.sin(j)*np.cos(k))
                yf.append(i*np.sin(j)*np.sin(k))
                zf.append(i*np.cos(j))
                
                T1 = (((Sn0*Rf**2)/(6*kf))*((1-(i/Rf)**2) + (3/10)*b*(1-(i/Rf)**4)) + ((Sn0*Rf**2)/(3*kc))*(1+(3/5)*b)*(1 - Rf/Rc) + Tout)
                Tf.append(T1)
                rfpoints = aux +1

            else:
                xc.append(i*np.sin(j)*np.cos(k))
                yc.append(i*np.sin(j)*np.sin(k))
                zc.append(i*np.cos(j))

                T2 = ((Sn0*Rf**2)/(3*kc))*(1+(3/5)*b)*(Rf/i - Rf/Rc) + Tout
                Tc.append(T2)
    aux = aux + 1


Tt = np.concatenate([Tf, Tc])
xt = np.concatenate([xf, xc])
yt = np.concatenate([yf, yc])
zt = np.concatenate([zf, zc])

T_tempo = []
T_tempo2 = []

# Resolução das equaçoes diferenciais por meio da funcao  solve_bvp do modulo scipy.integrate
def dTf_dr(rf, T):
    return np.array([(-np.exp(t) / kf) * (rf/3 + (b / Rf**2) * (rf ** 3) / 5)])

def bc(Ta, Tb):
   
    return np.array([Tb[0] - solution2.y[0][0]])

def dTc_dr(rc, T):
    return np.array([(-np.exp(t) / kc) *(1/3 + b / 5) * (Rf**3) /rc**2])

def bc2(Ta, Tb):
   
    return np.array([Tb[0] - Tout])


raiof = np.linspace(0, Rf, rfpoints)
T_initial = np.ones((1, raiof.size)) 

raioc = np.linspace(Rf, Rc, rpoints-rfpoints)
T_initial2 = np.ones((1, raioc.size)) 

# Resolve as equaçoes diferenciais de primeira ordem, com variavel depende r.
# Porem utiliza o termo de geracao de calor depedente do tempo por meio de uma exponencial "np.exp(t)".

for i in range(300):
    t = i/30
    num_repeats = len(phi)*len(theta)
    solution2 = solve_bvp(dTc_dr, bc2, raioc, T_initial2)

    #repeated_array2 = np.repeat(solution2.y[0], num_repeats)

    if len(solution2.y[0]) == len(raioc):
        repeated_array2 = np.repeat(solution2.y[0], num_repeats)
    elif len(solution2.y[0]) > len(raioc):
        diferenca = len(solution2.y[0]) - len(raioc)
        repeated_array2 = np.repeat(solution2.y[0][diferenca:], num_repeats)
    
    else:
        diferenca = -len(solution2.y[0]) + len(raioc)
        for i in range(diferenca):
            solucao_aux = np.insert(solution2.y[0], 0, solution2.y[0][0])
            repeated_array2 = np.repeat(solucao_aux, num_repeats)
    
    
    repeated_list2 = repeated_array2.tolist()
    T_tempo2.append(repeated_list2)


    solution = solve_bvp(dTf_dr, bc, raiof, T_initial)
    # Usar numpy.repeat para gerar a nova lista
    if len(solution.y[0]) == len(raiof):
        repeated_array = np.repeat(solution.y[0], num_repeats)
    elif len(solution.y[0]) > len(raiof):

        diferenca = len(solution.y[0]) - len(raiof)
        repeated_array = np.repeat(solution.y[0][diferenca:], num_repeats)

    else:

        diferenca = -len(solution.y[0]) + len(raiof)
        for i in range(diferenca):
            solucao_aux = np.insert(solution.y[0], 0, solution.y[0][0])
            repeated_array = np.repeat(solucao_aux, num_repeats)

        #repeated_array = np.repeat(np.insert(solution.y[0], 0, solution.y[0][0]), num_repeats)


     #Converter para lista, se necessário
    repeated_list = repeated_array.tolist()
    
    T_tempo.append(repeated_list)
        
    print("\nLen Tf:\n" + str(len(solution.y)) + "\n"+ "Len Tc:\n" + str(len(solution2.y)) + "\n")

    # print("\nRaio F:\n" + str(solution.x) + "\n"+ "Temperatura F:\n" + str(solution.y[0]) + "\n")
    # print("\nRaio C:\n" + str(solution2.x) + "\n"+ "Temperatura C:\n" + str(solution2.y[0]) + "\n")

#print("\nRaio F:\n" + str(solution.x) + "\n"+ "Temperatura F:\n" + str(solution.y[0]) + "\n")
#print("\nRaio C:\n" + str(solution2.x) + "\n"+ "Temperatura C:\n" + str(solution2.y[0]) + "\n")

# Supondo que T_tempo e T_tempo2 são listas
T_tempo_np = np.array(T_tempo[-1])  # Converter para NumPy array
T_tempo2_np = np.array(T_tempo2[0])  # Converter para NumPy array
# Encontrar o intervalo comum (mínimo e máximo globais)

min_value = min(T_tempo_np.min(), T_tempo2_np.min())
max_value = max(T_tempo_np.max(), T_tempo2_np.max())

# Criar normalização compartilhada
norm = Normalize(vmin=min_value, vmax=max_value)

# Plota os graficos 3D
scatter = ax.scatter(xf, yf, zf, c = T_tempo[0], cmap='coolwarm', norm=norm, s=1)
scatter2 = ax.scatter(xc, yc, zc, c = T_tempo2[0], cmap='coolwarm', norm=norm, s=1)

ax.set_xlim([-Rc, Rc])
ax.set_ylim([-Rc, Rc])
ax.set_zlim([-Rc, Rc])
ax.set_box_aspect([1, 1, 1])  # Manter proporções iguais para todos os eixos
ax.set_xticks([])             # Remover números do eixo X
ax.set_yticks([])             # Remover números do eixo Y
ax.set_zticks([])             # Remover números do eixo Z
ax.grid(False) 

plt.colorbar(scatter, ax=ax)

elev = 20  # Ângulo de elevação inicial

def updates(frame):
        
        #atualiza as cores no grafico 3D
        scatter.set_array(T_tempo[frame])
        scatter2.set_array(T_tempo2[frame])
        ax.view_init(elev, frame)  # Incremento do azimute para girar a esfera

        return scatter2,

 

ani = FuncAnimation(fig, updates, frames=len(T_tempo), interval=1)

plt.show()