# -*- coding: utf-8 -*-

# Pacotes utilizados
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# LAA 27/8/24

#**************************************************************************
#* Simulacao de um sistema exemplo, o atrator Corda


def dvSistema(x, ux, uy, t, par):
    #
    # function dvSistema(x, ux, uy, t, par)
    # simula como exemplo um sistema caotico que gera o atrator corda
    # x vetor de estado
    # como este sistema nao tem entradas: ux=uy=0
    # xd eh a derivada temporal de x, ou seja, eh o campo vetorial em x


    
    xd=np.zeros(3)
    
    # 1a Equacao diferencial 
    xd[0] = -x[1] - x[2] - par[0] * (x[0] - par[2]) 
     
    # 2a Equacao diferencial
    xd[1] = x[0] * x[1] - par[1] * x[0] * x[2] - x[1] + par[3]
    
    # 3a Equacao diferencial 
    xd[2] = par[1] * x[0] * x[1] + x[0] * x[2] - x[2]
    
    xdot=xd
    
    return xdot


#%%

# Funcao rkSistema
def rkSistema(x0, ux, uy, h, t, par):
    # function x=rkSistema(x0, ux, uy, h, t, par)
    #
    # Implementa a integracao numerica usando Runge-Kutta de 4a ordem
    # x0 eh o vetor de estado ANTES de chamar a funcao, i.e. a condicao inicial a cada passo de integracao 
    # ux e uy se forem diferentes de zero sao forcas externas como, por exemplo, entradas ou acoes de controle 
    # assume-se que as entradas nao mudam durante um passo de integracao, h
    # as equacoes de estado estao implementadas em dvSistema.m
    # h intervalo de integracao
    # t eh o tempo antes de chamar a funcao.
    
    # 1a avaliacao 
    xd = dvSistema(x0,ux,uy,t,par)
    savex0 = x0
    phi = xd
    x0 = savex0+0.5*h*xd
    
    # 2a avaliacao
    xd = dvSistema(x0,ux,uy,t+0.5*h,par)
    phi = phi+2*xd
    x0 = savex0+0.5*h*xd
    
    # 3a avaliacao
    xd = dvSistema(x0,ux,uy,t+0.5*h,par)
    phi = phi+2*xd
    x0 = savex0+h*xd
    
    # 4a avaliacao
    xd = dvSistema(x0,ux,uy,t+h,par)
    x = savex0+(phi+xd)*h/6
    
    return x


#%%
    

a = 0.258;
b = 4.033;
F = 8;
G = 1;
Kas = np.array([a, b, F, G])

# Condicoes Iniciais 
x0 = np.array([0.1, 0.1, 0.1])


t0 = 0 # tempo inicial
tf = 150 # tempo final
# intervalo de integracão
h = 0.01
t = np.arange(t0, tf, h) # vetor de tempo


# inicializacao
x = np.zeros((len(x0), len(t)))
u1 = np.zeros((len(t)))
u2 = np.zeros((len(t)))

    
for i in range(1, len(t)):
    x[:, i] = rkSistema(x[:, i-1], u1[i-1], u2[i-1], h, t[i-1], Kas)
        
#%% 

# gera graficos 
plt.figure(1, figsize=(10, 8.107))
plt.plot(t,x[0,:],color='black',linewidth=1.5)
plt.plot(t,x[1,:],color='blue',linewidth=1.5)
plt.plot(t,x[2,:],color='red',linewidth=1.5)
plt.xlabel('t', fontsize=17.5)
plt.ylabel('vaiaveis de estado', fontsize=17.5, labelpad=15)
plt.show()


fig = plt.figure(2, figsize=(10, 8.107))
ax = fig.add_subplot(111, projection='3d')
ax.plot(x[0,:],x[1,:],x[2,:])

# Add labels
ax.set_xlabel('X Axis')
ax.set_ylabel('Y Axis')
ax.set_zlabel('Z Axis')
plt.show()

