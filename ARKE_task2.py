
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib import cm

def model(T,t):

    T=T.reshape(-1,1)

    dx = 0.1  # spacing

    global n
    global Tb

    k = 2000
    rho = 3100
    Cp = 750
    alpha = k / (rho * Cp)

    #setting up A in A*T+B
    A = np.zeros([n*n, n*n])

    #setting up sparse matrix

    for i in range(0, n*n):
        for j in range(0, n*n):
            if i==j:
                A[i,j]=-4
            elif i+1==j or i-1==j:
                A[i, j]=1
            elif i+n==j or i-n==j:
                A[i,j]=1

    #changing values in sparse matrix
    for i in range(0, n*n):
        for j in range(0, n*n):
            if i % n==n-1 and j%n==0:
                A[i,j] = 0
            elif (i-1)%n==n-1 and (j-1)%n==n-2:
                A[i,j] = 0

    A1=A*alpha/(dx**2)

    b=np.zeros([n*n,1])
    for i in range (1,n*n+1):
        if i==1 or i==n or i==n*n or i==n**2-n+1:
            b[i-1]=2;
        elif i%n==1 or i%n==0 or i<n and i>1 or i>n*n-n and i<n*n:
            b[i-1]=1;
        else:
            b[i-1] = 0;
    B=alpha*Tb*b/(dx**2);
    A1_T=np.matmul(A1,T)
    dTdt=A1_T+B;

    dTdt=dTdt.reshape(1,-1)
    dTdt = dTdt.flatten()

    return dTdt

# Temp
T0=20
Tb=500

n=4 #number of nodes

start_time=0
final_time=200

for i in range (start_time,final_time,10):

    Initial_Temp=np.linspace(T0,T0,n*n)
    heating_time=np.linspace(0,i)
    T=odeint(model,Initial_Temp,heating_time)

    fig=plt.figure()
    ax=fig.add_subplot(111)

    vmin=0
    vmax=Tb

    x = np.arange(0, n)
    y = np.arange(0, n)
    T_sol = T[n, :]
    z = T_sol.reshape(n, n)
    levels = np.linspace(0, Tb, 10)
    plot=ax.contourf(x, y, z, levels)

    cbar = fig.colorbar(plot)
    cbar.set_label("Temperature / oC")

    plt.show()

    if i < final_time:
        plt.clf()