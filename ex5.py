#exercicio 5 - FINAL 
def gaussian(A, w, x0, p,xs):
    u = ((A*w)/w)*np.exp(-(xs-x0)**2/(2*w**2) + (1j)*p*(xs-x0)) 
    u = u / np.sqrt((np.abs(u)**2).sum())
    return u
def sol(A, w, x0, p,xs,t):
    u = ((A*w)/np.sqrt(w**2 + 2*(1j)*t))*np.exp(-(xs-x0-2*p*t)**2/(2*w**2 + 4*(1j)*t) + (1j)*p*(xs-x0) - (1j)*(p**2)*t) 
    u = u / np.sqrt((np.abs(u)**2).sum())
    return u
def probability(state):
    prob = np.abs(state)**2
    return prob / prob.sum()
import numpy as np
import matplotlib.pyplot as plt
import os, sys
import matplotlib
matplotlib.rc('font', size=18)
matplotlib.rc('font', family='Arial')
file=open('a.dat','w')
N = 2001 
dt = 5.e-5 
L = float(50) 
nsteps = 1000 
dx = L/(N-1)
nplot = 20 
alpha = ((1j)*dt)/(2*dx**2) 
A = np.zeros((N,N),dtype='complex')
B = np.zeros((N,N),dtype='complex')
#define matrices A, B and b array
for i in range(N):
    if i==0:
        A[i,:] = [1+2*alpha if j==0 else -alpha if j==1 or j==N-1 else 0 for j in range(N)]
        B[i,:] = [1-2*alpha if j==0 else alpha if j==1 or j==N-1 else 0 for j in range(N)]
    elif i==N-1:
        A[i,:] = [-alpha if j==0 or j==N-2 else 1+2.*alpha if j==N-1 else 0 for j in range(N)]
        B[i,:] = [alpha if j==0 or j==N-2 else 1-2.*alpha if j==N-1 else 0 for j in range(N)]
    else:
        A[i,:] = [-alpha if j==i-1 or j==i+1 else 1+2.*alpha if j==i else 0 for j in range(N)]
        B[i,:] = [alpha if j==i-1 or j==i+1 else 1-2.*alpha if j==i else 0 for j in range(N)]
x = np.linspace(0, L, N)
t= np.arange(0,dt*nsteps,dt)
u = gaussian(1.0, L/20.0, L/4.,2.0, x) 
bb = B.dot(u[:]) 
sol=np.asarray([sol(1.0, L/20.0, L/4.,2.0, x,xx) for xx in t])
inv_A = np.linalg.inv(A)
fig=plt.figure(figsize=(25.0,15.0))
ax1 = fig.add_subplot(121)
ax1.set_xlim(left=0,right=L)
ax1.set_ylim(bottom=-0.22,top=0.22)
ax1.tick_params(axis='x', labelsize=24)
ax1.tick_params(axis='y', labelsize=24)
ax1.set_title('(a) Simulação', size=24)
ax1.plot(x,u.real,linewidth=2)
ax1.plot(x,u.imag,linewidth=2)
ax1.plot(x,probability(u),linewidth=2)
ax2 = fig.add_subplot(122)
ax2.set_xlim(left=0,right=L)
ax2.set_ylim(bottom=-0.22,top=0.22)
ax2.tick_params(axis='x', labelsize=24)
ax2.tick_params(axis='y', labelsize=24)
ax2.set_title('(a) Solução exata', size=24)
ax2.plot(x,sol[0].real,linewidth=2)
ax2.plot(x,sol[0].imag,linewidth=2)
ax2.plot(x,probability(sol[0]),linewidth=2)
#ax1.legend(['Real','imag','prob'],prop={'size':10})
filename = 'foo' + str(c+1).zfill(3) + '.jpg';
plt.savefig(filename)
plt.clf()
tempo =0.
time=[]
vv=[]
c = 0
for j in range(nsteps):
    u[:] = inv_A @ bb
    v=0.0
    for k in range(len(u)):
        v = v + abs(u[k])
    tempo = tempo+dt
    vv.append(v)
    time.append(tempo)
    print(v)
    print(j)
    bb = B.dot(u[:]) 
    if(j%nplot==0):
        fig=plt.figure(figsize=(25.0,15.0))
        ax1 = fig.add_subplot(121)
        ax1.set_xlim(left=0,right=L)
        ax1.set_ylim(bottom=-0.22,top=0.22)
        ax1.tick_params(axis='x', labelsize=24)
        ax1.tick_params(axis='y', labelsize=24)
        ax1.set_title('(a) Simulação', size=24)
        ax1.plot(x,u.real,linewidth=2)
        ax1.plot(x,u.imag,linewidth=2)
        ax1.plot(x,probability(u),linewidth=2)
        ax2 = fig.add_subplot(122)
        ax2.set_xlim(left=0,right=L)
        ax2.set_ylim(bottom=-0.22,top=0.22)
        ax2.tick_params(axis='x', labelsize=24)
        ax2.tick_params(axis='y', labelsize=24)
        ax2.set_title('(a) Solução exata', size=24)
        ax2.plot(x,sol[j].real,linewidth=2)
        ax2.plot(x,sol[j].imag,linewidth=2)
        ax2.plot(x,probability(sol[j]),linewidth=2)
        ax1.legend(['Real','imag','prob'],prop={'size':10})
        ax2.legend(['Real','imag','prob'],prop={'size':10})
        filename = 'foo' + str(c+1).zfill(3) + '.jpg';
        plt.savefig(filename)
        plt.clf()
        c += 1

#os.system("ffmpeg -r 5 -y -i 'foo%03d.jpg' heat_equation.m4v")
#os.system("rm -f *.jpg")
