## algoritmo per risolvere numericamente problemi a piu gradi di liberta
import numpy as np
from math import cos
from scipy.linalg import eigh
from scipy.linalg import inv
from scipy.linalg import det
from matplotlib import pyplot as plt

#setup parameters (in SI units)
m = 10
M = 50 
r = 0.5
J0 = 0.5*M*(r**2)
F0 = 5
k = 1000
omega = 1
dof = 3
end_time = 60
time_step = 1.0e-3

#setup matrices
W = np.array([[M+(J0/(9*(r**2))),(-J0/(9*(r**2))),0],[((-J0)/(9*(r**2))),((3*m)+((J0)/(9*(r**2)))),0],[0,0,m]]) #mass matrix
print("This is the weight matrix:")
print(W)
K = np.array([[(41*k/9),(-8*k/9),(-8*k/3)],[(-8*k/9),(2*k/9),(2*k/3)],[(-8*k/3),(2*k/3),(5*k)]]) #stiffness matrix
print("This is the stiffness matrix:")
print(K)
I = np.identity(dof)
A = np.zeros((2*dof,2*dof)) 
B = np.zeros((2*dof,2*dof))
Y = np.zeros((2*dof,1))

A[0:3,0:3] = W #la sottomatrice (012)x(012) e M 
A[3:6,3:6] = I #la sottomatrice (345)x(345) e I
B[0:3,3:6] = K #la sottomatrice (012)x(345) e K
B[3:6,0:3] = -I #la sottomatrice (345)x(012) e -I

#find natural frequencies and the modal forms of the system
ddet_W = det(W)
evals, evecs = eigh(K,W)
print("Those are the modal forms: " + str(evecs))
print("Those are the eigenvalues: ")
print(evals)
nat_frequencies = np.sqrt(evals) #nat_frequencies is a real array because K and W are positive-definite
print("Those are the natural frequencies of the system: " + str(nat_frequencies))

#Define a runge-kutta integrator (order 4)
A_inv = inv(A)
force=[]
X1 = []
X2 = []
X3 = []
V1 = []
V2 = []
V3 = []

def F(t):
    F = np.zeros((2*dof,1))
    F[:,0] = F0*cos(omega*t)
    return F

def G(Y,t): 
  return A_inv.dot((F(t) - B.dot(Y)))

def RK4(Y, t, dt):
  k1 = G(Y,t)
  k2 = G(Y+0.5*k1*dt, t+0.5*dt)
  k3 = G(Y+0.5*k2*dt, t+0.5*dt)
  k4 = G(Y+k3*dt, t+dt)

  return dt * (k1 + 2*k2 + 2*k3 + k4) / 6

for t in np.arange(0, end_time, time_step):
    Y= Y + RK4(Y, t, time_step)
    force.append(F(t)[1])
    X1.append(Y[3])
    X2.append(Y[4])
    X3.append(Y[5])
    V1.append(Y[0])
    V2.append(Y[1])
    V3.append(Y[2])
    
#plot results
time = [round(t,5) for t in np.arange(0, end_time, time_step)]
plt.figure(dpi=600)
plt.plot(time,X1, color="dodgerblue", linewidth=0.8, markersize=7)
plt.plot(time,V1, color="firebrick", linewidth=0.8, markersize=7)
plt.xlabel('time [s]')
plt.title('Response 1')
plt.legend(['X1', 'V1'], loc='upper right')
plt.show()

plt.figure(dpi=600)
plt.plot(time,X2, color="dodgerblue", linewidth=0.8, markersize=7)
plt.plot(time,V2, color="firebrick", linewidth=0.8, markersize=7)
plt.xlabel('time [s]')
plt.title('Response 2')
plt.legend(['X2', 'V2'], loc='upper right')
plt.show()

plt.figure(dpi=600)
plt.plot(time,X3, color="dodgerblue", linewidth=0.8, markersize=7)
plt.plot(time,V3, color="firebrick", linewidth=0.8, markersize=7)
plt.xlabel('time [s]')
plt.title('Response 3')
plt.legend(['X3', 'V3'], loc='upper right')
plt.show()

plt.figure(dpi=600)
plt.plot(time,X1, linewidth=0.8, markersize=7)
plt.plot(time,X2, linewidth=0.8, markersize=7)
plt.plot(time,X3, linewidth=0.8, markersize=7)
plt.plot(time,V1, linewidth=0.8, markersize=7)
plt.plot(time,V2, linewidth=0.8, markersize=7)
plt.plot(time,V3, linewidth=0.8, markersize=7)
plt.xlabel('time [s]')
plt.legend(['X1', 'X2', 'X3', 'V1', 'V2', 'V3'], loc='upper right')
plt.show()

plt.figure(dpi=600)
plt.plot(X1,V1, color="firebrick", linewidth=0.8, markersize=7)
plt.title("Poincare map 1")
plt.show()

plt.figure(dpi=600)
plt.plot(X2,V2, color="firebrick", linewidth=0.8, markersize=7)
plt.title("Poincare map 2")
plt.show()

plt.figure(dpi=600)
plt.plot(X3,V3, color="firebrick", linewidth=0.8, markersize=7)
plt.title("Poincare map 3")
plt.show()

plt.figure(dpi=600)
plt.plot(time,force, linewidth=0.8, markersize=7)
plt.xlabel('time [s]')
plt.ylabel('F [N]')
plt.show()

    
#find the solution as a linear combination of the modal forms