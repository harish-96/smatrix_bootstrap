from cvxpy import *
import numpy as np

n = 50
an, bn, cn = Variable(n), Variable(n), Variable(n)

def sigma_1(s, an, bn, cn):
    ch = [s**n for n in range(n)]
    return sum((an+1j*bn)*ch)

def sigma_2(s, an, bn, cn):
    ch = [s**n for n in range(n)]
    return sum(cn*ch)

def sigma_3(s, an, bn, cn):
    ch = [s**n for n in range(n)]
    return sum((an-1j*bn)*ch)

def sigma_sing(s, an, bn, cn, D=26):
    return (D-2)*sigma_1(s,an,bn,cn) + sigma_2(s,an,bn,cn) + sigma_3(s,an,bn,cn)

def sigma_sym(s, an, bn, cn, D=26):
    return sigma_2(s,an,bn,cn) + sigma_3(s,an,bn,cn)

def sigma_anti(s, an, bn, cn, D=26):
    return sigma_2(s,an,bn,cn) - sigma_3(s,an,bn,cn)

objective = Minimize(an[1])
constraints = [1 >= abs(sigma_sing(-4*1j*(np.exp(1j*x)-1)/(np.exp(1j*x)+1), an, bn, cn)) for x in np.linspace(0,np.pi,100,endpoint=False)] \
               + [1 >= abs(sigma_sym(-4*1j*(np.exp(1j*x)-1)/(np.exp(1j*x)+1), an, bn, cn)) for x in np.linspace(0,np.pi,100,endpoint=False)] \
               + [1 >= abs(sigma_anti(-4*1j*(np.exp(1j*x)-1)/(np.exp(1j*x)+1), an, bn, cn)) for x in np.linspace(0,np.pi,100,endpoint=False)]
prob = Problem(objective, constraints)
prob.solve(verbose=True)
