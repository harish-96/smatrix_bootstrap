from cvxpy import *
import numpy as np

n = 50
D = 4
alpha2 = (D-26)/384/np.pi
beta3 = 0

c4 = Variable()
an, bn, cn = Variable(n-5), Variable(n-5), Variable(n-5)

a04 = [0,0,0,-alpha2,-beta3]
b04 = [0,0,0,0,0]
c03 = [0,1,0.25,0]


def sigma_1(s, an, bn, cn):
    sn = [s**n for n in range(5,n)]
    return sum((an+1j*bn)*sn) + np.dot(a04, np.array((1,s,s**2,s**3,s**4)))

def sigma_2(s, an, bn, cn):
    sn = [s**n for n in range(5,n)]
    return sum(cn*sn) + np.dot(c03, np.array((1,s,s**2,s**3))) + c4*s**4

def sigma_3(s, an, bn, cn):
    s = -s
    sn = [s**n for n in range(5,n)]
    return sum((an-1j*bn)*sn) + np.dot(a04, np.array((1,s,s**2,s**3,s**4)))

def sigma_sing(s, an, bn, cn, D=4):
    return (D-2)*sigma_1(s,an,bn,cn) + sigma_2(s,an,bn,cn) + sigma_3(s,an,bn,cn)

def sigma_sym(s, an, bn, cn, D=4):
    return sigma_2(s,an,bn,cn) + sigma_3(s,an,bn,cn)

def sigma_anti(s, an, bn, cn, D=4):
    return sigma_2(s,an,bn,cn) - sigma_3(s,an,bn,cn)

objective = Minimize(c4 - 2*beta3)
constraints = [1 >= abs(sigma_sing(-4*1j*(np.exp(1j*x)-1)/(np.exp(1j*x)+1), an, bn, cn, D)) for x in np.linspace(0,np.pi,100,endpoint=False)] \
               + [1 >= abs(sigma_sym(-4*1j*(np.exp(1j*x)-1)/(np.exp(1j*x)+1), an, bn, cn, D)) for x in np.linspace(0,np.pi,100,endpoint=False)] \
               + [1 >= abs(sigma_anti(-4*1j*(np.exp(1j*x)-1)/(np.exp(1j*x)+1), an, bn, cn, D)) for x in np.linspace(0,np.pi,100,endpoint=False)]
prob = Problem(objective, constraints)
prob.solve(verbose=True)
