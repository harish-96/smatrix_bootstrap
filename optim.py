from cvxpy import *
import numpy as np

n = 50 #Order of taylor expansion
D = 4 #Dimension of ambient spacetime
N_constr = 200 #Number of points on the circle at which constraints are applied

alpha2 = (D-26)/384/np.pi
beta3 = 0

an, bn, cn = Variable(n), Variable(n), Variable(n)
ns = np.arange(n)

#Coeffs in s plane in terms of coeffs in \chi plane
a0_t, b0_t = sum(an), sum(bn)

a1_t, b1_t = sum(-ns/2*bn), sum(ns/2*an)

a2_t, b2_t = sum(-ns**2/8*an), sum(-ns**2/8*bn)

a3_t, b3_t = sum(ns*(2*ns**2+1)/96*bn), sum(-ns*(2*ns**2+1)/96*an)

c0_t = sum(cn)

c1_t_by_i = sum(ns/2*cn) #c1_t/1j

c2_t = sum(-ns**2/8*cn)

c3_t_by_i = sum(-ns*(2*ns**2+1)/96*cn) #c3_t/1j

def sigma_1(chi, an, bn, cn):
    ch = chi**ns
    return sum((an+1j*bn)*ch)

def sigma_2(chi, an, bn, cn):
    ch = chi**ns
    return sum(cn*ch)

def sigma_3(chi, an, bn, cn):
    ch = chi**ns
    return sum((an-1j*bn)*ch)

def sigma_sing(chi, an, bn, cn, D=4):
    return (D-2)*sigma_1(chi,an,bn,cn) + sigma_2(chi,an,bn,cn) + sigma_3(chi,an,bn,cn)

def sigma_sym(chi, an, bn, cn, D=4):
    return sigma_2(chi,an,bn,cn) + sigma_3(chi,an,bn,cn)

def sigma_anti(chi, an, bn, cn, D=4):
    return sigma_2(chi,an,bn,cn) - sigma_3(chi,an,bn,cn)

alpha3 = c3_t_by_i - beta3
objective = Minimize(alpha3)

def cheby_pts(n):
    k = np.arange(n)
    return np.cos(np.pi*(2*k+1)/(2*n+2))
constr_pts = (cheby_pts(N_constr)+1)*np.pi/2
# constr_pts = np.linspace(0, np.pi, N_constr, endpoint=False)

constraints = [a0_t == 0, b0_t == 0, a1_t == 0, b1_t ==0, a2_t == 0, b2_t == -alpha2, a3_t == alpha2/4] \
               + [b3_t == -beta3, c0_t == 1, c1_t_by_i == 0.25, c2_t == -1/32] \
               + [(abs(sigma_sing(np.exp(1j*x), an, bn, cn, D)) <= 1) for x in constr_pts] \
               + [(abs(sigma_sym(np.exp(1j*x), an, bn, cn, D)) <= 1) for x in constr_pts] \
               + [(abs(sigma_anti(np.exp(1j*x), an, bn, cn, D)) <= 1) for x in constr_pts]

prob = Problem(objective, constraints)
prob.solve(solver=GUROBI, BarQCPConvTol=1e-50, verbose=True)
# prob.solve(solver=ECOS, verbose=True)
