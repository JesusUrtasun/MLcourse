# Gaussian Processes
# Check N3LO data

import pdb
import numpy as np
from scipy import special as sc
from matplotlib import pyplot as plt
from matplotlib import gridspec

# Setting parameters
np.random.seed(1)
ca = 3
cf = 4 / 3
nf = 5
alphas = 0.1126
# beta0 = (33 - 2 * nf) / 3
beta0 = (33 - 2 * nf) / (12 * np.pi)

##### Large N contributions #####

# Expressions (2.17) and Table 2
g01 = (2 / np.pi) * (3 * pow(np.euler_gamma, 2) + pow(np.pi, 2) + 11/4)
g02 = 40.10
g12 = 2 * ca / np.pi
g21 = 2 * ca * np.euler_gamma * np.pi
# Used for nnlo only
g01bar = 4.9374
g02bar = 10.92
g03bar = 0 # Unknown!!!

# Appendix (A.12) - NLO
b10 = 0
b11 = 2 * g12
# Appendix (A.13) - NNLO
b20 = (1 / pow(np.pi, 2)) * (( -(101 / 27) + (11 / 3) * sc.zeta(2) + (7 / 2) * sc.zeta(3)) * pow(ca, 2) + ((14 / 27) - (2 / 3) * sc.zeta(2)) * ca * nf)
b21 = (1 / pow(np.pi, 2)) * (( (67 / 9) - 2 * sc.zeta(2)) * pow(ca, 2) - (10 / 9) * ca * nf)
b22 = (1 / pow(np.pi, 2)) * ( -(11 / 3) * pow(ca, 2) + (2 / 3) * ca * nf)
# Appendix (A.14) - N3LO
b30 = (1 / pow(np.pi, 3)) * ((-(297029 / 23328) + (6139 / 324) * sc.zeta(2) + (2509 / 108) * sc.zeta(3) - (187 / 60) * pow(sc.zeta(2), 2) - (11 / 6) * sc.zeta(2) * sc.zeta(3) - 6 * sc.zeta(5)) * pow(ca, 3)
        + ((31313 / 11664) - (1837 / 324) * sc.zeta(2) - (155 / 36) * sc.zeta(3) + (23 / 30) * pow(sc.zeta(2), 2)) * pow(ca, 2) * nf
        + ((1711 / 864) - (1 / 2) * sc.zeta(2) - (19 / 18) * sc.zeta(3) - (1 / 5) * pow(sc.zeta(2), 2)) * ca * cf * nf
        + (-(58 / 729) + (10 / 27) * sc.zeta(2) + (5 / 27) * sc.zeta(3)) * ca * pow(nf, 3))
b31 = (1 / pow(np.pi, 3)) * (pow(ca, 3) + pow(ca, 2) * nf + ca * cf * nf + ca * pow(nf, 3))
b32 = (1 / pow(np.pi, 3)) * ((-(445 / 27) + (11 / 3) * sc.zeta(2)) * pow(ca, 3) +((289 / 54) - (2 / 3) * sc.zeta(2)) * pow(ca, 2) * nf + (1 / 2) * ca * cf * nf - (10 / 27) * ca * pow(nf, 3))
b33 = (1 / pow(np.pi, 3)) * ((121 / 27) * pow(ca, 3) - (44 / 27) * pow(ca, 2) * nf + (4 / 27) * ca * cf * nf)

# Expressions D (A.6a)
def D0(n):
        return - np.euler_gamma - sc.polygamma(0, n)

def D1(n):
        return (1 / 12) * (6 * pow(np.euler_gamma, 2) + pow(np.pi, 2) + 6 * sc.polygamma(0, n) * (2 * np.euler_gamma + sc.polygamma(0, n)) - 6 * sc.polygamma(1, n))

def D2(n):
        return - (1 / 6) * np.log(n) * (6 * pow(np.euler_gamma, 2) + pow(np.pi, 2) + 2 * np.log(n) * (3 * np.euler_gamma + np.log(n)))

def D3(n):
        return sc.polygamma(2, n)

# Expressions Dlog (A.6b)
def Dlog0(n):
        return - np.log(n)

def Dlog1(n):
        return (1 / 2) * np.log(n) * (2 * np.euler_gamma + np.log(n))

def Dlog2(n):
        return - (1 / 6) * np.log(n) * (6 * pow(np.euler_gamma, 2) + pow(np.pi, 2) + 2 * np.log(n) * (3 * np.euler_gamma + np.log(n)))

def Dlog3(n):
        return sc.polygamma(2, n)

# Expressions Dhat (A.6c)
def Dhat0(n):
        return - np.euler_gamma - sc.polygamma(0, n)

def Dhat1(n):
        return 0.5 * (pow(sc.polygamma(0, n), 2) + 2 * np.euler_gamma * sc.polygamma(0, n) + sc.zeta(2) + pow(np.euler_gamma, 2))

def Dhat2(n):
        return (1 / 12) * (- sc.polygamma(2, n) - 2 * (pow(np.pi, 2) * (np.euler_gamma + sc.polygamma(0, n)) + 2 * pow(np.euler_gamma + sc.polygamma(0, n), 3) + 4 * sc.zeta(3) ))

def Dhat3(n):
        return sc.polygamma(2, n)

# Prpare soft2 ()
def get_D_soft2(d_hat):

        return lambda n, d_hat = d_hat: 2 * d_hat(n) - 3 * d_hat(n + 1) + 2 * d_hat(n + 2)

def D_general(k, n, exp_type):
        
        # first key to a list of functions
        D_conf = {"soft1": [Dhat0, Dhat1, Dhat2],
                "soft2": [get_D_soft2(Dhat0), get_D_soft2(Dhat1), get_D_soft2(Dhat2)],
                "soft0": [D0, D1, D2]}
        
        # current conf is just a list of conf names
        current_conf = D_conf[exp_type]
        # Shitf argument for 2.36a
        if exp_type == "soft1":
                n = n + 1
        return current_conf[k](n)

# S(N) expressions in (2.36)
def s1(n, exp_type):
        return b10 * D_general(0, n, exp_type) + b11 * D_general(1, n, exp_type) # (2.36a)

def s2(n, exp_type):
        return b20 * D_general(0, n, exp_type) + b21 * D_general(1, n, exp_type) + b22 * D_general(2, n, exp_type) # (2.36b)

def s3(n, exp_type):
        return b30 * D_general(0, n, exp_type) + b31 * D_general(1, n, exp_type) + b32 * D_general(2, n, exp_type) + b33 * D_general(3, n, exp_type) # (2.36a)

##### Asymptotic expressions #####

# NLO large N
def nlo_large_n(n, exp_type):

        # Soft approximations (2.40)
        csoft = s1(n, exp_type) + g01bar

        return csoft

# NNLO large N
def nnlo_large_n(n, exp_type):

        # Soft approximations (2.40)
        csoft = (1 / 2) * pow(s1(n, exp_type), 2) + s2(n, exp_type) + g01bar * s1(n, exp_type) + g02bar

        return csoft

# N3LO large N
def n3lo_large_n(n, exp_type):

        # Soft approximations (2.40)
        # csoft = (1 / 6) * pow(s1(n, exp_type), 3) + s1(n, exp_type) * s2(n, exp_type + s3(n, exp_type)) + g01bar * ((1 / 2) * pow(s1(n, exp_type), 2) + s2(n)) + g02bar * s1(n) + g30bar
        csoft = (1 / 6) * pow(s1(n, exp_type), 3) + s1(n, exp_type) * s2(n, exp_type) + g01bar * ((1 / 2) * pow(s1(n, exp_type), 2) + s2(n, exp_type)) + g02bar * s1(n, exp_type) + g03bar

        return csoft