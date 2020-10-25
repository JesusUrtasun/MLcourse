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

##### Small N contributions #####

# Appendix (A.20)
c10 = 2.28
c20 = 4.12
c30 = 8.64
c11 = 5.66
c21 = 10.54

# Expressions (2.52)
e00 = (-11 * ca + 2 * nf * (2 * (cf / ca) - 1)) / (12 * np.pi)
e01 = ca / np.pi
e11 = ((13 * cf) / (18 * pow(np.pi, 2)) - (23 * ca) / (36 * pow(np.pi, 2))) * nf
e12 = 0
e22 = 1
e23 = np.nan # Unknown!!!

# Expressions in (2.51)
def gamma0(n):
        return  e01 / (n - 1) + e00

def gamma1(n):
        return e11 / (n - 1)
        # return e12 / pow((n - 1), 2) + e11 / (n - 1)

def gamma2(n):
        return e23 / pow((n - 1), 3) + e22 / pow((n - 1), 2)

##### Asymptotic expresssions #####

# NLO small N
def nlo_small_n_abf(n):

        # Small N asymptotics (2.53)
        cabf = 2 * c10 * gamma0(n)

        return  cabf

# NLO small N - subleading
def nlo_small_n_abf_subl(n):

        # Small N asymptotics (2.54)
        cabf_subl = nlo_small_n_abf(n) - 2 * nlo_small_n_abf(n + 1) + nlo_small_n_abf(n + 2)

        return  cabf_subl

# NNLO small N
def nnlo_small_n_abf(n):

        # Small N asymptotics (2.53)
        cabf = (2 * c20 + c11) * pow(gamma0(n), 2) - 2 * c20 * beta0 * gamma0(n) + 2 * c10 * gamma1(n)

        return  cabf

# NNLO small N - subleading
def nnlo_small_n_abf_subl(n):

        # Small N asymptotics (2.54)
        cabf_subl = nnlo_small_n_abf(n) - 2 * nnlo_small_n_abf(n + 1) + nnlo_small_n_abf(n + 2)

        return  cabf_subl

# N3LO small N
def n3lo_small_n_abf(n):

        # Small N asymptotics (2.53)
        cabf = (c30 + c21) * 2 * pow(gamma0(n), 3) - (3 * c30 + 2 * c21) * 2 * beta0 * pow(gamma0(n), 2) + 4 * c30 * pow(beta0, 2) * gamma0(n)
        + (2 * c20 + c11) * 2 * gamma0(n) * gamma1(n) - 4 * c20 * beta0 * gamma1(n) + 2 * c10 * gamma2(0)

        return  cabf

# N3LO small N - subleading
def n3lo_small_n_abf_subl(n):

        # Small N asymptotics (2.54)
        cabf_subl = n3lo_small_n_abf(n) - 2 * n3lo_small_n_abf(n + 1) + n3lo_small_n_abf(n + 2)

        return  cabf_subl