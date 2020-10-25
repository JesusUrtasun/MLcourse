# Gaussian Processes
# Check data

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
beta0 = (33 - 2 * nf) / 3

# Small N asymptotics
def sigma_small_n(x):

        # Small N contributions
        c10 = 2.28
        c20 = 4.12
        c30 = 8.64
        c11 = 5.66
        c21 = 10.54
        
        # Expressions (2.52)
        e00 = (-11 * ca + 2 * nf * (2 * (cf / ca) + 1)) / (12 * np.pi)
        e01 = ca / np.pi
        e11 = ((13 * cf) / (18 * np.pi) + (23 * ca) / (36 * np.pi)) * nf
        e12 = 0

        # Expressions (2.51)
        gamma0 = e01 / (x - 1) + e00

        # Small N asymptotics
        smallN = 2 * c10 * gamma0

        return  smallN

# Large N asymptotics
def sigma_large_n(x):
        
        # Expressions (...)
        g01 = (2 / np.pi) * (3 * pow(np.euler_gamma, 2) + pow(np.pi, 2) + 11/4)
        g02 = 40.10
        g12 = 2 * ca / np.pi
        g21 = 2 * ca * np.euler_gamma * np.pi

        # Expressions (2.27) and (2.28)
        Dlog = 0.5 * (pow(np.log(x), 2) + 2 * np.euler_gamma * np.log(x))
        Dhat = 0.5 * (pow(sc.polygamma(0, x), 2) + 2 * np.euler_gamma * sc.polygamma(0, x) + sc.zeta(2) + pow(sc.polygamma(0, x), 2))

        # Large N asymptotics (2.19)
        largeN = 4 * (ca / np.pi) * Dlog + g01 
        # largeN = 2 * g12 * Dhat + g01 + (-6 * (pow(np.euler_gamma, 2) / np.pi) - np.pi)

        return largeN

print("Gaussian Processes - check data")

# Loading data files
print("\nLoading datafiles")

# Known points in our n range
n_small = np.linspace(1.05, 1.5, 20)
n_large = np.linspace(2.5, 12, 20)
n = np.concatenate((n_small, n_large), axis = None).reshape(-1, 1)
# From Mathematica - I
n_small2 = np.loadtxt("data/mathematica/N_small_long.txt")
n_large2 = np.loadtxt("data/mathematica/N_large_long.txt")
n2 = np.concatenate((n_small2, n_large2), axis = None).reshape(-1, 1)
# From Mathematica - II
n_small3 = np.loadtxt("data/mathematica/N_small2_long.txt")
n_large3 = np.loadtxt("data/mathematica/N_large2_long.txt")
n3 = np.concatenate((n_small3, n_large3), axis = None).reshape(-1, 1)

# Known values of the cross section - NLO inclusive
sigma_n_small = sigma_small_n(n_small)
sigma_n_large = sigma_large_n(n_large)
sigma = np.concatenate((sigma_n_small, sigma_n_large), axis = None)
# From Mathematica - I
sigma_n_small2 = np.loadtxt("data/mathematica/sigma_N_small_long.txt")
sigma_n_large2 = np.loadtxt("data/mathematica/sigma_N_large_long.txt")
sigma2 = np.concatenate((sigma_n_small2, sigma_n_large2), axis = None)
# From Mathematica - II
sigma_n_small3 = np.loadtxt("data/mathematica/sigma_N_small2_long.txt")
sigma_n_large3 = np.loadtxt("data/mathematica/sigma_N_large2_long.txt")
sigma3 = np.concatenate((sigma_n_small3, sigma_n_large3), axis = None)

# Exact coefficient function - mesh the input space for the prediction and its MSE
n_exact = np.loadtxt("data/mathematica/N_exact_long.txt").reshape(-1, 1)
sigma_exact = np.loadtxt("data/mathematica/sigma_N_exact_long.txt")
# Exact coefficient function from ggHiggs https://www.ge.infn.it/~bonvini/higgs/
ggHiggs_data = np.loadtxt("data/ggHiggs_test_long.txt").T
n_ggHiggs = ggHiggs_data[0]
nlo_exact = ggHiggs_data[1]
nnlo_exact = ggHiggs_data[2]
n3lo_exact = ggHiggs_data[3]

print("n: {}".format(n.shape))
print("sigma: {}".format(sigma.shape))
print("n mesh: {}".format(n_exact.shape))
print("sigma exact: {}".format(sigma_exact.shape))

# Check different implementations of the inclusive cross section
def plot_comparison():
        plt.figure()
        plt.plot(n_ggHiggs, nlo_exact, "r", label = "NLO - ggHiggs")
        plt.plot(n_exact, sigma_exact, "b", label = "NLO - aprox Rgg")
        plt.ylim((0, 13))
        plt.ylim((0, 35))
        plt.xlabel(r"$N$")
        plt.ylabel(r"$C^{(1)}(N)$")
        plt.legend(loc = "upper right")
        plt.show()
plot_comparison()

# Check different asymptotics at small N
def plot_small_n():
        plt.figure()
        plt.plot(n_ggHiggs, nlo_exact, "r", label = "NLO C(N) ggHiggs")
        plt.plot(n_small, sigma_n_small, "r.", label = "Asymptotics - small N")
        plt.plot(n_small2, sigma_n_small2, "g.", label = "Asymptotics - (2.53)")
        plt.plot(n_small3, sigma_n_small3, "b.", label = "Asymptotics - (2.54)")
        plt.xlim((1, 2))
        plt.ylim((0, 200))
        plt.xlabel(r"$N$")
        plt.ylabel(r"$C(N)$")
        plt.legend(loc = "upper right")
        plt.show()
plot_small_n()

# Check different asymptotics at large N
def plot_large_n():
        plt.figure()
        plt.plot(n_ggHiggs, nlo_exact, "r", label = "NLO C(N) ggHiggs")
        plt.plot(n_large, sigma_n_large, "r.", label = "Asymptotics - Python")
        plt.plot(n_large2, sigma_n_large2, "g.", label = "Asymptotics - (2.19)")
        plt.plot(n_large3, sigma_n_large3, "b.", label = "Asymptotics - (2.28)")
        plt.xlim((1, 12))
        plt.ylim((0, 35))
        plt.xlabel(r"$N$")
        plt.ylabel(r"$C(N)$")
        plt.legend(loc = "upper right")
        plt.show()
plot_large_n()