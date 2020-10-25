# Gaussian Processes
# Check N3LO data

import pdb
import numpy as np
from scipy import special as sc
from matplotlib import pyplot as plt
from matplotlib import gridspec
# For computing asymptotics
import contributions_small_n
import contributions_large_n

# Setting parameters
np.random.seed(1)
ca = 3
cf = 4 / 3
nf = 5
alphas = 0.1126
beta0 = (33 - 2 * nf) / (12 * np.pi)

##### Small N contributions #####

# Appendix (A.20)
c10 = 2.28
c20 = 4.12
c30 = 8.64
c11 = 5.66
c21 = 10.54

# Expressions (2.52)
e00 = (- 11 * ca + 2 * nf * (2 * (cf / ca) - 1)) / (12 * np.pi)
e01 = ca / np.pi
e11 = ((13 * cf) / (18 * pow(np.pi, 2)) - (23 * ca) / (36 * pow(np.pi, 2))) * nf
e12 = 0
e22 = ((pow(ca, 3) * sc.zeta(3)) / (2 * pow(np.pi, 3))) + ((11 * pow(ca, 3) * sc.zeta(2)) / (12 * pow(np.pi, 3))) - ((395 * pow(ca, 3)) / (108 / pow(np.pi, 3))) + (((pow(ca, 2) * sc.zeta(2)) / (6 * pow(np.pi, 3))) - ((71 * pow(ca, 2)) / (108 * pow(np.pi, 3))) - ((ca * cf * sc.zeta(2)) / (3 * pow(np.pi, 3))) + ((71 * ca * cf) / (54 * pow(np.pi, 3)))) * nf
e23 = 0

# Expressions in (2.51)
def gamma0(n):
        return  e01 / (n - 1) + e00

def gamma1(n):
        # return e12 / pow((n - 1), 2) + e11 / (n - 1)
        return e11 / (n - 1)

def gamma2(n):
        # return e23 / pow((n - 1), 3) + e22 / pow((n - 1), 2)
        return e22 / pow((n - 1), 2)

##### Large N contributions #####

# Expressions (2.17) and Table 2
# g01 = (2 / np.pi) * (3 * pow(np.euler_gamma, 2) + pow(np.pi, 2) + 11/4)
g01 = 8.7153
g02 = 40.10
g12 = 2 * ca / np.pi
g21 = 4 * ca * np.euler_gamma / np.pi
# Used for nnlo only
r1 = 3.7779
r2 = 29.18
r3 = 114.7
g01bar = 4.9374
g02bar = 10.92
g03bar = 2.0 # updated version https://arxiv.org/pdf/1404.3204.pdf
# g03bar = -r3

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
        + (-(58 / 729) + (10 / 27) * sc.zeta(2) + (5 / 27) * sc.zeta(3)) * ca * pow(nf, 2))
# b30 = 0
b31 = (1 / pow(np.pi, 3)) * (((15503 / 648) - (188 / 9) * sc.zeta(2) - 11 * sc.zeta(3) + (11 / 5) * pow(sc.zeta(2), 2)) * pow(ca, 3)
        + (-(2051 / 324) + 6 * sc.zeta(2)) * pow(ca, 2) * nf + (-(55 / 24) + 2 * sc.zeta(3)) * ca * cf * nf
        + ((25 / 81) - (4 / 9) * sc.zeta(2)) * ca * pow(nf, 2))
# b31 = 0
b32 = (1 / pow(np.pi, 3)) * (-((445 / 27) + (11 / 3) * sc.zeta(2)) * pow(ca, 3) +((289 / 54) - (2 / 3) * sc.zeta(2)) * pow(ca, 2) * nf
        + (1 / 2) * ca * cf * nf - (10 / 27) * ca * pow(nf, 2))
# b32 = 0
b33 = (1 / pow(np.pi, 3)) * ((121 / 27) * pow(ca, 3) - (44 / 27) * pow(ca, 2) * nf + (4 / 27) * ca * pow(nf, 2))
b33 = 0

# Expressions D (A.6a)
def D0(n):
        return - np.euler_gamma - sc.polygamma(0, n)

def D1(n):
        return (1 / 12) * (6 * pow(np.euler_gamma, 2) + pow(np.pi, 2) + 6 * sc.polygamma(0, n) * (2 * np.euler_gamma + sc.polygamma(0, n)) - 6 * sc.polygamma(1, n))

def D2(n):
        return (1 / 6) * (- 2 * pow(np.euler_gamma, 3) - np.euler_gamma * pow(np.pi, 2) - sc.polygamma(0, n) * (6 * pow(np.euler_gamma, 2) + pow(np.pi, 2)
        + 2 * sc.polygamma(0, n) * (3 * np.euler_gamma + sc.polygamma(0, n))) + 6 * (np.euler_gamma + sc.polygamma(0, n)) * sc.polygamma(1, n) - 2 * sc.polygamma(2, n) - 4 * sc.zeta(3))

def D3(n):
        return (1 / 80) * (20 * pow(np.euler_gamma, 4) + 20 * pow(np.euler_gamma * np.pi, 2) + 3 * pow(np.pi, 4) - 20 * sc.polygamma(3, n)
        + 20 * (4 * np.euler_gamma * pow(sc.polygamma(0, n), 3) + pow(sc.polygamma(0, n), 4) + pow(sc.polygamma(0, n), 2) * (6 * pow(np.euler_gamma, 2) + pow(np.pi, 2) - 6 * sc.polygamma(1, n))
        - (6 * pow(np.euler_gamma, 2) + pow(np.pi, 2)) * sc.polygamma(1, n) + 3 * pow(sc.polygamma(1, n), 2) + 4 * np.euler_gamma * (sc.polygamma(2, n) + 2 * sc.zeta(3)) 
        + 2 * sc.polygamma(0, n) * (2 * pow(np.euler_gamma, 3) + np.euler_gamma * pow(np.pi, 2) - 6 * np.euler_gamma * sc.polygamma(1, n) + 2 * sc.polygamma(2, n) + 4 * sc.zeta(3))))

# Expressions Dlog (A.6b)
def Dlog0(n):
        return - np.log(n)

def Dlog1(n):
        return (1 / 2) * np.log(n) * (2 * np.euler_gamma + np.log(n))

def Dlog2(n):
        return - (1 / 6) * np.log(n) * (6 * pow(np.euler_gamma, 2) + pow(np.pi, 2) + 2 * np.log(n) * (3 * np.euler_gamma + np.log(n)))

def Dlog3(n):
        return (1 / 4) * np.log(n) * ((2 * np.euler_gamma + np.log(n)) * (2 * pow(np.euler_gamma, 2) + pow(np.pi, 2) + 2 * np.euler_gamma * np.log(n) + pow(np.log(n), 2)) + 8 * sc.zeta(3))

def Dlog_general(k, n):
        Dlog_conf = [Dlog0, Dlog1, Dlog2, Dlog3]
        current_conf = Dlog_conf[k]
        return current_conf(n)

# Expressions Dhat (A.6c)
def Dhat0(n):
        return - np.euler_gamma - sc.polygamma(0, n)

def Dhat1(n):
        return (1 / 12) * (pow(np.pi, 2) + 6 * pow(np.euler_gamma + sc.polygamma(0, n), 2))

def Dhat2(n):
        return (1 / 12) * (- sc.polygamma(2, n) - 2 * (pow(np.pi, 2) * (np.euler_gamma + sc.polygamma(0, n)) + 2 * pow(np.euler_gamma + sc.polygamma(0, n), 3) + 4 * sc.zeta(3)))

def Dhat3(n):
        return (1 / 80) * (20 * pow(np.euler_gamma, 4) + 20 * pow(np.euler_gamma * np.pi, 2) + 3 * pow(np.pi, 4)
        + 20 * ((6 * pow(np.euler_gamma, 2) + pow(np.pi, 2)) * pow(sc.polygamma(0, n), 2) + 4 * np.euler_gamma * pow(sc.polygamma(0, n), 3) + pow(sc.polygamma(0, n), 4)
        + np.euler_gamma * (sc.polygamma(2, n) + 8 * sc.zeta(3)) + sc.polygamma(0, n) * (4 * pow(np.euler_gamma, 3) + 2 * np.euler_gamma * pow(np.pi, 2) + sc.polygamma(2, n) + 8 * sc.zeta(3))))

# Prpare soft2 ()
def get_D_soft2(d_hat):

        return lambda n, d_hat = d_hat: 2 * d_hat(n) - 3 * d_hat(n + 1) + 2 * d_hat(n + 2)

def D_general(k, n, exp_type):
        
        # first key to a list of functions
        D_conf = {"soft1": [Dhat0, Dhat1, Dhat2, Dhat3],
                "soft2": [get_D_soft2(Dhat0), get_D_soft2(Dhat1), get_D_soft2(Dhat2), get_D_soft2(Dhat3)],
                "soft0": [D0, D1, D2, D3]}
        
        # current conf is just a list of conf names
        current_conf = D_conf[exp_type]
        # Shitf argument for 2.36a
        if exp_type == "soft1":
                n = n + 1
        return current_conf[k](n)

def check_A8():
        print("Check Dhat - D in (A.8)")
        for k in [0, 1, 2, 3]:
                for n in [2, 5, 20, 50, 100, 1000]:
                        diff = D_general(k, n, "soft1") - D_general(k, n, "soft0")
                        print("A8: Dhat{}({}) - D = {}".format(k, n, diff))

# check_A8()

def check_A9():
        print("Check Dhat - Dlog in (A.9)")
        rhside = [-np.euler_gamma, (1 / 12) * (6 * pow(np.euler_gamma, 2) + pow(np.pi, 2)),
                (1 / 6) * (-2 * pow(np.euler_gamma, 3) - np.euler_gamma * pow(np.pi, 2) - 4 * sc.zeta(3)),
                (1 / 4) * (pow(np.euler_gamma, 4) + pow(np.euler_gamma * np.pi, 2) + (3 / 20) * pow(np.pi, 4) + 8 * np.euler_gamma * sc.zeta(3))]
        for k in [0, 1, 2, 3]:
                for n in [2, 5, 20, 50, 100, 1000]:
                        lhs = D_general(k, n, "soft1") - Dlog_general(k, n)
                        rhs = rhside[k]
                        diff = lhs - rhs
                        print("A9: Dhat{}({}) - Dlog = {}".format(k, n, diff))

# check_A9()

# S(N) expressions in (2.36)
def s1(n, exp_type):
        return b10 * D_general(0, n, exp_type) + b11 * D_general(1, n, exp_type) # (2.36a)

def s2(n, exp_type):
        return b20 * D_general(0, n, exp_type) + b21 * D_general(1, n, exp_type) + b22 * D_general(2, n, exp_type) # (2.36b)

def s3(n, exp_type):
        return b30 * D_general(0, n, exp_type) + b31 * D_general(1, n, exp_type) + b32 * D_general(2, n, exp_type) + b33 * D_general(3, n, exp_type) # (2.36a)

##### Asymptotics at NLO #####

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

# NLO large N
def nlo_large_n(n, exp_type):

        # Soft approximations (2.41a)
        csoft = s1(n, exp_type) + g01bar

        return csoft

##### Asymptotics at NNLO #####

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

# NNLO large N
def nnlo_large_n(n, exp_type):

        # Soft approximations (2.41b)
        csoft = (1 / 2) * pow(s1(n, exp_type), 2) + s2(n, exp_type) + g01bar * s1(n, exp_type) + g02bar

        return csoft

##### Asymptotics at N3LO #####

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

# N3LO large N
def n3lo_large_n(n, exp_type):

        # Soft approximations (2.41c)
        csoft = (1 / 6) * pow(s1(n, exp_type), 3) + s1(n, exp_type) * s2(n, exp_type) + s3(n, exp_type) + g01bar * ((1 / 2) * pow(s1(n, exp_type), 2) + s2(n, exp_type)) + g02bar * s1(n, exp_type) + g03bar

        return csoft

print("Gaussian Processes - check data")

pdb.set_trace()

# Known points in our n range
n_small = np.linspace(0.02, 1.5, 20)
n_large = np.linspace(2.5, 12, 30)
n = np.concatenate((n_small, n_large), axis = None).reshape(-1, 1)
# Exact coefficient function from ggHiggs https://www.ge.infn.it/~bonvini/higgs/
ggHiggs_data = np.loadtxt("data/ggHiggs_mesh_correct.txt").T
n_ggHiggs = ggHiggs_data[0]
nlo_exact = ggHiggs_data[1]
nnlo_exact = ggHiggs_data[2]
n3lo_exact = ggHiggs_data[3]

# Check different orders of ggHiggs
def plot_ggHiggs():
        plt.figure()
        plt.plot(n_ggHiggs, nlo_exact, "r", label = "NLO C(N) ggHiggs")
        plt.plot(n_ggHiggs, nnlo_exact, "b", label = "NNLO C(N) ggHiggs")
        plt.plot(n_ggHiggs, n3lo_exact, "g", label = "N3LO C(N) ggHiggs")
        plt.xlim((1, 12))
        plt.ylim((0, 350))
        plt.xlabel(r"$N$")
        plt.ylabel(r"$C(N)$")
        plt.legend(loc = "upper right")
        plt.show()
# plot_ggHiggs()

# Check figure 1 NLO (left)
def plot_figure1_nlo():
        
        # Generate n range
        n = np.linspace(0, 13, 1000)
        
        # Plot exact and large N asymptotics
        plt.figure()
        plt.plot(n_ggHiggs, nlo_exact, "r", label = "exact")
        plt.plot(n, nlo_large_n(n, "soft0"), "r:", label = "soft0")
        plt.plot(n, nlo_large_n(n, "soft1"), "b:", label = "soft1")
        plt.plot(n, nlo_large_n(n, "soft2"), "g:", label = "soft2")
        plt.xlim((0, 13))
        plt.ylim((0, 35))
        plt.xlabel(r"$N$")
        plt.ylabel(r"$C^{(1)}(N)$")
        plt.legend(loc = "lower right")
        plt.show()

plot_figure1_nlo()

# Check figure 1 NNLO (right)
def plot_figure1_nnlo():

        # Generate n range
        n = np.linspace(0, 13, 1000)
        
        # Plot exact and large N asymptotics
        plt.figure()
        plt.plot(n_ggHiggs, nnlo_exact, "r", label = "exact")
        plt.plot(n, nnlo_large_n(n, "soft0"), "r:", label = "soft0")
        plt.plot(n, nnlo_large_n(n, "soft1"), "b:", label = "soft1")
        plt.plot(n, nnlo_large_n(n, "soft2"), "g:", label = "soft2")
        plt.xlim((0, 13))
        plt.ylim((0, 350))
        plt.xlabel(r"$N$")
        plt.ylabel(r"$C^{(2)}(N)$")
        plt.legend(loc = "lower right")
        plt.show()

plot_figure1_nnlo()

# Check large N N3LO
def plot_figure1_n3lo():

        # Generate n range
        n = np.linspace(0, 13, 1000)
        
        # Plot exact and large N asymptotics
        plt.figure()
        plt.plot(n_ggHiggs, n3lo_exact, "r", label = "exact")
        plt.plot(n, n3lo_large_n(n, "soft0"), "r:", label = "soft0")
        plt.plot(n, n3lo_large_n(n, "soft1"), "b:", label = "soft1")
        plt.plot(n, n3lo_large_n(n, "soft2"), "g:", label = "soft2")
        plt.xlim((0, 13))
        plt.ylim((0, 5000))
        plt.xlabel(r"$N$")
        plt.ylabel(r"$C^{(3)}(N)$")
        plt.legend(loc = "lower right")
        plt.show()

plot_figure1_n3lo()

# Check figure 2 NLO (left)
def plot_figure2_nlo():

        # Generate n range
        n = np.linspace(0, 13, 1000)

        # Plot exact an small N asymptotics
        plt.figure()
        plt.plot(n_ggHiggs, nlo_exact, "r", label = "exact")
        plt.plot(n, nlo_small_n_abf(n), "r:", label = "ABF (2.53)")
        plt.plot(n, nlo_small_n_abf_subl(n), "b:", label = "ABF-subl (2.54)")
        plt.hlines(0, n_ggHiggs[0], n_ggHiggs[-1])
        plt.xlim((1, 4))
        plt.ylim((-5, 30))
        plt.xlabel(r"$N$")
        plt.ylabel(r"$C^{(1)}(N)$")
        plt.legend(loc = "upper right")
        plt.show()

        # Plot ratio
        # plt.figure()
        # plt.title("Higgs NLO - figure 1")
        # plt.hlines(1, n_ggHiggs[0], n_ggHiggs[-1])
        # plt.plot(n_ggHiggs, nlo_small_n_abf(n_ggHiggs) / nlo_exact, "r:", label = "ABF (2.53)")
        # plt.plot(n_ggHiggs, nlo_small_n_abf_subl(n_ggHiggs) / nlo_exact, "b:", label = "ABF-subl (2.54)")
        # plt.xlim((1, 1.1))
        # plt.xlabel(r"$N$")
        # plt.ylabel(r"$C_{small N} / C$")
        # plt.legend(loc = "upper right")
        # plt.show()

plot_figure2_nlo()

# Check figure 2 NNLO (right)
def plot_figure2_nnlo():

        # Generate n range
        n = np.linspace(0, 13, 1000)

        # Plot exact an small N asymptotics        
        plt.figure()
        plt.plot(n_ggHiggs, nnlo_exact, "r", label = "exact")
        plt.plot(n_ggHiggs, nnlo_small_n_abf(n_ggHiggs), "r:", label = "ABF (2.53)")
        plt.plot(n_ggHiggs, nnlo_small_n_abf_subl(n_ggHiggs), "b:", label = "ABF-subl (2.54)")
        plt.hlines(0, n_ggHiggs[0], n_ggHiggs[-1])
        plt.xlim((1, 4))
        plt.ylim((-50, 300))
        plt.xlabel(r"$N$")
        plt.ylabel(r"$C^{(2)}(N)$")
        plt.legend(loc = "upper right")
        plt.show()

plot_figure2_nnlo()

# Check small N N3LO
def plot_figure2_n3lo():

        # Generate n range
        n = np.linspace(0, 13, 1000)

        # Plot exact an small N asymptotics        
        plt.figure()
        plt.plot(n_ggHiggs, n3lo_exact, "r", label = "exact")
        plt.plot(n_ggHiggs, n3lo_small_n_abf(n_ggHiggs), "r:", label = "ABF (2.53)")
        plt.plot(n_ggHiggs, n3lo_small_n_abf_subl(n_ggHiggs), "b:", label = "ABF-subl (2.54)")
        plt.hlines(0, n_ggHiggs[0], n_ggHiggs[-1])
        plt.xlim((1, 4))
        plt.ylim((-500, 3000))
        plt.xlabel(r"$N$")
        plt.ylabel(r"$C^{(3)}(N)$")
        plt.legend(loc = "upper right")
        plt.show()

plot_figure2_n3lo()

# Check figure 5
def plot_figure5():

        # Generate n range
        n = np.linspace(0, 13, 1000)

        # Plot exact an small N asymptotics        
        plt.figure()
        plt.plot(n_ggHiggs, [1 for i in range(len(n_ggHiggs))], "r", label = "LO")
        plt.plot(n_ggHiggs, 1 + alphas * nlo_exact, "b", label = "NLO")
        plt.plot(n_ggHiggs, 1 + alphas * nlo_exact + pow(alphas, 2) * nnlo_exact, "g", label = "NNLO")
        plt.plot(n_ggHiggs, 1 + alphas * nlo_exact + pow(alphas, 2) * nnlo_exact + pow(alphas, 3) * n3lo_exact, "y", label = "N3LO")
        plt.xlim((1, 5))
        plt.ylim((0, 7))
        plt.xlabel(r"$N$")
        plt.ylabel(r"$C^{(3)}(N)$")
        plt.legend(loc = "upper right")
        plt.show()

plot_figure5()