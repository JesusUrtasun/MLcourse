# Gaussian Processes
# Higgs pt distribution

"""Doc string - None"""
# print(__doc__)

import sys
import numpy as np
from scipy import special as sc
from scipy.integrate import quad
from matplotlib import pyplot as plt
from matplotlib import gridspec
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C
# Play with pdfs
import lhapdf  
import pdb


# Setting random seed
np.random.seed(1)

# True function to approximate
def Sigma_Exact(x):
        # Parameters of the Hypergeometric
        a, b, c = 0.5, 2.0, 2.5
        y = pow((np.sqrt(x + 1) - np.sqrt(x)), 4)
        return sc.beta(2, 0.5) * sc.hyp2f1(a, b, c, y) / x

def Sigma_Low_pt(x):
        # Asymptotic behavior of the FO cross section at low pt
        return  -0.5 * np.log(x) / x - (sc.polygamma(0, 2) + np.euler_gamma) / x

def Sigma_Large_pt(x):
        # Asymptotic behavior of the FO cross section at low pt
        return sc.beta(2, 0.5) / x

print("Gaussian Processes\nHiggs pt distribution")

# Load PDF from LHAPDF
print("\nLoading PDF set")
# 0 - central memeber (study documentation ...)
pdf = lhapdf.mkPDF("NNPDF31_nlo_as_0118", 0)
Q2 = 10
x_values = np.logspace(-3, 0)
pdf_values = []
for i in x_values:
        pdf_values.append(pdf.xfxQ2(i, Q2))

# Acces the pdfs in a varying range of x
g = []
u = []
d = []
dbar = []
ubar = []
for pdf_fix_x in pdf_values:
        g.append(pdf_fix_x[21])
        d.append(pdf_fix_x[1])
        dbar.append(pdf_fix_x[-1])
        u.append(pdf_fix_x[2])
        ubar.append(pdf_fix_x[-2])
# Build valence partons and turn to numpy array
g = np.array(g)
dval = np.array(d) - np.array(dbar)
uval = np.array(u) - np.array(ubar)

# Plot the function, the prediction and the 95% confidence interval based on the MSE
plt.figure()
plt.title("PDF Q2 = 10")
plt.plot(x_values, g/10, "r:", markersize = 10, label = "g")
plt.plot(x_values, dval, "b:", markersize = 10, label = "dval")
plt.plot(x_values, uval, "g:", markersize = 10, label = "uval")
plt.xscale("log")
plt.xlabel(r"$x$")
plt.yscale("log")
plt.ylabel(r"$xf(x)$")
plt.legend(loc = "upper right")
plt.show()

# Build hadronic observables
# 1. build luminosity
# 2. mellin transform
# 3. Convolution. Just element wise prod in mellin

# partonic cross sectio FO
pt = np.logspace(-5, 4, 1000).reshape(-1, 1)
sigma_part = Sigma_Exact(pt).reshape(1000)  

# Add normalization options: sigma log or normalized to the integral (default)
if len(sys.argv) > 1 and sys.argv[1] == "log":
        is_log_norm = True
else:
        is_log_norm = False

# Preparing data files
print("\nLoading datafiles")

# Known points in our pt range
pt_low = np.logspace(-5, -3, 10)
pt_large = np.logspace(2.5, 4, 10)
pt = np.concatenate((pt_low, pt_large), axis = None).reshape(-1, 1)
# Known values of the cross section
sigma_low_pt = Sigma_Low_pt(pt_low)
sigma_large_pt = Sigma_Large_pt(pt_large)
sigma = np.concatenate((sigma_low_pt, sigma_large_pt), axis = None)
dsigma = 0.0001 * sigma
print("pt: {}".format(pt.shape))
print("\nsigma: {}".format(sigma.shape))

# Mesh the input space for evaluations of the real functions, the prediction and its MSE
pt_exact = np.logspace(-5, 4, 1000).reshape(-1, 1)
sigma_exact = Sigma_Exact(pt_exact).reshape(1000)  
print("\npt mesh: {}".format(pt_exact.shape))
print("sigma exact: {}".format(sigma_exact.shape))

# Normalizing data files
print("\nNormalizing datafiles")
# Normalize by computing log of sigma
if is_log_norm:
        sigma = np.log(sigma)
        sigma_exact = np.log(sigma_exact)
# Normalize to the integral of the cross section
else:
        norm = np.sum(sigma_exact)
        sigma = sigma / norm
        dsigma = dsigma / norm
        sigma_exact = sigma_exact / norm
        print("pt: {}".format(pt.shape))
        print("sigma: {}".format(sigma.shape))
        print("pt mesh: {}".format(pt_exact.shape))
        print("sigma exact: {}".format(sigma_exact.shape))

# Instantiate a Gaussian Process model
print("\nInstantiate a Gaussian Process model")
kernel = C(1.0, (1e-3, 1e3)) * RBF(10, (1e-2, 1e2))
# kernel = Matern(10, (1e-2, 1e2))
gp = GaussianProcessRegressor(kernel = kernel, alpha = dsigma, n_restarts_optimizer = 30)

# Fit to data using the Maximum Likelihood Estimation of the parameters
print("Fitting Gaussian Process")
gp.fit(np.log(pt), sigma)

# Make the predictin on the meshed x-axis (ask for MSE as well)
print("\nMaking predictions")
sigma_pred, std_dev = gp.predict(np.log(pt_exact), return_std = True)
ratio_GP = abs(sigma_pred / sigma_exact)
print("sigma pred: {}".format(sigma_pred.shape))
print("sigma exact: {}".format(sigma_exact.shape))
print("ratio: {}".format(ratio_GP.shape))

# Plot the function, the prediction and the 95% confidence interval based on the MSE
plt.figure()
gs = gridspec.GridSpec(3, 1, height_ratios = [3.5, 1, 1])
plt.subplot(gs[0])
plt.title("Higgs pt distribution - LO in pt")
plt.plot(pt_exact, sigma_exact, "r:", markersize = 10, label = "FO theory")
# plt.plot(pt, sigma, "r.", markersize = 10, label = "Asymptotics")
plt.errorbar(pt, sigma, dsigma, fmt = "r.", markersize = 10, label = "Asymptotics")
plt.plot(pt_exact, sigma_pred, "b-", label = "GP Prediction")
plt.fill(np.concatenate([pt_exact, pt_exact[::-1]]),
        np.concatenate([sigma_pred - 1.9600 * std_dev, (sigma_pred + 1.9600 * std_dev)[::-1]]),
        alpha = 0.5, fc = "b", ec = "None", label = "95% confidence interval")
plt.xscale("log")
plt.xlabel(r"$p_{\perp}^{2} / m_{h}^{2}$")
if is_log_norm:
        # Log sigma values in the y axis
        plt.ylabel(r"$\frac{d\sigma}{dp_{\perp}}(N, p_{\perp})$")
else:
        # Sigma values normalized to the cross section
        plt.ylabel(r"$\frac{d\sigma}{dp_{\perp}}(N, p_{\perp})$")
plt.legend(loc = "upper right")

# Plot ratio of GP vs the exact theory
plt.subplot(gs[1])
plt.plot(pt_exact, ratio_GP, "g", markersize = 10)
plt.hlines(1, pt_exact[0], pt_exact[-1], color = "r")
plt.fill(np.concatenate([pt_exact, pt_exact[::-1]]),
        np.concatenate([ratio_GP - (1 / sigma_exact * 1.9600 * std_dev), (ratio_GP + (1 / sigma_pred) * 1.9600 * std_dev)[::-1]]),
        alpha = 0.5, fc = "b", ec = "None", label = "95% confidence interval")
plt.ylim((-20, 20))
plt.xscale("log")
plt.xlabel(r"$p_{\perp}^{2} / m_{h}^{2}$")
plt.ylabel(r"$\sigma_{GP} / \sigma_{theory}$")

# Plot ratio of GP vs the exact theory
plt.subplot(gs[2])
plt.plot(pt_exact, ratio_GP, "g", markersize = 10)
plt.hlines(1, pt_exact[0], pt_exact[-1], color = "r")
plt.fill(np.concatenate([pt_exact, pt_exact[::-1]]),
        np.concatenate([ratio_GP - (1 / sigma_exact * 1.9600 * std_dev), (ratio_GP + (1 / sigma_pred) * 1.9600 * std_dev)[::-1]]),
        alpha = 0.5, fc = "b", ec = "None", label = "95% confidence interval")
plt.xscale("log")
plt.xlabel(r"$p_{\perp}^{2} / m_{h}^{2}$")
plt.ylim((0, 2))
plt.ylabel(r"$\sigma_{GP} / \sigma_{theory}$")
plt.show()

# Checking integral
# print("Checking integral of GP prediction")
# def Sigma_Pred_Fun(pt):
#         return gp.predict(np.log(np.array([pt]).reshape(-1, 1)))[0]

# import ipdb
# ipdb.set_trace()

# print("Int sigma theory = {:.5f}".format( quad(Sigma_Exact, 1e-3, 10**(1.5))[0]) )
# print("Int sigma GP = {:.5f}".format( quad(Sigma_Pred_Fun , 1e-3, 10**(1.5))[0]) )