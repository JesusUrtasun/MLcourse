# Gaussian Processes
# Higgs inclusive

"""Doc string - None"""
# print(__doc__)

import pdb
import argparse
import numpy as np
from scipy import special as sc
from matplotlib import pyplot as plt
from matplotlib import gridspec
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel as C
from sklearn.gaussian_process.kernels import RBF, Matern, ExpSineSquared
from sklearn.gaussian_process.kernels import DotProduct, RationalQuadratic, PairwiseKernel
from sklearn.gaussian_process.kernels import Sum, Product, Exponentiation
from sklearn.model_selection import GridSearchCV, cross_validate

# Setting random seed
np.random.seed(1)

# Setting parameters
delta = 0.00001

# Setting normalization conditions
parser = argparse.ArgumentParser()
parser.add_argument("-o", "--order", help = "lo, nlo, nnlo, n3lo", required = True)
parser.add_argument("-r", "--restarts", help = "Number of restarts to ask the GP for", default = 100, type = int)
args = parser.parse_args()

print("Gaussian Processes - Higgs inclusive")

# Loading data files
print("\nLoading data files")
sigma_low_raw = np.loadtxt("data/sigma_low_tau.txt").T
sigma_large_raw = np.loadtxt("data/sigma_large_tau_test.txt").T
sigma_exact_raw = np.loadtxt("data/sigma_exact.txt").T
# Storing data
print("\nStoring data")
# Known points in our tau range
tau_low = sigma_low_raw[2]
tau_large = sigma_large_raw[2]
tau = np.concatenate((tau_low, tau_large), axis = None).reshape(-1, 1)
# Mesh the input space for evaluations of the real functions, the prediction and its MSE
tau_exact = sigma_exact_raw[2].reshape(-1, 1)
if args.order == "lo":
        # Known values of the cross section - LO
        sigma_tau_low = sigma_low_raw[3]
        sigma_tau_large = sigma_large_raw[3]
        sigma = np.concatenate((sigma_tau_low, sigma_tau_large), axis = None)
        dsigma = delta * sigma
        sigma_exact = sigma_exact_raw[3]
elif args.order == "nlo":
        # Known values of the cross section - NLO
        sigma_tau_low = sigma_low_raw[4]
        sigma_tau_large = sigma_large_raw[4]
        sigma = np.concatenate((sigma_tau_low, sigma_tau_large), axis = None)
        dsigma = delta * sigma
        # Mesh the input space for evaluations of the real functions, the prediction and its MSE
        sigma_exact = sigma_exact_raw[4]
elif args.order == "nnlo":
        # Known values of the cross section - NNLO
        sigma_tau_low = sigma_low_raw[5]
        sigma_tau_large = sigma_large_raw[5]
        sigma = np.concatenate((sigma_tau_low, sigma_tau_large), axis = None)
        dsigma = delta * sigma
        # Mesh the input space for evaluations of the real functions, the prediction and its MSE
        sigma_exact = sigma_exact_raw[5]
elif args.order == "n3lo":
        # Known values of the cross section - NNNLO
        sigma_tau_low = sigma_low_raw[6]
        sigma_tau_large = sigma_large_raw[6]
        sigma = np.concatenate((sigma_tau_low, sigma_tau_large), axis = None)
        dsigma = delta * sigma
        # Mesh the input space for evaluations of the real functions, the prediction and its MSE
        sigma_exact = sigma_exact_raw[6]
else:
        raise NotImplementedError(f"The order {args.order} is not yet implemented")

print("tau: {}".format(tau.shape))
print("sigma: {}".format(sigma.shape))
print("tau mesh: {}".format(tau_exact.shape))
print("sigma exact: {}".format(sigma_exact.shape))

# Normalizing data files
print("\nNormalizing data files")
# norm = np.sum(sigma_exact)
# sigma = sigma / norm
# dsigma = dsigma / norm
# sigma_exact = sigma_exact / norm
# print("tau: {}".format(tau.shape))
# print("sigma: {}".format(sigma.shape))
# print("tau mesh: {}".format(tau_exact.shape))
# print("sigma exact: {}".format(sigma_exact.shape))

print("\nPreparing for exhaustive GridSearch():")

# Set the parameters to hyperoptimize
kernel1 = DotProduct(sigma_0 = 1.0, sigma_0_bounds = (1e-5, 1e5))
kernel2 = ExpSineSquared(length_scale = 1.0, periodicity = 1.0, length_scale_bounds = (1e-5, 1e5))
kernel3 = Exponentiation(RBF(length_scale = 1.0, length_scale_bounds = (1e-2, 1e2)), 2)
kernel4 = Matern(length_scale = 1.0, length_scale_bounds = (1e-5, 1e5), nu = 1.5)
kernel5 = PairwiseKernel(gamma = 1.0, gamma_bounds = (1e-5, 1e5))
kernel6 = Product(RBF(1.0, (1e-5, 1e5)), Matern(1.0, (1e-5, 1e5), nu = 1.5))
kernel7 = RBF(length_scale = 1.0, length_scale_bounds = (1e-5, 1e5))
kernel8 = RationalQuadratic(length_scale = 1.0, alpha = 1.0, length_scale_bounds = (1e-5, 1e5), alpha_bounds = (1e-5, 1e5))
kernel9 = Sum(RBF(1.0, (1e-2, 1e2)), Matern(10, (1e-2, 1e2), nu = 1.5))
# List of hyperparameters given to the GridSearchCV()
tuned_parameters = [{"kernel": [kernel1, kernel2, kernel3, kernel4, kernel5, kernel6, kernel7, kernel8, kernel9]}]
kernel_names = ["DP", "ES", "Exp", "Mat", "PW", "Prod", "RBF", "RQ", "Sum"]

# Instance the GridSearchCV()
gp_hyper = GridSearchCV(GaussianProcessRegressor(), tuned_parameters)
gp_hyper.fit(tau, sigma)

# Compute means and standard deviations
print("\nChecking performance")
means = gp_hyper.cv_results_["mean_test_score"]
stds = gp_hyper.cv_results_["std_test_score"]
print("means: {}".format(means))
print("stds: {}".format(stds))

# Import cv_results into a Pandas Datarame
# (...)

# Check means and stds for each parameter configuration
print("\nSummary")
for mean, std, params in zip(means, stds, gp_hyper.cv_results_['params']):
        print("%0.3f (+/-%0.03f) for %r"
              % (mean, std * 2, params))

# Choose the best set of parameters for our model
print("\nBest parameters set found on development set:")
print(gp_hyper.best_params_)
# cv_results_["params"] contains a dictionary with the best parameters found
best_model_dict = gp_hyper.cv_results_['params'][gp_hyper.best_index_]
best_kernel = best_model_dict["kernel"]

# Define the kernel as product with the constant kernel
kernel = C(1.0, (1e-3, 1e3)) * best_kernel
best_model = GaussianProcessRegressor(kernel = best_kernel, alpha = 1e-4 * sigma, n_restarts_optimizer = args.restarts)

# Fit to data using the Maximum Likelihood Estimation of the parameters
print("\nFitting Gaussian Process: {}".format(best_kernel))
best_model.fit(tau, sigma)

# Make the predictin on the meshed x-axis (ask for MSE as well)
print("\nMaking predictions")
sigma_pred2, std_dev2 = best_model.predict(tau_exact, return_std = True)

# Computing chi2 estimator
num = pow(sigma_exact - sigma_pred2, 2)
# Correct negative variances
corrected_std = np.maximum(pow(std_dev2 , 2), 1e-10)
chi2 = np.sum(num / corrected_std) / len(sigma_pred2)

# Computing ratio GP vs exact
ratio_GP2 = abs(sigma_pred2 / sigma_exact)
print("sigma pred: {}".format(sigma_pred2.shape))
print("sigma exact: {}".format(sigma_exact.shape))
print("ratio: {}".format(ratio_GP2.shape))
print("chi2: {}".format(chi2))

# Plot the function, the prediction and the confidence intervals for the best model
print("\nBest model found for Kernel = {}".format(best_model_dict["kernel"]))
plt.figure()
gs = gridspec.GridSpec(2, 1, height_ratios = [4, 1])
plt.subplot(gs[0])
if args.order == "lo":
        plt.title("Higgs gg - LO inclusive - GridSearch")
elif args.order == "nlo":
        plt.title("Higgs gg - NLO inclusive - GridSearch")
elif args.order == "nnlo":
        plt.title("Higgs gg - NNLO inclusive - GridSearch")
elif args.order == "nnnlo":
        plt.title("Higgs gg - N3LO inclusive - GridSearch")
plt.plot(tau_exact, sigma_exact, "r:", markersize = 10, label = "FO theory")
plt.errorbar(tau, sigma, dsigma, fmt = "r.", markersize = 10, label = "Asymptotics")
plt.plot(tau_exact, sigma_pred2, "b-", label = "GP Prediction")
plt.fill(np.concatenate([tau_exact, tau_exact[::-1]]),
        np.concatenate([sigma_pred2 - 1.9600 * std_dev2,
                        (sigma_pred2 + 1.9600 * std_dev2)[::-1]]),
        alpha = 0.5, fc = "gray", ec = "None", label = "95% confidence interval")
plt.fill(np.concatenate([tau_exact, tau_exact[::-1]]),
        np.concatenate([sigma_pred2 - 1.0 * std_dev2,
                        (sigma_pred2 + 1.0 * std_dev2)[::-1]]),
        alpha = 0.5, fc = "b", ec = "None", label = "68% confidence interval")
plt.xscale("log")
plt.xlabel(r"$m_{H}^{2} / s$")
plt.ylabel(r"$\sigma(\tau) \quad (pb)$")
plt.legend(loc = "upper right")

# Plot ratio of GP vs the exact theory
plt.subplot(gs[1])
plt.plot(tau_exact, ratio_GP2, "g", markersize = 10)
plt.hlines(1, tau_exact[0], tau_exact[-1], color = "r")
plt.fill(np.concatenate([tau_exact, tau_exact[::-1]]),
        np.concatenate([ratio_GP2 - abs(1 / sigma_exact * 1.9600 * std_dev2),
                        (ratio_GP2 + abs(1 / sigma_exact * 1.9600 * std_dev2))[::-1]]),
        alpha = 0.5, fc = "b", ec = "None", label = "95% confidence interval")
plt.ylim((-20, 20))
plt.xscale("log")
plt.xlabel(r"$m_{H}^{2} / s$")
plt.ylabel(r"$\frac{\sigma_{GP}}{\sigma_{theory}}$")
plt.show()

print("End of exhaustive GridSearch")
pdb.set_trace()

print("\nLoop iteration over kernels\nChi2 as scores")

# Prepare list of kernels for the gaussian regression
# https://scikit-learn.org/stable/modules/classes.html#module-sklearn.gaussian_process
kernels = [kernel1, kernel2, kernel3, kernel4, kernel5, kernel6, kernel7, kernel8, kernel9]
# Prepare lists for storing the predictions
sigma_pred_list = []
ratio_list = []
std_dev_list = []
chi2_list = []

# Instantiate a Gaussian Process model
print("\nInstantiate a Gaussian Process model")

for kernel_input, name in zip(kernels, kernel_names):

        # Define each time the kernel as product with the constant kernel
        kernel = C(1.0, (1e-3, 1e3)) * kernel_input
        gp = GaussianProcessRegressor(kernel = kernel, alpha = 1e-4 * sigma, n_restarts_optimizer = args.restarts)

        # Fit to data using the Maximum Likelihood Estimation of the parameters
        print("\nFitting Gaussian Process: {}".format(name))
        gp.fit(tau, sigma)

        # Make the predictin on the meshed x-axis (ask for MSE as well)
        print("\nMaking predictions")
        sigma_pred, std_dev = gp.predict(tau_exact, return_std = True)
                
        # Computing ratio GP vs exact
        ratio_GP = abs(sigma_pred / sigma_exact)
        print("sigma pred: {}".format(sigma_pred.shape))
        print("sigma exact: {}".format(sigma_exact.shape))
        print("ratio: {}".format(ratio_GP.shape))

        # Computing chi2 estimator
        num = pow(sigma_exact - sigma_pred, 2)
        # Correct negative variances
        corrected_std = np.maximum(pow(std_dev , 2), 1e-10)
        chi2 = np.sum(num / corrected_std) / len(sigma_pred)

        # Store the predictions in the prepared lists
        sigma_pred_list.append(sigma_pred)
        ratio_list.append(ratio_GP)
        std_dev_list.append(std_dev)
        chi2_list.append(chi2)

# Computing chi2 for the different kernels
print("\nChecking performance of different kernels")
for chi2, name in zip(chi2_list, kernel_names):
        print("{}: chi2 = {}".format(name, chi2))

# Choose the model best chi2
min_chi2 = min(chi2_list)
i = 0
example = 0
for chi2 in chi2_list:
        if chi2 == min_chi2:
                example = i
        i = i+1

# Plot the function, the prediction and the confidence intervals for the best model
print("\nBest model found for\nExample = {}, Kernel = {}, chi2 = {}".format(example, kernel_names[example], chi2_list[example]))
plt.figure()
gs = gridspec.GridSpec(2, 1, height_ratios = [4, 1])
plt.subplot(gs[0])
if args.order == "lo":
        plt.title("Higgs gg - LO inclusive")
elif args.order == "nlo":
        plt.title("Higgs gg - NLO inclusive")
elif args.order == "nnlo":
        plt.title("Higgs gg - NNLO inclusive")
elif args.order == "nnnlo":
        plt.title("Higgs gg - N3LO inclusive")
plt.plot(tau_exact, sigma_exact, "r:", markersize = 10, label = "FO theory")
plt.errorbar(tau, sigma, dsigma, fmt = "r.", markersize = 10, label = "Asymptotics")
plt.plot(tau_exact, sigma_pred_list[example], "b-", label = "GP Prediction")
plt.fill(np.concatenate([tau_exact, tau_exact[::-1]]),
        np.concatenate([sigma_pred_list[example] - 1.9600 * std_dev_list[example],
                        (sigma_pred_list[example] + 1.9600 * std_dev_list[example])[::-1]]),
        alpha = 0.5, fc = "gray", ec = "None", label = "95% confidence interval")
plt.fill(np.concatenate([tau_exact, tau_exact[::-1]]),
        np.concatenate([sigma_pred_list[example] - 1.0 * std_dev_list[example],
                        (sigma_pred_list[example] + 1.0 * std_dev_list[example])[::-1]]),
        alpha = 0.5, fc = "b", ec = "None", label = "68% confidence interval")
plt.xscale("log")
plt.xlabel(r"$m_{H}^{2} / s$")
plt.ylabel(r"$\sigma(\tau) \quad (pb)$")
plt.legend(loc = "upper right")

# Plot ratio of GP vs the exact theory
plt.subplot(gs[1])
plt.plot(tau_exact, ratio_GP, "g", markersize = 10)
plt.hlines(1, tau_exact[0], tau_exact[-1], color = "r")
plt.fill(np.concatenate([tau_exact, tau_exact[::-1]]),
        np.concatenate([ratio_list[example] - abs(1 / sigma_exact * 1.9600 * std_dev_list[example]),
                        (ratio_list[example] + abs(1 / sigma_exact * 1.9600 * std_dev_list[example]))[::-1]]),
        alpha = 0.5, fc = "b", ec = "None", label = "95% confidence interval")
plt.ylim((-20, 20))
plt.xscale("log")
plt.xlabel(r"$m_{H}^{2} / s$")
plt.ylabel(r"$\frac{\sigma_{GP}}{\sigma_{theory}}$")
plt.show()