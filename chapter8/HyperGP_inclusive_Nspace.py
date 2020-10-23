# Gaussian Processes
# Higgs inclusive partonic

"""Doc string - None"""
# print(__doc__)

import pdb
import argparse
import numpy as np
from scipy import special as sc
from matplotlib import pyplot as plt
from matplotlib import gridspec
from sklearn import datasets
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel as C
from sklearn.gaussian_process.kernels import RBF, Matern, ExpSineSquared
from sklearn.gaussian_process.kernels import DotProduct, RationalQuadratic, PairwiseKernel
from sklearn.gaussian_process.kernels import Sum, Product, Exponentiation
from sklearn.model_selection import GridSearchCV, cross_validate
from HyperGP import HyperGP

from asymptotics import asymptotics

# Setting normalization conditions
parser = argparse.ArgumentParser()
parser.add_argument("-o", "--order", help = "Perturbative order: nlo, nnlo, n3lo", required = True, type = str)
parser.add_argument("-n", "--normalization", help = "Normalization", default = "True", type = str)
parser.add_argument("-cd", "--check_data", help = "Show data use for training and validate", default = "False", type = str)
parser.add_argument("-s", "--summary", help = "Show summary and performance", default = "False", type = str)
parser.add_argument("-r", "--restarts", help = "Number of restarts to ask the GP for", default = 0, type = int)
parser.add_argument("-vp", "--validation_points", help = "Number of points dedicated to validation", default = 6, type = int)
args = parser.parse_args()

print("Gaussian Processes - Higgs inclusive")

# Loading data files
print("\nLoading datafiles")

# Computed using ggHiggs https://www.ge.infn.it/~bonvini/higgs/
ggHiggs_data = np.loadtxt("data/ggHiggs_mesh_correct.txt").T
n_exact = ggHiggs_data[0].reshape(-1, 1)

# Known values of the cross section - NLO inclusive
if args.order == "nlo":
        # Exact coefficient function - mesh the input for the prediction and its MSE
        sigma_exact = ggHiggs_data[1]
        # Known points in our n range
        n_small_stop = 1.075
        n_large_stop = 2.1
        n_points = 25
        n_small = np.linspace(n_exact[0], n_small_stop, n_points)
        n_large = np.linspace(n_large_stop, n_exact[-1], 100)
        n = np.concatenate((n_small, n_large), axis = None).reshape(-1, 1)
        # Compute coefficient functions from asymptotics
        nlo_small_n_abf = asymptotics.nlo_small_n_abf(n_small)
        nlo_large_n_soft = asymptotics.nlo_large_n(n_large, "soft2")
        sigma = np.concatenate((nlo_small_n_abf, nlo_large_n_soft), axis = None)
if args.order == "nnlo":
        # Exact coefficient function - mesh the input for the prediction and its MSE
        sigma_exact = ggHiggs_data[2]
        # Known points in our n range
        n_small_stop = 1.075
        n_large_stop = 2.1
        n_points = 25
        n_small = np.linspace(n_exact[0], n_small_stop, n_points)
        n_large = np.linspace(n_large_stop, n_exact[-1], 100)
        n = np.concatenate((n_small, n_large), axis = None).reshape(-1, 1)
        # Compute coefficient functions from asymptotics
        nnlo_small_n_abf = asymptotics.nnlo_small_n_abf(n_small)
        nnlo_large_n_soft = asymptotics.nnlo_large_n(n_large, "soft2")
        sigma = np.concatenate((nnlo_small_n_abf, nnlo_large_n_soft), axis = None)
if args.order == "n3lo":
        # Exact coefficient function - mesh the input for the prediction and its MSE
        sigma_exact = ggHiggs_data[3]
        # Known points in our n range
        n_small_stop = 1.175
        n_large_stop = 2.5
        n_points = 25
        n_small = np.linspace(n_exact[0], n_small_stop, n_points)
        n_large = np.linspace(n_large_stop, n_exact[-1], 100)
        n = np.concatenate((n_small, n_large), axis = None).reshape(-1, 1)
        # Compute coefficient functions from asymptotics
        n3lo_small_n_abf = asymptotics.n3lo_small_n_abf(n_small)
        n3lo_large_n_soft = asymptotics.n3lo_large_n(n_large, "soft2")
        sigma = np.concatenate((n3lo_small_n_abf, n3lo_large_n_soft), axis = None)

if args.check_data == "True":

        # Print data structure
        print("n: {}".format(n.shape))
        print("sigma: {}".format(sigma.shape))
        print("n mesh: {}".format(n_exact.shape))
        print("sigma exact: {}".format(sigma_exact.shape))

        # Check asymptotics (trainig) and exact theory (validation)
        plt.figure()
        if args.order == "nlo":
                plt.plot(n_exact, sigma_exact, "r:", label = "NLO ggHiggs")
                plt.plot(n, sigma, "b.", label = "NLO asymptotics")
                plt.xlim((1, 4))
                plt.ylim((0, 50))
                plt.xlabel(r"$N$")
                plt.ylabel(r"$C^{(1)}(N)$")
                plt.legend(loc = "upper right", fontsize = 7)
                # plt.savefig("../../Update_Notes/2020_04_HyperGP_inclusive/plots/nlo_asymptotics.png")
        if args.order == "nnlo":
                plt.plot(n_exact, sigma_exact, "r:", label = "NNLO ggHiggs")
                plt.plot(n, sigma, "b.", label = "NNLO asymptotics")
                plt.xlim((1, 4))
                plt.ylim((0, 500))
                plt.xlabel(r"$N$")
                plt.ylabel(r"$C^{(2)}(N)$")
                plt.legend(loc = "upper right", fontsize = 7)
                # plt.savefig("../../Update_Notes/2020_04_HyperGP_inclusive/plots/nnlo_asymptotics.png")
        if args.order == "n3lo":
                plt.plot(n_exact, sigma_exact, "r:", label = "N3LO ggHiggs")
                plt.plot(n, sigma, "b.", label = "N3LO asymptotics")
                plt.xlim((1, 4))
                plt.ylim((0, 5000))
                plt.xlabel(r"$N$")
                plt.ylabel(r"$C^{(2)}(N)$")
                plt.legend(loc = "upper right", fontsize = 7)
                # plt.savefig("../../Update_Notes/2020_04_HyperGP_inclusive/plots/n3lo_asymptotics.png")
        plt.show()

        pdb.set_trace()

if args.normalization == "True":

        # Normalize data files
        print("\nNormalizing data files")
        # norm = np.sum(sigma_exact)    
        # sigma = sigma / norm
        # sigma_exact = sigma_exact / norm
        sigma = np.log10(sigma)
        sigma_exact = np.log10(sigma_exact)

# Calling HyperGP for hyperoptimization
print("\nPreparing hyperopt")
# https://scikit-learn.org/stable/modules/classes.html#module-sklearn.gaussian_process
tuned_parameters = [
        # {"kernel": ["dotprod"],
        # "length": [i/10 for i in range(1, 15)],
        # "sigma_0": [i/ 10 for i in range(5, 15)]
        # },
        # {"kernel": ["expsin"],0.0
        # "length": [i/10 for i in range(1, 15)],
        # "periodicity": [i/ 10 for i in range(5, 15)]
        # },
        # {"kernel": ["exp"],
        # "length": [i/10 for i in range(1, 15)]
        # },
        {"kernel": ["matern"],
        "length": [i/100 for i in range(40, 90)],
        "nu": [i/ 100 for i in range(60, 110)]
        },
        # {"kernel": ["pairwise"],
        # "length": [i/10 for i in range(5, 15)],
        # "gamma": [i/ 10 for i in range(15, 25)]
        # },
        # {"kernel": ["rbf"],
        # "length": [i/10 for i in range(5, 20)],
        # },
        # {"kernel": ["rquad"],
        # "length": [i/10 for i in range(5, 15)],
        # "alpha_rquad": [i/ 10 for i in range(5, 15)]
        # },
        # {"kernel": ["prod"],
        # "length": [i/10 for i in range(5, 15)],
        # "nu": [i/ 10 for i in range(5, 25)],
        # },
        # {"kernel": ["sum"],
        # "length": [0.5],
        # "nu": [1.6],
        # },
        {"kernel": ["prod2"],
        "length": [i/100 for i in range(100, 150)],
        "nu": [i/100 for i in range(150, 200)],
        "alpha_rquad": [i/100 for i in range(100, 150)]
        },
        {"kernel": ["sum2"],
        "length": [i/100 for i in range(100, 150)],
        "nu": [i/100 for i in range(200, 250)],
        "alpha_rquad": [i/100 for i in range(50, 100)]
        }
        ]

# Preparing k-folding
val_list = []
train_list = []
len_train = n.shape[0]

# Generate the validation sets
my_choice_n_small = [(len(n_small) - 1 - i) for i in range(int(args.validation_points))]
my_choice_n_large = [(len(n_small) + i) for i in range(int(args.validation_points))]
val_list.append(my_choice_n_small)
val_list.append(my_choice_n_large)

# Generate training set according to the validation sets
full = list(range(len_train))
full_set = set(full)
for i in val_list:
        val_set = set(i)
        train_set = full_set - val_set
        train_list.append(list(train_set))

# Instantiate the GridSearchCV with k-folding
print("Instantiate GridSearchCV")
gp_hyper = GridSearchCV(HyperGP(n_restarts_optimizer = args.restarts), tuned_parameters, cv = zip(train_list, val_list))
gp_hyper.fit(n, sigma)

# Compute means and standard deviations
all_params = gp_hyper.cv_results_['params']
means = gp_hyper.cv_results_["mean_test_score"]
stds = gp_hyper.cv_results_["std_test_score"]

if args.summary == "True":
        # Check scores for each parameter configuration
        # print("\nSummary")
        if args.order == "nlo":
                # f = open("summary/nlo_hyperopt2.txt", "a")
                # f.write("\nInterpolation N range: ({}, {}) - {} points".format(n_small_stop, n_large_stop, n_points))
                for mean, std, params in zip(means, stds, all_params):
                        print("%0.3f (+/-%0.03f); for %r" % (mean, std * 2, params))
                #         f.write("\n%0.3f (+/-%0.03f); for %r" % (mean, std * 2, params))
                # f.close()
        if args.order == "nnlo":
                # f = open("summary/nnlo_hyperopt2.txt", "a")
                # f.write("\nInterpolation N range: ({}, {}) - {} points".format(n_small_stop, n_large_stop, n_points))
                for mean, std, params in zip(means, stds, all_params):
                        print("%0.3f (+/-%0.03f); for %r" % (mean, std * 2, params))
                #         f.write("\n%0.3f (+/-%0.03f); for %r" % (mean, std * 2, params))
                # f.close()
        if args.order == "n3lo":
                f = open("summary/n3lo_hyperopt3.txt", "a")
                f.write("\nInterpolation N range: ({}, {}) - {} points".format(n_small_stop, n_large_stop, n_points))
                for mean, std, params in zip(means, stds, all_params):
                        print("%0.3f (+/-%0.03f); for %r" % (mean, std * 2, params))
                        f.write("\n%0.3f (+/-%0.03f); for %r" % (mean, std * 2, params))
                f.close()

pdb.set_trace()

# Extract the best set of parameters for our model
print("\nBest parameters set:")
print(gp_hyper.best_params_)
best_model_dict = gp_hyper.cv_results_['params'][gp_hyper.best_index_]

# Instantiate a Gaussian Process with the best parameter set
print("\nInstantiate a Gaussian Process model")
# Since we are using the parameters through the HyperGP class
# we need to also instantiate the GaussianProcessRegressor through the same class
best_model = HyperGP(n_restarts_optimizer = args.restarts)
best_model.set_params(**best_model_dict)

# Fit to data using the Maximum Likelihood Estimation of the parameters
print("\nFitting Gaussian Process")
best_model.fit(n, sigma)

# Make the prediction on the meshed x-axis (ask for MSE as well)
print("\nMaking predictions")
sigma_pred, std_dev = best_model.predict(n_exact, return_std = True)

# Compute chi2 estimator
num = pow(sigma_exact - sigma_pred, 2)
chi2 = np.sum(num / std_dev) / len(sigma_pred)
# Define loss of the best model after training
score = chi2
print("chi2 = {}".format(chi2))
if args.order == "nlo":
        f = open("summary/nlo_summary.txt", "a")
        f.write("\nInterpolation N range: ({}, {}) - {} points".format(n_small_stop, n_large_stop, n_points))
        f.write("\nBest dict: {} chi2 = {:.5f}".format(str(best_model_dict), chi2))
        f.close()
if args.order == "nnlo":
        f = open("summary/nnlo_summary.txt", "a")
        f.write("\nInterpolation N range: ({}, {}) - {} points".format(n_small_stop, n_large_stop, n_points))
        f.write("\nBest dict: {} chi2 = {:.5f}".format(str(best_model_dict), chi2))
        f.close()
if args.order == "n3lo":
        f = open("summary/n3lo_summary.txt", "a")
        f.write("\nInterpolation N range: ({}, {}) - {} points".format(n_small_stop, n_large_stop, n_points))
        f.write("\nBest dict: {} chi2 = {:.5f}".format(str(best_model_dict), chi2))
        f.close()

# Computing ratio GP vs exact
ratio_GP = abs(sigma_pred / sigma_exact)
print("sigma pred: {}".format(sigma_pred.shape))
print("sigma exact: {}".format(sigma_exact.shape))
print("ratio: {}".format(ratio_GP.shape))

# Process data and predictions for plot
sigma = pow(10, sigma)
sigma_pred = pow(10, sigma_pred)
sigma_exact = pow(10, sigma_exact)

# pdb.set_trace()

# Plot the function, the prediction and the confidence intervals for the best model
plt.figure()
gs = gridspec.GridSpec(2, 1, height_ratios = [4, 1])
plt.subplot(gs[0])
plt.plot(n_exact, sigma_exact, "r:", markersize = 10, label = "exact")
plt.plot(n, sigma, "r.", markersize = 10, label = "Asymptotics")
plt.plot(n_exact, sigma_pred, "b-", markersize = 10, label = "GP prediction")
plt.fill(np.concatenate([n_exact, n_exact[::-1]]),
        np.concatenate([sigma_pred - 1.9600 * std_dev,
                        (sigma_pred + 1.9600 * std_dev)[::-1]]),
        alpha = 0.5, fc = "gray", ec = "None", label = "95% confidence interval")
plt.fill(np.concatenate([n_exact, n_exact[::-1]]),
        np.concatenate([sigma_pred - 1.0 * std_dev,
                        (sigma_pred + 1.0 * std_dev)[::-1]]),
        alpha = 0.5, fc = "b", ec = "None", label = "68% confidence interval")
if args.order == "nlo":
        plt.xlim((1, 4))
        plt.ylim((0, 60))
        plt.xlabel(r"$N$")
        plt.ylabel(r"$C^{(1)}(N)$")
if args.order == "nnlo":
        plt.xlim((1, 4))
        plt.ylim((0, 600))
        plt.xlabel(r"$N$")
        plt.ylabel(r"$C^{(2)}(N)$")
if args.order == "n3lo":
        plt.xlim((1, 4))
        plt.ylim((0, 5000))
        plt.xlabel(r"$N$")
        plt.ylabel(r"$C^{(3)}(N)$")
plt.legend(loc = "upper right", fontsize = 6)

# Plot ratio of GP vs the exact theory
plt.subplot(gs[1])
plt.plot(n_exact, ratio_GP, "g", markersize = 10)
plt.hlines(1, n_exact[0], n_exact[-1], color = "r")
plt.fill(np.concatenate([n_exact, n_exact[::-1]]),
        np.concatenate([ratio_GP - abs(1 / sigma_exact * 1.9600 * std_dev),
                        (ratio_GP + abs(1 / sigma_exact * 1.9600 * std_dev))[::-1]]),
        alpha = 0.5, fc = "b", ec = "None", label = "95% confidence interval")
plt.xlim((1, 4))
plt.ylim((0, 2))
if args.order == "nlo":
        plt.xlabel(r"$N$")
        plt.ylabel(r"$\frac{C^{(1)}(N)_{GP}}{C^{(1)}(N)_{exact}}$")
        plt.savefig("../../Update_Notes/2020_04_HyperGP_inclusive/plots/nlo_test/nlo_test_{}_chi2:{}.png".format(str(best_model_dict), chi2))
if args.order == "nnlo":
        plt.xlabel(r"$N$")
        plt.ylabel(r"$\frac{C^{(2)}(N)_{GP}}{C^{(2)}(N)_{exact}}$")
        plt.savefig("../../Update_Notes/2020_04_HyperGP_inclusive/plots/nnlo_test/nnlo_test{}_chi2:{}.png".format(str(best_model_dict), chi2))
if args.order == "n3lo":
        plt.xlabel(r"$N$")
        plt.ylabel(r"$\frac{C^{(3)}(N)_{GP}}{C^{(3)}(N)_{exact}}$")
        plt.savefig("../../Update_Notes/2020_04_HyperGP_inclusive/plots/n3lo_test/n3lo_test{}_chi2:{}.png".format(str(best_model_dict), chi2))
plt.show()

# Prepare violin plot for hyperopt
# https://github.com/NNPDF/nnpdf/blob/master/n3fit/src/n3fit/hyper_optimization/plotting.py