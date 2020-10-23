# Gaussian Processes
# Higgs pt distribution

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

# Setting parameters
np.random.seed(1)
n = 2

# True function to approximate
def sigma_exact(x):
        # Parameters of the Hypergeometric
        a, b, c = 0.5, 2.0, 2.5
        y = pow((np.sqrt(x + 1) - np.sqrt(x)), 4)
        return sc.beta(n, 0.5) * sc.hyp2f1(a, b, c, y) / x

def sigma_small_pt(x):
        # Asymptotic behavior of the FO cross section at low pt
        return  -0.5 * np.log(x) / x - (sc.polygamma(0, n) + np.euler_gamma) / x

def sigma_large_pt(x):
        # Asymptotic behavior of the FO cross section at low pt
        return sc.beta(n, 0.5) / x

# Setting normalization conditions
parser = argparse.ArgumentParser()
parser.add_argument("-n", "--normalization", help = "Normalization", default = "True", type = str)
parser.add_argument("-r", "--restarts", help = "Number of restarts to ask the GP for", default = 100, type = int)
parser.add_argument("-s", "--summary", help = "Show summary and performance", default = "False", type = str)
parser.add_argument("-k", "--k_folding_splits", help = "Number of splits done in the k-folding", default = 5, type = int)
parser.add_argument("-vp", "--validation_points", help = "Number of points dedicated to validation", default = 4, type = int)
args = parser.parse_args()

print("Gaussian Processes - Higgs pt distribution")

# Known points in our pt range
pt_small = np.logspace(-5, -3, 10)
pt_large = np.logspace(0.5, 1.5, 10)
pt = np.concatenate((pt_small, pt_large), axis = None).reshape(-1, 1)
# Known values of the cross section1
sigma_small_pt = sigma_small_pt(pt_small)
sigma_large_pt = sigma_large_pt(pt_large)
sigma = np.concatenate((sigma_small_pt, sigma_large_pt), axis = None)
# Mesh the input space for evaluations of the real functions, the prediction and its MSE
pt_exact = np.logspace(-5, 1.5, 1000).reshape(-1, 1)
sigma_exact = sigma_exact(pt_exact).reshape(1000)
print("pt: {}".format(pt.shape))
print("sigma: {}".format(sigma.shape))
print("pt mesh: {}".format(pt_exact.shape))
print("sigma exact: {}".format(sigma_exact.shape))

if args.normalization == "True":
        # Normalizing data files
        print("\nNormalizing data files")
        norm = np.sum(sigma_exact)
        sigma = sigma / norm
        sigma_exact = sigma_exact / norm
        # sigma = np.log10(sigma)
        # sigma_exact = np.log10(sigma_exact)

# Calling HyperGP for hyperoptimization
print("\nPreparing hyperopt:")
# Set list of parameters to be hyperoptimized
# https://scikit-learn.org/stable/modules/classes.html#module-sklearn.gaussian_process
tuned_parameters = [ {"kernel": ["dotprod", "exp", "matern", "expsin", "pairwise", "prod", "rbf", "rquad", "sum"],
                "length": [0.9, 1, 1.1], 
                "length_bounds": [(1e-5, 1e5)]} ]

# Preparing k-folding
print("\nPreparing k-folding:")
# Split the training set in k different folds
# k-1 will be used for training and remaining one for validation
val_list = []
train_list = []
len_train = pt.shape[0]
# Generate the validation sets
for i in range(args.k_folding_splits):
        # Randomly select 4 elements of the training set
        random_choice = list(np.random.choice(range(len_train), args.validation_points, replace = False))
        val_list.append(random_choice)
# Generate training set according to the validation sets
full = list(range(len_train))
full_set = set(full)
for i in val_list:
        val_set = set(i)
        train_set = full_set - val_set
        train_list.append(list(train_set))

# Instantiate the GridSearchCV with k-folding
print("Instantiate GridSearchCV")
gp_hyper = GridSearchCV(HyperGP(n_restarts_optimizer = args.restarts), tuned_parameters)
gp_hyper.fit(np.log10(pt), sigma)

# Compute means and standard deviations
all_params = gp_hyper.cv_results_['params']
means = gp_hyper.cv_results_["mean_test_score"]
stds = gp_hyper.cv_results_["std_test_score"]

if args.summary == "True":
        # Check scores for each parameter configuration
        print("\nSummary")
        for mean, std, params in zip(means, stds, all_params):
                print("%0.3f (+/-%0.03f); for %r"
                % (mean, std * 2, params))

# Extract the best set of parameters for our model
print("\nBest parameters set:")
print(gp_hyper.best_params_)
best_model_dict = gp_hyper.cv_results_['params'][gp_hyper.best_index_]

# Instantiate a Gaussian Process with the best parameter set
print("\nInstantiate a Gaussian Process model")
# Since we are using the parameters through the HyperGP class
# the GaussianProcessRegressor needs to be instantiated through the same class
best_model = HyperGP(n_restarts_optimizer = args.restarts)
best_model.set_params(**best_model_dict)

# Fit to data using the Maximum Likelihood Estimation of the parameters
print("\nFitting Gaussian Process")
best_model.fit(np.log10(pt), sigma)

# Make the prediction on the meshed x-axis (ask for MSE as well)
print("\nMaking predictions")
sigma_pred, std_dev = best_model.predict(np.log10(pt_exact), return_std = True)

# Computing ratio GP vs exact
ratio_GP = abs(sigma_pred / sigma_exact)
print("sigma pred: {}".format(sigma_pred.shape))
print("sigma exact: {}".format(sigma_exact.shape))
print("ratio: {}".format(ratio_GP.shape))

pdb.set_trace()

# Plot the function, the prediction and the confidence intervals for the best model
plt.figure()
gs = gridspec.GridSpec(2, 1, height_ratios = [4, 1])
plt.subplot(gs[0])
plt.title("Higgs pt distribution - LO in pt")
# plt.text(6, 0.005, "chi2 = {}\nk-folding score = {}".format(round(chi2, 5), round(score_mean, 5)), fontsize = 6)
plt.plot(pt_exact, sigma_exact, "r:", markersize = 10, label = "FO theory")
plt.plot(pt, sigma, "r.", markersize = 10, label = "Asymptotics")
plt.plot(pt_exact, sigma_pred, "b-", label = "GP Prediction")
plt.fill(np.concatenate([pt_exact, pt_exact[::-1]]),
        np.concatenate([sigma_pred - 1.9600 * std_dev,
                        (sigma_pred + 1.9600 * std_dev)[::-1]]),
        alpha = 0.5, fc = "gray", ec = "None", label = "95% confidence interval")
plt.fill(np.concatenate([pt_exact, pt_exact[::-1]]),
        np.concatenate([sigma_pred - 1.0 * std_dev,
                        (sigma_pred + 1.0 * std_dev)[::-1]]),
        alpha = 0.5, fc = "b", ec = "None", label = "68% confidence interval")
plt.xscale("log")
plt.xlabel(r"$p_{\perp}^{2} / m_{h}^{2}$")
plt.ylabel(r"$\frac{d\sigma}{dp_{\perp}}(N, p_{\perp})$")
plt.legend(loc = "upper right", fontsize = 6)

# Plot ratio of GP vs the exact theory
plt.subplot(gs[1])
plt.plot(pt_exact, ratio_GP, "g", markersize = 10)
plt.hlines(1, pt_exact[0], pt_exact[-1], color = "r")
plt.fill(np.concatenate([pt_exact, pt_exact[::-1]]),
        np.concatenate([ratio_GP - abs(1 / sigma_exact * 1.9600 * std_dev),
                        (ratio_GP + abs(1 / sigma_exact * 1.9600 * std_dev))[::-1]]),
        alpha = 0.5, fc = "b", ec = "None", label = "95% confidence interval")
plt.ylim((-20, 20))
plt.xscale("log")
plt.xlabel(r"$p_{\perp}^{2} / m_{h}^{2}$")
plt.ylabel(r"$\frac{\sigma_{GP}}{\sigma_{theory}}$")
plt.show()