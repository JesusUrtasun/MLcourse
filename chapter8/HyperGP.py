import pdb
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel as C
from sklearn.gaussian_process.kernels import RBF, Matern, ExpSineSquared
from sklearn.gaussian_process.kernels import DotProduct, RationalQuadratic, PairwiseKernel
from sklearn.gaussian_process.kernels import Sum, Product, Exponentiation

# List of params to hyperoptimize over not included by default in the GPR
EXTRA_PARAMS = ["n_small_stop", "n_large_stop", "length", "length_bounds", "sigma_0", "periodicity", "nu", "gamma", "alpha_rquad"]

# Inherited class from GaussianProcessRegressor
class HyperGP(GaussianProcessRegressor):
    """
        This class inherits from GaussianProcessRegressor and it is used in order
        to hyperoptimize the parameters of the kernels and the GPR at the same time
    """

    # Init method
    def __init__(self, n_small_stop = None, n_large_stop = None, kernel = None, length = None, length_bounds = None,
                sigma_0 = None, periodicity = None, nu = None, gamma = None, alpha_rquad = None, **kwargs):

        # Set parameters as class attributes - equivalent to doing self.param = None for each param
        for param in EXTRA_PARAMS:
            setattr(self, param, None)

        # Not needed in principle. Some problem is going on with super()
        # Kernel and length will be given by the tuned_parameters list
        super().__init__(kernel = kernel, **kwargs)

    # Override the get_params method called at the beggining of each fit
    def get_params(self, *args, **kwargs):
        
        # Call get_param() of the abstract class
        gp_params = super().get_params(*args, **kwargs)
        # Add length as an accepted parameter
        for param in EXTRA_PARAMS:
            if param not in gp_params:
                gp_params[param] = getattr(self, param)

        return gp_params
        
    # Define score used for the grid search
    # Arguments x, y are already the validation lists x_val, y_val
    def score(self, x, y):

        # Make the prediction on the meshed x-axis (ask for MSE as well)
        y_pred, std_dev = super().predict(x, return_std = True)

         # Correct negative variances
        corrected_std = np.maximum(pow(std_dev, 2), 1e-15)
    
        # Compute chi2 estimator
        num = pow(y - y_pred, 2)
        chi2 = np.sum(num / corrected_std) / len(y_pred)

        # Define score
        score = 1 - chi2

        # Trying statistical distances
        # kl_score = np.sum( y_pred * np.log10(y / y_pred) )
        # hell_score = 0.5 * np.sum(pow((y - y_pred), 2))

        return score

    # Override the set_params method called during of each fit
    # Kernels will be always given in the list, but set a default for length
    def set_params(self, n_small_stop = None,  n_large_stop = None, kernel = None, length = 1.0, length_bounds = (1e-5, 1e5),
                    sigma_0 = 1.0, periodicity = 1.0, gamma = 1.0, alpha_rquad = 1.0, nu = 1.0, **kwargs):

        # Define the interpolation region in N
        self.n_small_stop = n_small_stop
        self.n_large_stop = n_large_stop
        self.kernel_name = kernel

        # Kernels taken from tuned_parameters list
        if kernel == "dotprod":
            new_kernel = C(length, (1e-3, 1e3)) * DotProduct(sigma_0 = sigma_0, sigma_0_bounds = (1e-5, 1e5))
        elif kernel == "expsin":
            new_kernel = C(1.0, (1e-3, 1e3)) * ExpSineSquared(length_scale = length, periodicity = periodicity, length_scale_bounds = length_bounds)
        elif kernel == "exp":
            new_kernel = C(1.0, (1e-3, 1e3)) * Exponentiation(RBF(length_scale = length, length_scale_bounds = length_bounds), 2)
        elif kernel == "matern":
            new_kernel = C(1.0, (1e-3, 1e3)) * Matern(length_scale = length, length_scale_bounds = length_bounds, nu = nu)
        elif kernel == "pairwise":
            new_kernel = C(1.0, (1e-3, 1e3)) * PairwiseKernel(gamma = gamma, gamma_bounds = (1e-5, 1e5))
        elif kernel == "rbf":
            new_kernel = C(1.0, (1e-3, 1e3)) * RBF(length_scale = length, length_scale_bounds = length_bounds)
        elif kernel == "rquad":
            new_kernel = C(1.0, (1e-3, 1e3)) * RationalQuadratic(length_scale = length, alpha = alpha_rquad, length_scale_bounds = length_bounds, alpha_bounds = (1e-5, 1e5))
        # Combinations of basic kernels
        elif kernel == "prod":
            new_kernel = C(1.0, (1e-3, 1e3)) * Product(RBF(length, length_bounds), Matern(length, length_bounds, nu = nu))
        elif kernel == "sum":
            new_kernel = C(1.0, (1e-3, 1e3)) * Sum(RBF(length, length_bounds), Matern(length, length_bounds, nu = nu))
        elif kernel == "prod2":
            new_kernel = C(1.0, (1e-3, 1e3)) * Product(RationalQuadratic(length_scale = length, alpha = alpha_rquad, length_scale_bounds = length_bounds, alpha_bounds = (1e-5, 1e5)),
                                                Matern(length, length_bounds, nu = nu))
        elif kernel == "sum2":
            new_kernel = C(1.0, (1e-3, 1e3)) * Sum(RationalQuadratic(length_scale = length, alpha = alpha_rquad, length_scale_bounds = length_bounds, alpha_bounds = (1e-5, 1e5)),
                                                Matern(length, length_bounds, nu = nu))
        else:
            raise NotImplementedError(f"Kernel not implemented {kernel}")

        # Call the super() class (GaussianProcessRegressor)
        # and give manually the parameteres taken from tuned_parameters list
        result = super().set_params(kernel = new_kernel, **kwargs)

        return result