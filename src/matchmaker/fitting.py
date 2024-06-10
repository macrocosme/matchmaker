import numpy as np
from scipy.optimize import curve_fit

# Powerlaw
powerlaw = lambda x, A, k: A * x**k
def fit_power_law(x, y, y_err=None, absolute_sigma=True):
    """
    Fits a power law function to the given data points.

    Parameters:
    - x: array-like
        The x-coordinates of the data points.
    - y: array-like
        The y-coordinates of the data points.
    - y_err: array-like, optional
        The uncertainties in the y-coordinates. If provided, the fit will take into account the errors.
    - absolute_sigma: bool, optional
        If True, the uncertainties in y are absolute. If False, they are relative.

    Returns:
    - a: float
        The fitted parameter 'a' in the power law function.
    - k: float
        The fitted parameter 'k' in the power law function.
    - d_a: float
        The uncertainty in the fitted parameter 'a'.
    - d_k: float
        The uncertainty in the fitted parameter 'k'.
    """
    from scipy.optimize import curve_fit
    if y_err is not None:
        a_fit, cov = curve_fit(powerlaw, x, y, sigma=y_err, absolute_sigma=absolute_sigma)
    else:
        a_fit, cov = curve_fit(powerlaw, x, y)
    a = a_fit[0]
    k = a_fit[1]
    d_a = np.sqrt(cov[0][0])
    d_k = np.sqrt(cov[1][1])

    return a, k, d_a, d_k

# Linear
linear = lambda x, a, b: a*x+b
linear_no_b = lambda x, a: a*x
def fit_linear_a_and_b(x, y, theta):
    """
    Fits a linear function to the given data points (x, y) starting from index theta.

    Parameters:
    x (array-like): The x-coordinates of the data points.
    y (array-like): The y-coordinates of the data points.
    theta (int): The starting index for fitting the linear function.

    Returns:
    tuple: A tuple containing the fitted parameters a and b, where a is the slope and b is the intercept of the linear function.
    """
    a_fit, cov = curve_fit(linear, x[theta:], y[theta:])
    a, b = a_fit[0], a_fit[1]
    return a, b

def fit_linear_a_only(x, y, theta):
    """
    Fits a linear function with only the slope parameter 'a' to the given data.

    Parameters:
    x (array-like): The x-values of the data points.
    y (array-like): The y-values of the data points.
    theta (int): The index from which to start fitting the data.

    Returns:
    float: The fitted value of the slope parameter 'a'.
    """
    a_fit, cov = curve_fit(linear_no_b, x[theta:], y[theta:])
    a = a_fit[0]
    return a

# Rayleigh
rayleigh = lambda x, sigma: (x/sigma**2) * np.exp(-(x**2 / (2*sigma**2)))
def fit_rayleigh(x, y, theta=None):
    """
    Fits a Rayleigh distribution to the given data points.

    Parameters:
    - x: The x-values of the data points.
    - y: The y-values of the data points.
    - theta: The optional index to limit the data points used for fitting.

    Returns:
    - sigma_fit: The fitted value of the Rayleigh distribution parameter sigma.
    """
    sigma_fit, cov = curve_fit(rayleigh, x if theta is None else x[:theta], y if theta is None else y[:theta])

    return sigma_fit[0]
