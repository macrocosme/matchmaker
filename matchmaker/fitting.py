import numpy as np
from scipy.optimize import curve_fit

# Powerlaw
powerlaw = lambda x, A, k: A * x**k
def fit_power_law(x, y, y_err=None, absolute_sigma=True):
    """Fit power law to data

    Parameters
    ----------
    x:
    y:
    y_err:
    absolute_sigma:

    Returns
    -------
        a:
        k:
        d_a:
        d_k:
    """
    from scipy.optimize import curve_fit
    if y_err is not None:
        a_fit,cov=curve_fit(powerlaw, x, y, sigma=y_err, absolute_sigma=absolute_sigma)
    else:
        a_fit,cov=curve_fit(powerlaw, x, y)
    a = a_fit[0]
    k = a_fit[1]
    d_a = np.sqrt(cov[0][0])
    d_k = np.sqrt(cov[1][1])

    return a, k, d_a, d_k

# Linear
linear = lambda x, a, b: a*x+b
linear_no_b = lambda x, a: a*x
def fit_linear_a_and_b(x:np.array, y:np.array, theta:int):
    """Fit linear relation, returning a and b

    Parameters
    ----------
    x:numpy.array
    y:numpy.array
    theta:int

    Returns
    -------
    a, b
    """
    a_fit,cov=curve_fit(linear, x[theta:], y[theta:])
    a, b = a_fit[0], a_fit[1]
    return a, b

def fit_linear_a_only(x:np.array, y:np.array, theta:int):
    """Fit linear relation, returning a only

    TODO: merge code with fit_linear_a_and_b

    Parameters
    ----------
    x:numpy.array
    y:numpy.array
    theta:int

    Returns
    -------
    a
    """
    a_fit,cov=curve_fit(linear_no_b, x[theta:], y[theta:])
    a = a_fit[0]
    return a

# Rayleigh
rayleigh = lambda x, sigma: (x/sigma**2) * np.exp(-(x**2 / (2*sigma**2)))
def fit_rayleigh(x:np.array, y:np.array, theta:int=None):
    """Fit Rayleigh distribution

    Parameters
    ----------
    x:numpy.array
    y:numpy.array
    theta:int

    Returns
    -------
    exponent TODO:double check if this is indeed the exponent
    """
    sigma_fit,cov=curve_fit(rayleigh, x if theta is None else x[:theta], y if theta is None else y[:theta])

    return sigma_fit[0]
