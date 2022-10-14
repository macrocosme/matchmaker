import numpy as np
from scipy.optimize import curve_fit

# Powerlaw
powerlaw = lambda x, A, k: A * x**k
def fit_power_law(x, y, y_err=None, absolute_sigma=True):
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
def fit_linear_a_and_b(x, y, theta):
    a_fit,cov=curve_fit(linear, x[theta:], y[theta:])
    a, b = a_fit[0], a_fit[1]
    return a, b

def fit_linear_a_only(x, y, theta):
    a_fit,cov=curve_fit(linear_no_b, x[theta:], y[theta:])
    a = a_fit[0]
    return a

# Rayleigh
rayleigh = lambda x, sigma: (x/sigma**2) * np.exp(-(x**2 / (2*sigma**2)))
def fit_rayleigh(x, y, theta=None):
    # lmfit:
    # rmodel = Model(rayleigh)
    # params = rmodel.make_params(sigma=10)
    # out = rmodel.fit(y[:40]-np.min(y), params, x=x[:40])
    # _y = y * rayleigh(x, out.params['sigma'].value)
    sigma_fit,cov=curve_fit(rayleigh, x if theta is None else x[:theta], y if theta is None else y[:theta])

    return sigma_fit[0]
