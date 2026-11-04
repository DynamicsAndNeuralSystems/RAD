import warnings

import numpy as np


def RAD(x, tau=1, doCentre=True):
    """
    Harris, Gollo, & Fulcher's rescaled auto-density (RAD) noise-insensitive
    metric for inferring the distance to criticality.

    Parameters
    ----------
    x: array
        A time-series input vector

    tau: integer
        The embedding and differencing delay in units of the timestep

    doCentre : boolean
        Whether to centre the time series at its median and take absolute
        values (use this unless the time series represents radial values)

    Returns
    -------
    The RAD feature value, or NaN (with a warning) when either half of the
    median split has fewer than two values or no spread, as happens for short
    or heavily tied (e.g., integer-valued) time series.
    """

    # Before the arguments were unified across languages, the signature was
    # RAD(x, centre=False, tau=1): catch old-style positional calls
    if isinstance(tau, (bool, np.bool_)):
        raise TypeError("RAD(x, tau=1, doCentre=True): tau must be an integer "
                        "(the argument order changed from RAD(x, centre, tau))")
    if not float(tau).is_integer() or tau < 1:
        raise ValueError("tau must be a positive integer")
    tau = int(tau)
    if doCentre not in (True, False):
        raise TypeError("RAD(x, tau=1, doCentre=True): doCentre must be True or False "
                        "(the argument order changed from RAD(x, centre, tau))")

    # ensure that x is a 1-D numpy array
    x = np.asarray(x, dtype=float)
    if x.ndim != 1:
        if x.ndim == 2 and 1 in x.shape:
            x = x.ravel()
        else:
            raise ValueError("x must be a vector")

    # if specified: centre the time series and take the absolute value
    if doCentre:
        x = x - np.median(x)
        x = np.abs(x)

    # Delay embed at interval tau
    y = x[tau:]
    x = x[:-tau]

    # Median split
    subMedians = x < np.median(x)
    superMedianValues = x[~subMedians]
    subMedianValues = x[subMedians]
    if len(superMedianValues) < 2 or len(subMedianValues) < 2:
        warnings.warn("RAD: a half of the median split has fewer than two values; returning NaN",
                      RuntimeWarning, stacklevel=2)
        return np.nan
    superMedianSD = np.std(superMedianValues, ddof=1)
    subMedianSD = np.std(subMedianValues, ddof=1)
    if superMedianSD == 0 or subMedianSD == 0:
        warnings.warn("RAD: a half of the median split has no spread (tied values); returning NaN",
                      RuntimeWarning, stacklevel=2)
        return np.nan

    # Properties of the auto-density
    sigma_dx = np.std(y - x, ddof=1)
    densityDifference = (1/superMedianSD) - (1/subMedianSD)

    # return RAD
    return sigma_dx * densityDifference
