"""Module for statistics.

Written by Jesse Bloom."""


import re
import math
import os
import string
import random
import transcendental


def AllNone(xlist):
    """Returns True if all entries in a list are 'None'; 'False' otherwise."""
    for x in xlist:
        if x != None:
            return False
    else:
        return True



def Median(numlist):
    """Returns the median of a list of numbers.
        
    If any entries of the list are 'None' or '-', they are removed first."""
    assert isinstance(numlist, list)
    xlist = []  # make a copy of the list
    for x in numlist:
        if isinstance(x, (int, float)):
            xlist.append(x)
        elif x in [None, '-']:
            pass
        else:
            raise ValueError, "Invalid value of %r in list." % x
    if len(xlist) == 0:
        raise ValueError, "Empty list."
    xlist.sort()
    n = len(xlist)
    if n % 2 == 0: # even length list, average two middle entries
        med = xlist[n / 2] + xlist[n / 2 - 1]
        med = float(med) / 2.0
        return med
    else: # odd length list, get middle entry
        return xlist[n / 2]



def Mean(numlist):
    """Returns the mean of a list of numbers.
        
    If any entries of the list are 'None' or '-', they are removed first."""
    mean = 0.0
    n = 0
    for x in numlist:
        if x in [None, '-']:
            continue
        if not isinstance(x, (int, float)):
            raise ValueError, "Invalid entry of %r." % x
        mean += x
        n += 1
    assert n == len(numlist) - numlist.count(None) - numlist.count('-')
    if n <= 0:
        raise ValueError, "Empty list."
    return mean / float(n)



def Variance(numlist, unknownmean = False):
    """Returns the variance of a list of numbers.
    Call is 'v = Variance(numlist, [unknownmean = False])'
    'numlist' is a list of numbers.  If any entries of the list are
        'None' or '-', they are removed first.
    'unknownmean' is an optional switch, set to False by default, 
        that specifies that the mean is also being estimated
        from the sample.  In this case, the denominator is 
        N - 1 rather than N.
    'v' is returned as the variance."""
    sum = sum2 = 0.0
    n = 0
    for x in numlist:
        if x in [None, '-']:
            continue
        assert isinstance(x, int) or isinstance(x, float)
        sum += x
        sum2 += x * x
        n += 1
    if unknownmean:
        var = (1.0 / (n - 1)) * (sum2 - n * (sum / n)**2)
    else:
    	mean2 = sum2 / float(n)
       	mean = sum / float(n)
        var = mean2 - mean * mean
    assert var >= -0.00000001
    return var 



def StandardError(numlist):
    """Returns the standard error of a list of numbers.

    If any entries of the list are 'None' or '-', they are
        removed first.  This is just the standard deviation
        divided by the square root of the number of points.
    """
    return StandardDeviation(numlist) / math.sqrt(len(numlist))



def StandardDeviation(numlist):
    """Returns the standard deviation of a list of numbers.

    If any entries of the list are 'None' or '-', they are removed first."""
    return math.sqrt(Variance(numlist))



def PearsonCorrelation(xdata, ydata):
    """Computes the Pearson linear correlation between two data sets.
   
    Call is '(r, p, n) = PearsonCorrelation(xdata, ydata)'
    The input data is given in the two lists 'xdata' and 'ydata' which should
        be of the same length.  If entry i of either list is 'None', this
        entry is disregarded in both lists.
    Returns Pearson's correlation coefficient, the two-tailed P-value, 
        and the number of data points as a tuple '(r, p, n)'."""
    if len(xdata) != len(ydata):
        raise ValueError, "Data sets are of different lengths."
    xlist = []
    ylist = []
    xylist = []
    for i in range(len(xdata)):
        if xdata[i] != None and ydata[i] != None:
            xlist.append(xdata[i])
            ylist.append(ydata[i])
            xylist.append(xdata[i] * ydata[i])
    n = len(xlist)
    if n <= 1:
        raise ValueError, "One or less data points: %r." % n
    xmean = Mean(xlist)
    ymean = Mean(ylist)
    xymean = Mean(xylist)
    xsd = StandardDeviation(xlist)
    ysd = StandardDeviation(ylist)
    try:
        r = float(xymean - xmean * ymean) / (xsd * ysd)
    except ZeroDivisionError:
        raise ValueError, "Standard deviation is zero."
    assert -1.0000000001 <= r < 1.0000000001
    z = math.fabs(r) * math.sqrt(n) / math.sqrt(2.0)
    p = Complementary_Error_Function(z)
    assert 0 <= p <= 1.0
    return (r, p, n)



def Complementary_Error_Function(z):
    """Calculates the error function of z.

    The complementary error function of z is defined as:
        erfc(z) = 2 / sqrt(pi) * integral(e^(t^2) dt) where the integral 
        is from z to infinity.
    Can be used to calculate cumulative normal probabilities: given a 
        distribution with mean m and standard deviation s,
        the probability of observing x > m  when x > 0 is: 
        P = 0.5 * erfc((x - m) / (s * sqrt(2)))
    Calculated according to Chebyshev fitting given by Numerical Recipes 
        in C, page 220-221."""
    x = math.fabs(z)
    t = 1.0 / (1.0 + 0.5 * x)
    ans = t * math.exp(-x * x - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))))
    if z > 0.0:
        return ans
    else:
        return 2.0 - ans



def BetaDistribution(x, alpha, beta, a=0.0, b=1.0):
    """Returns the value of the beta distribution probability density function.

    The beta distribution probability density function is defined over
        the interval a <= x <= b as
        f(x) = (x - a)^(alpha - 1) (b - x)^(beta - 1) / (B(alpha, beta) (b - a)^(alpha + beta - 1))
        where B(alpha, beta) is the beta function.
    'x' is the number for which we want to compute the value of the beta distribution probability
        density function.  It must satisfy a <= x <= b.
    'alpha' specifies the alpha parameter (also sometimes written as p), must be nonnegative
    'beta' specifies the beta parameter (also sometimes written as q), must be nonnegative
    'a' specifies the lower limit of the range of the beta distribution (0 by default)
    'b' specifies the upper limit of the range of the beta distribution (1 by default).
    """
    assert alpha >= 0 and beta >= 0
    betafunc = transcendental.beta(alpha, beta)
    return (x - a)**(alpha - 1) * (b - x)**(beta - 1) / (betafunc * (b - a)**(alpha + beta - 1))


def DerivativeBetaDistribution(x, alpha, beta, a=0.0, b=1.0):
    """The derivative of the beta distribution probability density function with respect to x.
    
    The beta distribution probability density function is defined over
        the interval a <= x <= b as, and is evaluated by 'BetaDistribution' as
        f(x) = (x - a)^(alpha - 1) (b - x)^(beta - 1) / (B(alpha, beta) (b - a)^(alpha + beta - 1))
        where B(alpha, beta) is the beta function.
    'x' is the number for which we want to compute the value of the beta distribution probability
        density function.  It must satisfy a <= x <= b.
    'alpha' specifies the alpha parameter (also sometimes written as p), must be nonnegative
    'beta' specifies the beta parameter (also sometimes written as q), must be nonnegative
    'a' specifies the lower limit of the range of the beta distribution (0 by default)
    'b' specifies the upper limit of the range of the beta distribution (1 by default).

    This function returns the derivative of the beta distribution at the value x with
        respect to x.

    >>> (alpha, beta) = (1.4, 1.6)
    >>> (a, b) = (-20.0, 20.0)
    >>> x = 5.6
    >>> Bx = BetaDistribution(x, alpha, beta, a, b)
    >>> dBx = DerivativeBetaDistribution(x, alpha, beta, a, b)
    >>> dx = 0.01
    >>> Bx2 = BetaDistribution(x + dx, alpha, beta, a, b)
    >>> 1000 * abs(Bx2 - Bx - dx * dBx) < abs(Bx2 - Bx)
    True
    >>> dBx2 = Bx * ((alpha - 1) / (x - a) - (beta - 1) / (b - x))
    >>> abs((dBx - dBx2) / dBx) < 0.000000001
    True
    """
    assert a <= x <= b and a < b
    assert alpha >= 0 and beta >= 0
    betafunc = transcendental.beta(alpha, beta)
    return ((alpha - 1) * (x - a)**(alpha - 2) * (b - x)**(beta - 1) - (beta - 1) * (b - x)**(beta - 2) * (x - a)**(alpha - 1)) / (betafunc * (b - a)**(alpha + beta - 1))


def ParameterizeBetaDistribution(center, mode_mean, parameter_sum, a=0.0, b=1.0):
    """Determines the alpha and beta values for a beta distribution.

    This method is used to find the alpha and beta values to parameterize a beta
        distribution like the one defined by 'BetaDistribution' so that it is
        centered at a given value.
    'center' is a number specifying where the beta distribution should be centered.
        It must satisfy a <= center <= b.
    'mode_mean' is a string specifying how the beta distribution is centered.
        If the value is 'MODE' then the mode (peak) of the distribution is at
        'center'.  If the value is 'MEAN' then the mean (weighted average) of
        the distribution is at 'center'.
    'parameter_sum' is a number > 2 that specifies the sum of the alpha and
        beta parameters.
    'a' gives the minimum defined range for the beta distribution; 0 by default.
    'b' gives the maximum defined range for the beta distribution; 1 by default.
    The returned value is the 2-tuple '(alpha, beta)' where these are the
        beta distribution parameters that give a beta distribution that
        satisfies the specified criteria.
    """
    assert a <= center <= b and a < b
    assert parameter_sum > 2
    if mode_mean == 'MEAN':
        # formula is mean = a + (b - a) alpha / (alpha + beta)
        alpha = (center - a) / float(b - a) * parameter_sum
    elif mode_mean == 'MODE':
        # formula is mode = a + (b - a) (alpha - 1) / (alpha + beta - 2)
        alpha = (center - a) / float(b - a) * (parameter_sum - 2.0) + 1
    else:
        raise ValueError, "Undefined value for 'mode_mean': %s." % (str(mode_mean))
    beta = parameter_sum - alpha
    assert alpha >= 0 and beta >= 0
    return (alpha, beta)



def FitExponential(x, y):
    """Fits an exponential curve to a set of data points using least squares method.

    This method fits an exponential curve of the form y = a * exp(b * x) to
        a set of data.  The fitting is based on that suggested by MathWorld
        (http://mathworld.wolfram.com/LeastSquaresFittingExponential.html)
        with the points weighted equally.  That is, the method minimizes the
        function:
            sum_i yi (log yi - log a - b * xi)**2
    The calling variables are the lists x and y, which should be of the same 
        length.  Elements x[i] and y[i] are paired data points (i.e.
        y[i] = a * exp(b * x[i]).  Each of these points should be a number.
        The y values must be > 0.
    The returned variables are the values a and b that minimize the weighted
        least squares error in the equation y = a * exp(b * x).
    """
    if len(x) != len(y):
        raise ValueError, "x and y of different lengths."
    x2_y = y_log_y = x_y = x_y_log_y = y_sum = 0.0
    for (xi, yi) in zip(x, y):
        if yi <= 0.0:
            raise ValueError, "yi is not > 0."
        log_yi = math.log(yi)
        xi_yi = xi * yi
        x2_y += xi * xi_yi
        y_log_y += yi * log_yi
        x_y += xi_yi
        x_y_log_y += xi_yi * log_yi
        y_sum += yi
    log_a = (x2_y * y_log_y - x_y * x_y_log_y) / (y_sum * x2_y - x_y**2)
    b = (y_sum * x_y_log_y - x_y * y_log_y) / (y_sum * x2_y - x_y**2)
    a = math.exp(log_a)
    return (a, b)



def Normalize(p):
    """Normalizes a vector so that its entries sum one.

    On call, p should be a list of numbers.  The returned
        value is a new list with the entries normalized (multiplied
        by a scalar) so that they sum to one.
    """
    sum = 0.0
    for pi in p:
        sum += pi
    if sum == 0.0:
        raise ValueError, "Cannot normalize since entries sum to zero."
    return [pi / sum for pi in p]



def TTest(x1, x2):
    """Do two distributions have different means?

    Call is: '(t, p) = TTest(x1, x2)'
    This method uses Student's t test, allowing for unequal variances,
        to test if two distributions have different means.  It is 
        based on the 'tutest' method given by Numerical Recipes in
        C (chapter 14.2).
    'x1' and 'x2' are lists of numbers drawn from the two different
        distributions.
    Returns the 2-tuple '(t, p)' where 't' is the Student's t statistic,
        and p is the probability that the difference in means would be at least
        this large if the samples were drawn from distributions
        with the same mean."""
    (n1, n2) = (float(len(x1)), float(len(x2)))
    (x1avg, x2avg) = (Mean(x1), Mean(x2))
    (x1var, x2var) = (Variance(x1, True), Variance(x2, True))
    t = (x1avg - x2avg) / math.sqrt(x1var / n1 + x2var / n2)
    df = (x1var / n1 + x2var / n2)**2 / ((x1var / n1)**2 / (n1 - 1) + ((x2var / n2)**2 / (n2 - 1)))
    # evaluate the incomplete beta function
    p = transcendental.incbet(0.5 * df, 0.5, df / (df + t**2))
    return (t, p)


def GeometricStatistics(xlist):
    """Computes the geometric mean, geometric standard deviation factor, geometric standard error factor.

    Takes as input a list of N numbers.
    This method returns a 3-tuple (gm, gsd, gse) where:
    'gm' is the geometric mean, which is the product of all of the numbers taken to the Nth root.
    'gsd' is the standard deviation of the geometric mean which is a factor giving the spread
        of the numbers.  This is equal to exp(sqrt((sum_i(ln x_i - gm)**2)/N))
    'gse' is the standard error of the geometric mean, which is exp(ln(gsd / sqrt(N))).

    >>> xlist = [1.1, 0.9, 1.2]
    >>> (gm, gsd, gse) = GeometricStatistics(xlist)
    >>> print round(gm, 3)
    1.059
    >>> print round(gsd, 3)
    1.128
    >>> print round(gse, 3)
    1.072
    """
    n = len(xlist)
    prod = 1
    for x in xlist:
        assert x > 0, "Trying to compute the geometric mean with a number not > 0."
        prod *= x
    gm = prod**(1.0 / n)
    log_gm = math.log(gm)
    squarediffsum = 0.0
    for x in xlist:
        squarediffsum += (math.log(x) - log_gm)**2
    log_gsd = math.sqrt(squarediffsum / n)
    gsd = math.exp(log_gsd)
    gse = math.exp(log_gsd / math.sqrt(n))
    return (gm, gsd, gse)


if __name__ == "__main__": # run doctest
    import doctest
    doctest.testmod()

