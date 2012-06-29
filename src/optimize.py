"""Module for function optimization.

Written by Jesse Bloom, 2004-2009."""


import math


class OptimizationError(Exception):
    """Error in optimization routine."""


def ConjugateGradient(x0, f, f_and_df, itmax=1000, tolerance=1.0e-10):
    """Performs a conjugate gradient minimization of a function.

    Takes as input a function to minimize 'f' and a function that gives
        the derivative of 'f', 'df'.  These functions take as arguments
        a list of numbers of the same length as the list 'x0'. 
        returns a scalar which is the value of the function at 'x'. 
        'df(x)' returns a list which gives the gradient of 'f' and 'x'.
    'x0' is the list giving the initial guess for the parameters of the
        function be optimized.
    'f' is a function that takes as its single calling variable the
        list 'x' of size x0 and returns a number giving the value f(x).  This function
        seeks the value of x that minimizes f.
    'f_and_df' is a function that takes as its single calling variable the
        list 'x', and returns a 2-tuple of the form '(fx, dfx)'.  'fx' is the
        value of the function 'f' evaluated at 'x', while 'dfx' is a list of the 
        same length as 'x0'.  'dfx[i]' is the derivative of 'f' with respect to
        element 'i' of 'x'.
    'itmax' is the maximum number of iterations to perform.  If convergence
        is not achieved in this many iterations, an OptimizationError exception is raised.
    'tolerance' is the convergence tolerance for the minimization.

    The returned variable is the 3-tuple '(x, fmin, its)' where 'x' is
        the value of the calling variables that minimize 'f', 'fmin'
        is the value of the function at 'f', and 'its' is the number of iterations
        needed to minimize the functiion.

    The method is the Polak-Ribiere conjugate gradient formula,
        taken from Numerical Recipes, where it is given as the function
        'frprmn'.
        
    If this method fails to converge, either in the conjugate gradient search itself
        or in the line search, an OptimizationError exception is raised.
    """
    epsilon = 1.0e-10 # in case function converges exactly to zero
    # Error check on input variables
    assert isinstance(x0, list)
    assert isinstance(itmax, int) and itmax >= 1
    assert isinstance(tolerance, (int, float))
    # make 'x' a copy of 'x0'
    x = list(x0)
    n = len(x) # length of the vector
    # Compute initial values of function, derivative
    (fx, dfx) = f_and_df(x)
    assert isinstance(fx, (int, float))
    assert isinstance(dfx, list) and len(dfx) == n
    # Compute initial values for the lists 'g' and 'h', make 'dfx' negative
    g = []
    h = []
    for j in range(n): # loop over vector components
        dfx[j] = -dfx[j]
        g.append(dfx[j])
        h.append(dfx[j])
    # begin iteractions
    for its in range(itmax):
        (fmin, x) = LineMinimization(x, f, dfx) 
        # see if we have converged
        if (2.0 * math.fabs(fmin - fx)) <= tolerance * (math.fabs(fmin) + math.fabs(fx) + epsilon):
            return (x, fmin, its) # return
        # we have not converged, continue
        (fx, dfx) = f_and_df(x) # new value of f(x), df(x)
        gg = dgg = 0.0
        for j in range(n): # loop over vector components
            gg += g[j]**2
            dgg += (dfx[j] + g[j]) * dfx[j]
        if gg == 0.0: # gradient is exactly zero, return
            return (x, fmin, its) # return
        gam = dgg / gg
        for j in range(n):
            g[j] = -dfx[j]
            dfx[j] = h[j] = g[j] + gam * h[j]
    # if we are, we did too many iterations
    raise OptimizationError("Failed to converge in %r iterations." % itmax)



def LineMinimization(x0, f, dd, epsilon=1.0e-7, threshold=1.0e-8, maxiterations=100):
    """Performs a line search.

    'f' is a function that takes as an argument a list of the same length
        as the list 'x0'.
    'dd' is a vector of the same length as 'x0' specifying the direction of
        the line search.
    'fmin' is the value of 'f' at the minimum along the line, and
        the returned value 'x' is the point at which this minimum lies.
    'epsilon' is the initial step size for the line search.
    'threshold' is the convergence criterium for the line search.  Both 
        the 'x' points and the bracketing 'f' values must be within
        the threshold distance
    'maxiterations' is the maximum number of iterations that are performed
        before an OptimizationError exception is raised.

    Performs a line search beginning from 'x0' in the direction of 'dd'.
    
    Will raise an OptimizationError exception if it fails to converge.
    """
    assert isinstance(x0, list)
    # get initial points along line
    xn = list(x0)
    xnminus = list(x0)
    xnplus = list(x0)
    # Get initial values of f and check that 'f' is valid
    fn = fnminus = fnplus = f(xn) # value of 'f' at 'xn'
    assert isinstance(fn, (int, float))
    # propagate the line search by bounding the minimum between 'xn' and 'xnplus'
    for n in range(maxiterations):
        xnminus = list(xn)
        xn = list(xnplus)
        for i in range(len(xn)):
            xnplus[i] = xn[i] + dd[i] * epsilon
        epsilon *= 2.0 # double step size
        fnminus = fn
        fn = fnplus
        fnplus = f(xnplus)
        if fnplus > fn:
            break # terminate loop
    else:
        raise OptimizationError("Line search did not bound minimum in %r iterations." % maxiterations)
    # Now use the triples to find the minimum
    for n in range(maxiterations):
        xmidminus = [(xn[i] + xnminus[i]) / 2.0 for i in range(len(xn))]
        xmidplus = [(xn[i] + xnplus[i]) / 2.0 for i in range(len(xn))]
        fmidminus = f(xmidminus)
        fmidplus = f(xmidplus)
        if (fn < fmidminus) and (fn < fmidplus):
            # minimum between 'xmidminus' and 'xmidplus'
            xnminus = list(xmidminus)
            xnplus = list(xmidplus)
            fnminus = fmidminus
            fnplus = fmidplus
        elif (fmidminus <= fn):
            # minimum between 'xnminus' and 'xn'
            xnplus = list(xn)
            fnplus = fn
            xn = list(xmidminus)
            fn = fmidminus
        elif (fmidplus <= fn):
            # minimum between 'xn' and 'xnplus'
            xnminus = list(xn)
            fnminus = fn
            xn = list(xmidplus)
            fn = fmidplus
        else:
            raise ValueError("Minimum not bracketed: fnminus = %f, fmidminus = %f, fn = %f, fmidplus = %f, fnplus = %f" % (fnminus, fmidminus, fn, fmidplus, fnplus))
        # compute 'fdiff', the spread from 'fmidminus' to 'fmidplus'
        fdiff = math.fabs(fmidplus - fmidminus)
        if fdiff < threshold: # we may have converged
            xdist2 = 0.0
            for i in range(len(xn)):
                xdist2 += (xnplus[i] - xnminus[i])**2
            if xdist2 < threshold**2:
                return (fn, xn)
    else:
        raise OptimizationError("Line search did not find minimum in %r iterations." % maxiterations)
