"""Contains functions for analyzing Markov matrices and their stationary distributions.

It uses matrices in the numpy.ndarray format.

References for the material in this module:

"The analysis of panel data under a Markov assumption." J. D. Kalbfleisch and J. F. Lawless,
Journal of the American Statistical Association, 80:863-871 (1985).

"The role of the group generalized inverse in the theory of finite Markov Chains."
CD Meyer, SIAM Review, 17:443-464 (1975).

"Using the QR factorization and group inversion to compute, differentiate, and
estimate the sensitivity of stationary probabilities for Markov chains."
GH Golub and CD Meyer, SIAM J. Alg. Disc. Meth., 7:273-281 (1986).

Written by Jesse Bloom, 2009."""


import numpy


def DerivMatrixExponential(dG, alpha, S, Sinv, D):
    """Computes the derivative of an exponential of a Markov matrix.

    Computes dM/dy where M = exp(alpha * G) and y is a parameter on which the entries of G depend.
    'G' is defined by G = P - I where P is the one-step left stochastic matrix
        of a Markov process.  So each column of G should sum to unity.  Note that
        'G' is not actually a calling variable to this method, although it can
        in principle be reconstituted from S, Sinv, and D.
    'dG' is the element-by-element derivative of of 'G' with respect to
        some parameter y.
    'alpha' is a number by which 'G' is multiplied in the exponential.
    'S' is the matrix with columns equal to the right eigenvectors of G,
        'Sinv' is the inverse of 'S', and 'D' is a one-dimensional
        array of the eigenvalues of 'G' in the same order the eigenvectors
        are listed in 'S'.  Therefore, G = numpy.dot(numpy.dot(S, numpy.diagflat(D)), Sinv)
    All matrices should be in numpy.ndarray format.
    The returned matrix dM is the specified derivative.
    Reference is: 
        "The analysis of panel data under a Markov assumption." J. D. Kalbfleisch and J. F. Lawless,
        Journal of the American Statistical Association, 80:863-871 (1985).

    >>> G = numpy.array([[-1.  ,  0.5 ,  0.5 ,  0.25],
    ...            [ 0.5 , -1.  ,  0.25,  0.25],
    ...            [ 0.5 ,  0.5 , -1.  ,  0.25],
    ...            [ 0.  ,  0.  ,  0.25, -0.75]])
    >>> dG = numpy.array([[-0.5 , 0. , 0. , 0.],
    ...             [ 0.5 , -1.0 , 0. , 0.],
    ...             [ 0. ,  0. , 0. , 0.],
    ...             [ 0. ,  1.0 , 0. , 0.]])
    >>> alpha = 1.7
    >>> (D, S) = numpy.linalg.eig(G)
    >>> Sinv = numpy.linalg.inv(S)
    >>> numpy.allclose(G, numpy.dot(numpy.dot(S, numpy.diagflat(D)), Sinv))
    True
    >>> dM = DerivMatrixExponential(dG, alpha, S, Sinv, D)
    >>> dy = 0.001
    >>> M = numpy.dot(numpy.dot(S, numpy.diagflat(numpy.exp(alpha * D))), Sinv)
    >>> Gdy = G + dy * dG
    >>> (Ddy, Sdy) = numpy.linalg.eig(Gdy)
    >>> Sinvdy = numpy.linalg.inv(Sdy)
    >>> Mdy = numpy.dot(numpy.dot(Sdy, numpy.diagflat(numpy.exp(alpha * Ddy))), Sinvdy)
    >>> not numpy.allclose(Mdy, M)
    True
    >>> numpy.allclose(Mdy, M + dy * dM) 
    True
    """
    (n1, n2) = dG.shape
    assert n1 == n2, "dG is not a square matrix."
    n = n1
    assert S.shape == (n, n), 'S does not have the correct dimensions.'
    assert Sinv.shape == (n, n), 'S does not have the correct dimensions.'
    assert D.shape == (n, ), 'D does not have the correct dimensions.'
    assert isinstance(alpha, (int, float)) or alpha.shape == ()
    B = numpy.dot(numpy.dot(Sinv, dG), S)
    expalphaD = numpy.exp(alpha * D)
    V = numpy.ndarray((n, n))
    for x in range(n):
        for y in range(n):
            if x != y:
                V[x, y] = B[x, y] * (expalphaD[x] - expalphaD[y]) / (D[x] - D[y])
            else:
                V[x, y] = B[x, x] * alpha * expalphaD[x]
    return numpy.dot(numpy.dot(S, V), Sinv)



def StationaryDistributionAndDerivative(P, dP, pi_Agi=None, left_stochastic=False):
    """Computes the stationary distribution and its derivative of a Markov Chain.

    Given a Markov chain with stationary distribution pi satisfying pi P = pi,
        this function calculates the derivative dpi of pi with respect to some
        dummy variable t, where dP gives the element-by-element derivatives of
        P with respect to this dummy variable.
    'P' is a square matrix encoded as a numpy.ndarray that gives the (right stochastic)
        one-step transition matrix of a Markov process, such that x_i P = x_(i + 1).
        This means that each row of 'P' should sum to one, and all entries should 
        be >= 1.  If 'pi_Agi' is specified, then 'P' is actually not needed,
        and so can be set to 'None'.
    'dP' is the element-by-element derivative of 'P' with respect to some dummy variable.
        So for example, if we just want the derivative of the stationary distribution
        with respect to the first element of 'P', then 'dP' is all zeros except for
        this element, which is one.
    'pi_Agi' is set to 'None' by default.  However, if the stationary distribution 'pi'
        and the group inverse 'Agi' of the Markov chain are already known, then they
        can be specified in this variable to reduce the computational time.  Specifically,
        'Agi' is the group inverse of the matrix A = I - P, where I is the identity
        matrix.  If 'pi' and 'Agi' are already known, then 'pi_Agi'
        should be set equal to the 2-tuple '(pi, Agi)'.  In general, we can compute
        pi_Agi = MarkovGroupInverse(I - P) where 'MarkovGroupInverse' is another
        function in this module.
    'left_stochastic' is a switch specifying that we are working with a left stochastic
        matrix, such that P x_i = x_(i + 1).  Unless this option is set to 'True',
        we assume that we instead have a right stochastic matrix.  If 'left_stochastic'
        is True, then if 'pi_Agi' is specified, it should also have been computed
        with 'left_stochastic' as True, so that 
        'pi_Agi = MarkovGroupInverse(I - P, left_stochastic=True)'
    The returned variable is the 2-tuple '(pi, dpi)' where 'pi' is the stationary
        distribution and 'dpi' is its derivative.  These are both numpy.ndarray
        one-dimensional arrays.

    The method follows the prescription given in the following reference (Theorem 3.2):
        "Using the QR factorization and group inversion to compute, differentiate,
        and estimate the sensitivity of stationary probabilities for Markov chains."
        GH Golub and CD Meyer, SIAM J. Alg. Disc. Meth., 7:273-281 (1986).

    >>> P = numpy.array([[ 0.  , 0.5 , 0.5 , 0.  ],
    ...                    [0.5 ,  0.  , 0.5 , 0.  ],
    ...                    [0.5 , 0.25,  0.  , 0.25],
    ...                    [0.25, 0.25, 0.25,  0.25]])
    >>> dP = numpy.array([[ -1.  , 1 , 0. , 0.  ],
    ...                    [0. ,  0.  , 0. , 0.  ],
    ...                    [0. , 0.,  0.  , 0.],
    ...                    [0., 0., 0.,  0.]])
    >>> (pi, dpi) = StationaryDistributionAndDerivative(P, dP)
    >>> numpy.allclose(pi, numpy.dot(pi, P))
    True
    >>> print numpy.array2string(pi, precision=8)
    [ 0.31578947  0.26315789  0.31578947  0.10526316]
    >>> dt = 0.001
    >>> (pidt, x) = StationaryDistributionAndDerivative(P + dt * dP, dP)
    >>> numpy.allclose(pidt, (pi + dpi * dt)) and not numpy.allclose(pi, pidt)
    True
    >>> (n, n) = P.shape
    >>> I = numpy.identity(n)
    >>> (pi2, dpi2) = StationaryDistributionAndDerivative(P, dP, pi_Agi=MarkovGroupInverse(I - P))
    >>> numpy.allclose(pi, pi2) and numpy.allclose(dpi, dpi2)
    True
    >>> (Pt, dPt) = (P.transpose(), dP.transpose())
    >>> (pi3, dpi3) = StationaryDistributionAndDerivative(Pt, dPt, left_stochastic=True)
    >>> numpy.allclose(pi, pi3) and numpy.allclose(dpi, dpi3)
    True
    >>> numpy.allclose(pi, numpy.dot(Pt, pi))
    True
    >>> pi_Agi = MarkovGroupInverse(I - Pt, left_stochastic=True)
    >>> (pi4, dpi4) = StationaryDistributionAndDerivative(Pt, dPt,
    ...      pi_Agi=pi_Agi, left_stochastic=True)
    >>> numpy.allclose(pi, pi4) and numpy.allclose(dpi, dpi4)
    True
    >>> (pi5, dpi5) = StationaryDistributionAndDerivative(None, dPt,
    ...      pi_Agi=pi_Agi, left_stochastic=True)
    >>> numpy.allclose(pi, pi5) and numpy.allclose(dpi, dpi5)
    True
    """
    tol = 1.0e-10 # tolerance for numerical error
    assert isinstance(dP, numpy.ndarray), "'dP' is not a numpy.ndarray."
    (n1, n2) = dP.shape
    if n1 != n2:
        raise ValueError("'dP' is not a square matrix, but has dimension %d X %d." % (n1, n2))
    n = n1 
    assert (P == None and pi_Agi != None) or (isinstance(P, numpy.ndarray) and P.shape == (n, n)), "'P' is not valid."
    if left_stochastic:
        if P != None:
            P = P.transpose()
        dP = dP.transpose()
    if __debug__ and P != None: # check to make sure rows sum to one
        for x in numpy.add.reduce(P.transpose()):
            if abs(1.0 - x) > tol:
                raise ValueError("'A' contains a row that does not come close to summing to one.")
    if pi_Agi:
        (pi, Agi) = pi_Agi
        assert isinstance(pi, numpy.ndarray) and len(pi) == n, "'pi' is not a valid stationary distribution."
        assert isinstance(Agi, numpy.ndarray) and Agi.shape == (n, n), "'Agi' is not a valid group inverse."
        if left_stochastic:
            Agi = Agi.transpose()
    else:
        (pi, Agi) = MarkovGroupInverse(numpy.identity(n) - P)
    dpi = numpy.dot(pi, numpy.dot(dP, Agi))
    return (pi, dpi)



def MarkovGroupInverse(A, left_stochastic=False):
    """Computes the stationary distribution and group inverse of a Markov Chain.

    'A' is a square matrix encoded as a numpy.ndarray that is related to the transition
        matrix of a finite and homogenous Markov chain.  Specifically, let 'T' be the
        right stochastic Markov one-step transition matrix, such that x_i T = x_(i + 1).  
        Then 'A' is defined by A = I - T, where I is the identity matrix.  This means
        that each row of 'A' should sum to zero.
    'left_stochastic' is a switch specifying that we are working with a left stochastic
        matrix, such that T x_i = x_(i + 1).  In this case, the calculation is performed
        so that the returned arrays correspond to the stationary distributions and group
        inverses of a left stochastic matrix.
    This function returns the 2-tuple '(w, Agi)'.  Here 'w' is a one-dimensional
        numpy.ndarray of length n (where 'A' is n X n) giving the stationary distribution
            of the Markov chain with transition matrix T.  Here 'Agi' is the group inverse
            of the matrix 'A'.

    The method follows the prescription given in section 5 of:
        "The role of the group generalized inverse in the theory of finite Markov Chains."
        CD Meyer, SIAM Review, 17:443-464 (1975).

    >>> A = numpy.array([[ 1.  , -0.5 , -0.5 , -0.  ],
    ...                    [-0.5 ,  1.  , -0.5 , -0.  ],
    ...                    [-0.5 , -0.25,  1.  , -0.25],
    ...                    [-0.25, -0.25, -0.25,  0.75]])
    >>> (w, Agi) = MarkovGroupInverse(A)
    >>> print numpy.array2string(w, precision=8)
    [ 0.31578947  0.26315789  0.31578947  0.10526316]
    >>> print numpy.array2string(Agi, precision=8)
    [[ 0.48938135 -0.11265005 -0.17728532 -0.19944598]
     [-0.17728532  0.55401662 -0.17728532 -0.19944598]
     [-0.21237304 -0.25300092  0.45429363  0.01108033]
     [-0.38781163 -0.28808864 -0.38781163  1.06371191]]
    >>> At = A.transpose()
    >>> (wt, Atgi) = MarkovGroupInverse(At, left_stochastic=True)
    >>> numpy.allclose(w, wt)
    True
    >>> (n, n) = A.shape
    >>> I = numpy.identity(n)
    >>> T = I - A
    >>> numpy.allclose(w, numpy.dot(w, T))
    True
    >>> numpy.allclose(w, numpy.dot(T.transpose(), w))
    True
    >>> numpy.allclose(-MarkovGroupInverse(A)[1], MarkovGroupInverse(-A)[1])
    True
    """
    tol = 1.0e-10 # tolerance for numerical error
    assert isinstance(A, numpy.ndarray), "'A' is not a numpy.ndarray."
    (n1, n2) = A.shape
    if n1 != n2:
        raise ValueError("'A' is not a square matrix, but has dimension %d X %d." % (n1, n2))
    n = n1 
    if left_stochastic:
        A = A.transpose()
    if __debug__: # check to make sure rows sum to zero
        for x in numpy.add.reduce(A.transpose()):
            if abs(x) > tol:
                raise ValueError("'A' contains a row that does not come close to summing to zero.")
    # U, c, d, and alpha are as defined in Theorem 5.2 of the reference
    U = A[0 : n - 1, 0 : n - 1]
    alpha = A[n - 1, n - 1]
    c = A[0 : n -1, n - 1]
    d = A[n - 1, 0 : n - 1]
    Uinv = numpy.linalg.inv(U) # compute the inverse of U, which Thereom 5.1 guarantees will exist
    j = numpy.ones(n - 1) # vector of ones
    # h, delta, beta, and F are computed as defined in Theorem 5.2
    h = numpy.dot(d, Uinv)
    beta = 1.0 - numpy.dot(h, j)
    # w is the stationary distribution of the Markov Chain.
    w = numpy.resize(-h, n) # this and the next two lines calculate w according to Theorem 5.3
    w[n - 1] = 1.0
    w /= beta
    # compute W according to Theorem 2.3
    W = numpy.resize(w, (n, n))
    AAgi = numpy.identity(n) - W # from definition of W given on page 444
    # compute Agi according to (5.5) under Theorem 5.4
    X = numpy.zeros((n, n)) # array shown in (5.5) after next line
    X[0 : n - 1, 0 : n - 1] = Uinv
    Agi = numpy.dot(AAgi, numpy.dot(X, AAgi))
    if left_stochastic:
        Agi = Agi.transpose()
    return (w, Agi)



if __name__ == "__main__": # run doctest
    import doctest
    doctest.testmod()
