"""This module contains amino acid hydrophobicity scales.

Written by Jesse Bloom, 2009."""



def KyteDoolittle(aa):
    """Returns the hydrophobicity of an amino acid based on the Kyte Doolittle scale.

    'aa' is the one letter code for an amino acid.
    The returned value is a number giving the hydrophobicity, as defined by Kyte
        and Doolittle in: 
            J. Kyte & R. F. Doolittle: 
            "A simple method for displaying the hydropathic character of a protein." 
            J Mol Biol, 157, 105-132
    More positive values indicate higher hydrophobicity, while more negative values
        indicate lower hydrophobicity.

    >>> print round(KyteDoolittle('F'), 2)
    2.8

    >>> print round(KyteDoolittle('N'), 2)
    -3.5

    >>> print round(KyteDoolittle('X'), 2)
    Traceback (most recent call last):
        ...
    ValueError: Invalid amino acid code of X.
    """
    aa = aa.upper()
    d = {'A':1.8, 'C':2.5, 'D':-3.5, 'E':-3.5, 'F':2.8, 'G':-0.4,
         'H':-3.2, 'I':4.5, 'K':-3.9, 'L':3.8, 'M':1.9, 'N':-3.5,
         'P':-1.6, 'Q':-3.5, 'R':-4.5, 'S':-0.8, 'T':-0.7, 
         'V':4.2, 'W':-0.9, 'Y':-1.3}
    try:
        return d[aa]
    except KeyError:
        raise ValueError("Invalid amino acid code of %s." % aa)


# Test with doctest
if __name__ == '__main__':
    import doctest
    doctest.testmod()
