"""Module for Bayesian inferences of ddG values from sequences.

This module performs tasks related to inferring the ddG values for single
mutations to proteins from sequence and phylogenetic information.

Underlying assumption in the code is that all ddG values are independent
and additive.

Written by Jesse Bloom, 2007-2009."""



import sys
import math
import random
import copy
import cPickle 
import shelve
import os
import copy
import numpy
from pips import tree
from pips import stats
from pips import clog_likelihood
from pips import cddg_inference
from pips import hydrophobicity
from pips import markovmatrix
from pips import optimize



def AminoAcids():
    """Returns a list of the one-letter codes of all amino acids.
    
    Importantly, these are ordered so that i = AminoAcidIndexDict()[AminoAcids()[i]]."""
    return ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']



def AminoAcidIndexDict():
    """Returns a dictionary that maps each amino acid to an integer code.

    The returned value is a dictionary with keys equal to all of the amino acid
        codes, and values being integers.  For example, 'AminoAcidIndexDict()["A"] = 0',
        'AminoAcidIndexDict()["C"] = 1', etc.  The integers range from 0 to 19,
        and are according to alphabetical order of the amino acids codes.
    Importantly, satisfies i = AminoAcidIndexDict()[AminoAcids()[i]].
    """
    return {'A':0, 'C':1, 'D':2, 'E':3, 'F':4, 'G':5, 'H':6, 'I':7, 'K':8, 'L':9, 'M':10, 'N':11, 'P':12, 'Q':13, 'R':14, 'S':15, 'T':16, 'V':17, 'W':18, 'Y':19} 



def NAminoAcids():
    """Returns the number of amino acids."""
    return 20



def RandomSequence(length):
    """Returns a random amino acid sequence of the specified length.

    On input, 'length' is an integer specifying the length of the sequence.
    At each position, each amino acid is equally likely to appear."""
    aas = AminoAcids()
    return ''.join([random.choice(aas) for i in range(length)])



class DDGSet(object):
    """Class for representing the set of ddG values for single protein mutations.
   
    The ddG values are assumed to be for free energies of folding.  That is, positive
        values represent destabilizing mutations.  The units are arbitrary.
    More detailed information is in the doc string for the __init__ method.
    """
    
    def __init__(self, protseq, treedata, c_set, initial_values, ddgs_prior, mutrate, gamma=0.8, beta=0.265, ddg_range=20, use_c_extensions=True, underflow=5, rounding=0.001, runtestcode=False, ddgs_prior_boundfrac=0.95):
        """Creates the object storing the ddG values for single mutations to a protein.

        'protseq' is a string specifying the protein sequence for which we are trying to
            estimate ddG values.  We use this initial sequence to specify the priors (the
            reference amino acids are those found in this sequence) in addition to determining
            the length of the protein for which are finding ddG values.

        'treedata' is a 'tree.Tree' object that holds the phylogenetic tree structure
            and also all of the protein sequences.  That is, a sequence must be 
            specified for each node of the tree.  These sequences must be composed
            entirely of one-letter amino acid codes or the - symbol for gaps.  These
            sequences must all be of the same length, which must equal 'ddgset.Length()'.

        'c_set' specifies the probabilities for transitions among the various amino acids in
            the absence of any selection due to stability (i.e. just taking into account
            mutations and the underlying structure of the genetic code.  It can be set to
            any value that is a valid calling argument to the 'CreateCSet' function.  Currently,
            these options are:
                * 'ALL_EQUAL' meaning that a nucleotide mutation to any codon is equally likely
                    to yield a new codon for any of the 20 amino acids.
                * 'CODON_EQUAL' meaning that each amino acid 'y' is considered to be equally
                    likely to occupy any of its constituent codons.  Each of the nine nucleotide
                    mutations to each of these codons are then considered to be equally likely
                    when computing the probability that a nucleotide mutation to the codon causes
                    a transition to some other amino acid 'x'.
                * '("TRANSITION_TRANSVERSION_RATIO", r)' is similar to 'CODON_EQUAL' in that each
                    amino acid is considered equally likely to occupy any of its codons, but
                    now nucleotide mutations are considered to have a transition transversion
                    bias of 'r'.  If 'r' is 0.5, then there is no bias.  In most real biological
                    data, transitions are more likely than transversions, so you may want to set
                    'r' to a larger value.  For example, for influenza 'r' is thought to be 
                    about 5.

        'initial_values' specifies the initial values for the ddG values.  It is a tuple of some
            sort.  If any ddG value is greater than 'ddg_range' it is set to 'ddg_range', and if
            it is less than '- ddg_range' they are set to '- ddg_range'.
            There are several forms for this tuple:
            * '("CONSTANT", ddg_value)' indicates that all of the mutations to 'protseq' have the
                same value.  All mutations to this reference sequence have a value of 'ddg_value'.
            * '("GAUSSIAN", ddg_mean, ddg_sd)' indicates that all of the mutations to  'protseq'
                follow a Gaussian distribution.  Mutations to 'protseq' have ddG values that
                are randomly drawn from a Gaussian distribution with mean 'ddg_mean' and standard
                deviation 'ddg_sd'.  
            * '("SPECIFIED", ddgs, mean_noise, sd_noise)' means that initial
                values for all of the ddG values are given.  There is the option of specifying
                some Gaussian noise on these values; if you don't want to specify any noise then
                just set the noise parameters to zero.  In this case, 'ddgs' should be given as 
                a dictionary that is keyed by by integers i representing the residue numbers
                (0 <= i < 'len(protseq)').  'ddgs[i]' is the 2-tuple '(wt, iddgs)'.  'wt' 
                should equal 'protseq[i]'.  'iddgs[mut]' is defined
                for 'mut' equal to all 20 amino acid codes except for 'wt'.  'iddgs[mut]' should
                give the ddG value for mutating residue 'i' from 'wt' to 'mut'.  

        'ddgs_prior' specifies the prior over the ddg values.  Essentially, this is used to specify
            a prior over each ddG individual value.  The overall prior for the DDGSet is the
            product of the priors over the individual amino acids.  In all cases, 'ddgs_prior'
            should be a tuple.  Note that the ddG values over which the priors are calculated are
            those for mutating 'protseq'.  'ddgs_prior' should be a 3-tuple of the following
            form: ("BETA", parameter_sum, centers).  This specifies the shapes and
            centers of the beta distributions used as the priors.  The prior value at
            each ddG value is the beta distribution evaluated at this value.  
            'parameter_sum' is a number > 2.  It is equal to the sum of the two shape parameters
                alpha and beta, which are themselves determined by the value of 'centers'.
                In general, smaller values of 'parameter_sum' (closer to 2) give
                broader distributions, while larger values of 'parameter_sum' give
                more sharply peaked distributions.
            'centers' specifies how we choose the values used to center the beta distribution.
                If 'centers' is a number, then all of the ddG values (for mutating 'protseq')
                are centered around this number.  If 'centers' is the string 'INITIAL' then the
                prior for each ddG value is centered around the initial ddG value that is
                set by the option 'initial_values'.  In general, you would probably want to
                use 'INITIAL' if you have specified specific ddG values based on some sort
                of biophysical measurement or computational prediction.  If 'centers'
                is the 3-tuple ('KYTE_DOOLITTLE_HYDROPHOBICITY', scale, shift), then the prior
                for a mutation is centered on the absolute value of the difference in
                hydrophobicities multiplied by the number 'scale', then added to the number
                'shift'.  The mode of the beta distribution is equal to the center for
                that mutation.  Note that any center that is outside or too close to the
                limits of the ddG range is adjusted up/down as determined by 
                ddgs_prior_boundfrac.

        'mutrate' specifies how many mutations to a codon occur before an actual amino acid
            substitution occurs, on average.  A good value might be something like 8.0.  It should
            always be a number >= one.

        'gamma' and 'beta' give parameters that help determine exactly how the ddG values
            influence the substitution probabilities.  Specifically, the substitution 
            probabilities are equal to 
                0.5 - 0.5 * tanh(beta * ddg - 0.5 * ln[gamma / (1 - gamma)]).
            'gamma' therefore gives the substitution probability for a mutation with 
            ddG = 0, 'beta' gives the sensitivity of the substitution probabilities
            to variation in the ddG values.  By default, gamma = 0.8 and beta = 0.265.
            If you change these parameters, you may want to change them in concert
            with each other and with 'ddg_range' to make sure that the effects
            of ddG on substitution probabilities remains sensible.  They should both
            always be greater than zero.

        'ddg_range' specifies the range of the ddG values.  They can have any integer value
            ddg such that -ddg_range <= ddg <= ddg_range.  So 'ddg_range' must be a
            number >= 1.  By default, 'ddg_range' is set to 20.  If you change it, you
            may want to change 'gamma' and 'beta' also.

        'use_c_extensions' specify that we perform certain calculations using C extension
            code.  This will in general be much faster.  By default, this is 'True'.

        'underflow' specifies how we correct to attempt to avoid numerical underflow.
            It is an integer argument > 0.  It specifies that every 'underflow' nodes,
            we scale the conditional probabilities up to that point in the tree
            by the maximum conditional probability, as described on pg 426 of
            Z Yang, J Mol Evol, 51:423-432 (2000).  In general, smaller values of
            underflow will lead to more frequent scaling to avoid underflow.  A value
            of 1 means that scaling is performed at every node.  A value greater than
            twice the number of nodes in the tree means that no scaling is ever
            performed.  By default, 'underflow' is set to 5, which should be more
            than sufficient to prevent underflow.

        'rounding' specifies whether we round branch lengths in the tree.  Such rounding
            can substantially accelerate the calculation of the mutation probabilities
            if many branches are rounded to have
            the same length.  If 'rounding' is set to None, then no rounding is performed.
            Otherwise, 'rounding' should be set to a number > 0.  Each branch length is
            then rounded to the nearest multiple of 'rounding'.  For example, if 'rounding'
            is set to 0.01 (a reasonable value), branches of lengths 0.009 and 0.011 would
            both be rounded to have the same branch length of 0.01.  One concern is that
            due to rounding, branch lengths could be rounded down to zero.  Therefore,
            any branch length that would be zero after rounding (including any that
            are already zero) are instead set to be equal to 'rounding / 2.0'.

        'runtestcode' is a Boolean switch that specifies that we run testing on the
            DDGSet object to make sure all calculations are working correctly.
            This is quite time consuming, will alter the initial ddG values, and is
            probably not something that you want to do all the time.  However, it is
            a good check after installing the progrm on a new operating system or
            after modifying the source code.

        'ddgs_prior_boundfrac' is a number specifying how close the center of the 
            prior probability distributions over the ddG values can be to the bounds
            of the allowable ddG values, as specified by 'ddg_range'.  The centers
            of all prior probability distributions is constrained to satisfy:
            -ddg_range * ddgs_prior_boundfrac <= center < ddg_range * ddgs_prior_boundfrac
            Therefore, we must always have 0 < ddgs_prior_boundfrac < 1.
            Any center of a prior that does not meet this range is adjusted up/down to
            be at the limit of this range.
            This constraint is necessary to ensure that the beta distribution is always peaked
            internal to the ddg_range.
        """
        self._n_aa = NAminoAcids()
        self._n_aa2 = self._n_aa**2
        self._amino_acids = AminoAcids()
        self._aa_index = AminoAcidIndexDict()
        self._reverse_aa_index = dict([(ires, res) for (res, ires) in self._aa_index.iteritems()])
        self._use_c_extensions = use_c_extensions
        assert isinstance(self._use_c_extensions, bool)
        self._mutrate = mutrate
        assert isinstance(self._mutrate, (int, float)) and self._mutrate >= 0.
        assert isinstance(protseq, str) and len(protseq) >= 1
        if __debug__:
            for r in protseq:
                if r not in self._amino_acids:
                    raise ValueError, "Invalid amino acid of %s." % r
        self._protseq = protseq
        self._length = len(protseq)
        assert isinstance(beta, (int, float)) and beta > 0
        self._beta = float(beta)
        assert isinstance(gamma, (int, float)) and gamma > 0
        self._gamma = float(gamma)
        assert isinstance(ddg_range, (int, float)) and ddg_range > 0
        self._ddg_range = float(ddg_range)
        # store _c_set, giving the codon transition probabilities
        c_set = CreateCSet(c_set)
        assert isinstance(c_set, dict)
        if __debug__: # do some error checking on c_set
            for aa_2 in self._amino_acids:
                sum = 0.0
                for aa_1 in self._amino_acids:
                    c = c_set[(aa_1, aa_2)]
                    assert isinstance(c, (int, float)) and 0 <= c <= 1
                    sum += c
                assert 0.9999 < sum < 1.0001
        self._c_set = c_set
        # Now we create the list _indexed_c_set.  This list is of length self._n_aa2.  The element
        # yi * n_aa + xi contains the same value as self._c_set[(x, y)] where xi = AminoAcidIndexDict()[x].
        self._indexed_c_set = [0.0] * self._n_aa2
        index = 0
        for y in self._amino_acids:
            for x in self._amino_acids:
                self._indexed_c_set[index] = float(self._c_set[(x, y)])
                index += 1
        #
        # The ddG values are stored in self._ddgs, which is a list of length self._length.  Element
        # ddgs[ires][aa] (0 <= ires < self._length) gives the ddG value for mutating residue ires from
        # protseq[ires] to aa.  All entries are floats.
        assert isinstance(initial_values, tuple) and len(initial_values) > 0
        self._ddgs = []
        for ires in range(self._length):
            wt = protseq[ires]
            iddgs = {}
            if initial_values[0] == 'CONSTANT':
                (initialize_type, ddg_value) = initial_values
                for aa in self._amino_acids:
                    if aa != wt:
                        iddgs[aa] = ddg_value
            elif initial_values[0] == 'GAUSSIAN':
                (initialize_type, ddg_mean, ddg_sd) = initial_values
                for aa in self._amino_acids:
                    if aa != wt:
                        iddgs[aa] = random.gauss(ddg_mean, ddg_sd)
            elif initial_values[0] == 'SPECIFIED':
                (initialize_type, ddgs, mean_noise, sd_noise) = initial_values
                assert len(ddgs) == self._length and sd_noise >= 0
                iddgs = dict(ddgs[ires][1])
                if ddgs[ires][0] != wt:
                    raise ValueError("Invalid 'wt' residue.")
                if mean_noise or sd_noise: # add gaussian noise
                    for mut in self._amino_acids:
                        if mut != wt:
                            iddgs[mut] += random.gauss(mean_noise, sd_noise)
            else:
                raise ValueError, "Invalid initial values of %s." % str(initial_values)
            assert len(iddgs) == 19
            self._ddgs.append(iddgs)
        # Now make sure all values are floats within ddg_range
        for ires in range(self._length):
            wt = protseq[ires]
            iresdict = self._ddgs[ires]
            for aa in self._amino_acids:
                if aa != wt:
                    ddg = float(iresdict[aa])
                    ddg = min(ddg, self._ddg_range - 1e-10)
                    ddg = max(ddg, -self._ddg_range + 1e-10)
                    iresdict[aa] = ddg
        #
        # Information about the priors is stored in the variable self._priors.
        # This is a dictionary keyed by residue number with entries that are 
        # dictionaries in turn keyed by residue for all non-wildtype residues
        # at that position.  Each entry is in turn the 2-tuple '(alpha, beta)'
        # giving the beta distribution alpha and beta parameters for that prior.
        assert 0 < ddgs_prior_boundfrac < 1
        assert isinstance(ddgs_prior, tuple), str(ddgs_prior)
        self._priors = {}
        (prior_type, parameter_sum, centers) = ddgs_prior
        if prior_type != 'BETA':
            raise ValueError('Invalid ddgs_prior.')
        assert parameter_sum > 2.0
        if isinstance(centers, (int, float)):
            if not (-ddg_range < centers < ddg_range):
                raise ValueError("Center for prior is not inside ddg range.")
            newcenter = min(ddg_range * ddgs_prior_boundfrac, centers)
            newcenter = max(-ddg_range * ddgs_prior_boundfrac, newcenter)
            if newcenter != centers:
                print "WARNING: prior center of %f has been adjusted to %f." % (centers, newcenter)
                sys.stderr.write("WARNING: prior center of %f has been adjusted to %f." % (centers, newcenter))
            (alpha, beta) = stats.ParameterizeBetaDistribution(newcenter, 'MODE', parameter_sum, -ddg_range, ddg_range)
            for ires in range(self._length):
                wt = protseq[ires]
                self._priors[ires] = dict([(aa, (alpha, beta)) for aa in self._amino_acids if aa != wt])
        elif centers == 'INITIAL':
            for ires in range(self._length):
                wt = protseq[ires]
                iresdict = {}
                for aa in self._amino_acids:
                    if aa == wt:
                        continue
                    ddg = self._ddgs[ires][aa]
                    newcenter = min(ddg_range * ddgs_prior_boundfrac, ddg)
                    newcenter = max(-ddg_range * ddgs_prior_boundfrac, newcenter)
                    if newcenter != ddg:
                        print "WARNING: prior center of %f has been adjusted to %f." % (ddg, newcenter)
                    (alpha, beta) = stats.ParameterizeBetaDistribution(newcenter, 'MODE', parameter_sum, -ddg_range, ddg_range)
                    iresdict[aa] = (alpha, beta)
                self._priors[ires] = iresdict                        
        elif isinstance(centers, tuple) and len(centers) == 3 and centers[0] == 'KYTE_DOOLITTLE_HYDROPHOBICITY': 
            (scale, shift) = (centers[1], centers[2])
            assert isinstance(scale, (int, float)) and isinstance(shift, (int, float))
            for ires in range(self._length):
                wt = protseq[ires]
                wt_kdh = hydrophobicity.KyteDoolittle(wt)
                iresdict = {}
                for aa in self._amino_acids:
                    if aa == wt:
                        continue
                    aa_kdh = hydrophobicity.KyteDoolittle(aa)
                    center = abs(wt_kdh - aa_kdh) * scale + shift
                    newcenter = min(ddg_range * ddgs_prior_boundfrac, center)
                    newcenter = max(-ddg_range * ddgs_prior_boundfrac, newcenter)
                    if newcenter != center:
                        print "WARNING: prior center of %f has been adjusted to %f." % (center, newcenter)
                    (alpha, beta) = stats.ParameterizeBetaDistribution(newcenter, 'MODE', parameter_sum, -ddg_range, ddg_range)
                    iresdict[aa] = (alpha, beta)
                self._priors[ires] = iresdict
        else:
            raise ValueError("Invalid value of 'centers': %s" % (str(centers)))
        # create variables for evaluating the log likelihood
        # DEFINITIONS OF PRIVATE OBJECT VARIABLES FOR LOG LIKELIHOOD CALCULATION (all preceded by self._)
        # rounding is the value of the calling parameter of the same name.
        # underflow is the value of the calling parameter of the same name.
        # n_internal is the integer number of internal nodes in the tree (usual index is i).
        # n_tips is the integer number of tips in the tree (usual index is t).
        # n_aa is the integer number of amino acids (usual index is x or y).  Each amino
        #    acid type is encoded by an integer x, 0 <= x < n_aa.  The mapping from amino
        #    acid code to integer is given by AminoAcidIndexDict().
        # tips is an integer list of length * n_tips, with tips[r * n_tips + t] (0 <= t < n_tips,
        #    0 <= r < length) giving the identity of amino acid r at tip t as x, where 0 <= x < n_aa.
        #    If an amino acid is a gap, it is encoded as -1.
        # descendents and descendent_is_tip are both integer lists of length 2 * n_internal. 
        #    They give information about the descendents of internal node i, where 0 <= i < n_aa.
        #    If descendent_is_tip[2 * i] is 1, then the right descendent of i is a tip node.  In
        #    this case, the identity of the amino acid at residue r of this descendent is given by
        #    element r * n_tips + descendents[2 * i] of tips (so 0 <= descendents[2 * i] < n_tips).
        #    If descendent_is_tip[2 * i] is 0, then the right descendent of i is not a tip node.  In
        #    this case, descendents[2 * i] gives the number j of the internal node (0 <= j < n_internal)
        #    that is the right descendent of node j.  Crucially, the nodes are numbered so that 
        #    j will always be less than i.  The information about the left descendent is stored
        #    similarly, but in element 2 * i + 1 of descendents and descendent_is_tip.  Because of
        #    the ordering where j < i, this guarantees that the last entry (i = n_aa - 1) corresponds
        #    to the root node of the tree.
        # p is a double (float) list of length n_internal * n_aa.  It is recomputed for each
        #    residue each time the function is called.  p[i * n_aa + y] is the cumulative 
        #    probability that internal node i is residue y given the information in the subtree
        #    below i.  Because of the ordering of the internal node numbering, the entries for
        #    i = n_aa - 1 correspond to the root node of the tree.
        # tips_ut is a double (float) list of length n_tips.  tips_ut[t] (0 <= t < n_tips) is
        #    the length of the branch leading to tip t.  Note that these values are rounded according
        #    to the value of rounding.
        # internals_ut is a double (float) list of length n_internal.  internals_ut[i]
        #    (0 <= i < n_internal) is the length of the branch leading to internal node i.
        # unique_uts is a double (float) list of all of the unique ut values stored in tips_ut and internals_ut.
        #    So it will be of a size of <= n_tips * n_internal, with the size being < if there
        #    are any duplicated uts.
        # n_unique is the size of unique_uts, and so gives the number of unique branch lengths.
        # unique_uts_index is a list of integers, with the list size being n_tips + n_internal.
        #    For any tip node t (0 <= t < n_tips), unique_uts_index[t] gives the index in unique_uts
        #    corresponding to the entry with value ut equal to the branch length for this tip node.  
        #    For any internal node i (0 <= i < n_internal), unique_uts_index[n_tips + i] gives the 
        #    the index in unique_uts corresponding to the entry value value ut equal to the branch
        #    length for this internal node.
        # unique_uts_index_n_aa2 is just like unique_uts_index, except that every entry is multiplied by
        #    n_aa**2.
        # mlist_only_need is a list of integers.  It is used when calling the DDGSet._MList
        #    function, and specifies that we not need to compute certain values in mlist.
        #    It is of size length * n_unique.  For any tip node t with a unique ut value (branch
        #    length) that is not shared by any other node, mlist_only_need[r * n_unique + unique_uts_index[t]],
        #    is the value of tips[r * n_tips + t], giving the value of residue r at tip t.  All other entries
        #    in mlist_only_need are equal to -2, indicating that we need to consider all possible amino acids.
        if not (rounding == None or (isinstance(rounding, (int, float)) and rounding > 0)):
            raise ValueError, "Invalid value of 'rounding'."
        self._rounding = rounding
        if not (isinstance(underflow, int) and underflow >= 0):
            raise ValueError, "Invalid value of 'underflow'."
        self._underflow = underflow
        assert isinstance(treedata, tree.Tree)
        self._tree = copy.deepcopy(treedata) # make a copy of the tree since we are renumbering it
        self._root = self._tree.GetRoot()
        # This next method numbers the tip and internal nodes in a way that satisfies the constraints
        # on internal node numbers.
        (self._n_tips, self._n_internal) = tree.RecursivelySetTipAndInternalNumbers(self._root, 0, 0)
        assert self._root.number == self._n_internal - 1
        tip_list = []
        internal_list = []
        tree.ListsOfNodes(self._root, tip_list, internal_list)
        assert self._n_tips == len(tip_list) and self._n_internal == len(internal_list)
        for tip in tip_list:
            if len(tip.seq) != self._length:
                raise ValueError("Tip sequence is of length %d, not %d." % (len(tip.seq), self._length))
        # Begin constructing the private object variables
        self._tips = [0] * self._length * self._n_tips
        for tip in tip_list: # encode sequence information in self._tips
            for r in range(self._length):
                try:
                    self._tips[r * self._n_tips + tip.number] = self._aa_index[tip.seq[r]]
                except KeyError:
                    if tip.seq[r] == '-':
                        self._tips[r * self._n_tips + tip.number] = -1
                    else:
                        raise
        self._descendents = [0] * 2 * self._n_internal
        self._descendent_is_tip = [0] * 2 * self._n_internal
        for internal in internal_list: # encode information about internal node descendents
            self._descendents[2 * internal.number] = internal.rightdescendent.number
            self._descendents[2 * internal.number + 1] = internal.leftdescendent.number
            if internal.rightdescendent.tip:
                self._descendent_is_tip[2 * internal.number] = 1
            if internal.leftdescendent.tip:
                self._descendent_is_tip[2 * internal.number + 1] = 1 
        # perform branch length computations
        self._tips_ut = [1.0] * self._n_tips
        self._internals_ut = [1.0] * self._n_internal
        if rounding == None:
            # store exact branch lengths
            for tip in tip_list: # encode information about branch lengths is self._tips_ut
                assert tip.ancestorbranch != None
                self._tips_ut[tip.number] = tip.ancestorbranch
            for internal in internal_list: # encode information about branch lengths in self._internals_ut
                assert (internal.ancestorbranch != None) or internal == self._root
                if internal == self._root:
                    self._internals_ut[internal.number] = 1.0 # value for root branch is irrelevant, but needs to be made a number for subsequent code to function
                else:
                    self._internals_ut[internal.number] = internal.ancestorbranch
        else:
            # store rounded branch lengths
            # the rounded value of ut is: round(ut / rounding) * rounding
            for tip in tip_list: # encode information about branch lengths is self._tips_ut
                assert tip.ancestorbranch != None
                self._tips_ut[tip.number] = max(rounding / 2.0, round(tip.ancestorbranch / rounding) * rounding)
            for internal in internal_list: # encode information about branch lengths in self._internals_ut
                assert (internal.ancestorbranch != None) or internal == self._root
                if internal == self._root:
                    self._internals_ut[internal.number] = 1.0 # value for root branch is irrelevant, but needs to be made a number for subsequent code to function
                else:
                    self._internals_ut[internal.number] = max(rounding / 2.0, round(internal.ancestorbranch / rounding) * rounding)
        self._unique_uts = []
        self._unique_uts_index = []
        for t in range(self._n_tips):
            ut = self._tips_ut[t]
            try:
                # there is already a value for this ut stored
                unique_index = self._unique_uts.index(ut)
                self._unique_uts_index.append(unique_index)
            except ValueError:
                # this ut is unique, no value yet stored
                self._unique_uts_index.append(len(self._unique_uts))
                self._unique_uts.append(ut)
        for i in range(self._n_internal):
            ut = self._internals_ut[i]
            try:
                # there is already a value for this ut stored
                unique_index = self._unique_uts.index(ut)
                self._unique_uts_index.append(unique_index)
            except ValueError:
                # this ut is unique, no value yet stored
                self._unique_uts_index.append(len(self._unique_uts))
                self._unique_uts.append(ut)
        self._n_unique = len(self._unique_uts)
        unique_tip_branch_length = [True] * self._n_tips # entry t is true iff t's branch length is unique
        for t in range(self._n_tips):
            ut = self._tips_ut[t]
            unique_index = self._unique_uts.index(ut)
            if self._unique_uts_index.count(unique_index) > 1:
                unique_tip_branch_length[t] = False
        n_aa2 = self._n_aa * self._n_aa
        self._unique_uts_index_n_aa2 = [index * n_aa2 for index in self._unique_uts_index]
        # fill self._mlist_only_need
        self._mlist_only_need = [-2] * (self._length * self._n_unique)
        for r in range(self._length):
            r_n_unique = r * self._n_unique
            r_n_tips = r * self._n_tips
            for t in range(self._n_tips):
                if unique_tip_branch_length[t]: # branch length is unique
                    self._mlist_only_need[r_n_unique + self._unique_uts_index[t]] = self._tips[r_n_tips + t]
        if self._use_c_extensions: # perform the initialization needed for the C extension 
            clog_likelihood.Initialize(self._length, self._n_aa, self._n_tips, self._n_internal, self._tips, self._descendents, self._descendent_is_tip, self._n_unique, self._unique_uts_index_n_aa2, self._underflow)
        self._unique_uts = [ut * self._mutrate for ut in self._unique_uts] # scale by mutrate
        #
        # Now compute internal properties:
        #
        # self._indexed_fs[r] is a list of length n_aa * n_aa.  Element yi * n_aa + xi contains  
        # the fxyr value for mutating residue r (0 <= r < length) from yi to xi, where
        # yi = aa_index[y] is the integer code for residue y, and xi is the same for residue x.  
        # Elements corresponding to xi = yi are zero.
        #
        # For each residue r, self._grs[r] is a 3-tuple of (diag_Gr, P, and P_inv) such that 
        # diag_Gr is diagonal elements and Gr = P numpy.diagflat(diag_Gr), P_inv.  Entries 
        # are 'None' if they are not yet computed
        #
        # For each residue r, self._group_inverses[r] is a numpy.ndarray of n_aa X n_aa which
        # is the group inverse of the negative of the Gr matrix.
        #
        # For each residue r, self._dgrs[r] is a list of length n_aa.  Element zi gives the
        # derivative of the Gr matrix for residue r with respect to the ddG for mutating
        # that residue from the wildtype residue at that position to residue z, where
        # zi = AminoAcidIndexDict()[z].  If z is equal to the wildtype residue, self._dgrs[r][zi]
        # is None.
        #
        # self._pi[n_aa * r + yi] gives the equilibrium probability that
        # residue r (0 <= r < self._length) is residue y where yi = AminoAcidIndexDict()[y] and
        # n_aa = NAminoAcids(). 
        #
        # self._dpi[n_aa2 * r + zi * n_aa + yi] gives the derivative of self._pi[n_aa * r + yi] with
        # respect to the ddG for mutating residue r from the wildtype residue to residue z.
        # Entries for zi equal to the wildtype residue are undefined.
        #
        # self._residue_log_likelihoods[r] is the log likelihood for residue r
        #
        # self._d_residue_log_likelihoods[r][zi] is the derivative of the log likelihood 
        # for residue r with respect to the ddG of mutating from the wildtype residue to residue
        # z, with entries for z equal to the wildtype residue undefined, 
        # where zi = AminoAcidIndexDict()[z]
        #
        # self._residue_log_priorprob[r] is the log of the prior probability of residue r
        # having the current ddG value.
        #
        # self._d_residue_log_priorprob[r][zi] is the derivative of self._residue_log_priorprob[r]
        # with respect to the ddG for mutating from the wildtype residue to residue z, with entries
        # for z equal to the wildtype residue undefined, where zi = AminoAcidIndexDict()[z]
        #
        # self._residue_log_posterior[r] is the log of the posterior probability of residue r
        # having the current ddG value.
        #
        # self._d_residue_log_posterior[r][zi] is the derivative of self._residue_log_posterior[r]
        # with respect to the ddG for mutating from the wildtype residue to residue z, with entries
        # for z equal to the wildtype residue undefined, where zi = AminoAcidIndexDict()[z]
        #
        # self._values_current[r] is True if self._indexed_fs, self._grs, self._group_inverses,
        # self._pi, self_residue_log_likelihoods, and self._resiude_log_priorprob are current
        # for residue r, and False otherwise.
        #
        # self._derivatives_current[r] is True if self._dpi, self._dgrs, 
        # self._d_residue_log_likelihoods, and self._d_residue_log_priorprob are current 
        # for residue 'r', and False otherwise.  If self._derivatives_current[r] is True,
        # then self._values_current[r] is also True (the converse statement is not true).
        self._indexed_fs = [None] * self._length
        self._grs = [None] * self._length
        self._group_inverses = [None] * self._length
        self._dgrs = []
        for r in range(self._length):
            self._dgrs.append(self._n_aa * [None])
        self._pi = [None] * self._length * self._n_aa
        self._dpi = [None] * self._length * self._n_aa2 
        self._residue_log_likelihoods = [None] * self._length
        self._d_residue_log_likelihoods = [None] * self._length
        self._residue_log_priorprob = [None] * self._length
        self._d_residue_log_priorprob = []
        for r in range(self._length):
            self._d_residue_log_priorprob.append(self._n_aa * [None])
        self._residue_log_posterior = [None] * self._length
        self._d_residue_log_posterior = []
        for r in range(self._length):
            self._d_residue_log_posterior.append(self._n_aa * [None])
        self._values_current = [False] * self._length
        self._derivatives_current = [False] * self._length
        for r in range(self._length):
            self._ComputeInternals(r, compute_derivatives=True) # computes these internal properties
        # test the code
        if runtestcode:
            self._Test(print_output=True)

    def LogProbabilities(self, r):
        """Returns the log posterior probability, log likelihood, and log prior probability.

        The returned log probabilities are for the ddG values currently held in the DDGSet.
        'r' specifies the residue for which we want the probabilities.  If it is an integer
            giving a residue number (0 <= r < DDGSet.Length()) then we return the log
            probabilities just for the ddG values associated with that residue.
            If 'r' is instead equal to the string 'ALL' then we return the overall log
            probabilities for all of the residues (the sums of the log probabilities for
            the individual residues).

        The returned variable is the tuple '(logposterior, loglikelihood, logprior)'

        If the relevant probabilities are already computed, they are not re-computed.
            However, if they still need to be computed, this method will recompute them.
        """
        if isinstance(r, int):
            assert 0 <= r < self._length
            self._ComputeInternals(r, compute_derivatives=False)
            return (self._residue_log_posterior[r], self._residue_log_likelihoods[r], self._residue_log_priorprob[r])
        elif r == 'ALL':
            logposterior = loglikelihood = logprior = 0
            for r in range(self._length):
                self._ComputeInternals(r, compute_derivatives=False)
                logposterior += self._residue_log_posterior[r]
                loglikelihood += self._residue_log_likelihoods[r]
                logprior += self._residue_log_priorprob[r]
            return (logposterior, loglikelihood, logprior)
        else:
            raise ValueError("Invalid value of r: %r" % r)

    def _Test(self, d_ddg=0.002, print_output=False):
        """Performs some testing of internal class methods.

        The computation of derivatives is tested numerically using 'd_ddg'.
            The derivatives must accurately predict (using numpy.allclose) the
            effect of changing a ddG value by d_ddg.  If d_ddg is too small there
            is no sensitivity since the original and new values will be equal
            (this causes an exception to be raised).  Conversely, if it is too
            large then the derivatives will not accurately predict the changes
            even if they are calculated accurately if the second derivative is
            large.  The current value of d_ddg appears to be a reasonable compromise
            between these competing concerns for testing purposes.
        'print_output' specifies that we print output summaries to standard output."""
        assert isinstance(print_output, bool)
        assert isinstance(d_ddg, (int, float))    
        random.seed(1) # seed the random number generator so that this test is consistent
        n_aa = self._n_aa
        for r in range(self._length): # mutate ddG values for each residue
            nequal_m = ntotal_m = 0 # makes sure that most of the times, new and old M matrices are not equal
            ntotal_ll = nflagged_ll = 0 # make sure that most of the log likelihoods pass the allclose tests
            for z in self._amino_acids:
                if z == self._protseq[r]:
                    continue # this is reference amino acid
                zi = self._aa_index[z]
                if print_output:
                    print "Testing for changes in ddG for mutating residue %d to %s..." % (r, z)
                old_ddg = self._ddgs[r][z]
                (old_gr_diag, old_p, old_p_inv) = self._grs[r]
                if self._use_c_extensions:
                    old_gr = numpy.dot(old_p.reshape((n_aa, n_aa)), numpy.dot(numpy.diagflat(old_gr_diag), old_p_inv.reshape((n_aa, n_aa))))
                else:
                    old_gr = numpy.dot(old_p, numpy.dot(numpy.diagflat(old_gr_diag), old_p_inv))
                old_dgr = self._dgrs[r][zi]
                if self._use_c_extensions:
                    old_dgr = old_dgr.reshape((n_aa, n_aa))
                old_residue_log_likelihood = self._residue_log_likelihoods[r]
                old_d_residue_log_likelihoods = copy.copy(self._d_residue_log_likelihoods[r])
                old_mlist = self._MList(r)
                old_dmlist = self._dMList(r)
                old_pi = self._pi[self._n_aa * r : self._n_aa * (r + 1)]
                old_dpi = self._dpi[self._n_aa2 * r + self._n_aa * zi : self._n_aa2 * r + self._n_aa * (zi + 1)]
                old_residue_log_priorprob = self._residue_log_priorprob[r]
                old_d_residue_log_priorprob = copy.copy(self._d_residue_log_priorprob[r])
                self._ddgs[r][z] += d_ddg
                self._values_current[r] = False
                self._derivatives_current[r] = False
                old_residue_log_posterior = self._residue_log_posterior[r]
                old_d_residue_log_posterior = copy.copy(self._d_residue_log_posterior[r])
                # first test that we compute the same values with and without compute_derivatives
                self._ComputeInternals(r, compute_derivatives=False)
                noderiv_residue_log_likelihood = self._residue_log_likelihoods[r]
                self._ComputeInternals(r, compute_derivatives=True)
                if noderiv_residue_log_likelihood != self._residue_log_likelihoods[r]:
                    raise ValueError("Different log likelihoods are computed depending on whether or not we are computing the derivatives!")
                #
                # First test that we are computing the derivatives of Gr correctly:
                (gr_diag, p, p_inv) = self._grs[r]
                if self._use_c_extensions:
                    gr = numpy.dot(p.reshape((n_aa, n_aa)), numpy.dot(numpy.diagflat(gr_diag), p_inv.reshape((n_aa, n_aa))))
                else:
                    gr = numpy.dot(p, numpy.dot(numpy.diagflat(gr_diag), p_inv))
                if numpy.allclose(gr, old_gr):
                    raise ValueError("New and old Gr matrices are nearly equal.")
                if not numpy.allclose(gr, old_gr + d_ddg * old_dgr):
                    print "new Gr - old Gr ", gr - old_gr
                    print "new Gr - (old Gr + d_ddG * old dGr)", gr - (old_gr + d_ddg * old_dgr)
                    raise ValueError("New Gr matrix is not nearly equal to old one plus derivative.")
                #
                # Now test that we are computing the derivatives of the M values correctly
                mlist = self._MList(r)
                # mlist entries are indexed by ut_index * n_aa2 + xi * n_aa + yi
                # dmlist entries are indexed by ut_index * n_aa3 + zi * n_aa2 + xi * n_aa + yi
                n_aa = self._n_aa
                n_aa2 = self._n_aa2
                n_aa3 = n_aa * n_aa2
                for iut in range(len(self._unique_uts)):
                    xi = self._mlist_only_need[r * self._n_unique + iut]
                    if xi == -1: 
                        continue # transition probabilities not defined here
                    elif xi == -2:
                        # consider all possible transitions
                        index = iut * n_aa2
                        index2 = iut * n_aa3 + zi * n_aa2
                        old_m_ut = old_mlist[index : index + n_aa2]
                        old_dm_ut = old_dmlist[index2 : index2 + n_aa2]
                        m_ut = mlist[index : index + n_aa2]
                    else:
                        # only transitions to xi are considered
                        index = iut * n_aa2 + xi * n_aa
                        index2 = iut * n_aa3 + zi * n_aa2 + xi * n_aa
                        old_m_ut = old_mlist[index : index + n_aa]
                        old_dm_ut = old_dmlist[index2 : index2 + n_aa]
                        m_ut = mlist[index : index + n_aa]
                    ntotal_m += 1
                    if numpy.allclose(m_ut, old_m_ut):
                        nequal_m += 1
                    if not numpy.allclose(m_ut, old_m_ut + d_ddg * old_dm_ut):
                        raise ValueError("New M matrix is not nearly equal to the old one plus derivative for ut = %f." % self._unique_uts[iut])
                #
                # Now test that we are computing the derivatives of the pi values correctly
                pi = self._pi[self._n_aa * r : self._n_aa * (r + 1)]
                if numpy.allclose(pi, old_pi):
                    raise ValueError("New pi values are nearly equal to the old ones.")
                old_plus_derivative = [old_pi[j] + d_ddg * old_dpi[j] for j in range(self._n_aa)]
                if not numpy.allclose(pi, old_plus_derivative):
                    raise ValueError("New pi values are not nearly equal to the old ones plus derivative.")
                # 
                # Now test that we are computing the derivatives of the log likelihoods correctly
                ntotal_ll += 1
                if numpy.allclose(old_residue_log_likelihood, self._residue_log_likelihoods[r], rtol=1e-7) or not numpy.allclose(old_residue_log_likelihood + d_ddg * old_d_residue_log_likelihoods[zi], self._residue_log_likelihoods[r], rtol=1e-7):
                    nflagged_ll += 1
                    old_new_diff = abs(old_residue_log_likelihood - self._residue_log_likelihoods[r])
                    oldplusderiv_new_diff = abs(old_residue_log_likelihood + d_ddg * old_d_residue_log_likelihoods[zi] - self._residue_log_likelihoods[r])
                    if (oldplusderiv_new_diff > old_new_diff / 50. and oldplusderiv_new_diff >= 2e-6):
                        print "old = %f, old + derivative = %f, new = %f, old/new diff = %f, oldplusderiv/new diff %f" % (old_residue_log_likelihood, old_residue_log_likelihood + d_ddg * old_d_residue_log_likelihoods[zi], self._residue_log_likelihoods[r], old_new_diff, oldplusderiv_new_diff)
                        raise ValueError("The old log likelihood plus the derivative is not much closer to the new log likelihood than the old log likelihood without the derivative.")
                #
                # Now test that we are computing the log prior probability derivatives correctly
                if numpy.allclose(old_residue_log_priorprob, self._residue_log_priorprob[r], rtol=1e-9) or not numpy.allclose(old_residue_log_priorprob + d_ddg * old_d_residue_log_priorprob[zi], self._residue_log_priorprob[r], rtol=1e-9):
                    if abs(old_residue_log_priorprob + d_ddg * old_d_residue_log_priorprob[zi] - self._residue_log_priorprob[r]) > 0.01 * abs(old_residue_log_priorprob - self._residue_log_priorprob[r]) and abs(old_residue_log_priorprob - self._residue_log_priorprob[r]) > 1e-7:
                        print "old =", old_residue_log_priorprob, "new =", self._residue_log_priorprob[r], "old + derivative =", old_residue_log_priorprob + d_ddg * old_d_residue_log_priorprob[zi]
                        raise ValueError("Problem with prior probability derivatives.")
                #
                # Now test that we are computing the log posterior probability derivatives correctly
                old_new_diff = abs(old_residue_log_posterior - self._residue_log_posterior[r])
                oldplusderiv_new_diff = abs(old_residue_log_posterior + d_ddg * old_d_residue_log_posterior[zi] - self._residue_log_posterior[r])
                if numpy.allclose(old_new_diff, 0.0, rtol=1e-7) or not numpy.allclose(oldplusderiv_new_diff, 0, rtol=1e-7):
                    if (oldplusderiv_new_diff > old_new_diff / 50. and oldplusderiv_new_diff >= 2e-6):
                        print "old/new diff = %f, oldplusderiv/new diff %f" % (old_new_diff, oldplusderiv_new_diff)
                        raise ValueError("The old log posterior plus the derivative is not much closer to the new log posterior than the old log posterior without the derivative.")
                #
                # Print some output
                if print_output:
                    print "Completed testing for changes for this ddG."
            if nequal_m / float(ntotal_m) > 0.1:
                raise ValueError("More than 10% of the time the new and old M matrices were nearly equal.")
            if nflagged_ll / float(ntotal_ll) > 0.5:
                raise ValueError("Flagging more than 50% of the log likelihood values.")

    def __setstate__(self, d):
        """Sets the state of an object.

        This method is explicitly defined in case a DDGSet object is being
            unpickled from a serialized instance.  It ensures that the
            'Initialize' method is called for calculating log likelihoods
            with the C extensions.
        """
        self.__dict__ = d.copy()
        if self._use_c_extensions and not clog_likelihood.IsInitialized():
            clog_likelihood.Initialize(self._length, self._n_aa, self._n_tips, self._n_internal, self._tips, self._descendents, self._descendent_is_tip, self._n_unique, self._unique_uts_index_n_aa2, self._underflow)

    def Sequence(self):
        """Returns the reference sequence for this ddG set.

        The returned value is the protein sequence that defines the reference amino
            acids, and was set during object initialization."""
        return self._protseq

    def DDGDict(self):
        """Returns the ddG values in a dictionary.
       
        These are the ddG values for mutating the sequence that was initialized
            with this DDGSet object (the one returned by Sequence()).
        The returned variable specifies a set of ddG values in the form of a
            dictionary 'ddgs'.  The keys of this dictionary are integers specifying
            the residue number (0 <= r < length).  'ddgs[r]' is the 2-tuple 
            '(wt, iddgs)' where 'wt' is the identity of residue 'r' and
            'iddgs' is a dictionary with 'iddgs[mut]' giving the ddG value
            for mutating this residue from 'wt' to 'mut'.
        """
        ddgs = {}
        for r in range(self._length):
            wt = self._protseq[r]
            iddgs = self._ddgs[r].copy()
            ddgs[r] = (wt, iddgs)
        return ddgs

    def Length(self):
        """Returns the length of the protein for which the ddG values are specified."""
        return self._length

    def DDGRange(self):
        """Returns the allowable range of ddG values.

        The returned variable is an integer giving the value of ddg_range
            that was set upon initialization.  All ddG values are 
            constrained by -DDGRange() <= ddg <= DDGRange()
        """
        return self._ddg_range

    def _MList(self, residue_to_compute):
        """Returns a variable 'mlist' holding mutation probabilities for a residue.

        This list holds the "M" values.  Mxy(ut) gives the probability that after
            an average of ut mutations per codon have occurred, an amino acid
            that was originally y will now be x.
        Here is a brief description of 'mlist'.  It is a numpy ndarray
            of type 'float_' (C type double) that is one dimensional and of length
            n_aa2 * _n_unique where n_aa2 is the square of the number of
            amino acids.  The element of this array corresponding to
            the probability of mutating residue r 'residue_to_compute'
            from y (one letter amino acid code) to x (one letter amino acid code)
            in the time ut = _unique_uts[ut_index] is indexed in this list by
            ut_index * n_aa2 + xi * n_aa + yi
            where xi = AminoAcidIndexDict()[x] and yi = AminoAcidIndexDict()[y].
        The calling variable is 'residue_to_compute'.  This should be
            an integer specifying a residue number (0 <= residue_to_compute < Length()),
            specifying for which residue we recompute the entries in 'mlist'.

        This method performs no hashing of computed values.  So calling
            the method twice leads to recomputing the values entirely anew.
        """
        uts = self._unique_uts
        only_need = self._mlist_only_need
        assert 0 <= residue_to_compute < self._length and isinstance(residue_to_compute, int)
        assert isinstance(uts, list) and len(uts) >= 1
        assert isinstance(only_need, list) and len(only_need) == self._length * len(uts)
        n_aa = self._n_aa
        n_aa2 = self._n_aa2
        mlist = numpy.ndarray(n_aa2 * self._n_unique, dtype = 'float_')
        if self._use_c_extensions: # perform computation in C
            # use the c extension to fill '_mlist' with the appropriate values.
            # Note that is extremely important that the arrays in self._grs be C style contiguous.
            # The arrays should either be of float or complex type.
            # This should have been assured when they are generated.  If you get an undiagnosable
            # errors, this is probably the problem
            return_value = cddg_inference.MList(uts, only_need, n_aa, self._length, self._grs, mlist, residue_to_compute)
            return mlist
        # if we made it here, perform computation in python
        only_need_index = residue_to_compute * self._n_unique
        mlist_index = 0
        (gr_diag, p, p_inv) = self._grs[residue_to_compute]
        for ut in uts:
            assert ut >= 0 and isinstance(ut, (int, float)), "Invalid ut of %s" % (str(ut))
            xi = only_need[only_need_index]
            if xi == -1:
                mlist_index += n_aa2
            else:
                # compute all entries of mlist for this value of ut
                gr_diag_ut = numpy.exp(ut * gr_diag) # exponentiate the diagonal matrix
                gr_ut = numpy.dot(p * gr_diag_ut, p_inv)
                # convert to list with proper indexing by xi, yi
                for mi in numpy.abs(gr_ut.reshape(-1).real):
                    mlist[mlist_index] = mi
                    mlist_index += 1
            only_need_index += 1
        return mlist

    def _dMList(self, residue_to_compute):
        """Returns 'dmlist' holding derivatives of mutation probabilities for a residue.

        This list holds the "dM/d(ddG)" values.  Mxy(ut) gives the probability that after
            an average of ut mutations per codon have occurred, an amino acid
            that was originally y will now be x.  dMxy(ut) / d(ddG) is the
            derivative with respect to the ddG value.
        Here is a brief description of 'mlist'.  It is a numpy ndarray
            of type 'float_' (C type double) that is one dimensional and of length
            n_aa3 * _n_unique where n_aa3 is the cube of the number of
            amino acids.  The element of this array corresponding to derivative of
            the probability of mutating residue r = 'residue_to_compute'
            from y (one letter amino acid code) to x (one letter amino acid code)
            in the time ut = _unique_uts[ut_index] with respect to the ddG value from
            the wildtype residue to residue z is indexed in this list by
            ut_index * n_aa3 + zi * n_aa2 + xi * n_aa + yi
            where xi = AminoAcidIndexDict()[x], yi = AminoAcidIndexDict()[y]
            and zi = AminoAcidIndexDict()[z].  Values for z equal to the wildtype
            residue are undefined (may have arbitrary values).
        The calling variable is 'residue_to_compute'.  This should be
            an integer specifying a residue number (0 <= residue_to_compute < Length()),
            specifying for which residue we recompute the entries in 'mlist'.

        Currently, this method without using the C extension is EXTREMELY inefficient,
            since the python method called to compute the derivatives performs no
            hashing of quantities that are computed repeatedly.

        This method performs no hashing of computed values.  So calling
            the method twice leads to recomputing the values entirely anew.
        """
        uts = self._unique_uts
        only_need = self._mlist_only_need
        assert 0 <= residue_to_compute < self._length and isinstance(residue_to_compute, int)
        assert isinstance(uts, list) and len(uts) >= 1
        assert isinstance(only_need, list) and len(only_need) == self._length * len(uts)
        iwt = self._aa_index[self._protseq[residue_to_compute]] # wildtype residue index
        n_aa = self._n_aa
        n_aa2 = self._n_aa2
        n_aa3 = n_aa * n_aa2
        dmlist = numpy.zeros(n_aa3 * self._n_unique, dtype = 'float_')
        (gr_diag, p, p_inv) = self._grs[residue_to_compute]
        if self._use_c_extensions: # perform computation in C
            # use the c extension to fill '_dmlist' with the appropriate values.
            # Note that is extremely important that the arrays in self._grs be C style contiguous.
            # The arrays should either be of float or complex type.
            # This should have been assured when they are generated.  If you get an undiagnosable
            # errors, this is probably the problem
            brs = []
            p = p.reshape((n_aa, n_aa))
            p_inv = p_inv.reshape((n_aa, n_aa))
            for zi in range(n_aa):
                if zi == iwt:
                    brs.append(None)
                else:
                    dgr = self._dgrs[residue_to_compute][zi].reshape((n_aa, n_aa))
                    brs.append(numpy.dot(p_inv, numpy.dot(dgr, p)).ravel())
            return_value = cddg_inference.dMList(uts, only_need, n_aa, self._length, self._grs, dmlist, residue_to_compute, iwt, brs)
            return dmlist
        # if we made it here, perform computation in python
        assert not self._use_c_extensions
        only_need_index = residue_to_compute * self._n_unique
        dmlist_index = 0
        for ut in uts:
            assert ut >= 0 and isinstance(ut, (int, float)), "Invalid ut of %s" % (str(ut))
            xi = only_need[only_need_index]
            if xi == -1:
                dmlist_index += n_aa3
            else:
                for zi in range(n_aa):
                    # compute derivative with respect to ddG from wildtype to zi
                    if zi == iwt:
                        # derivative not defined when zi is the wildtype residue
                        dmlist_index += n_aa2
                    else:
                        dgr = self._dgrs[residue_to_compute][zi]
                        if self._use_c_extensions:
                            dgr = dgr.reshape((n_aa, n_aa))
                        dm = markovmatrix.DerivMatrixExponential(dgr, ut, p, p_inv, gr_diag)
                        # convert to list with proper indexing by xi, yi
                        for dmi in dm.reshape(-1).real:
                            dmlist[dmlist_index] = dmi
                            dmlist_index += 1
            only_need_index += 1
        return dmlist

    def C(self, x, y):
        """Returns the mutation probability.

        'x' and 'y' are one letter amino acid codes.
        This method returns the probability that a random single nucleotide mutation to a codon
            for amino acid y yields the new amino acid x.
        """
        assert x in self._amino_acids and y in self._amino_acids
        return self._c_set[(x, y)]

    def MutRate(self):
        """Returns the mutation rate u.

        This is the value that was set using 'mutrate' on initialization.
            It represents the number of mutations to a codon before an 
            amino acid substitution occurs, on average.
        """
        return self._mutrate

    def MaximizePosterior(self, nrandomstarts, printprogress=False, randomstartchoice=(0.0, 0.2)):
        """Adjusts the ddG values to maximize the posterior probability.

        The method adjust the ddG values of the DDGSet in order to find
            a maximum in the posterior probability of the phylogenetic
            tree data given the priors for the ddG values.  The method
            proceeds residue by residue.  For each residue, it performs a conjugate
            gradient optimization to maximize the posterior.  There is of course no
            guarantee that the optimum found will be the global maximum.
        'nrandomstarts' specifies that we try the posterior maximization from
            several different starting ddG values.  For each residue, we always try
            starting from the initial ddG values.  If 'nrandomstarts' is zero,
            then we simply perform this one maximization from the initial ddG values
            for each residue.  If 'nrandomstarts' is an integer greater than zero,
            we then try 'nrandomstarts' additional maximizations starting from
            other ddG values chosen according to 'randomstartchoice'.  For each
            residue, we keep the ddG values corresponding to the highest local
            maximum achieved.
        'printprogress' is a Boolean switch indicating whether we print progress
            during the maximization.  If it is True, we print an update after 
            every residue maximization.  If it is False (default value), we
            print no such update.
        'randomstartchoice' indicates how we choose the random starting ddG values
            for each of the 'nrandomstarts' choices of random starting ddG values
            for each residue.  'randomstartchoice' should be a 2-tuple of the form
            '(meanscale, sdscale)'.  Each ddG value is chosen at random from a Gaussian
            distribution with a mean of 'meanscale * ddgset.DDGRange()' and a standard
            deviation of 'meansd * ddgset.DDGRange()'.  Any ddG values generated that
            are outside of 'DDGRange()' are rounded up or down to be in this range.
            By default, 'randomstartchoice' is (0.0, 0.5).

        At the conclusion of this method, the ddG values are set to the values that give
            the maximum posterior probability.
        """
        ddg_range_tol = 1e-8 # we let ddG values get no closer to this than ddG range limits
        ddgtol = 0.25 # tolerance for considering ddG values roughly equal
        maxadditional = nrandomstarts + 2 # maximum number of times we will try restarting if optimization fails
        (ddgmin, ddgmax) = (-self._ddg_range + ddg_range_tol, self._ddg_range - ddg_range_tol)
        (meanscale, sdscale) = randomstartchoice
        if printprogress:
            print "Performining maximization of posterior probability with respect to ddG values..."
        for r in range(self._length):
            nonwt_aas = [aa for aa in self._amino_acids if aa != self._protseq[r]]
            nonwt_aa_indices = [self._aa_index[aa] for aa in nonwt_aas]
            if printprogress:
                print "\nPerforming maximization for residue %d of %d..." % (r, self._length)
                sys.stdout.flush()

            def F(ddgs):
                """Takes as input list of ddG values and returns negative residue log posterior.
                
                'ddgs'[i] is the ddG for mutating from wildtype residue to nonwt_aa_indices[i]"""
                assert len(ddgs) == len(nonwt_aas)
                for (aa, ddg) in zip(nonwt_aas, ddgs):
                    ddg_original = ddg
                    ddg = max(ddg, ddgmin)
                    ddg = min(ddg, ddgmax)
                    self._ddgs[r][aa] = ddg
                self._values_current[r] = False
                self._derivatives_current[r] = False
                self._ComputeInternals(r, compute_derivatives=False)
                return -self._residue_log_posterior[r]

            def F_and_dF(ddgs):
                """Takes as input list of ddGs, returns negative residue log posterior and derivatives.
                
                'ddgs'[i] is the ddG for mutating from wildtype residue to nonwt_aa_indices[i].
                Returned value is 2-tuple.  First element is negative log posterior.  Second element is
                    list with element 'i' giving the derivative of negative log posterior with
                    respect to 'ddgs[i]'"""
                assert len(ddgs) == len(nonwt_aas)
                for (aa, ddg) in zip(nonwt_aas, ddgs):
                    ddg_original = ddg
                    ddg = max(ddg, ddgmin)
                    ddg = min(ddg, ddgmax)
                    self._ddgs[r][aa] = ddg
                self._values_current[r] = False
                self._derivatives_current[r] = False
                self._ComputeInternals(r, compute_derivatives=True)
                return (-self._residue_log_posterior[r], [-self._d_residue_log_posterior[r][zi] for zi in nonwt_aa_indices])

            best_log_posterior = None
            irandomstart = nadditional = 0
            while irandomstart <= nrandomstarts + nadditional:
                if irandomstart == 0:
                    if printprogress:
                        self._ComputeInternals(r, compute_derivatives=True)
                        print "\tBeginning maximization from initial ddG values (initial log posterior = %g)..." % self._residue_log_posterior[r]
                else:
                    for aa in nonwt_aas:
                        rand = random.gauss(meanscale * self._ddg_range, sdscale * self._ddg_range)
                        rand = max(rand, ddgmin)
                        rand = min(rand, ddgmax)
                        self._ddgs[r][aa] = rand
                    if printprogress:
                        self._ComputeInternals(r, compute_derivatives=True)
                        print "\tBeginning maximization from randomized ddG values (start %d, initial log posterior = %g)..." % (irandomstart, self._residue_log_posterior[r])
                    self._values_current[r] = False
                    self._derivatives_current[r] = False
                ddgs0 = [self._ddgs[r][aa] for aa in nonwt_aas]
                try:
                    (ddgs, logposterior, iterations) = optimize.ConjugateGradient(ddgs0, F, F_and_dF)
                    logposterior = -logposterior
                    assert abs(self._residue_log_posterior[r] - logposterior) < 1.0e-8
                    if printprogress:
                        print "\tReached maximum log posterior of %g in %d iterations of conjugate gradient maximization." % (logposterior, iterations)
                        if best_log_posterior != None:
                            if numpy.allclose(ddgs, best_ddgs, atol=ddgtol, rtol=0):
                                print "\tddG values are equal to those from previous best maximum within a tolerance of %g." % ddgtol
                            else:
                                print "\tddG values differ from those from previous best maximum within a tolerance of %g." % ddgtol
                    if best_log_posterior == None or logposterior > best_log_posterior:
                        best_ddgs = ddgs
                        best_log_posterior = logposterior
                except optimize.OptimizationError:
                    exception_info = str(sys.exc_info()[1])
                    if nadditional >= maxadditional:
                        raise ValueError("Maximization is failing to converge for residue %d; we have performed %d tries with new starting values." % (r, nadditional))
                    if printprogress:
                        print "\tWARNING: Maximization failed, probably due to attempts by conjugate gradient search to move outside of ddG range.  The error message was %s.  Trying again with new initial values." % exception_info
                    nadditional += 1
                irandomstart += 1
            for (aa, ddg) in zip(nonwt_aas, best_ddgs):
                self._ddgs[r][aa] = ddg
            self._values_current[r] = False
            self._derivatives_current[r] = False
            if nrandomstarts > 0 and printprogress:
                print "\tHighest log posterior achieved was %g; saving these ddG values for residue %d." % (best_log_posterior, r)

    def _ComputeFs(self, r):
        """Computes the f values (fixation probabilities) based on current ddG values.

        The equation is:
            f = 1/2 - (1/2) tanh(beta * ddg - 0.5 ln(gamma/(1 - gamma)))

        The returned variable is an indexed list giving the f values.
        """
        assert 0 <= r < self._length
        # We now compute the list indexed_fs containing the fxyr values for residue r.  
        # This list is of length n_aa * n_aa.  The element yi * n_aa + xi contains the fxyr value
        # for mutating residue r (0 <= r < length) from yi to xi, where yi = aa_index[y] is
        # the integer code for residue y, and xi is the same for residue x.  Elements
        # corresponding to xi = yi are zero.
        fs = []
        wt = self._protseq[r] # this is the wildtype, or reference, amino acid for this residue
        ddgs_r = self._ddgs[r]
        amino_acids = self._amino_acids
        beta = self._beta
        half_log_gamma = 0.5 * math.log((self._gamma / (1.0 - self._gamma)))
        for y in amino_acids:
            try: 
                ref_to_y = ddgs_r[y]
            except KeyError:
                assert y == wt
                ref_to_y = 0.0 # y is the reference amino acid
            for x in amino_acids:
                if x != y: # do the computation in this case, otherwise value stays at 0
                    try:
                        ref_to_x = ddgs_r[x]
                    except KeyError:
                        assert x == wt
                        ref_to_x = 0.0 # x is the reference amino acid
                    fs.append(0.5 - 0.5 * math.tanh(beta * (ref_to_x - ref_to_y) - half_log_gamma))
                else:
                    fs.append(0.0)
        assert len(fs) == self._n_aa2
        return fs

    def _ComputeGr(self, r):
        """Computes the value of Gr for residue 'r'.  Also computes the group inverse of -Gr.

        The f values (self._indexed_fs, computed by self._ComputeFs) must be current.

        Returned variable is the tuple:
            ((gr_diag, p, p_inv), w, groupinverse)
        where w is the stationary distribution.
        """
        assert 0 <= r < self._length
        # The cddg_inference method for this task is currently obsolete, so computation is
        # performed in pure Python code.
        if False: # we currently skip this option
            pass 
        #if self._use_c_extensions:
        #    gr = cddg_inference.ComputeGr(self._n_aa, indexed_ddgs, self._p0sums, self._indexed_c_set)
        else:
            n_aa = self._n_aa
            index = 0
            gr = numpy.zeros(shape=(n_aa, n_aa), dtype = 'float_')
            indexed_fs = self._indexed_fs[r]
            indexed_c_set = self._indexed_c_set
            for yi in range(n_aa):
                gryy = 0.0
                for xi in range(n_aa):
                    if xi == yi:
                        index += 1
                        continue
                    frxy = indexed_fs[index]
                    grxy = frxy * indexed_c_set[index]
                    index += 1
                    gryy -= grxy
                    gr[(xi, yi)] = grxy
                gr[(yi, yi)] = gryy
        # we must diagonalize this matrix
        (gr_diag, p) = numpy.linalg.eig(gr)
        p_inv = numpy.linalg.inv(p)
        # now compute the group inverse of -Gr
        (w, groupinverse) = markovmatrix.MarkovGroupInverse(-gr, left_stochastic=True)
        # we now have gr = numpy.dot(p, numpy.dot(numpy.diagflat(gr_diag), p_inv))
        if self._use_c_extensions:
            # It is very important that we call ravel() if we are using the C extensions
            return ((gr_diag.ravel(), p.ravel(), p_inv.ravel()), w, groupinverse)
        else:
            return ((gr_diag, p, p_inv), w, groupinverse)

    def _Compute_dGr(self, r, zi):
        """Computes the derivative of Gr with respect to the ddG from the wildtype residue to zi.

        The f values in self._indexed_fs, computed by self._ComputeFs, must be current.
        
        The computation performed by this method can probably be accelerated dramatically
            by noting that dGr matrices for different values of zi differ only in the signs
            of some of the elements, that some computations are repeated, and that
            some multiplications might be performed in numpy."""
        assert 0 <= r < self._length
        assert 0 <= zi < self._n_aa
        assert self._protseq[r] != self._amino_acids[zi]
        n_aa = self._n_aa
        index = 0
        dgr = numpy.zeros(shape=(n_aa, n_aa), dtype='float_')
        indexed_fs = self._indexed_fs[r]
        indexed_c_set = self._indexed_c_set
        twobeta = 2.0 * self._beta
        for yi in range(n_aa):
            for xi in range(n_aa):
                frxy = indexed_fs[index]
                cxy = indexed_c_set[index]
                if zi == xi != yi:
                    dgr[(xi, yi)] = cxy * twobeta * frxy * (frxy - 1.0)
                elif zi == yi != xi:
                    dgr[(xi, yi)] = -cxy * twobeta * frxy * (frxy - 1.0)
                elif xi == yi != zi:
                    index2 = yi * n_aa + zi
                    czy = indexed_c_set[index2]
                    frzy = indexed_fs[index2]
                    dgr[(xi, yi)] = -czy * twobeta * frzy * (frzy - 1.0)
                elif xi == yi == zi:
                    index2 = yi * n_aa
                    for wi in range(n_aa):
                        if wi != zi:
                            cwy = indexed_c_set[index2]
                            frwy = indexed_fs[index2]
                            dgr[(xi, yi)] += cwy * twobeta * frwy * (frwy - 1.0)
                        index2 += 1
                index += 1
        if self._use_c_extensions:
            return dgr.ravel()
        else:
            return dgr

    def _GetPiValues(self, r):
        """Returns a list giving the pi values for residue 'r'.

        It is assumed that the current value of self._grs is already
            computed for r.  The returned variable is a list of
            real numbers giving the eigenvector corresponding to
            the eigenvalue of gr closest to zero.   If this
            eigenvector is complex, real components are taken.
            It is normalized to sum to zero.  So entry yi gives
            the equilibrium probability that residue r is amino
            acid yi.
        """
        n_aa = self._n_aa
        assert 0 <= r < self._length
        (gr_diag, p, p_inv) = self._grs[r]
        # One of the eigenvalues in gr_diag should be zero (or within numerical error of it).
        # We need to find the index of this eigenvalue as zero_eig_index, since it corresponds
        # to the principal eigenvector of gr + identity_matrix.
        zero_eig = abs(gr_diag[0])
        zero_eig_index = 0
        index = 1
        for eig in gr_diag[1 : ]:
            if abs(eig) < zero_eig:
                zero_eig_index = index
                zero_eig = abs(eig)
            index += 1
        if zero_eig > 1.0e-11:
            raise ValueError, "Eigenvalue %d is not close to zero in %s." % (zero_eig_index, str(gr_diag))
        if self._use_c_extensions:
            pir = p[zero_eig_index::n_aa].copy() # p has been raveled
        else:
            pir = p.transpose()[zero_eig_index].copy() # p has not been raveled
        pir = pir.real # make pir real
        pir /= numpy.sum(pir) # normalize pir
        pir = numpy.abs(pir) # make all entries positive
        assert 1.0e-12 > abs(1.0 - numpy.sum(pir)), "Values are not all of the same sign: %s" % (str(pir))
        return pir.tolist() # return as list

    def _ComputeInternals(self, r, compute_derivatives):
        """Computes internal variables used in calculations for a single residue.

        'r' specifies the residue for which we are computing the internal variables.
            The internal variables are recomputed based on the current ddG values
            for that residue, as stored in self._ddgs.

        'compute_derivatives' is a Boolean switch specifying whether or not we
            also compute the derivatives (we do if it is True, do not if it is
            False).
       
        This method needs to be run every time the internal ddG values are adjusted.

        After this method has run, the following variables are always guaranteed to be computed
            for the residue in question:            
                self._values_current
                self._derivatives_current
                self._indexed_fs
                self._grs
                self._group_inverses
                self._pi
                self._residue_log_likelihoods
                self._residue_log_priorprob
                self._residue_log_posterior
            Additionally, if the method is run with 'compute_derivatives' equal to True,
            the following variables are also guaranteed to be computed:
                self._dgrs
                self._dpi
                self._d_residue_log_likelihoods
                self._d_residue_log_priorprob
                self._d_residue_log_posterior

        After completion of this method, 'self._values_current[r]' is True.  If the
            method is called with 'compute_derivatives' True, then 'self._derivatives_current[r]'
            is also True, otherwise it is 'False'.  If this method is called with
            'compute_derivatives' False, and 'self._values_current[r]' is already True,
            then no computation is performed.  If this method is called with 'compute_derivatives'
            True, and self._derivatives_current' is already True, then no computation is
            performed.

        No economization of computation occurs if this method is called with 'compute_derivatives'
            True after earlier calling it with 'compute_derivatives' False.  Therefore, if you
            are going to eventualy call the method with 'compute_derivatives' True before
            adjusting the ddG values, then you should not both first calling it with
            'compute_derivatives' False.  On the other hand, if you are never going to call it
            with 'compute_derivatives' True for the current ddG values, then you save computation
            by calling it with 'compute_derivatives' False.
        """
        assert 0 <= r < self._length
        assert isinstance(compute_derivatives, bool)
        #
        # first check if needed values are current
        if self._derivatives_current[r] and not self._values_current[r]:
            raise ValueError("Problem with tracking of computations.")
        if compute_derivatives and self._derivatives_current[r]:
            return
        elif self._values_current[r] and not compute_derivatives:
            return
        #
        # self._values_current[r] is True if self._indexed_fs, self._grs, self._group_inverses,
        # self._pi, self._residue_log_likelihoods, and self._residue_log_priorprob are current
        # for residue 'r', and False otherwise
        self._values_current[r] = True
        #
        # self._derivatives_current[r] is True if self._dpi, self._dgrs, and 
        # self._d_residue_log_likelihoods are current for residue 'r', and False otherwise.
        self._derivatives_current[r] = compute_derivatives
        #
        # self._indexed_fs[r] is a list of length n_aa * n_aa.  Element yi * n_aa + xi contains  
        # the fxyr value for mutating residue r (0 <= r < length) from yi to xi, where
        # yi = aa_index[y] is the integer code for residue y, and xi is the same for residue x.  
        # Elements corresponding to xi = yi are zero.
        self._indexed_fs[r] = self._ComputeFs(r)
        #
        # self._grs[r] holds the matrix Gr for mutations to residue r.  This is stored
        # as a 3-tuple: (Gr_diag, P, P_inv) where Gr = P numpy.diagflat(Gr_diag) P_inv and
        # Gr_diag is the diagonal of the diagonal matrix.  Gr[(x, y)] is for mutating residue r
        # from y to x. 
        (gr, w, groupinverse) = self._ComputeGr(r)
        self._grs[r] = gr
        self._group_inverses[r] = groupinverse
        #
        # For each residue r, self._dgrs[r] is a list of length n_aa.  Element zi gives the
        # derivative of the Gr matrix for residue r with respect to the ddG for mutating
        # that residue from the wildtype residue at that position to residue z, where
        # zi = AminoAcidIndexDict()[z].  If z is equal to the wildtype residue, self._dgrs[r][zi]
        # is None.
        zi = 0
        for z in self._amino_acids:
            if z != self._protseq[r]:
                if compute_derivatives:
                    dgr = self._Compute_dGr(r, zi)
                    self._dgrs[r][zi] = dgr
                else:
                    self._dgrs[r][zi] = None
            zi += 1
        #
        # self._pi[n_aa * r + yi] gives the equilibrium probability that residue
        # r is y where yi = AminoAcidIndexDict()[y] and
        # n_aa = NAminoAcids(). 
        pir = self._GetPiValues(r)
        assert numpy.allclose(pir, w, rtol=1e-2, atol=1e-3), "pi = %s\nw = %s" % (str(pir), str(w))
        self._pi[r * self._n_aa : (r + 1) * self._n_aa] = pir
        #
        # self._dpi[n_aa2 * r + zi * n_aa + yi] gives the derivative of self._pi[n_aa * r + yi] with
        # respect to the ddG for mutating residue r from the wildtype residue to residue z.
        # Entries for zi equal to the wildtype residue are undefined.
        zi = 0
        for z in self._amino_acids:
            if z != self._protseq[r]:
                index = self._n_aa2 * r + zi * self._n_aa
                if compute_derivatives:
                    dgr = self._dgrs[r][zi]
                    if self._use_c_extensions:
                        dgr = dgr.reshape((self._n_aa, self._n_aa))
                    (pi, dpi) = markovmatrix.StationaryDistributionAndDerivative(None, dgr, pi_Agi=(w, groupinverse), left_stochastic=True)
                    self._dpi[index : index + self._n_aa] = dpi.tolist()
                else:
                    self._dpi[index : index + self._n_aa] = [None] * self._n_aa
            zi += 1
        #
        # Now that we have computed all pi and gr values, compute the residue log likelihoods and
        # their derivatives
        if compute_derivatives:
            (self._residue_log_likelihoods[r], self._d_residue_log_likelihoods[r]) = self._ResidueLogLikelihoodAndDeriv(self._MList(r), self._dMList(r), r)
        else:
            self._residue_log_likelihoods[r] = self._ResidueLogLikelihood(self._MList(r), r)
            self._d_residue_log_likelihoods[r] = None
        #
        # self._residue_log_priorprob[r] is the log of the prior probability of residue r
        # having the current ddG value.
        #
        # self._d_residue_log_priorprob[r][zi] is the derivative of self._residue_log_priorprob[r]
        # with respect to the ddG for mutating from the wildtype residue to residue z, with entries
        # for z equal to the wildtype residue undefined, where zi = AminoAcidIndexDict()[z]
        zi = 0
        self._residue_log_priorprob[r] = 0.0
        g = self._ddg_range
        for z in self._amino_acids:
            if z != self._protseq[r]:
                ddg = self._ddgs[r][z]
                (alpha, beta) = self._priors[r][z]
                self._residue_log_priorprob[r] += math.log(stats.BetaDistribution(ddg, alpha, beta, -g, g))
                if compute_derivatives:
                    self._d_residue_log_priorprob[r][zi] = (alpha - 1.0) / (ddg + g) - (beta - 1.0) / (g - ddg)
                else:
                    self._d_residue_log_priorprob[r][zi] = None
            zi += 1
        #
        # self._residue_log_posterior[r] is the log of the posterior probability of residue r
        # having the current ddG value.
        #
        # self._d_residue_log_posterior[r][zi] is the derivative of self._residue_log_posterior[r]
        # with respect to the ddG for mutating from the wildtype residue to residue z, with entries
        # for z equal to the wildtype residue undefined, where zi = AminoAcidIndexDict()[z]
        self._residue_log_posterior[r] = self._residue_log_likelihoods[r] + self._residue_log_priorprob[r]
        if compute_derivatives:
            zi = 0
            for z in self._amino_acids:
                if z != self._protseq[r]:
                    self._d_residue_log_posterior[r][zi] = self._d_residue_log_likelihoods[r][zi] + self._d_residue_log_priorprob[r][zi]
                zi += 1


    def _ResidueLogLikelihood(self, mlist, r):
        """Computes the log likelihood for a specific residue r.

        It is assumed that self._pi and self._grs are already computed for this residue
        This method does not use recursion, but instead relies on precomputed lists
            that are designed in such a way that the computation can be performed without
            recursion.  The basic idea is that the nodes have been numbered in such a way
            that any internal node that is a descendent of some other internal node has
            a smaller number.  The method then progresses up the set of internal node
            numbers, computing the likelihood of the sub-trees below the nodes.  By the
            time it reaches the final internal node (the root node), it is has computed
            the likelihood of the whole tree.
        The returned variable is the number 'rloglikelihood'.
            'rloglikelihood' is a number giving the log likelihood for this residue.
        """
        assert isinstance(r, int) and 0 <= r < self._length
        if self._use_c_extensions: 
            rloglikelihood = clog_likelihood.ResidueLogLikelihood(mlist, self._pi, r)
        else:
            # Perform calculation in pure python code.
            # p is a list of length n_internal * n_aa.  It is recomputed for each
            #    residue each time the function is called.  p[i * n_aa + y] is the cumulative 
            #    probability that internal node i is residue y given the information in the subtree
            #    below i.  Because of the ordering of the internal node numbering, the entries for
            #    i = n_aa - 1 correspond to the root node of the tree.
            #
            # mlist entries are indexed by ut_index * n_aa2 + xi * n_aa + yi
            n_tips = self._n_tips
            n_aa = self._n_aa
            n_aa2 = self._n_aa**2
            p = [1.0] * self._n_internal * n_aa
            unique_uts_index_n_aa2 = self._unique_uts_index_n_aa2
            underflow = self._underflow
            r_log_scaling_sum = 0.0  # total number we add to residue log likelihood to correct for scaling
            r_times_n_tips = r * n_tips
            for i in range(self._n_internal): # loop over internal nodes
                right_index = 2 * i
                left_index = right_index + 1
                i_times_n_aa = i * n_aa
                i_times_n_aa2 = i * n_aa2
                if self._descendent_is_tip[right_index] and self._descendent_is_tip[left_index]:
                    # both descendents are tips
                    tr = self._descendents[right_index]
                    tl = self._descendents[left_index]
                    xr = self._tips[r_times_n_tips + tr]
                    xl = self._tips[r_times_n_tips + tl]
                    if xr == -1: # a gap on right descendent
                        if xl == -1: # both left and right are gaps
                            for y in range(0, n_aa):
                                p[i_times_n_aa + y] = 1.0
                        else: # gap is only on right descendent
                            index = unique_uts_index_n_aa2[tl] + xl * n_aa
                            for y in range(0, n_aa):
                                p[i_times_n_aa + y] = mlist[index]
                                index += 1
                    elif xl == -1: # a gap only on left descendent
                        index = unique_uts_index_n_aa2[tr] + xr * n_aa
                        for y in range(0, n_aa):
                            p[i_times_n_aa + y] = mlist[index]
                            index += 1
                    else: # no gaps on either descendent
                        indexr = unique_uts_index_n_aa2[tr] + xr * n_aa
                        indexl = unique_uts_index_n_aa2[tl] + xl * n_aa
                        for y in range(0, n_aa):
                            (mr, ml) = (mlist[indexr], mlist[indexl])
                            p[i_times_n_aa + y] = mr * ml
                            indexr += 1
                            indexl += 1
                elif self._descendent_is_tip[right_index] and not self._descendent_is_tip[left_index]:
                    # right descendent is tip, left descendent is not
                    tr = self._descendents[right_index]
                    xr = self._tips[r_times_n_tips + tr]
                    jl = self._descendents[left_index]
                    jl_times_n_aa = jl * n_aa
                    jl_times_n_aa2 = jl_times_n_aa * n_aa
                    if xr == -1: # right descendent is a gap
                        indexl = unique_uts_index_n_aa2[n_tips + jl]
                        for y in range(0, n_aa):
                            indexly = indexl + y
                            piy = 0.0
                            for x in range(0, n_aa):
                                mly = mlist[indexly]
                                pjlx = p[jl_times_n_aa + x]
                                piy += mly * pjlx
                                indexly += n_aa
                            p[i_times_n_aa + y] = piy
                    else: # not a gap 
                        indexl = unique_uts_index_n_aa2[n_tips + jl]
                        indexr = unique_uts_index_n_aa2[tr] + xr * n_aa
                        for y in range(0, n_aa):
                            indexly = indexl + y
                            piy = 0.0
                            dp_index = i_times_n_aa2 + y * n_aa
                            for x in range(0, n_aa):
                                mly = mlist[indexly]
                                pjlx = p[jl_times_n_aa + x]
                                piy += mly * pjlx
                                indexly += n_aa
                            mr = mlist[indexr]
                            p[i_times_n_aa + y] = piy * mr
                            indexr += 1
                elif not self._descendent_is_tip[right_index] and self._descendent_is_tip[left_index]:
                    # left descendent is tip, right descendent is not
                    jr = self._descendents[right_index]
                    jr_times_n_aa = jr * n_aa
                    tl = self._descendents[left_index]
                    xl = self._tips[r_times_n_tips + tl]
                    if xl == -1: # left descendent is a gap
                        indexr = unique_uts_index_n_aa2[n_tips + jr]
                        for y in range(0, n_aa):
                            indexry = indexr + y
                            piy = 0.0
                            for x in range(0, n_aa):
                                mry = mlist[indexry]
                                pjrx = p[jr_times_n_aa + x]
                                piy += mry * pjrx
                                indexry += n_aa
                            p[i_times_n_aa + y] = piy
                    else: # not a gap 
                        indexr = unique_uts_index_n_aa2[n_tips + jr]
                        indexl = unique_uts_index_n_aa2[tl] + xl * n_aa
                        for y in range(0, n_aa):
                            indexry = indexr + y
                            piy = 0.0
                            for x in range(0, n_aa):
                                mry = mlist[indexry]
                                pjrx = p[jr_times_n_aa + x]
                                piy += mry * pjrx
                                indexry += n_aa
                            ml = mlist[indexl]
                            p[i_times_n_aa + y] = piy * ml
                            indexl += 1
                else: 
                    # neither descendent is a tip; both are internal nodes
                    jr = self._descendents[right_index]
                    jr_times_n_aa = jr * n_aa
                    jl = self._descendents[left_index]
                    jl_times_n_aa = jl * n_aa
                    indexr = unique_uts_index_n_aa2[n_tips + jr]
                    indexl = unique_uts_index_n_aa2[n_tips + jl]
                    for y in range(0, n_aa):
                        piyr = 0.0
                        piyl = 0.0
                        indexry = indexr + y
                        indexly = indexl + y
                        for x in range(0, n_aa):
                            pjrx = p[jr_times_n_aa + x]
                            mry = mlist[indexry]
                            piyr += mry * pjrx
                            pjlx = p[jl_times_n_aa + x]
                            mly = mlist[indexly]
                            piyl += mly * pjlx
                            indexry += n_aa
                            indexly += n_aa
                        p[i_times_n_aa + y] = piyr * piyl
                if not ((i + 1) % underflow): # scale probabilities for this node to avoid underflow
                    # We determine the maximum conditional probability of the subtree below node i
                    # given that i is some particular amino acid y.  We them divide all of the
                    # conditional probabilities for node i being y by this maximum, and save the
                    # logarithm of this scaling factor to be added to the final log likelihood.
                    cond_prob_list = [p[i_times_n_aa + y] for y in range(0, n_aa)]
                    scale = max(cond_prob_list) # the scaling factor
                    if scale != 0: # don't scale if the scaling factor is zero
                        for y in range(0, n_aa):
                            p[i_times_n_aa + y] /= scale # scale each conditional probability
                        r_log_scaling_sum += math.log(scale) # add the log of this scaling factor
            # we have now computed p for the root node.  Compute overall likelihood.
            root = self._n_internal - 1 # index of root node
            pr = 0.0
            dpr = [0.0] * n_aa
            r_times_n_aa = r * n_aa
            root_times_n_aa = root * n_aa
            for y in range(0, n_aa):
                piry = self._pi[r_times_n_aa + y]
                pry = p[root_times_n_aa + y]
                pr += pry * piry
            if pr == 0:
                raise ValueError, "Underflow occurred in the computation for residue %d." % r
            rloglikelihood = math.log(pr) # logarithm of scaled probability
            # Now correct for the scaling
            rloglikelihood += r_log_scaling_sum
        if not (-1e300 < rloglikelihood < 0):
            # print some diagnostic information before raising an exception
            sys.stderr.write("The value of grs is being serialized to '_debug_gr.pickle'.\n")
            sys.stderr.write("Here are the values of pi for this residue:\n")
            for y in range(0, self._n_aa):
                sys.stderr.write("\t%f\n" % self._pi[r * self._n_aa + y])
            cPickle.dump(self._grs[r], open('_debug_gr.pickle', 'w'))
            raise ValueError, "Invalid log likelihood of %s for residue %d." % (str(rloglikelihood), r)
        return rloglikelihood # return log likelihood

    def _ResidueLogLikelihoodAndDeriv(self, mlist, dmlist, r):
        """Computes the log likelihood for a specific residue r, and derivatives with respect to ddGs.

        It is assumed that self._pi, self._dpi, and self._grs are already computed for this residue
        This method does not use recursion, but instead relies on precomputed lists
            that are designed in such a way that the computation can be performed without
            recursion.  The basic idea is that the nodes have been numbered in such a way
            that any internal node that is a descendent of some other internal node has
            a smaller number.  The method then progresses up the set of internal node
            numbers, computing the likelihood of the sub-trees below the nodes.  By the
            time it reaches the final internal node (the root node), it is has computed
            the likelihood of the whole tree.
        The returned variable is the 2-tuple '(rloglikelihood, drloglikelihood)'.
            'rloglikelihood' is a number giving the log likelihood for this residue.
            'drloglikelihood' is a list, with 'drloglikelihood[zi]' the 
                derivative of this log likelihood with respect
                to the ddG for mutating from the wildtype residue to residue
                zi = AminoAcidIndex()[z].  Entries for zi corresponding to the wildtype
                residue are undefined.
        """
        assert isinstance(r, int) and 0 <= r < self._length
        if self._use_c_extensions: 
            (rloglikelihood, drloglikelihood) = clog_likelihood.ResidueLogLikelihoodAndDeriv(mlist, dmlist, self._pi, self._dpi, r, self._aa_index[self._protseq[r]])
        else:
            # Perform calculation in pure python code.
            # p is a list of length n_internal * n_aa.  It is recomputed for each
            #    residue each time the function is called.  p[i * n_aa + y] is the cumulative 
            #    probability that internal node i is residue y given the information in the subtree
            #    below i.  Because of the ordering of the internal node numbering, the entries for
            #    i = n_aa - 1 correspond to the root node of the tree.
            # dp is a list of length n_internal * n_aa2.  dp[i * n_aa2 + y * n_aa + z] is the 
            #    derivative of p[i * n_aa + y] (the cumulative probability that internal node i
            #    is residue y given the information in the subtree below i) with respect to the
            #    ddG for mutating the residue from its wildtype identity to z, with entries
            #    for z equal to the wildtype identity being undefined.  
            #
            # mlist entries are indexed by ut_index * n_aa2 + xi * n_aa + yi
            # dmlist entries are indexed by ut_index * n_aa3 + zi * n_aa2 + xi * n_aa + yi
            n_tips = self._n_tips
            n_aa = self._n_aa
            n_aa2 = self._n_aa**2
            p = [1.0] * self._n_internal * n_aa
            dp = [0.0] * self._n_internal * n_aa2
            unique_uts_index_n_aa2 = self._unique_uts_index_n_aa2
            underflow = self._underflow
            r_log_scaling_sum = 0.0  # total number we add to residue log likelihood to correct for scaling
            r_times_n_tips = r * n_tips
            nonwt_aa_indices = [z for z in range(n_aa) if z != self._aa_index[self._protseq[r]]]
            for i in range(self._n_internal): # loop over internal nodes
                right_index = 2 * i
                left_index = right_index + 1
                i_times_n_aa = i * n_aa
                i_times_n_aa2 = i * n_aa2
                if self._descendent_is_tip[right_index] and self._descendent_is_tip[left_index]:
                    # both descendents are tips
                    tr = self._descendents[right_index]
                    tl = self._descendents[left_index]
                    xr = self._tips[r_times_n_tips + tr]
                    xl = self._tips[r_times_n_tips + tl]
                    if xr == -1: # a gap on right descendent
                        if xl == -1: # both left and right are gaps
                            for y in range(0, n_aa):
                                p[i_times_n_aa + y] = 1.0
                                dp_index = i_times_n_aa2 + y * n_aa
                                for z in nonwt_aa_indices:
                                    dp[dp_index + z] = 0.0
                        else: # gap is only on right descendent
                            index = unique_uts_index_n_aa2[tl] + xl * n_aa
                            dm_index = n_aa * unique_uts_index_n_aa2[tl] + xl * n_aa
                            for y in range(0, n_aa):
                                p[i_times_n_aa + y] = mlist[index]
                                dp_index = i_times_n_aa2 + y * n_aa
                                for z in nonwt_aa_indices:
                                    dp[dp_index + z] = dmlist[dm_index + n_aa2 * z]
                                index += 1
                                dm_index += 1
                    elif xl == -1: # a gap only on left descendent
                        index = unique_uts_index_n_aa2[tr] + xr * n_aa
                        dm_index = n_aa * unique_uts_index_n_aa2[tr] + xr * n_aa
                        for y in range(0, n_aa):
                            p[i_times_n_aa + y] = mlist[index]
                            dp_index = i_times_n_aa2 + y * n_aa
                            for z in nonwt_aa_indices:
                                dp[dp_index + z] = dmlist[dm_index + n_aa2 * z]
                            index += 1
                            dm_index += 1
                    else: # no gaps on either descendent
                        indexr = unique_uts_index_n_aa2[tr] + xr * n_aa
                        indexl = unique_uts_index_n_aa2[tl] + xl * n_aa
                        dm_indexr = n_aa * unique_uts_index_n_aa2[tr] + xr * n_aa
                        dm_indexl = n_aa * unique_uts_index_n_aa2[tl] + xl * n_aa
                        for y in range(0, n_aa):
                            (mr, ml) = (mlist[indexr], mlist[indexl])
                            p[i_times_n_aa + y] = mr * ml
                            dp_index = i_times_n_aa2 + y * n_aa
                            for z in nonwt_aa_indices:
                                z_n_aa2 = z * n_aa2
                                dp[dp_index + z] = dmlist[dm_indexr + z_n_aa2] * ml + mr * dmlist[dm_indexl + z_n_aa2]
                            indexr += 1
                            indexl += 1
                            dm_indexr += 1
                            dm_indexl += 1
                elif self._descendent_is_tip[right_index] and not self._descendent_is_tip[left_index]:
                    # right descendent is tip, left descendent is not
                    tr = self._descendents[right_index]
                    xr = self._tips[r_times_n_tips + tr]
                    jl = self._descendents[left_index]
                    jl_times_n_aa = jl * n_aa
                    jl_times_n_aa2 = jl_times_n_aa * n_aa
                    if xr == -1: # right descendent is a gap
                        indexl = unique_uts_index_n_aa2[n_tips + jl]
                        dm_indexl = n_aa * indexl
                        for y in range(0, n_aa):
                            indexly = indexl + y
                            dm_indexly = dm_indexl + y
                            piy = 0.0
                            dp_index = i_times_n_aa2 + y * n_aa
                            for x in range(0, n_aa):
                                mly = mlist[indexly]
                                pjlx = p[jl_times_n_aa + x]
                                piy += mly * pjlx
                                dpjlx_index = jl_times_n_aa2 + x * n_aa
                                for z in nonwt_aa_indices:
                                    dp[dp_index + z] += dmlist[dm_indexly + z * n_aa2] * pjlx + mly * dp[dpjlx_index + z]
                                indexly += n_aa
                                dm_indexly += n_aa
                            p[i_times_n_aa + y] = piy
                    else: # not a gap 
                        indexl = unique_uts_index_n_aa2[n_tips + jl]
                        indexr = unique_uts_index_n_aa2[tr] + xr * n_aa
                        dm_indexl = n_aa * indexl
                        dm_indexr = unique_uts_index_n_aa2[tr] * n_aa + xr * n_aa
                        for y in range(0, n_aa):
                            indexly = indexl + y
                            dm_indexly = dm_indexl + y
                            piy = 0.0
                            dp_index = i_times_n_aa2 + y * n_aa
                            for x in range(0, n_aa):
                                mly = mlist[indexly]
                                pjlx = p[jl_times_n_aa + x]
                                piy += mly * pjlx
                                dpjlx_index = jl_times_n_aa2 + x * n_aa
                                for z in nonwt_aa_indices:
                                    dp[dp_index + z] += dmlist[dm_indexly + z * n_aa2] * pjlx + mly * dp[dpjlx_index + z]
                                indexly += n_aa
                                dm_indexly += n_aa
                            mr = mlist[indexr]
                            p[i_times_n_aa + y] = piy * mr
                            for z in nonwt_aa_indices:
                                dp[dp_index + z] *= mr
                                dp[dp_index + z] += dmlist[dm_indexr + z * n_aa2] * piy
                            indexr += 1
                            dm_indexr += 1
                elif not self._descendent_is_tip[right_index] and self._descendent_is_tip[left_index]:
                    # left descendent is tip, right descendent is not
                    jr = self._descendents[right_index]
                    jr_times_n_aa = jr * n_aa
                    jr_times_n_aa2 = jr_times_n_aa * n_aa
                    tl = self._descendents[left_index]
                    xl = self._tips[r_times_n_tips + tl]
                    if xl == -1: # left descendent is a gap
                        indexr = unique_uts_index_n_aa2[n_tips + jr]
                        dm_indexr = n_aa * indexr
                        for y in range(0, n_aa):
                            indexry = indexr + y
                            dm_indexry = dm_indexr + y
                            piy = 0.0
                            dp_index = i_times_n_aa2 + y * n_aa
                            for x in range(0, n_aa):
                                mry = mlist[indexry]
                                pjrx = p[jr_times_n_aa + x]
                                piy += mry * pjrx
                                dpjrx_index = jr_times_n_aa2 + x * n_aa
                                for z in nonwt_aa_indices:
                                    dp[dp_index + z] += dmlist[dm_indexry + z * n_aa2] * pjrx + mry * dp[dpjrx_index + z]
                                indexry += n_aa
                                dm_indexry += n_aa
                            p[i_times_n_aa + y] = piy
                    else: # not a gap 
                        indexr = unique_uts_index_n_aa2[n_tips + jr]
                        indexl = unique_uts_index_n_aa2[tl] + xl * n_aa
                        dm_indexr = n_aa * indexr
                        dm_indexl = unique_uts_index_n_aa2[tl] * n_aa + xl * n_aa
                        for y in range(0, n_aa):
                            indexry = indexr + y
                            dm_indexry = dm_indexr + y
                            piy = 0.0
                            dp_index = i_times_n_aa2 + y * n_aa
                            for x in range(0, n_aa):
                                mry = mlist[indexry]
                                pjrx = p[jr_times_n_aa + x]
                                piy += mry * pjrx
                                dpjrx_index = jr_times_n_aa2 + x * n_aa
                                for z in nonwt_aa_indices:
                                    dp[dp_index + z] += dmlist[dm_indexry + z * n_aa2] * pjrx + mry * dp[dpjrx_index + z]
                                indexry += n_aa
                                dm_indexry += n_aa
                            ml = mlist[indexl]
                            p[i_times_n_aa + y] = piy * ml
                            for z in nonwt_aa_indices:
                                dp[dp_index + z] *= ml
                                dp[dp_index + z] += dmlist[dm_indexl + z * n_aa2] * piy
                            indexl += 1
                            dm_indexl += 1
                else: 
                    # neither descendent is a tip; both are internal nodes
                    jr = self._descendents[right_index]
                    jr_times_n_aa = jr * n_aa
                    jr_times_n_aa2 = jr_times_n_aa * n_aa
                    jl = self._descendents[left_index]
                    jl_times_n_aa = jl * n_aa
                    jl_times_n_aa2 = jl_times_n_aa * n_aa
                    indexr = unique_uts_index_n_aa2[n_tips + jr]
                    indexl = unique_uts_index_n_aa2[n_tips + jl]
                    dm_indexr = n_aa * indexr
                    dm_indexl = n_aa * indexl
                    for y in range(0, n_aa):
                        dp_index = i_times_n_aa2 + y * n_aa
                        piyr = 0.0
                        piyl = 0.0
                        indexry = indexr + y
                        indexly = indexl + y
                        dm_indexry = dm_indexr + y
                        dm_indexly = dm_indexl + y
                        dpiyr = [0.0] * n_aa
                        dpiyl = [0.0] * n_aa
                        for x in range(0, n_aa):
                            pjrx = p[jr_times_n_aa + x]
                            mry = mlist[indexry]
                            piyr += mry * pjrx
                            pjlx = p[jl_times_n_aa + x]
                            mly = mlist[indexly]
                            piyl += mly * pjlx
                            dpjrx_index = jr_times_n_aa2 + x * n_aa
                            dpjlx_index = jl_times_n_aa2 + x * n_aa
                            for z in nonwt_aa_indices:
                                z_n_aa2 = z * n_aa2
                                dpiyr[z] += dmlist[dm_indexry + z_n_aa2] * pjrx + mry * dp[dpjrx_index + z]
                                dpiyl[z] += dmlist[dm_indexly + z_n_aa2] * pjlx + mly * dp[dpjlx_index + z]
                            indexry += n_aa
                            indexly += n_aa
                            dm_indexry += n_aa
                            dm_indexly += n_aa
                        p[i_times_n_aa + y] = piyr * piyl
                        for z in nonwt_aa_indices:
                            dp[dp_index + z] = dpiyr[z] * piyl + piyr * dpiyl[z]
                if not ((i + 1) % underflow): # scale probabilities for this node to avoid underflow
                    # We determine the maximum conditional probability of the subtree below node i
                    # given that i is some particular amino acid y.  We them divide all of the
                    # conditional probabilities for node i being y by this maximum, and save the
                    # logarithm of this scaling factor to be added to the final log likelihood.
                    cond_prob_list = [p[i_times_n_aa + y] for y in range(0, n_aa)]
                    scale = max(cond_prob_list) # the scaling factor
                    if scale != 0: # don't scale if the scaling factor is zero
                        for y in range(0, n_aa):
                            p[i_times_n_aa + y] /= scale # scale each conditional probability
                            dp_index = i_times_n_aa2 + y * n_aa
                            for z in nonwt_aa_indices:
                                dp[dp_index + z] /= scale
                        r_log_scaling_sum += math.log(scale) # add the log of this scaling factor
            # we have now computed p for the root node.  Compute overall likelihood.
            root = self._n_internal - 1 # index of root node
            pr = 0.0
            dpr = [0.0] * n_aa
            r_times_n_aa = r * n_aa
            root_times_n_aa = root * n_aa
            root_times_n_aa2 = root * n_aa2
            for y in range(0, n_aa):
                piry = self._pi[r_times_n_aa + y]
                pry = p[root_times_n_aa + y]
                pr += pry * piry
                dp_index = root_times_n_aa2 + n_aa * y
                dpi_index = n_aa2 * r + y
                for z in nonwt_aa_indices:
                    dpr[z] += self._dpi[dpi_index + z * n_aa] * pry + piry * dp[dp_index + z]
            if pr == 0:
                raise ValueError, "Underflow occurred in the computation for residue %d." % r
            rloglikelihood = math.log(pr) # logarithm of scaled probability
            drloglikelihood = [dprz / pr for dprz in dpr] # convert to derivatives of logarithms
            # Now correct for the scaling
            rloglikelihood += r_log_scaling_sum
        if not (-1e300 < rloglikelihood < 0):
            # print some diagnostic information before raising an exception
            sys.stderr.write("The value of grs is being serialized to '_debug_gr.pickle'.\n")
            sys.stderr.write("Here are the values of pi for this residue:\n")
            for y in range(0, self._n_aa):
                sys.stderr.write("\t%f\n" % self._pi[r * self._n_aa + y])
            cPickle.dump(self._grs[r], open('_debug_gr.pickle', 'w'))
            raise ValueError, "Invalid log likelihood of %s for residue %d." % (str(rloglikelihood), r)
        return (rloglikelihood, drloglikelihood) # return log likelihood



def RescaleDDGs(ddgs, value, parameter, recenter=None, min_max=None):
    """Rescales a set of ddG values.

    On input, 'ddgs' specifies a set of ddG values in the form of a dictionary.
        The keys of this dictionary are integers specifying the residue number.
        'ddgs[i]' is the 2-tuple '(wt, iddgs)' where 'wt' is the identity of 
        residue 'i' and 'iddgs' is a dictionary with 'iddgs[mut]' giving the 
        ddG value for mutating this residue from 'wt' to 'mut'.
    'value' specifies the value that is used to rescale the ddG values.
        They are rescaled so that the descriptive parameter provided by
        'parameter' has a value equal to 'value'.
    'parameter' specifies how the ddG values are rescaled.  Possible values:
        'STANDARD_DEVIATION' means that the ddG values are rescaled so that
            their standard deviation is equal to 'value'.
        'FRAC_NEGATIVE' means that the ddG values are rescaled so that a
            fraction equal to 'value' of them (0 < 'value' < 1) are less
            than zero.  Note that if 'value' is > 0.5, then the median
            ddG value must be negative; and if value < 0.5, then the median
            ddG value must be positive.  An exception is raised if this condition
            is not satisfied.  This option can only be used if 'recenter' is
            not equal to 'None', since it is not possible to change the fraction
            negative without fixing the center.
        '10TH_TO_90TH' means that the ddG values are rescaled
            so that the difference between the 10th and 90th percentile ddG
            values is equal to 'value'.
    'recenter' is an optional parameter that gives the possibility of shifting
        the mean value of the ddG values.  By default, this parameter is 'None',
        which means that the mean of the ddG values is scaled by the same amount
        as all of the ddG values.  If, then the ddG values are all shifted so
        that their mean after rescaling is equal to the value of 'recenter'.
    'min_max' is an optional argument that specifies that we adjust all of the 
        ddG values to fall within some specified range.  IF THIS OPTION IS
        USED, THE SPREAD AND CENTERS SPECIFIED BY THE OTHER PARAMETERS
        MAY NOT BE EXACTLY ACHIEVED.  'min_max' by default is None, meaning
        that we place no restrictions on the range of the ddG values.  If
        'min_max' is set to another value, then it should be the 2-tuple
        '(ddgmin, ddgmax)' where 'ddgmin' and 'ddgmax' are both numbers
        with 'ddgmin < ddgmax'.  In this case, after performing all of
        the recentering and rescaling specified by the other parameters,
        we change any ddG value less than 'ddgmin' to be 'ddgmin', and
        any ddG value greater than 'ddgmax' to be 'ddgmax'.  Because
        this adjusting is done after the recentering and rescaling,
        the ddG values may not end up exactly matching the values
        specified for rescaling and recentering.
    The returned value is a new dictionary 'rescaled_ddgs' in the same format
        as 'ddgs' containing the rescaled ddG values."""
    assert value > 0
    if recenter == None:
        recentered_ddgs = ddgs # no recentering
    elif isinstance(recenter, (int, float)):
        meanddg = stats.Mean(DDGsToList(ddgs))
        shift = recenter - meanddg
        recentered_ddgs = {}
        for (i, (wt, iddgs)) in ddgs.iteritems():
            recentered_iddgs = {}
            for (mut, ddg) in iddgs.iteritems():
                recentered_iddgs[mut] = ddg + shift
            recentered_ddgs[i] = (wt, recentered_iddgs)
    else:
        raise ValueError, "Invalid value of 'recenter': %s" % (str(recenter))
    ddgslist = DDGsToList(recentered_ddgs)
    if parameter == 'STANDARD_DEVIATION':
        scale = float(value) / stats.StandardDeviation(ddgslist)
    elif parameter == 'FRAC_NEGATIVE':
        assert 0 < value < 1
        if recenter == None:
            raise ValueError("Cannot adjust fraction negative without fixing the center.")
        ddgslist.sort()
        median = stats.Median(ddgslist)
        if (value < 0.5 and median < 0) or (value > 0.5 and median > 0) or (value == 0.5 and median != 0):
            raise ValueError("Cannot rescale ddGs with this fraction negative (%f) and median (%f)." % (value, median))
        ddgtomakezero = ddgslist[int(value * len(ddgslist))] # this ddG value should be set to zero
        scale = float(recenter) / (recenter - ddgtomakezero)
        assert scale > 0
    elif parameter == '10TH_TO_90TH':
        ddgslist.sort()
        dist = ddgslist[int(0.9 * len(ddgslist))] - ddgslist[int(0.1 * len(ddgslist))]
        scale = float(value) / dist
    else:
        raise ValueError, "Invalid value of parameter: %s" % parameter
    rescaled_ddgs = {}
    for (i, (wt, iddgs)) in recentered_ddgs.iteritems():
        rescaled_iddgs = {}
        for (mut, ddg) in iddgs.iteritems():
            if recenter != None:
                rescaled_iddgs[mut] = (ddg - recenter) * scale + recenter
            else:
                rescaled_iddgs[mut] = ddg * scale
        rescaled_ddgs[i] = (wt, rescaled_iddgs)
    assert recenter == None or abs(recenter - stats.Mean(DDGsToList(rescaled_ddgs))) < 1e-10
    if min_max != None:
        (ddgmin, ddgmax) = min_max
        assert ddgmin < ddgmax and isinstance(ddgmin, (int, float)) and isinstance(ddgmax, (int, float))
        for (i, (wt, iddgs)) in rescaled_ddgs.iteritems():
            for (mut, ddg) in iddgs.iteritems():
                ddg = max(ddg, ddgmin)
                ddg = min(ddg, ddgmax)
                iddgs[mut] = ddg
    return rescaled_ddgs



def DDGsToList(ddgs):
    """Converts a dictionary specifying ddG values to a list.

    On input, 'ddgs' specifies a set of ddG values in the form of a dictionary.
        The keys of this dictionary are integers specifying the residue number.
        'ddgs[i]' is the 2-tuple '(wt, iddgs)' where 'wt' is the identity of 
        residue 'i' and 'iddgs' is a dictionary with 'iddgs[mut]' giving the 
        ddG value for mutating this residue from 'wt' to 'mut'.
    The returned value is a list of all of these ddG values.
    """
    ddglist = []
    for (i, (wt, iddgs)) in ddgs.iteritems():
        for (mut, ddg) in iddgs.iteritems():
            ddglist.append(ddg)
    return ddglist


def SortedDDGList(ddgs):
    """Gives a sorted list of the ddG values for all mutations to a protein.

    On input, 'ddgs' specifies a set of ddG values in the form of a dictionary.
        The keys of this dictionary are integers specifying the residue number.
        'ddgs[i]' is the 2-tuple '(wt, iddgs)' where 'wt' is the identity of
        residue 'i' and 'iddgs' is a dictionary with 'iddgs[mut]' giving the
        ddG value for mutating this residue from 'wt' to 'mut'.
    The returned value is a sorted list of the 2-tuple '(ddg, mut)' where
        'ddg' is the ddG value for the mutation, and mut is the string
        (i.e. M1A) giving the mutation identities.  They are sorted with
        the most negative (most stabilizing) ddG first.
    """
    sorted_ddgs = []
    for (i, (wt, iddgs)) in ddgs.iteritems():
        for (mut, ddg) in iddgs.iteritems():
            m = "%s%d%s" % (wt, i, mut)
            sorted_ddgs.append((ddg, m))
    sorted_ddgs.sort()
    return sorted_ddgs


def GetDDGs(ddg_dict, mutation_list):
    """Returns a list of ddG values corresponding to certain mutations.
    
    'ddgs_dict' should specify a set of ddG values in the form of a dictionary.
        The keys of this dictionary are integers specifying the residue number.
        'ddgs[i]' is the 2-tuple '(wt, iddgs)' where 'wt' is the identity of 
        residue 'i' and 'iddgs' is a dictionary with 'iddgs[mut]' giving the 
        ddG value for mutating this residue from 'wt' to 'mut'.
    'mutation_list' is a list of mutations for which we want to extract
        the ddG values.  Each entry should specify a mutation in the
        form (wt, r, mut) where 'wt' is the identity of the wildtype
        residue, r is the residue number, and 'mut' is the mutant residue.
        Alternatively, 'mutation_list' can be set to 'ALL' to get all ddG
        values in 'ddg_dict'.  In that case, the ddG values are returned
        for all mutations listed in 'ddg_dict' (the result is the same
        as would be obtained by calling 'DDGsToList' on 'ddg_dict'.
        Alternatively, 'mutation_list' can be a list of numbers.  Each of
        these numbers should then specify a residue number in 'ddg_dict';
        the returned value is a list of the ddGs for all mutations to 
        the residue numbers listed in 'mutation_list'.  If a residue does not
        have ddG values in ddg_dict, nothing is added to the list for
        this residue.
    The function returns a list of the ddG value for each of the mutations
        specified in 'mutation_list'.  The ddG values in the returned list
        are in the same order as the mutations are listed in 'mutation_list'
        (unless 'mutation_list' is 'ALL', in which case the order
        of the ddG values is arbitrary.
    """
    if mutation_list == 'ALL':
        return DDGsToList(ddg_dict)
    ddglist = []
    assert isinstance(mutation_list, list) 
    assert isinstance(ddg_dict, dict)
    if not mutation_list:
        raise ValueError('mutation_list is empty.')
    if isinstance(mutation_list[0], int):
        # we have a list of residue numbers
        for r in mutation_list:
            if r in ddg_dict:
                rddgs = ddg_dict[r][1]
                ddglist += rddgs.values()
        return ddglist
    # if we made it here, we have a list of specific mutations
    for mut_tuple in mutation_list:
        try:
            (wt, r, mut) = mut_tuple
        except ValueError:
            raise ValueError, "Invalid mutation tuple in 'mutation_list': %s." % str(mut_tuple)
        assert isinstance(r, int)
        try:
            (dict_wt, rddgs) = ddg_dict[r]
        except KeyError:
            raise ValueError, "Residue %d cannot be found in 'ddg_dict'." % r
        if dict_wt != wt:
            raise ValueError, "'ddg_dict' and 'mutation_list' define different wildtype residues (%s and %s) for residue %d." % (dict_wt, wt, r)
        try:
            ddg = rddgs[mut]
        except KeyError:
            raise ValueError, "Cannot find ddG for mutation to %s" % mut
        ddglist.append(ddg)
    return ddglist



def ZippedDDGLists(*args):
    """Returns lists of mutations and their ddG values.

    This function takes as input one or more dictionaries specifying ddG values
        for mutations to a protein sequence.  For mutations that have ddG values
        specified in all dictionaries, the function returns a list of these
        mutations as well as the corresponding ddG values.
    Each argument must be a dictionary specifying ddG values.  The keys of the 
        dictionaries are integers 'r' specifying the residue number.  The values
        are 2-tuples of the form '(wt, iddgs)' where 'wt' is the identity
        of the residue, and 'iddgs[mut]' gives the ddG value for mutating this
        this residue from 'wt' to 'mut'.  Each dictionary must be for
        the same sequence (same residue numbers and wildtype identities).   If any
        ddG value is 'None', then it is considered to be undefined.
    This function returns lists of ddG values for mutations that have defined ddG
        values in all dictionaries.  The lists are returned in the form of a tuple:
            '(mutlist, ddgs1list, ddgs2list, ...)'
        All of these lists are of the same length.  'mutlist' lists all of the mutations
        that are present in all ddG dictionaries in the standard form (for example,
        Y96A for a mutation of residue 96 from Tyr to Ala).  'ddgs1list' is a list of
        the ddG values for the first ddG dictionary argument, 'ddgs2list' for the second,
        etc.
    """
    ndicts = len(args)
    if not ndicts:
        raise ValueError("No argument dictionaries.")
    protlength = len(args[0])
    for d in args[1 : ]:
        if protlength != len(d):
            raise ValueError("Dictionaries not all of the same length.")
    amino_acids = AminoAcids()
    return_tuple = tuple([[] for i in range(ndicts + 1)])
    for r in range(protlength):
        wt = args[0][r][0]
        wts = [args[i][r][0] for i in range(ndicts)]
        if len([res for res in wts if res == wt]) != ndicts:
            raise ValueError("Not all dictionaries specify the same wildtype residue for %d." % r)
        for mut in amino_acids:
            if mut == wt:
                continue
            rddgs = [args[i][r][1][mut] for i in range(ndicts)]
            if None in rddgs:
                continue
            return_tuple[0].append("%s%d%s" % (wt, r, mut))
            for i in range(ndicts):
                return_tuple[i + 1].append(rddgs[i])
    return tuple(return_tuple)



def CorrelateDDGs(ddgs1, ddgs2, ddgs_ref = None):
    """Computes the correlation between two sets of ddG values.

    Takes as input two dictionaries, 'ddgs1' and 'ddgs2', specifying
        two sets of ddG values.  The keys of the dictionaries are
        integers 'r' specifying the residue number.  The values are
        2-tuples of the form '(wt, iddgs)' where 'wt' is the 
        identity of the residue, and 'iddgs[mut]' gives
        the ddG value for mutating this residue from 'wt' to 'mut'.
        The two sets of ddG values must be for the same sequence
        (same residue numbers and wildtype identities.  
    If any of the ddG values in ddgs1 and ddgs2 are stored as 'None', then they are
        disregarded in computing the correlation.
    'ddgs_ref' is an optional argument.  By default, it is 'None' meaning
        that nothing is done.  If it has another value, it should be set
        to a dictionary of ddG values of the same format as 'ddgs1' and 'ddgs2'.
        In this case, for each ddG value, if the value is 'None' in 'ddgs_ref',
        then this ddG is disregarded in computing the correlation between
        'ddgs1' and 'ddgs2', even if the the ddG is not 'None' in both
        of these dictionaries.
    The returned value is the tuple '(r, p, n)'.  'r' is the Pearson
        correlation coefficient, 'p' is the two-tailed P value, and
        'n' is the number of data points.
    """
    assert len(ddgs1) == len(ddgs2)
    assert ddgs_ref == None or len(ddgs_ref) == len(ddgs1)
    l1 = []
    l2 = []
    for r in ddgs1.iterkeys():
        assert r in ddgs2
        assert ddgs1[r][0] == ddgs2[r][0] and (ddgs_ref == None or ddgs_ref[r][0] == ddgs1[r][0])
        for mut in AminoAcids():
            if mut == ddgs1[r][0]:
                continue
            ddg1 = ddgs1[r][1][mut]
            ddg2 = ddgs2[r][1][mut]
            if ddgs_ref != None:
                ref = ddgs_ref[r][1][mut]
                if ref == None:
                    continue
            if ddg1 != None and ddg2 != None:
                l1.append(ddgs1[r][1][mut])
                l2.append(ddgs2[r][1][mut])
    return stats.PearsonCorrelation(l1, l2)



def CreateCSet(method):
    """Creates a set of "C" values specifying mutation probabilities among amino acids.

    The returned variable is the dictionary 'c_set'.  'c_set[(x, y)]' specifies the
        probability that a random nucleotide mutation to a codon for amino acid y
        yields amino acid x, where y and x are the one letter amino acid codes.  There
        are several options for creating this set, based on the value of 'method':
        * If 'method == "ALL_EQUAL"', then 'c_set[(x, y)] = 1 / 20. for all 'x'
            and 'y'.
        * If 'method == "CODON_EQUAL", then 'c_set[(x, y)]' is constructed by assuming
            that 'y' is equally likely to be any of its constituent codons.  Each of the
            nine single nucleotide mutatiions to each of these codons are then considered
            to be equally likely, and the resulting probability of transitioning to 'x'
            is computed.  If a mutation to a codon results in a stop codon, then the
            codon stays at its original identity.
        * If 'method' is the 2-tuple '(TRANSITION_TRANSVERSION_RATIO', r)' then
            'c_set' is constructed based on the assumption that all codons for
            an amino acid are equally likely, but that mutations occur with a bias
            in the transitions versus transversions.  'r' gives the transition/transversion
            ratio that quantifies this bias.  If 'r' is 0.5, then there is no bias.
            For real data, 'r' will probably tend to be larger than 0.5 since in most
            cases transitions are more likely than transversions.  For example,
            for influenza, 'r' is probably about 5.0.
    """
    #-----------------------------
    def NDifferences(s1, s2):
        """Returns the number of differences between two strings of the same length."""
        assert len(s1) == len(s2)
        n = 0
        for i in range(len(s1)):
            if s1[i] != s2[i]:
                n += 1
        return n
    #-----------------------------
    genetic_code = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L',
            'CTA':'L', 'CTG':'L', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'GTT':'V',
            'GTC':'V', 'GTA':'V', 'GTG':'V', 'TCT':'S', 'TCC':'S', 'TCA':'S',
            'TCG':'S', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'ACT':'T',
            'ACC':'T', 'ACA':'T', 'ACG':'T', 'GCT':'A', 'GCC':'A', 'GCA':'A',
            'GCG':'A', 'TAT':'Y', 'TAC':'Y', 'TAA':'STOP', 'TAG':'STOP',
            'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'AAT':'N', 'AAC':'N',
            'AAA':'K', 'AAG':'K', 'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
            'TGT':'C', 'TGC':'C', 'TGA':'STOP', 'TGG':'W', 'CGT':'R',
            'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGT':'S', 'AGC':'S', 'AGA':'R',
            'AGG':'R', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}
    transitions = {'A':'G', 'G':'A', 'T':'C', 'C':'T'}
    bases = ['A', 'T', 'G', 'C']
    c_set = {}
    if method == 'ALL_EQUAL':
        for x in AminoAcids():
            for y in AminoAcids():
                c_set[(x, y)] = 1 / 20.
    elif method == 'CODON_EQUAL':
    
        for y in AminoAcids():
            for x in AminoAcids():
                c_set[(x, y)] = 0.0
            ycodons = [codon_aa[0] for codon_aa in genetic_code.iteritems() if codon_aa[1] == y]
            nycodons = len(ycodons)
            assert nycodons
            for ycodon in ycodons:
                for x in AminoAcids() + ['STOP']:
                    xcodons = [codon_aa[0] for codon_aa in genetic_code.iteritems() if codon_aa[1] == x]
                    for xcodon in xcodons:
                        if NDifferences(xcodon, ycodon) == 1:
                            if x == 'STOP':
                                c_set[(y, y)] += 1 / 9. / nycodons
                            else:
                                c_set[(x, y)] += 1 / 9. / nycodons
    elif isinstance(method, tuple) and len(method) == 2 and method[0] == 'TRANSITION_TRANSVERSION_RATIO':
        r = float(method[1]) # transition/transversion ratio
        assert r > 0
        pr_transition = 1.0 / (3.0 + 3.0 / r) # probability that a random mutation to a codon is a specific transition
        pr_transversion = 1.0 / (6.0 * r + 6.0) # probability that a random mutation to a codon is a specific transversion
        assert 0.99999 < 3 * pr_transition + 6 * pr_transversion < 1.000001
        for y in AminoAcids():
            for x in AminoAcids():
                c_set[(x, y)] = 0.0
            ycodons = [codon_aa[0] for codon_aa in genetic_code.iteritems() if codon_aa[1] == y]
            nycodons = len(ycodons)
            for ycodon in ycodons:
                for i_nt in range(len(ycodon)):
                    nt = ycodon[i_nt]
                    for mut_nt in bases:
                        if mut_nt != nt:
                            xcodon = list(ycodon)
                            xcodon[i_nt] = mut_nt
                            xcodon = ''.join(xcodon)
                            x = genetic_code[xcodon]
                            if mut_nt == transitions[nt]:
                                # mutations is a transition
                                p_mut = pr_transition
                            else:
                                # mutation is a transversion
                                p_mut = pr_transversion
                            if x == 'STOP':
                                c_set[(y, y)] += 1.0 / nycodons * p_mut
                            else:
                                c_set[(x, y)] += 1.0 / nycodons * p_mut
    else:
        raise ValueError, "Invalid 'method' of %s." % (str(method))
    if __debug__: # check that sums over x all yield one
        for y in AminoAcids():
            sum = 0.0
            for x in AminoAcids():
                sum += c_set[(x, y)]
            if not (0.99999 < sum < 1.00001):
                raise ValueError, sum
    return c_set



def RestrictDDGsToPartner(ddgs, ddgs_partner):
    """Eliminates all ddG values from one set that are not also defined in another set.

    'ddgs' and 'ddgs_partner' are both dictionaries specifying ddG values for
        mutating a sequence.  'ddgs[r]' gives the ddG values for mutating
        residue 'r' of a sequence (0 <= r < length).  'ddgs[r]' is the
        2-tuple '(wt, iddgs)'.  'wt' is the identity of residue 'r'.
        'iddgs[mut]' gives the ddG value for mutating this residue to
        'mut' (mut != wt) as a number.  However, if the ddG is not
        defined then 'iddgs[mut]' is 'None'.
    The returned value is a copy of 'ddgs' in which any mutations which have a
        ddG value of 'None' in 'ddgs_partner' have been set to 'None'.
    """
    assert len(ddgs) == len(ddgs_partner)
    new_ddgs = {}
    for r in range(len(ddgs)):
        wt = ddgs[r][0]
        assert wt == ddgs_partner[r][0]
        rddgs = ddgs[r][1]
        rddgs_partner = ddgs_partner[r][1]
        newrddgs = {}
        for (mut, ddg) in rddgs.iteritems():
            if rddgs_partner[mut] != None:
                newrddgs[mut] = ddg
            else:
                newrddgs[mut] = None
        new_ddgs[r] = (wt, newrddgs)
    return new_ddgs



def ConsensusDDGs(protseq, sequences, pseudocounts, threshold=1):
    """Estimates ddG values based on "consensus energy" idea.

    This method estimates the ddG values of various mutations based on the
        "consensus energy" idea, in which the ddG of mutating a residue from
        from y to x is estimated to be proportional to -log(nx / ny), where nx
        is the number of appearances of x in the sequence alignment, and ny is 
        the number of appearances of y in the sequence alignment.
    'protseq' is the wildtype protein sequence for which we are estimating the 
        ddG values.
    'sequences' is a list composed of protein sequences.  These sequences should all be string/lists
        of size len(protseq), with each entry either a gap ('-') or a one-letter amino acid 
        code.  Alternatively, sequences can be a list of 2-tuples.  In that case, each 2-tuple
        is assumed to represent '(head, seq)' where 'head' is a sequence name or header that
        is ignored, and 'seq' is the sequence used to calculate the consensus ddGs.
    'pseudocounts' is an integer >= 0.  To compute the consensus ddGs, we count the
        number of appearances of amino acid x at a position among all of the sequences
        in 'sequences'.  To the count for each amino acid, we then add 'pseudocounts'.
    'threshold' is an optional argument specifying that we only estimate ddG values for
        residues that occur at least 'threshold' times once the counts in 'sequences'
        have the pseudocounts added.  By default, it is one.  You can set it be some
        other integer >= 1.
    The returned variable is a dictionary 'ddgs' that holds the ddG values.  The keys in
        this dictionary are the residue numbers 'r' (0 <= r < len(protseq)).  
        'ddgs[r]' is the 2-tuple '(wt, iddgs)'.  'wt' is the wildtype amino acid
        at position 'r', which is just 'protseq[r]'.  'iddgs' is another dictionary
        that gives the estimated ddG value for mutating 'r' from 'wt' to any of the
        other possible amino acids: 'iddgs[mut]' is the ddG value for mutating
        'r' from 'wt' to 'mut'.  If there are not enough counts to estimate
        this ddG value (the counts are less than 'threshold'), then the corresponding
        entry in 'iddgs' is 'None'.  This will be the case if either the number of occurrences 
        of the wildtype or of the mutant amino acid is less than threshold.
    """
    amino_acids = AminoAcids()
    assert isinstance(pseudocounts, int) and pseudocounts >= 0
    assert isinstance(threshold, int) and threshold >= 1
    counts = AAFrequencies(sequences, pseudocounts, raw_counts = True)
    assert len(counts) == len(protseq)
    ddgs = {}
    for r in range(len(protseq)):
        wt = protseq[r]
        nwt = float(counts[r][wt])
        iddgs = {}
        if nwt < threshold: # no occurrences of wildtype, so all ddGs are 'None' for this residue
            for mut in amino_acids:
                if mut == wt:
                    continue
                iddgs[mut] = None
        else: # estimate ddGs for this residue
            for mut in amino_acids:
                if mut == wt:
                    continue
                nmut = counts[r][mut]
                if nmut < threshold:
                    iddgs[mut] = None
                else:
                    iddgs[mut] = -math.log(nmut / nwt)
        ddgs[r] = (wt, iddgs)
    return ddgs


def AAFrequencies(sequences, pseudocounts, raw_counts=False):
    """Function for computing frequencies of amino acids at each position in sequence alignment.

    'sequences' is a list composed of sequences.  These sequences should all be strings/lists
        of the same length, with each entry either a gap ('-') or a one-letter amino acid
        code.  Alternatively, 'sequences' can be a list of 2-tuples of the form '(head, seq)'.
        In this case, 'head' is assumed to a be a header and is ignored, and 'seq' is used
        to calculate the amino acid frequencies.
    'pseudocounts' is an integer specifying the number of pseudocounts added to the tally.
    'raw_counts' is an optional boolean switch.  By default, it is 'False', meaning that
        we compute the frequencies of amino acids at each position.  If 'raw_counts' is
        set to 'True', then instead of computing frequencies, we return the raw number of
        counts of the amino acid (plus the pseudocounts).
    The returned variable is the list 'pi_set', which is the same length as the length
        of the sequences in 'sequences'.  'pi_set[r]' is a dictionary for 0 <= r < length.
        'pi_set[r][aa]' is defined for aa in 'AminoAcids()'.  'pi_set[r][aa]' gives
        the frequency of amino acid 'aa' at residue 'r'.  This frequency is based on the
        counts of 'aa' in the alignment in 'sequences' plus 'pseudocounts', divided by
        the total number of sequences and pseudocounts.  Gaps are ignored in the count.
        'pi_set[r]' sums to one when 'aa' is taken over all values."""
    amino_acids = AminoAcids()
    if not (isinstance(sequences, list) and len(sequences) > 0):
        raise ValueError, "Invalid value of 'sequences'."
    if isinstance(sequences[0], tuple) and len(sequences[0]) == 2:
        with_headers = True
        length = len(sequences[0][1])
    else:
        with_headers = False
        length = len(sequences[0])
    assert isinstance(pseudocounts, int) and pseudocounts >= 0
    pi_set = []
    for r in range(length):
        pi_set.append({})
        n = 0
        for aa in amino_acids:
            pi_set[r][aa] = pseudocounts
            n += pseudocounts
        for headseq in sequences:
            if with_headers:
                seq = headseq[1]
            else:
                seq = headseq
            aa = seq[r]
            if aa == '-':
                pass
            elif aa in amino_acids:
                pi_set[r][aa] += 1
                n += 1
            else:
                raise ValueError, "Unrecognized amino acid of %s." % aa
        if not n:
            raise ValueError, "No counts for residue %d." % r
        if not raw_counts:
            for aa in amino_acids:
                pi_set[r][aa] /= float(n)
    return pi_set



def FillEmptyDDGs(protseq, ddgs, fill_value):
    """Sets all missing ddG values to be some specified value.

    On call, 'protseq' should be a string representing a protein sequence.
    'ddgs' should be a ddG dictionary specifying some ddG values for this protein
        sequence.  'ddgs' is keyed by integers; it can have an entry for each residue
        number r where 0 <= r < len(protseq).  If there is no key for r, then the ddG values
        for mutating this residue are considered to be undefined.  Otherwise, 'ddgs[r]'
        should be the 2-tuple '(wt, rddgs)'.  'wt' should be the wildtype amino acid at
        residue r, equal to 'protseq[r]'.  'rddgs' is a dictionary keyed by all amino
        acids (one-letter symbols) except for 'wt'.  'rddgs[mut]' gives the ddG for mutating
        residue r from wt to mut.  Alternatively, 'rddgs[mut]' can be 'None', meaning that 
        the ddG value is undefined for this mutation.
    'fill_value' specifies the value that every undefined ddG value is set to.  It can
        be a number, in which case all undefined ddG values are set to this number.  It
        can be 'None', in which case all undefined ddG values are set to be 'None'.  Or
        it can be the string 'MEAN', in which case all undefined ddG values are set to the 
        mean of all of the defined ddG values in 'ddgs'.
    This function returns a new ddG dictionary called 'filled_ddgs'.  'filled_ddgs[r]' is
        defined for all 'r' in 'protseq', and all of the previously undefined ddG values
        have been set to 'fill_value'.
    """
    amino_acids = AminoAcids()
    filled_ddgs = {}
    if isinstance(fill_value, (int, float)):
        pass # valid value for 'fill_value'
    elif fill_value == None:
        pass # valid value for 'fill_value'
    elif fill_value == 'MEAN':
        ddgslist = DDGsToList(ddgs)
        fill_value = stats.Mean(ddgslist) # fill value is the mean ddG
    else:
        raise ValueError("Invalid value for fill_value: %s" % str(fill_value))
    for r in range(len(protseq)):
        wt = protseq[r]
        if r in ddgs:
            if wt != ddgs[r][0]:
                raise ValueError("'ddgs' and 'protseq' are not compatible at %d: the former has %s, while the latter has %s." % (r, wt, ddgs[r][0]))
            rddgs = ddgs[r][1].copy()
            for mut in amino_acids:
                if mut != wt and rddgs[mut] == None:
                    rddgs[mut] = fill_value
            filled_ddgs[r] = (wt, rddgs)
        else:
            rddgs = dict([(mut, fill_value) for mut in amino_acids if mut != wt])
            filled_ddgs[r] = (wt, rddgs)
    return filled_ddgs



def WriteDDGs(ddgs, filename, datetime):
    """Writes a set of ddG values to a text file.

    'ddgs' is a dictionary containing ddG values.  It is keyed by
        residue numbers 'r'.  For each residue 'r', 'ddgs[r]' is
        the 2-tuple (wt, ddgdict) where 'wt' is the wildtype residue
        at position 'r' and 'ddgdict[mut]' gives the ddG value for
        mutating residue 'r' from 'wt' to 'mut' as a number.  It can
        also be None if no ddG is defined for this mutation.
    'filename' is a string giving the name of the file to which we write
        the ddG values.  If the file already exists, it is overwritten.
        If it does not yet exist, it is created.  The ddGs are written
        in a format that should be self-explanatory -- they can be read
        back out of 'filename' using the method 'ReadDDGs'.
    'datetime' is a string that specifies the date and time at which the
        ddG values were stored.  This is of use if you want to keep track
        of when these values were obtained.  Typically, you
        can get the current 'datetime' by calling 'time.asctime()'.
    """
    f = open(filename, 'w')
    f.write("%s\n" % datetime)
    f.write("MUTATION\tDDG\n")
    for (i, (wt, iddgs)) in ddgs.iteritems():
        for (mut, ddg) in iddgs.iteritems():
            if ddg == None:
                f.write("%s%d%s\tNone\n" % (wt, i, mut))
            else:
                f.write("%s%d%s\t%f\n" % (wt, i, mut, ddg))
    f.close()



def ReadDDGs(filename):
    """Reads a set of ddG values from an existing text file.

    The file should have been created using the method 'WriteDDGs',
        and should have the name 'filename'.
    The returned value is a 2-tuple '(datetime, ddgs)' that specifies the ddG 
        values in 'filename' and the date/time that these ddG values were stored
        in the file.  The rationale behind the date/time is that it might be useful
        to know when the ddG values were obtained.
        'datetime' is a string that gives the date and time.  'ddgs' is a
        dictionary in the format described for 'WriteDDGs'.
    """
    if not os.path.isfile(filename):
        raise ValueError, "Cannot find the file %s which is supposed to contain the ddGs." % filename
    lines = open(filename).readlines()
    datetime = lines[0].strip()
    ddgs = {}
    for line in lines[2 : ]:
        (mutation, ddg) = line.split()
        mutation = mutation.strip()
        if ddg == 'None':
            ddg = None
        else:
            ddg = float(ddg)
        (wt, i, mut) = (mutation[0], int(mutation[1 : -1]), mutation[-1])
        if i not in ddgs:
            ddgs[i] = (wt, {})
        if wt != ddgs[i][0]:
            raise ValueError, "Inconsistent identity for residue %d in %s." % (i, filename)
        if mut in ddgs[i][1]:
            raise ValueError, "Duplicate ddG values for mutations %s%d%s in %s." % (wt, i, mut, filename)
        ddgs[i][1][mut] = ddg
    return (datetime, ddgs)
