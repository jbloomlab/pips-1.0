#!python

"""Script to peform a PIPS analysis to predict the ddG values from a protein phylogeny.

Written by Jesse Bloom, 2009.

This script expects to find an input file with the name PIPS_ANALYSIS.IN
in the current directory.  Using the configuration specified by this file,
the script creates two or more output files.  The first of these files,
PIPS_ANALYSIS.LOG is a log of the analysis progress.  The remaining files
contain the PIPS predictions for different priors, and have the 
names PIPS_PREDICTIONS_PRIORS-XXX.TXT where XXX is a particular prior
description that is specified in the input file (under the PRIORS key).  If any of these 
output files already exists, the program asks for user confirmation before continuing
and overriding the values.

All residues are numbered sequentially according to the protein sequence
being analyzed, with the first residue being number 1.

The input file, PIPS_ANALYSIS.IN should have the format as detailed below.
Anything come after a "#" character is considered a comment, and has no effect.
Format of PIPS_ANALYSIS.IN file:

-----------
# Anything coming after a pound sign is a comment.
# Typically there may be one or more comment lines explaining the
# subject matter of the input file.
#
# Each line in this file follows a key-value format.  The first word
# in the line gives the parameter that is being specified, and the subsequent
# words give the value of this parameter. 
#
# The PROTSEQ key specifies the sequence of the protein to be analyzed.  Typically, this
# will be a FASTA file in which the sequence of interest is the first entry, as shown below.
PROTSEQ protseq.fasta
#
# The ALIGNMENT key specifies the multiple sequence alignment of the homologous proteins.
# This alignment should be given by a FASTA file.
# The first protein sequence in this alignment must match that specified by PROTSEQ, with
# the possible addition of some gap characters.  Once the alignment has been read, positions
# corresponding to gaps in this sequence are stripped away from all sequences.  In addition,
# the names of the sequences must match those contained in the phylogenetic tree specified
# by the tree parameter.  The name of a sequence is considered to consist of all characters
# from the beginning of the header until the first colon.  For example, if the sequences are
# named "SEQ1" and "SEQ2", the file might look like this:
# >SEQ1: this sequence matches that given by PROTSEQ
# MT-DYST-FFT
# >SEQ2: second sequence
# MTFEFSTN-FS
ALIGNMENT renamed_alignment.fasta
#
# The TREE key specifies the phylogenetic tree giving the evolutionary relationships of
# the sequences in the multiple sequence alignment.  Typically, this tree is specified
# in a text file in the Newick format (http://evolution.gs.washington.edu/phylip/newick_doc.html).
# Branch lengths must be specified for all branches; the tree must be rooted; the sequence
# names must match those of the multiple sequence alignment specified by ALIGNMENT; there must
# be no spaces or quotations in the sequence names.
TREE tree.newick
#
# The PRIORS key specifies the centers for the prior probability distributions over the ddG
# values.  Each PIPS analysis can be performed with multiple different choices for the prior
# centers, and a different output file with the different predictions is output for each choice
# of the priors.  The priors themselves are beta distributions with the sum of the shape parameters
# specified by PRIORS_BETA_SUM.  The general format for a specification of the prior centers is
#    PRIORS prior_name prior_type
# Here prior_name indicates a string (no spaces allowed!) giving a short description of
# the prior centers.  The output file will have the name PIPS_PREDICTIONS_PRIORS-prior_name.TXT.
# So for example, if prior_name is set to REGULARIZING then the output file for these predictions
# will be PIPS_PREDICTIONS_PRIORS-REGULARIZING.TXT.
# The entry at prior_type indicates how the prior centers are chosen.  If you want all of the 
# priors to be centered at the same value, then prior_type should simply be that number.  A typical
# value would be 5, which specifies regularizing priors peaked at 5.
# If you want the priors to be proportional to the magnitude of the difference in the
# amino acid hydrophobicities (as determined by the Kyte-Doolittle scale), then the value
# of prior_type should be: KYTE_DOOLITTLE_HYDROPHOBICITY 1 0
# where the last two numbers specify respectively a scaling factor by which the magnitude of 
# the hydrophobicity difference is multiplied, and a factor that is then added to this
# rescaled value.  Typically values would simply be 1 and 0 as shown above.
# If you want to specify a specific prior estimate for each ddG value (for example, these estimates
# might be derived from some physicochemical modeling program), then you need to specify a file
# containing these prior estimates, as in:
#  FILE ddg_priors.txt rescale
# where rescale specifies how (if at all) the specified ddG values are rescaled.  If rescale
# is equal to the string NO_RESCALE then no rescaling is done.  Otherwise, rescale
# should be a set of four numbers:
#   10TH_TO_90TH RECENTER MIN MAX DEFAULT
# The ddG values are rescaled so that the difference between the 10th and 90th percentile
# values is equal to 10TH_TO_90TH, the mean value is equal to RECENTER, the minimum ddG value
# is equal to MIN, and the maximum ddG value is equal to MAX.  The ddG values are first rescaled
# according to 10TH_TO_90TH, then recentered according to RECENTER, and then have the extreme
# values adjusted according to MIN and MAX.  Finally, any mutations for which no ddG value is
# specified are set to the value of DEFAULT
# The file (ddg_priors.txt) itself should be a text file of the format that is output by this
# script as well as pips_run_foldx.py and pips_run_cupsat.py.  
PRIORS REGULARIZING 5
PRIORS HYDROPHOBICITY KYTE_DOOLITTLE_HYDROPHOBICITY 1 0
PRIORS CUPSAT FILE CUPSAT_PREDICTIONS.TXT 10 5 -5 15 5
PRIORS FOLDX FILE FOLDX_PREDICTIONS.TXT 10 5 -5 15 5
#
# PRIORS_BETA_SUM is a number > 2 specifying the sum of the beta distribution parameters.  The
# typical value is 3, as in:
PRIORS_BETA_SUM 3
# Larger values lead to "tighter" priors while smaller values lead to "looser" priors.
#
# The MUTPROBS key specifies how the mutation probabilities are set.  If you want it to be
# equally probable for any amino acid to mutate to any other, use the value ALL_EQUAL, as in
# MUTPROBS ALL_EQUAL
# If you want to assume that each amino acid is equally likely to be in any of its codons,
# and that these codons are equally likely to mutate to any other one-mutation neighbor codon,
# the use the value CODON_EQUAL, as in
# MUTPROBS CODON_EQUAL
# If you want to assume that each amino acid is equally likely to be in any of its codons,
# and that these codons mutate to adjacent codons according to a particular nucleotide
# transition-transversion bias, then use the value TRANSITION_TRANSVERSION_RATIO r (where r
# is the ratio), as in
MUTPROBS TRANSITION_TRANSVERSION_RATIO 4
#
# The MUTRATE key specifies the "mutation rate," which is meaningful with respect to the branch
# lengths of the phylogenetic tree.  Essentially, it gives the number of actual mutations that
# are assumed to have occurred for each unit change along the phylogenetic tree.  You might
# think of this as the number of mutations that occur before the occurrence of an evolutionarily
# "acceptable" one.  Reasonable values would probably be in the range from 5 to 20.
MUTRATE 10
#
# The NRANDOMSTARTS key specifies how many times we perform the likelihood maximization for each
# residue (we take the best value).  Typically, just one is adequate since the ddG values usually
# converge to essentially the same maximum, but you can try values larger than one just to be
# safe and increase the chances of hitting the global maximum.  Increasing the value will of course
# increase the programs runtime.
NRANDOMSTARTS 5
----------

"""

import sys
import os 
import time
import random
import math
import traceback
import pips.version
import pips.align
import pips.ddg_inference
import pips.fasta
import pips.tree
import pips.stats


class FileWriteExceptHook(object):
    """Custom exception handler that also prints exception information to file-like object.

    Objects are initialized by calling with a writeable file-like object.  
    The method 'ExceptHook' then duplicates exception information to this object.
    """

    def __init__(self, f):
        """f is the writeable file-like object to which the exception information is written."""
        self.f = f

    def ExceptHook(self, typ, value, tb):
        """Writes the exception information to 'f' before writing to standard error."""
        lines = traceback.format_exception(typ, value, tb)
        for line in lines:
            self.f.write(line)
        sys.__excepthook__(typ, value, tb)


def QueryFileOverwrite(filename, filedescription):
    """Asks the user whether they want to overwrite an existing file.

    This function will ask the user whether they want to overwrite an existing
        file.  It returns True if the user selects to overwrite the existing file,
        and False if they select not to overwrite the file.  The query is done by
        printing a message to standard output and waiting for the response
        to be entered at standard input.
    'filename' should give the name of the file that we might overwrite.
    'filedescription' should be a short string describing the file.
    """
    ans = None
    while ans not in ['Y', 'y', 'N', 'n']:
        ans = raw_input("The %s file (%s) already exists.  Should we overwrite this existing file [Y/N]? " % (filedescription, filename))
        ans = ans.strip()
        if ans in ['N', 'n']:
            return False
        elif ans in ['Y', 'y']:
            return True


def main():
    """The body of the script, performs the ddG analysis."""
    # Seed the random number generator to make output reproduceable
    random.seed(1)

    # Define input/output files
    firstresnum = 1 # number assigned to the first residue
    logfilename = 'PIPS_ANALYSIS.LOG'
    predictionsfileprefixsuffix = 'PIPS_PREDICTIONS_PRIORS-%s.TXT'
    infilename = 'PIPS_ANALYSIS.IN'
    if os.path.isfile(logfilename):
        if QueryFileOverwrite(logfilename, 'log'):
            print "Deleting existing %s file." % logfilename
            os.remove(logfilename)
        else:
            print "Canceling program."
            sys.exit()
    logfile = open(logfilename, 'w')

    # Make sure exception information is written to logfile
    exceptionwriter = FileWriteExceptHook(logfile)
    sys.excepthook = exceptionwriter.ExceptHook

    # print introductory information
    logfile.write("Beginning execution of pips_analysis.py at %s.\n" % time.asctime())
    logfile.write("Using version %s of PIPS, written by Jesse Bloom.\n" % pips.version.version)
    logfile.write("Using input/output files in the current directory, which is %s\n\n" % os.getcwd())
    logfile.write("Reading input from %s.\n" % infilename)

    # Read input data into the dictionary input_dict, which is keyed by input keywords and has
    # the input data as values
    try:
        input = open(infilename).read()
    except IOError:
        if not os.path.isfile(infilename):
            raise IOError("Cannot find pips_analysis.py input file, %s." % infilename)
        else:
            raise
    logfile.write("Contents of input file are listed below:\n\n--------\n%s\n-------\n\n" % input)
    input = [line for line in input.split('\n') if line and (not line.isspace()) and line[0] != '#']
    input_dict = {'PROTSEQ':None, 'ALIGNMENT':None, 'TREE':None, 'PRIORS':None, 
                  'PRIORS_BETA_SUM':None, 'MUTPROBS':None, 'MUTRATE':None, 'NRANDOMSTARTS':None}
    for line in input:
        entries = line.split()
        if entries[0] not in input_dict:
            raise ValueError("Input file contains invalid key of %s in line %s." % (entries[0], line))
        elif input_dict[entries[0]] != None and entries[0] != 'PRIORS':
            raise ValueError("Input file contains duplicate entry for key %s in line %s." % (entries[0], line))
        elif entries[0] == 'PROTSEQ':
            if len(entries) != 2:
                raise ValueError("Improper specification of PROTSEQ.  The correct format is:\nPROTSEQ protseq.fasta\nThe input file contained the following specification:\n%s" % line)
            if not os.path.isfile(entries[1]):
                raise IOError("Cannot find the specified PROTSEQ file %s." % entries[1])
            logfile.write("\nReading PROTSEQ protein sequence from %s...\n" % entries[1])
            protseq = pips.fasta.Read(entries[1])
            logfile.write("This file contains %d sequence(s).  The first one is being used.\n" % len(protseq))
            (protseqheader, protseq) = protseq[0]
            logfile.write("This protein sequence is as follows:\n")
            pips.fasta.Write([(protseqheader, protseq)], logfile, writable_file=True)
            input_dict['PROTSEQ'] = protseq
        elif entries[0] == 'ALIGNMENT':
            if len(entries) != 2:
                raise ValueError("Improper specification of ALIGNMENT.  The correct format is:\nALIGNMENT renamed_alignment.fasta\nThe input file contained the following specification:\n%s" % line)
            if not os.path.isfile(entries[1]):
                raise IOError("Cannot find the specified ALIGNMENT file %s." % entries[1])
            logfile.write("\nReading ALIGNMENT aligned protein sequences from %s...\n" % entries[1])
            alignment = pips.fasta.Read(entries[1])
            logfile.write("This file contains %d sequences.\n" % len(alignment))
            logfile.write("The aligned sequences are as follows:\n")
            pips.fasta.Write(alignment, logfile, writable_file=True)
            input_dict['ALIGNMENT'] = alignment
        elif entries[0] == 'TREE':
            if len(entries) != 2:
                raise ValueError("Improper specification of TREE.  The correct format is:\nTREE tree.newick\nThe input file contained the following specification:\n%s" % line)
            if not os.path.isfile(entries[1]):
                raise IOError("Cannot find the specified TREE file %s." % entries[1])
            logfile.write("\nReading TREE phylogenetic tree from %s...\n" % entries[1])
            newick_tree = open(entries[1]).read()
            newick_tree = newick_tree.replace('\n', '').replace('\t', '').replace(' ', '')
            logfile.write("The phylogenetic tree is as follows:\n%s\n" % newick_tree)
            input_dict['TREE'] = newick_tree
        elif entries[0] == 'PRIORS':
            if len(entries) <= 2:
                raise ValueError("PRIORS must contain at least three entries in line:\n%s" % line)
            prior_name = entries[1]
            predictionsfilename = predictionsfileprefixsuffix % prior_name
            logfile.write("\nReading from PRIORS the prior centers for the predictions that will be written to file %s...\n" % predictionsfilename)
            if os.path.isfile(predictionsfilename):
                if QueryFileOverwrite(predictionsfilename, 'predictions'):
                    logfile.write("Deleting existing %s file." % predictionsfilename)
                    print "Deleting existing %s file." % predictionsfilename
                    os.remove(predictionsfilename)
                else:
                    print "We will continue without performing a new PIPS calculation for %s." % predictionsfilename
                    logfile.write("The file %s already exists, and the user does not want to overwrite this file.  We will continue, but without performing a new calculation for this file.  New calculations may still be performed for other PRIORS values.\n" % predictionsfilename)
                    continue
            if len(entries) == 3:
                try:
                    prior_centers = float(entries[2])
                except ValueError:
                    raise ValueError("Cannot converts priors center to number in line: %s" % line)
                logfile.write("The prior estimates for all ddG values in this file will be peaked at %f.\n" % prior_centers)
            elif len(entries) == 5 and entries[2] == 'KYTE_DOOLITTLE_HYDROPHOBICITY':
                try:
                    (priors_scale, priors_shift) = (float(entries[3]), float(entries[4]))
                except ValueError:
                    raise ValueError("Cannot convert PRIORS scale and shift values to numbers in line: %s" % line)
                logfile.write("\nThe prior estimates for ddG values in this file will be equal to the magnitude of the difference in residue hydrophobicity (Kyte-Doolittle scale) multiplied by %f and then added to %f.\n" % (priors_scale, priors_shift))
                prior_centers = ('KYTE_DOOLITTLE_HYDROPHOBICITY', priors_scale, priors_shift)
            elif len(entries) == 9 and entries[2] == 'FILE':
                priors_filename = entries[3]
                if not os.path.isfile(priors_filename):
                    raise IOError("Cannot find the file %s specified by the PRIORS keyword in line: %s" % (priors_filename, line))
                try:
                    (prior_10_90, prior_recenter, prior_min, prior_max, prior_default) = (float(entries[4]), float(entries[5]), float(entries[6]), float(entries[7]), float(entries[8]))
                except ValueError:
                    raise ValueError("Cannot convert all five numbers specified in line: %s" % line)
                logfile.write("\nThe prior estimates for the ddG values in this file will be given the values specified in the file %s.  All values will be rescaled so that the difference between the 10th and 90th percentile is %f, then recentered so that the average value is %f, and then truncated so that the minimum and maximum values are %f and %f, respectively.  Values not specified in this file will be given the default value of %f.  Residues are assumed to be numbered sequentially with the first residue assigned the number %d.  Here are the contents of the file specifying the ddG values:\n\n-------\n%s\n-------\n" % (priors_filename, prior_10_90, prior_recenter, prior_min, prior_max, prior_default, firstresnum, open(priors_filename).read()))
                try:
                    (initial_ddgs_description, initial_ddgs) = pips.ddg_inference.ReadDDGs(priors_filename)
                except:
                    logfile.write("\nERROR ENCOUNTERED READING ddG VALUES FROM %s." % priors_filename)
                    sys.stderr.write("\nERROR ENCOUNTERED READING ddG VALUES FROM %s." % priors_filename)
                    raise
                logfile.write("Now rescaling these ddG values so that the difference between the 10th and 90th percentile is %f, then so that the average is %f, and then truncating extreme values to minimum/maximum of %f and %f.\n" % (prior_10_90, prior_recenter, prior_min, prior_max))
                initial_ddgs = pips.ddg_inference.RescaleDDGs(initial_ddgs, prior_10_90, '10TH_TO_90TH', recenter=prior_recenter, min_max=(prior_min, prior_max))
                prior_centers = ('INITIAL', initial_ddgs, prior_default)
            else:
                raise ValueError("Cannot parse the PRIORS value specified in line: %s" % line)
            if input_dict['PRIORS'] == None:
                input_dict['PRIORS'] = {}
            if predictionsfilename in input_dict['PRIORS']:
                raise ValueError("Duplicate predictionsfilename specified by PRIORS in line:\n%s" % line)
            input_dict['PRIORS'][predictionsfilename] = prior_centers
        elif entries[0] == 'PRIORS_BETA_SUM':
            if len(entries) != 2:
                raise ValueError("PRIORS_BETA_SUM has additional entries; correct format is: PRIORS_BETA_SUM 3")
            try:
                priors_beta_sum = float(entries[1])
            except ValueError:
                raise ValueError("Cannot convert beta distribution parameter sum to number in line: %s" % line)
            if priors_beta_sum <= 2:
                raise ValueError("Beta distribution parameter sum does not exceed two in line: %s" % line)
            input_dict['PRIORS_BETA_SUM'] = priors_beta_sum
            logfile.write('\nThe PRIORS_BETA_SUM sum of the beta distribution for the prior probability estimates is %f.' % priors_beta_sum)
        elif entries[0] == 'MUTPROBS':
            if len(entries) == 2 and entries[1] == 'CODON_EQUAL':
                logfile.write("\nMUTPROBS has the value of CODON_EQUAL, specifying that any codon is considered equally likely to mutate to any other codon.\n")
                input_dict['MUTPROBS'] = 'CODON_EQUAL'
            elif len(entries) == 2 and entries[1] == 'ALL_EQUAL':
                logfile.write("\nMUTPROBS has the value of ALL_EQUAL, specifying that any amino is considered equally likely to mutate to any other amino acid.\n")
                input_dict['MUTPROBS'] = 'ALL_EQUAL'
            elif len(entries) == 3 and entries[1] == 'TRANSITION_TRANSVERSION_RATIO':
                try:
                    r = float(entries[2])
                except ValueError:
                    raise ValueError("Cannot convert specified TRANSITION_TRANSVERSION_RATIO to a number in line: %s" % line)
                logfile.write("\nMUTPROBS has a value of TRANSITION_TRANSVERSION_RATIO with a ratio of %f, meaning that codons mutate with a transition / transversion ratio of %f.\n" % (r, r))
                input_dict['MUTPROBS'] = ('TRANSITION_TRANSVERSION_RATIO', r)
            else:
                raise ValueError("Invalid MUTPROBS value specified in line: %s" % line)
        elif entries[0] == 'MUTRATE':
            if len(entries) == 2:
                try:
                    mutrate = float(entries[1])
                except ValueError:
                    raise ValueError("Cannot convert specified MUTRATE to a number in line: %s" % line)
                logfile.write("\nMUTRATE specifies a mutation rate of %f mutations per unit branch length.\n" % mutrate)
                input_dict['MUTRATE'] = mutrate
            else:
                raise ValueError("Invalid MUTRATE value specified in line: %s" % line)
        elif entries[0] == 'NRANDOMSTARTS':
            if len(entries) == 2:
                try:
                    nstarts = int(entries[1])
                except ValueError:
                    raise ValueError("Cannot convert specified NRANDOMSTARTS to an integer in line: %s" % line)
                logfile.write("\nNRANDOMSTARTS specifies that we perform %d different starts for the residue likelihood maximization.\n" % nstarts)
                input_dict['NRANDOMSTARTS'] = nstarts - 1 # one start is assumed by program
            else:
                raise ValueError("Invalid NRANDOMSTARTS value specified in line: %s" % line)
        else:
            raise IOError("Problem parsing input file keys.")
    logfile.flush()

    # We have parsed the entire input file.  Make sure we have a value for each key:
    for (inputkey, inputvalue) in input_dict.iteritems():
        if inputvalue == None:
            raise ValueError("Failed to read an input value for the required key %s." % inputkey)

    # Construct the phylogenetic tree object as 'tree'
    logfile.write("\nConstructing the phylogenetic tree...\n")
    logfile.flush()
    input_dict['ALIGNMENT'] = pips.fasta.UnknownsToGaps(input_dict['ALIGNMENT']) # replace unknown amino acids with gaps
    input_dict['ALIGNMENT'] = pips.align.StripGapsToFirstSequence(input_dict['ALIGNMENT'])
    if input_dict['ALIGNMENT'][0][1] != input_dict['PROTSEQ']:
        raise ValueError("First protein sequence in ALIGNMENT does not match the sequence of PROTSEQ.")
    try:
        tree = pips.tree.Tree(input_dict['TREE'], tipnames_sequences=input_dict['ALIGNMENT'])
    except:
        logfile.write("\nERROR constructing the phylogenetic tree.  Most likely cause is either incorrect newick string for TREE, or ALIGNMENT failing to contain unique sequence names that match those in TREE.\n")
        sys.stderr.write("\nERROR constructing the phylogenetic tree.  Most likely cause is either incorrect newick string for TREE, or ALIGNMENT failing to contain unique sequence names that match those in TREE.\n")
        raise

    # construct the DDGSet object as 'ddgset', and perform maximization for each priors choice
    logfile.write("\nWe will make PIPS predictions using %d different sets of priors.\n" % (len(input_dict['PRIORS'])))
    iprior = 1
    for (predictionsfilename, priors) in input_dict['PRIORS'].iteritems():
        logfile.write("\nConstructing the ddG set for the set of priors that will have its predictions written to %s (this is prior set %d of %d)...\n" % (predictionsfilename, iprior, len(input_dict['PRIORS'])))
        iprior += 1
        logfile.flush()
        if isinstance(priors, tuple) and priors[0] == 'INITIAL':
            (initial_ddgs, initial_defaultvalue) = (priors[1], priors[2])
            renumbered_initial_ddgs = {} # renumbered so first residue is zero rather than firstresnum
            for (r, (rwt, rddgs)) in initial_ddgs.iteritems():
                renumbered_initial_ddgs[r - firstresnum] = (rwt, rddgs)
            try:
                renumbered_initial_ddgs = pips.ddg_inference.FillEmptyDDGs(input_dict['PROTSEQ'], renumbered_initial_ddgs, initial_defaultvalue)
            except:
                logfile.write("\nERROR setting up the ddG values specified by PRIORS.  Most likely cause is that the ddG values specified by PRIORS are for a sequence that does not match PROTSEQ.\n")
                sys.stderr.write("\nERROR setting up the ddG values specified by PRIORS.  Most likely cause is that the ddG values specified by PRIORS are for a sequence that does not match PROTSEQ.\n")
                raise
            initial_values = ('SPECIFIED', renumbered_initial_ddgs, 0, 0)
            priors = priors[0]
        else:  # set all initial ddG values to 5
            initial_values = ('CONSTANT', 5.0)
        # now create the ddG set
        ddgset = pips.ddg_inference.DDGSet(protseq=input_dict['PROTSEQ'], 
                    treedata=tree, 
                    c_set=input_dict['MUTPROBS'], 
                    initial_values=initial_values, 
                    ddgs_prior=('BETA', input_dict['PRIORS_BETA_SUM'], priors), 
                    mutrate=input_dict['MUTRATE'])
        # Now perform the maximization to obtain the ddG estimates
        logfile.write("\nNow beginning the maximization of posterior probability to estimate the ddG values from the protein phylogeny at %s...\n" % time.asctime())
        logfile.flush()
        try:
            old_stdout = sys.stdout # change standard output to logfile for printing of progress
            sys.stdout = logfile
            start_time = time.time()
            ddgset.MaximizePosterior(nrandomstarts=input_dict['NRANDOMSTARTS'], printprogress=True)
            end_time = time.time()
        finally:
            sys.stdout = old_stdout
        logfile.write("\nCompleted maximization of posterior probability at %s.  Total elapsed time of computation was %d seconds.\n" % (time.asctime(), end_time - start_time))
        # write the estimate ddG values to predictionsfilename
        logfile.write("\nThe PIPS estimated ddG values are now being written to the file %s...\n" % predictionsfilename)
        datetime = 'PIPS estimated ddG values computed using PIPS version %s at %s in the directory %s.' % (pips.version.version, time.asctime(), os.getcwd())
        pips_ddgs = ddgset.DDGDict()
        logfile.write("These PIPS ddGs are numbered sequentially with the first residue having number %d.\n" % firstresnum)
        renumbered_pips_ddgs = {} # renumber so first residue is firstresnum rather than zero
        for (r, (rwt, rddgs)) in pips_ddgs.iteritems():
            renumbered_pips_ddgs[r + firstresnum] = (rwt, rddgs)
        pips.ddg_inference.WriteDDGs(renumbered_pips_ddgs, predictionsfilename, datetime)
        logfile.write("These PIPS estimated ddG values are as shown below:\n\n--------\n%s\n--------\n" % open(predictionsfilename).read())

    # Program is complete.
    logfile.write("\nCompleted execution of pips_analysis.py at %s.\n" % time.asctime())
    logfile.close()


# run the script
main()
