#!python

"""Script to analyze the predicted ddG values for selected mutations.

Written by Jesse Bloom, 2009.

This script expects to find an input file with the name PIPS_CORRELATE_SELECTED_MUTATIONS.IN
in the current directory.  Using the configuration specified by this file, the script
creates two output file, plus possibly some graphs.  The first output file,
PIPS_CORRELATE_SELECTED_MUTATIONS.LOG is a log of the scripts progress.  The second
output file, PIPS_CORRELATE_SELECTED_MUTATIONS.TXT contains the results of
the analysis.  

Essentially, this script allows the user to specify certain mutations with known ddG values
It then examines the predicted ddG values for each of these mutations, and determines
how the predicted ddG values correlate with the known ones.

All residues are numbered sequentially according to the protein sequence
being analyzed, with the first residue being number 1.

The input file, PIPS_CORRELATE_SELECTED_MUTATIONS.IN should have the format as detailed below.
Anything come after a "#" character is considered a comment, and has no effect:

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
# The SELECTED_MUTATION key specifies a selected mutation with a known ddG value
# These mutations are indicated in the standard form (i.e. M2A for mutation of residue 2
# from Met to Ala).  The numbers for these mutations are based on sequential numbering
# starting with one for the first residue of PROTSEQ.  The mutations are followed by
# a numbering giving the ddG value.  One or more mutations can be
# specified by listing the key repeatedly with a single mutation following, as in:
#  SELECTED_MUTATION M2A -0.3
#  SELECTED_MUTATION F25I 5.7
# Alternatively, these mutations can be listed in a file.  In this case, each line of the file
# should list a separate mutation and associated ddG in the form "M2A -0.3", etc.  
# Lines of the file that begin with #
# and anything after a # character are considered to be comments, and are ignored.  Values
# from a file are specified by using the second key FILE and then the file name, as in:
SELECTED_MUTATION FILE known_ddG_values.txt
#
# The PREDICTION_FILE key specifies the names of one or more files containing the
# ddG predictions that we wish to gather for the selected mutations.  These files should
# be of the format that is output by PIPS_ANALYSIS.PY.  Each key should specify a string
# (no spaces) giving the name of the prediction method, plus the name of the prediction file.
# Spaces in the name of the prediction file should be represented as underscores; these are
# subsequently changed to spaces in the output.
PREDICTION_FILE CUPSAT CUPSAT_PREDICTIONS.TXT
PREDICTION_FILE PIPS_regularizing_priors PIPS_PREDICTIONS_PRIORS-REGULARIZING.TXT
#
----------

"""

import sys
import os 
import time
import random
import math
import traceback
import pips.version
import pips.ddg_inference
import pips.fasta
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
    logfilename = 'PIPS_CORRELATE_SELECTED_MUTATIONS.LOG' 
    outputfilename = 'PIPS_CORRELATE_SELECTED_MUTATIONS.TXT'
    infilename = 'PIPS_CORRELATE_SELECTED_MUTATIONS.IN'
    for (filename, filedescription) in [(logfilename, 'log'), (outputfilename, 'output')]:
        if os.path.isfile(filename):
            if QueryFileOverwrite(filename, filedescription):
                print "Deleting existing %s file." % filename
                os.remove(filename)
            else:
                print "Canceling program."
                sys.exit()
    logfile = open(logfilename, 'w')

    # Make sure exception information is written to logfile
    exceptionwriter = FileWriteExceptHook(logfile)
    sys.excepthook = exceptionwriter.ExceptHook

    # print introductory information
    logfile.write("Beginning execution of pips_correlate_selected_mutations.py at %s.\n" % time.asctime())
    logfile.write("Using version %s of PIPS, written by Jesse Bloom.\n" % pips.version.version)
    logfile.write("Using input/output files in the current directory, which is %s\n\n" % os.getcwd())
    logfile.write("Reading input from %s.\n" % infilename)

    # Read input data into the dictionary input_dict, which is keyed by input keywords and has
    # the input data as values
    try:
        input = open(infilename).read()
    except IOError:
        if not os.path.isfile(infilename):
            raise IOError("Cannot find input file, %s." % infilename)
        else:
            raise
    logfile.write("Contents of input file are listed below:\n\n--------\n%s\n-------\n\n" % input)
    input = [line for line in input.split('\n') if line and (not line.isspace()) and line[0] != '#']
    input_dict = {'PROTSEQ':None, 'SELECTED_MUTATION':None, 'PREDICTION_FILE':None}
    for line in input:
        entries = line.split()
        if entries[0] not in input_dict:
            raise ValueError("Input file contains invalid key of %s in line %s." % (entries[0], line))
        elif input_dict[entries[0]] != None and entries[0] not in ['SELECTED_MUTATION', 'PREDICTION_FILE']:
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
        elif entries[0] == 'SELECTED_MUTATION':
            if input_dict['SELECTED_MUTATION'] == None:
                input_dict['SELECTED_MUTATION'] = []
            if len(entries) == 3 and entries[1] != 'FILE':
                try:
                    (wt, r, mut, ddg) = (entries[1][0], int(entries[1][1 : -1]), entries[1][-1], float(entries[2]))
                except ValueError:
                    raise ValueError("Problem with SELECTED_MUTATION specified in line:\n%s" % line)
                logfile.write("\nThe SELECTED_MUTATION key specifies that the ddG is known to be %f for mutation %s%d%s.\n" % (ddg, wt, r, mut))
                input_dict['SELECTED_MUTATION'].append(((wt, r, mut), ddg))
            elif len(entries) == 3 and entries[1] == 'FILE':
                filename = entries[2]
                logfile.write("\nThe SELECTED_MUTATION key specifies that known ddG values for some mutations are listed in the file %s.\n" % filename)
                if not os.path.isfile(filename):
                    raise ValueError("Cannot find the SELECTED_MUTATION FILE %s." % filename)
                logfile.write("The contents of SELECTED_MUTATION FILE %s are listed below:\n\n*******\n%s******\n\n" % (filename, open(filename).read()))
                for line in open(filename).readlines():
                    if '#' in line:
                        line = line[ : line.index('#')]
                    if not line or line.isspace():
                        continue # comment or empty line
                    line = line.strip()
                    line_entries = line.split()
                    if len(line_entries) != 2:
                        raise ValueError("Not two entries (mutation, then ddG) in SELECTED_MUTATION FILE %s line:\n%s" % (filename, line))
                    try:
                        (wt, r, mut, ddg) = (line_entries[0][0], int(line_entries[0][1 : -1]), line_entries[0][-1], float(line_entries[1]))
                    except ValueError:
                        raise ValueError("Problem with SELECTED_MUTATION FILE %s on the following line:\n%s" % (filename, line))
                    logfile.write("The SELECTED_MUTATION FILE %s specifies that the ddG is known to be %f for mutation %s%d%s.\n" % (filename, ddg, wt, r, mut))
                    input_dict['SELECTED_MUTATION'].append(((wt, r, mut), ddg))
            else:
                raise ValueError("Invalid value for SELECTED_MUTATION in line:\n%s" % line)
        elif entries[0] == 'PREDICTION_FILE':
            if input_dict['PREDICTION_FILE'] == None:
                input_dict['PREDICTION_FILE'] = []
            if len(entries) != 3:
                raise ValueError("Invalid value for PREDICTION_FILE (not 3 entries) in line:\n%s" % line)
            prediction_name = entries[1].replace('_', ' ') # replace underscores with spaces
            predictionfile = entries[2]
            logfile.write("\nThe PREDICTION_FILE key specifies will included the ddG predictions from %s found in the file %s.\n" % (prediction_name, predictionfile))
            if not os.path.isfile(predictionfile):
                raise ValueError("Cannot find the PREDICTION_FILE %s" % predictionfile)
            logfile.write("Here are the contents of PREDICTION_FILE %s:\n***********\n%s\n*********\n" % (predictionfile, open(predictionfile).read()))
            (datetime, ddgs) = pips.ddg_inference.ReadDDGs(predictionfile)
            input_dict['PREDICTION_FILE'].append((prediction_name, ddgs))
        else:
            raise IOError("Problem parsing input file keys.")
    logfile.flush()

    # We have parsed the entire input file.  Make sure we have a value for each key:
    for (inputkey, inputvalue) in input_dict.iteritems():
        if inputvalue == None:
            raise ValueError("Failed to read an input value for the required key %s." % inputkey)

    # Check to make sure that all of the selected mutations are valid
    logfile.write('\nWe will look at the ddG values for %d different selected mutations, looking at the predictions made by %d different methods.\n' % (len(input_dict['SELECTED_MUTATION']), len(input_dict['PREDICTION_FILE'])))
    outputfile = open(outputfilename, 'w')
    outputfile.write("Analysis of selected mutations performed by pips_correlate_selected_mutations.py script in directory %s at %s.\n\n" % (os.getcwd(), time.asctime()))
    predicted_ddgs = {}
    for (predictiontype, ddgs) in input_dict['PREDICTION_FILE']:
        if predictiontype in predicted_ddgs: 
            raise ValueError("Duplicate PREDICTION_FILE of %s." % predictiontype)
        predicted_ddgs[predictiontype] = []
    actual_ddgs = []
    for ((wt, r, mut), ddg) in input_dict['SELECTED_MUTATION']:
        actual_ddgs.append(ddg)
        if input_dict['PROTSEQ'][r - 1] != wt:
            raise ValueError("The mutation list specifies that residue %d is %s, but PROTSEQ gives it as %s." % (r, wt, input_dict['PROTSEQ'][r - 1]))
        logfile.write('\nThe known ddG values for mutation %s%d%s is %f; the predicted values are:\n' % (wt, r, mut, ddg))
        outputfile.write('\nThe known ddG values for mutation %s%d%s is %f; the predicted values are:\n' % (wt, r, mut, ddg))
        for (predictiontype, ddgs) in input_dict['PREDICTION_FILE']:
            if r not in ddgs:
                logfile.write('\tPrediction method %s has no predicted ddG for this mutation.\n' % predictiontype)
                outputfile.write('\tPrediction method %s has no predicted ddG for this mutation.\n' % predictiontype)
                predicted_ddgs[predictiontype].append(None)
            else:
                if ddgs[r][0] != wt:
                    raise ValueError("Predicted ddGs from method %s do not have correct wildtype residue at position %r." % (predictiontype, r))
                predicted_ddg = ddgs[r][1][mut]
                predicted_ddgs[predictiontype].append(predicted_ddg)
                logfile.write("\tPrediction method %s predicts a ddG of %.2f.\n" % (predictiontype, predicted_ddg))
                outputfile.write("\tPrediction method %s predicts a ddG of %.2f.\n" % (predictiontype, predicted_ddg))
    logfile.write("\nHere are the Pearson correlation coefficients (R) and associated two-tailed P-values (P) and number of data points (N) for the correlations between the actual ddG values and those predicted by each method.\n")
    outputfile.write("\nHere are the Pearson correlation coefficients (R) and associated two-tailed P-values (P) and number of data points (N) for the correlations between the actual ddG values and those predicted by each method.\n")
    for (predictiontype, ddgs) in input_dict['PREDICTION_FILE']:
        (r, p, n) = pips.stats.PearsonCorrelation(actual_ddgs, predicted_ddgs[predictiontype])
        logfile.write("\tFor prediction method %s, R = %f, P = %f, and N = %d.\n" % (predictiontype, r, p, n))
        outputfile.write("\tFor prediction method %s, R = %f, P = %f, and N = %d.\n" % (predictiontype, r, p, n))
    outputfile.close()

    # Program is complete.
    logfile.write("\nCompleted execution of pips_correlate_selected_mutations.py at %s.\n" % time.asctime())
    logfile.close()


# run the script
main()
