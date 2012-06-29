#!python

"""Script to make consensus approach ddG predictions.

Written by Jesse Bloom, 2009.

This script calculates ddG values based on the "consensus" approach, which is
the idea that the frequency of an amino acid in the multiple sequence alignment
has a Boltzmann-like relationship to its ddG value.  Specifically, let N_i
be the number of times amino acid i occurs in a multiple sequence alignment
of N proteins.  Then the stability change for mutating from residue i to j
(ddG_i->j) is estimated to be:
    ddG_i->j = ln([N_j + X] / [N_i + X])
where X is a pseudocount that is added to each count.  If there are no
pseudocounts (X = 0), then some ddG values will be undefined.

This script expects to find an input file with the name PIPS_CONSENSUS.IN
in the current directory.  Using the configuration specified by this file,
the script creates two output files.  The first of these files,
PIPS_CONSENSUS.LOG is a log of the script's progress.  The second file,
CONSENSUS_PREDICTIONS.TXT gives the consensus predictions for the ddG values.  If either
of these output files already exists, the program asks for user confirmation before
continuing and overriding the values.

This script uses a numbering scheme where all residues are numbered sequentially 
according to the protein sequence being analyzed, with the first residue
being number 1.

The input file, PIPS_CONSENSUS.IN should have the format as detailed below.
Anything come after a "#" character is considered a comment, and has no effect.
Format of PIPS_CONSENSUS.IN file:

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
# The PSEUDOCOUNTS key specifies the number of pseudocounts that are added to each amino
# acid.  This should be an integer >= 0.  Note that if it set to zero, some ddG values
# will be undefined.
PSEUDOCOUNTS 1
#
# The ALIGNMENT key specifies the multiple sequence alignment of the homologous proteins.
# This alignment should be given by a FASTA file.
# The first protein sequence in this alignment must match that specified by PROTSEQ, with
# the possible addition of some gap characters.  Once the alignment has been read, positions
# corresponding to gaps in this sequence are stripped away from all sequences.  
# For example, if the sequences are
# named "SEQ1" and "SEQ2", the file might look like this:
#   >SEQ1: this sequence matches that given by PROTSEQ
#   MT-DYST-FFT
#   >SEQ2: second sequence
#   MTFEFSTN-FS
ALIGNMENT renamed_alignment.fasta
#
----------
"""


import sys
import os 
import time
import traceback
import pips.version
import pips.fasta
import pips.align
import pips.ddg_inference


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
    """The main body of the script."""

    # Define input/output files
    logfilename = 'PIPS_CONSENSUS.LOG' 
    infilename = 'PIPS_CONSENSUS.IN'
    predictionsfilename = 'CONSENSUS_PREDICTIONS.TXT'
    if os.path.isfile(logfilename):
        if QueryFileOverwrite(logfilename, 'log'):
            print "Deleting existing %s file." % (logfilename)
            if os.path.isfile(logfilename):
                os.remove(logfilename)
        else:
            print "Canceling program."
            sys.exit()
    logfile = open(logfilename, 'w')
    if os.path.isfile(predictionsfilename):
        if QueryFileOverwrite(predictionsfilename, 'predictions'):
            print "Deleting existing %s file." % (predictionsfilename)
            if os.path.isfile(predictionsfilename):
                os.remove(predictionsfilename)
        else:
            print "Canceling program."
            sys.exit()

    # Make sure exception information is written to logfile
    exceptionwriter = FileWriteExceptHook(logfile)
    sys.excepthook = exceptionwriter.ExceptHook

    # print introductory information
    logfile.write("Beginning execution of pips_consensus.py at %s.\n" % time.asctime())
    logfile.write("Using version %s of PIPS, written by Jesse Bloom.\n" % pips.version.version)
    logfile.write("Using input/output files in the current directory, which is %s\n\n" % os.getcwd())
    logfile.write("Reading input from %s.\n" % infilename)

    # Read input data into the dictionary input_dict, which is keyed by input keywords and has
    # the input data as values
    try:
        input = open(infilename).read()
    except IOError:
        if not os.path.isfile(infilename):
            raise IOError("Cannot find pips_consensus.py input file, %s." % infilename)
        else:
            raise
    logfile.write("Contents of input file are listed below:\n\n********\n%s\n******\n\n" % input)
    input = [line for line in input.split('\n') if line and (not line.isspace()) and line[0] != '#']
    input_dict = {'PROTSEQ':None, 'PSEUDOCOUNTS':None, 'ALIGNMENT':None}
    iline = 0
    while iline < len(input):
        line = input[iline]
        entries = line.split()
        if entries[0] not in input_dict:
            raise ValueError("Input file contains invalid key of %s in line %s." % (entries[0], line))
        elif input_dict[entries[0]] != None:
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
            input_dict['PROTSEQ'] = (protseqheader, protseq)
        elif entries[0] == 'PSEUDOCOUNTS':
            if len(entries) != 2:
                raise ValueError("Invalid specification of PSEUDOCOUNTS on line:\n%s" % line)
            pseudocounts = entries[1]
            try:
                pseudocounts = int(pseudocounts)
            except ValueError:
                raise ValueError("Cannot convert PSEUDOCOUNTS to integer in line:\n%s" % line)
            if pseudocounts < 0:
                raise ValueError("PSEUOCOUNTS must be >= 0 in line:\n%s" % line)
            logfile.write("\nRead PSEUDOCOUNTS value of %d.  This many pseudocounts will be added to the number of occurrences for each amino acid at each position.\n" % pseudocounts)
            input_dict['PSEUDOCOUNTS'] = pseudocounts
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
        else:
            raise IOError("Problem parsing input file keys.")
        iline += 1
    logfile.flush()

    # We have parsed the entire input file.  Make sure we have a value for each key:
    for (inputkey, inputvalue) in input_dict.iteritems():
        if inputvalue == None:
            raise ValueError("Failed to read an input value for the required key %s." % inputkey)

    # Comput the consensus ddG values
    logfile.write("\nComputing the consensus ddG values using %d pseudocounts...\n" % input_dict['PSEUDOCOUNTS'])
    logfile.flush()
    input_dict['ALIGNMENT'] = pips.fasta.UnknownsToGaps(input_dict['ALIGNMENT']) # replace unknown amino acids with gaps
    input_dict['ALIGNMENT'] = pips.align.StripGapsToFirstSequence(input_dict['ALIGNMENT'])
    if input_dict['ALIGNMENT'][0][1] != input_dict['PROTSEQ'][1]:
        raise ValueError("First protein sequence in ALIGNMENT does not match the sequence of PROTSEQ.")
    ddgs = pips.ddg_inference.ConsensusDDGs(input_dict['PROTSEQ'][1], input_dict['ALIGNMENT'], input_dict['PSEUDOCOUNTS'])
    logfile.write("Completed computing the consensus ddGs.\n")


    # Renumber the ddGs so that the first residue is number 1
    firstresnum = 1
    logfile.write("\nThe ddG values are numbered according to the scheme in which the first residue of PROTSEQ has number %d.\n" % firstresnum)
    renumbered_ddgs = {}
    for (r, (rwt, rddgs)) in ddgs.iteritems():
        renumbered_ddgs[r + firstresnum] = (rwt, rddgs)

    # Now write the ddGs to the predictions file
    logfile.write("\nThe consensus estimated ddG values are now being written to the file %s...\n" % predictionsfilename)
    datetime = 'Consensus ddGs computed at %s in the directory %s, utilized PIPS version %s.' % (time.asctime(), os.getcwd(), pips.version.version)
    pips.ddg_inference.WriteDDGs(renumbered_ddgs, predictionsfilename, datetime)
    logfile.write("These consensus ddG values are as shown below:\n\n**********\n%s\n*******\n" % open(predictionsfilename).read())

    # Program is complete.
    logfile.write("\nCompleted execution of pips_consensus.py at %s.\n" % time.asctime())
    logfile.close()


# run the script
main()
