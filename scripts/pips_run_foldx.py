#!python

"""Script to run FOLDX program for PIPS analyses of ddG values for protein mutations.

Written by Jesse Bloom, 2009.

This script is designed to run the FoldX program for predicting the stability changes
of proteins upon mutations.  It has been tested only with FoldX 3.0 Beta3 for
Mac OS X, and so is only guaranteed to work with that version on that operating
system.

FoldX is a force field that is available for download at http://foldx.crg.es.
The references for FoldX are:
   Francois Stricher, Tom Lenaerts, Joost Schymkowitz, Frederic Rousseau and 
       Luis Serrano (2008). FoldX 3.0.  In preparation
   Schymkowitz J., Borg J., Stricher F., Nys R., Rousseau F., Serrano L. (2005) 
       The FoldX web server: an online force field. 
       Nucleic Acids Research, vol 33, pW382-8.
    Schymkowitz J. W., Rousseau F., Martins I. C., Ferkinghoff-Borg J., 
       Stricher F., Serrano L. (2005) 
       Prediction of water and metal binding sites and their affinities by 
       using the Fold-X force field. 
       Proc Natl Acad Sci USA, vol 102, p 10147-52.
    Guerois R., Nielsen J. E., Serrano L. (2002) 
       Predicting changes in the stability of proteins and protein 
       complexes: a study of more than 1000 mutations. 
       J Mol Biol, vol 320, p369-87

FoldX requires a PDB structure to provide structural information on which to base
the predictions.  If a PDB structure for the exact protein in question is not available,
this script provides mechanisms to make predictions using the structure of a closely
related protein.  The resulting ddG values are then transferred to the protein of
interest by pairwise aligning the target protein sequence and the one in the PDB structure.
The ddG values for residues that mismatch are computed based on the assumption that
ddG values are additive.

This script expects to find an input file with the name PIPS_RUN_FOLDX.IN
in the current directory.  Using the configuration specified by this file,
the script creates two output files.  The first of these files,
PIPS_RUN_FOLDX.LOG is a log of the script's progress.  The second file,
FOLDX_PREDICTIONS.TXT gives the FOLDX predictions for the ddG values.  If either
of these output files already exists, the program asks for user confirmation before
continuing and overriding the values.

This script uses a numbering scheme where all residues are numbered sequentially 
according to the protein sequence being analyzed, with the first residue
being number 1.

The input file, PIPS_RUN_FOLDX.IN should have the format as detailed below.
Anything come after a "#" character is considered a comment, and has no effect.
Format of PIPS_RUN_FOLDX.IN file:

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
# The PDB_FILE key specifies the name of the PDB file containing the structure of the
# protein of interest.  Further information about this PDB structure is supplied
# the keys PDB_SEQS_FILE, PDB_CHAIN, and PDB_RES_MAPPING.
PDB_FILE pdb_structure.pdb
#
# The PDB_SEQS_FILE key specifies the name of a FASTA file containing the sequences
# of the chains in the PDB file.  This should be a FASTA file containing one or more
# sequences.  The header for each sequence should consist ONLY of the chain identifier.
# For example, if the PDB structure has two chains, A and B, the FASTA file might be:
#     >A
#     MTRSQFV
#     >B
#     MSRWNNQF
# If the chain has no ID (which sometimes happens when there is only one chain), then
# the header should be empty, as in:
#     >
#     MTRSQFV
# It is important that PDB_SEQS_FILE lists all of the chains in the PDB file even if they
# are not all being used in the FOLDX calculation, because the script uses this file
# to determine how many chains are in the PDB structure.
PDB_SEQS_FILE pdb_seqs.fasta
#
# The PDB_CHAIN key gives the chain ID of the chain that we are analyzing.  This ID
# must match that given in the header of PDB_SEQS_FILE.  If there is only one chain and
# it has no ID, then the chain ID should be set to the string NONE.
PDB_CHAIN A 
# 
# The PDB_RES_MAPPING specifies the relationship between the residue numbers in the PDB file
# and the residues in the PDB chain sequence given in PDB_SEQS_FILE.  In the simplest case,
# the residues in the chain sequence and in the PDB file exactly correspond, with the first
# residue in the chain being residue 1 in the PDB file, the second being residue 2, and
# so on.  In other words, the residues in the PDB file are numbered sequentially with
# no offset.  In this case, PDB_RES_MAPPING would have the following value:
# PDB_RES_MAPPING OFFSET 0
# It is also possible that the residues in the PDB file are numbered sequentially, but
# that there is an offset.  If the first residue in the PDB chain is numbered as X + 1 in
# the PDB file, the second residue in the PDB chain is numbered as X + 2, etc, then
# this is called an offset of X (where X is some integer).  In this case, PDB_RES_MAPPING 
# would have the following value:
# PDB_RES_MAPPING OFFSET X
# If either of these OFFSET options are used, then there must be matching residues in
# the PDB file up through the last residue of the PDB chain sequence (i.e. the PDB chain
# sequence specified in PDB_SEQS_FILE cannot contain additional C-terminal residues
# that are not present in the PDB structure file).
# Finally, there is the possibility that the residues in the PDB file are numbered 
# non-sequentially (or that we only want to look at certain residues).  In this case,
# we explicitly list the number mappings.  For example, say that the sequence of the
# chain in which we are interested is MTRSQFV.  Say that in the PDB file the first
# residue that appears is T which is numbered 4, then an R numbered 5, then an S numbered
# 7, a Q numbered 8, no F, and then a V numbered 10.  In this case, PDB_RES_MAPPING could be
# set as shown below:
# PDB_RES_MAPPING START_RESIDUE_LIST
# PDB_RES_MAPPING 2 4 A
# PDB_RES_MAPPING 3 5 A
# PDB_RES_MAPPING 4 7 A
# PDB_RES_MAPPING 5 8 A
# PDB_RES_MAPPING 7 10 A
# PDB_RES_MAPPING END_RESIDUE_LIST
# In other words, every line starts with PDB_RES_MAPPING, and each line between the 
# START_RESIDUE_LIST and END_RESIDUE_LIST keys gives the number of the residue in
# sequentially numbered sequence followed by the corresponding number in the PDB file. 
# The residue numbers are then followed by the chain ID for the residue in the 
# PDB structure, which in this case does not have to match the chain in PDB_SEQS_FILE.
# So if you wanted to analyze a protein that spans two protein sequence chains, you would
# concatenate the chains into a single sequence in PDB_SEQS_FILE, and then map them to
# the two different chains in the PDB structure using PDB_RES_MAPPING.
# Any residues that do not appear in the list are not analyzed.
PDB_RES_MAPPING OFFSET 0
#
#
# The ALIGNMENT_PROG key is used to specify information about a sequence alignment
# program.  This program is used to align the protein sequence specified by PROTSEQ
# with the PDB protein sequence specified by PDB_SEQS_FILE.  There should be two
# strings following the ALIGNMENT_PROG key, specifying the alignment program and a path
# to the alignment executable.  The two supported alignment programs are
# MUSCLE (http://www.drive5.com/muscle/index.htm) and PROBCONS
# (http://probcons.stanford.edu/).  This script has been tested and confirmed to work
# with version 3.6 of MUSCLE and version 1.12 of PROBCONS.  PROBCONS is the slower
# program, but in general is considered to be more accurate, so should probably be
# preferred unless there are time constraints due to the size of the alignment.
# The specified path should contain executables with the name of either 'muscle'
# or 'probcons'.  So valid entries would be:
# ALIGNMENT_PROG MUSCLE /Users/bloom/muscle3.6/
ALIGNMENT_PROG PROBCONS /Users/bloom/probcons/
#
# The FOLDX_PATH key is used to specify the path to the directory containing
# the FoldX executable (which has a name specified by FOLDX_PROG).  This
# directory must also contain the rotabase.txt file.
FOLDX_PATH /Users/bloom/FoldX/
#
# The FOLDX_PROG key specifies the name of the FoldX executable found in
# FOLDX_PATH
FOLDX_PROG FoldX.mac
----------
"""


import sys
import os 
import time
import traceback
import pips.version
import pips.align
import pips.fasta
import pips.ddg_inference
import pips.foldx


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


class Aligner(object):
    """Class defines objects for performing protein alignments.

    The method is initialized by calling with two arguments, the program type
        and the program directory.  Currently, the program type can be either
        'MUSCLE' or 'PROBCONS', while the program directory gives the directory
        where the executable resides.
    Each object defines a method 'Align' which takes as an argument a list
        of (header, sequence) 2-tuples, and returns a list 'aligned_headers_seqs'
        in which the sequences in the input list have had gaps inserted to align 
        them.
    """
    def __init__(self, program_type, program_directory):
        """Specifies type of program and the directory in which program is found."""
        self._program_type = program_type
        self._program_directory = program_directory
    def Align(self, headers_seqs):
        """Returns copy of 'headers_seqs' in which the sequence are aligned."""
        aligned_headers_seqs = pips.align.Align(headers_seqs, self._program_directory, self._program_type)
        return aligned_headers_seqs


def main():
    """The body of the script, performs alignment and builds phylogenetic trees."""

    # Define input/output files
    logfilename = 'PIPS_RUN_FOLDX.LOG' 
    infilename = 'PIPS_RUN_FOLDX.IN'
    predictionsfilename = 'FOLDX_PREDICTIONS.TXT'
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
    logfile.write("Beginning execution of pips_run_foldx.py at %s.\n" % time.asctime())
    logfile.write("Using version %s of PIPS, written by Jesse Bloom.\n" % pips.version.version)
    logfile.write("Using input/output files in the current directory, which is %s\n\n" % os.getcwd())
    logfile.write("Reading input from %s.\n" % infilename)

    # Read input data into the dictionary input_dict, which is keyed by input keywords and has
    # the input data as values
    try:
        input = open(infilename).read()
    except IOError:
        if not os.path.isfile(infilename):
            raise IOError("Cannot find pips_run_foldx.py input file, %s." % infilename)
        else:
            raise
    logfile.write("Contents of input file are listed below:\n\n********\n%s\n******\n\n" % input)
    input = [line for line in input.split('\n') if line and (not line.isspace()) and line[0] != '#']
    input_dict = {'PROTSEQ':None, 'PDB_FILE':None, 'PDB_SEQS_FILE':None, 
                  'PDB_CHAIN':None, 'PDB_RES_MAPPING':None, 'ALIGNMENT_PROG':None,
                  'FOLDX_PATH':None, 'FOLDX_PROG':None}
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
            input_dict['PROTSEQ'] = (protseqheader, protseq.upper())
        elif entries[0] == 'PDB_FILE':
            if len(entries) != 2:
                raise ValueError("Invalid specification of PDB_FILE on line:\n%s" % line)
            pdbfile = entries[1]
            if not os.path.isfile(pdbfile):
                raise IOError("Cannot find the PDB file of %s specified by PDB_FILE on line:\n%s" % (pdbfile, line))
            logfile.write("\nThe PDB_FILE key specifies that we will use the PDB file %s.\nHere are the contents of this file:\n\n***********\n%s\n**********\n\n" % (pdbfile, open(pdbfile).read()))
            input_dict['PDB_FILE'] = pdbfile
        elif entries[0] == 'PDB_SEQS_FILE':
            if len(entries) != 2:
                raise ValueError("Invalid specification of PDB_SEQS_FILE on line:\n%s" % line)
            if not os.path.isfile(entries[1]):
                raise IOError("Cannot find the PDB sequence file of %s specified by PDB_SEQS_FILE on line:\n%s" % (entries[1], line))
            logfile.write("\nThe PDB_SEQS_FILE key specifies that we will read the sequences of the proteins in the PDB file from the FASTA file %s.  Here are the contents of this file:\n\n********\n%s\n\n" % (entries[1], open(entries[1]).read()))
            input_dict['PDB_SEQS_FILE'] = pips.fasta.Read(entries[1])
            logfile.write("Read the sequences for %d sequences from the PDB_SEQS_FILE of %s.\n" % (len(input_dict['PDB_SEQS_FILE']), entries[1]))
        elif entries[0] == 'PDB_CHAIN':
            if len(entries) != 2:
                raise ValueError("Invalid specification of PDB_CHAIN on line:\n%s" % line)
            input_dict['PDB_CHAIN'] = entries[1]
            if input_dict['PDB_CHAIN'] == 'NONE':
                logfile.write("\nThe PDB_CHAIN key has a value of NONE, indicating that the PDB structure contains only one chain without a chain identifier.  It is this chain that will be used for the FOLDX analysis.\n")
            else:
                logfile.write("\nThe PDB_CHAIN key specifies that we will use chain %s from the PDB file.\n" % input_dict['PDB_CHAIN'])
        elif entries[0] == 'PDB_RES_MAPPING':
            if entries[1] == 'OFFSET' and len(entries) == 3:
                try:
                    offset = int(entries[2])
                except ValueError:
                    raise ValueError("OFFSET specified for PDB_RES_MAPPING is not an integer in line:\n%s" % line)
                logfile.write("\nThe PDB_RES_MAPPING key specifies that residues in the PDB file are expected to match those specified for the appropriate chain in the PDB file, with an offset of %d (the first residue in the PDB chain is numbered as %d in the PDB file).\n" % (offset, 1 + offset))
                input_dict['PDB_RES_MAPPING'] = offset
            elif entries[1] == 'START_RESIDUE_LIST' and len(entries) == 2:
                logfile.write("\nReading a residue mapping list for PDB_RES_MAPPING.\n")
                iline += 1
                line = input[iline]
                entries = line.split()
                mapping = {}
                while not (entries[0] == 'PDB_RES_MAPPING' and entries[1] == 'END_RESIDUE_LIST' and len(entries) == 2):
                    if entries[0] != 'PDB_RES_MAPPING' or len(entries) != 4:
                        raise ValueError("PDB_RES_MAPPING specified START_RESIDUE_LIST, but an invalid residue list entry has been encountered before the END_RESIDUE_LIST in the following line:\n%s" % line)
                    try:
                        rseq = int(entries[1])
                    except ValueError:
                        raise ValueError("Non-integer residue numbers specified in the following PDB_RES_MAPPING line:\n%s" % line)
                    rpdb = entries[2]
                    chain = entries[3]
                    if rseq in mapping:
                        raise ValueError("Multiple mappings specified for PROTSEQ residue %d PDB_RES_MAPPING in the line:\n%s" % (rseq, line))
                    mapping[rseq] = (chain, rpdb)
                    logfile.write("Residue %d of PROTSEQ is mapped to residue %s, chain %s of the PDB file.\n" % (rseq, rpdb, chain))
                    iline += 1
                    line = input[iline]
                    entries = line.split()
                logfile.write("Completed reading the residue mapping list PDB_RES_MAPPING.\n")
                input_dict['PDB_RES_MAPPING'] = mapping
            else:
                raise ValueError("Invalid specification PDB_RES_MAPPING on line:\n%s" % line)
        elif entries[0] == 'ALIGNMENT_PROG':
            logfile.write("\nReading specifications for the alignment program from ALIGNMENT_PROG.\n")
            if len(entries) != 3:
                raise ValueError("ALIGNMENT_PROG does not include three entries in the line:\n%s" % line)
            if entries[1] == 'MUSCLE':
                logfile.write("MUSCLE key specifies that the alignments will be done using MUSCLE.\n")
            elif entries[1] == 'PROBCONS':
                logfile.write("PROBCON key specifies that the alignments will be done using PROBCONS.\n")
            else:
                raise ValueError("Invalid alignment program specification in line:\n%s" % line)
            if os.path.isdir(entries[2]):
                logfile.write("The alignment program resides in directory %s\n" % entries[2])
            else:
                raise ValueError("Cannot find alignment program directory of %s" % entries[2])
            aligner = Aligner(entries[1], entries[2])
            input_dict['ALIGNMENT_PROG'] = aligner
        elif entries[0] == 'FOLDX_PATH':
            if len(entries) != 2:
                raise ValueError("Invalid number of entries in line specifying FOLDX_PATH:\n%s" % line)
            if not os.path.isdir(entries[1]):
                raise ValueError("Cannot find the FOLDX_PATH directory of %s specified in the following line:\n%s" % (entries[1], line))
            logfile.write("\nFOLDX_PATH specifies that the FoldX executable (specified by FOLDX_PROG) and the rotabase.txt file can be found in the directory %s\n" % entries[1])
            input_dict['FOLDX_PATH'] = entries[1]
        elif entries[0] == 'FOLDX_PROG':
            if len(entries) != 2:
                raise ValueError("Invalid number of entries in line specifying FOLDX_PROG:\n%s" % line)
            logfile.write("\nFOLD_PROG specifies that the FoldX executable (located in the directory specified by FOLDX_PATH) has the name %s" % entries[1])
            input_dict['FOLDX_PROG'] = entries[1]
        else:
            raise IOError("Problem parsing input file keys.")
        iline += 1
    logfile.flush()

    # We have parsed the entire input file.  Make sure we have a value for each key:
    for (inputkey, inputvalue) in input_dict.iteritems():
        if inputvalue == None:
            raise ValueError("Failed to read an input value for the required key %s." % inputkey)

    # Determine the specific sequence of the PDB protein, save as 'pdbseq'
    chain = input_dict['PDB_CHAIN']
    if chain == 'NONE':
        chain = ''
    pdbseq = [(head, seq) for (head, seq) in input_dict['PDB_SEQS_FILE'] if head.strip() == chain]
    if not chain:
        chain_id_for_foldx = None
    else:
        chain_id_for_foldx = input_dict['PDB_CHAIN']
    if len(pdbseq) != 1:
        raise ValueError("Expected to find one sequence in PDB_SEQS_FILE that had the name specified by PDB_CHAIN, but instead found %d." % len(pdbseq))
    pdbseq = pdbseq[0]
    logfile.write("\nThe PDB sequence in PDB_SEQS_FILE for PDB_CHAIN (chain %s) is shown below:\n>%s\n%s\n" % (input_dict['PDB_CHAIN'], pdbseq[0], pdbseq[1]))

    # Align 'PROTSEQ' with 'pdbseq'
    logfile.write("\nNow aligning PROTSEQ with this PDB sequence...\n")
    logfile.flush()
    alignment = aligner.Align([input_dict['PROTSEQ'], pdbseq])
    logfile.write("Completed the alignment of these two sequences.  The alignment is shown below:\n")
    pips.fasta.Write(alignment, logfile, writable_file=True)
    logfile.flush()
    alignment = (alignment[0][1], alignment[1][1])

    # Now build the 'pdb_matches' object to map the residue numbers to PDB numbers/chains
    if isinstance(input_dict['PDB_RES_MAPPING'], int):
        # protein and PDB sequences match, with an offset
        offset = input_dict['PDB_RES_MAPPING']
        mapping = dict([(r, r + 1 + offset) for r in range(len(pdbseq[1]))])
        pdb_matches = [(input_dict['PDB_FILE'], chain_id_for_foldx, alignment, mapping)]
    else:
        chain_dict = {}
        for (rseq, (chain, rpdb)) in input_dict['PDB_RES_MAPPING'].iteritems():
            if chain not in chain_dict:
                chain_dict[chain] = (input_dict['PDB_FILE'], chain, alignment, {})
            if rseq - 1 in chain_dict[chain][3]:
                raise ValueError("Duplicate PDB_RES_MAPPING entry specified for residue %d." % rseq)
            chain_dict[chain][3][rseq - 1] = rpdb
        pdb_matches = chain_dict.values()

    # Now make the FOLDX predictions using the webserver
    hashed_ddgs = 'FOLDX_HASHED_DDGS_FOR_%s.TXT' % input_dict['PDB_FILE']
    logfile.write("\nNow running the FOLDX predictions...\n")
    logfile.write("Hashed FOLDX ddG values will be read/written from file %s if the program terminates or has terminated prematurely.\n" % hashed_ddgs)
    logfile.flush()
    try:
        original_stdout = sys.stdout
        sys.stdout = logfile # assign standard output from FOLDX predictions to logfile
        foldx_ddgs = pips.foldx.DDGsAlignedToSeq(input_dict['PROTSEQ'][1], pdb_matches, input_dict['FOLDX_PATH'], foldxprog=input_dict['FOLDX_PROG'], hashed_ddgs=hashed_ddgs)
    finally:
        sys.stdout = original_stdout
    logfile.write("Completed making the FOLDX predictions.\n")

    # Renumber the FOLDX ddGs so that the first residue is number 1
    firstresnum = 1
    logfile.write("\nThe ddG values are numbered according to the scheme in which the first residue of PROTSEQ has number %d.\n" % firstresnum)
    renumbered_foldx_ddgs = {}
    for (r, (rwt, rddgs)) in foldx_ddgs.iteritems():
        renumbered_foldx_ddgs[r + firstresnum] = (rwt, rddgs)

    # Now write the FOLDX ddGs to the predictions file
    logfile.write("\nThe FOLDX estimated ddG values are now being written to the file %s...\n" % predictionsfilename)
    datetime = 'FOLDX ddGs computed from PDB file %s at %s in the directory %s, utilized PIPS version %s.' % (input_dict['PDB_FILE'], time.asctime(), os.getcwd(), pips.version.version)
    pips.ddg_inference.WriteDDGs(renumbered_foldx_ddgs, predictionsfilename, datetime)
    logfile.write("These FOLDX computed ddG values are as shown below:\n\n**********\n%s\n*******\n" % open(predictionsfilename).read())

    # Program is complete.
    logfile.write("\nCompleted execution of pips_run_foldx.py at %s.\n" % time.asctime())
    logfile.close()


# run the script
main()
