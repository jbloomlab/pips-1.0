#!python

"""Script to build a protein sequence alignment and phylogenetic tree for PIPS.

Written by Jesse Bloom, 2009.

This script can be used to create the sequence alignment and phylogenetic tree
that are subsequently used by 'pips_build_tree_and_alignment.py'.  In general, this script
can start with a list of putatively homologous sequences, as for example might be
returned by a BLAST search.  It can then be used to pare down these sequences
based on homology requirements, align the sequences, and use them to build
a phylogenetic tree.

This script expects to find an input file named PIPS_BUILD_TREE_AND_ALIGNMENT.IN
in the current directory.  Using the configuration specified by this file,
the script creates several output files.  The first of these files,
PIPS_BUILD_TREE_AND_ALIGNMENT.LOG is a log of the program progress.  The other
output files contain the constructed phylogenetic tree, the sequence alignment,
and the names that the names that the sequences are assigned in the
phylogenetic tree.  These output files have file names as specified in the input file.
If any of these output files already exist, the program asks for user confirmation
before overwriting them.

This program uses either the MUSCLE or PROBCONS programs for the sequence alignments
-- they must be present with executable paths specified in the input file.  This
program uses PHYLIP to build the phylogenetic tree -- this program must also be
present with executable paths specified in the input file.

The input file, PIPS_BUILD_TREE_AND_ALIGNMENT.IN should have the format detailed below.
Anything come after a "#" character is considered a comment, and has no effect.

Format of PIPS_BUILD_TREE_AND_ALIGNMENT.IN file:

***********
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
# The OUTPUT_ALIGNMENT key specifies the name of the file to which the aligned protein
# sequences are written.  This key should therefore be followed by a string giving
# a valid file name.  If the file already exists, the user will be prompted for confirmation
# before it is overwritten.  The aligned sequences are written to this file in FASTA format,
# with the first sequence in the file being the one specified by PROTSEQ.  The sequences
# are aligned by inserting gap ("-") characters as appropriate.  The sequences that are 
# aligned include PROTSEQ and all unique sequences specified by HOMOLOGOUS_SEQS
# that satisfy the criteria established by HOMOLOGY_THRESHOLD.  The headers
# of these sequences are replaced with short names, that are related to the original headers
# as described in the file specified by OUTPUT_SEQNAMES.  These names are the ones used
# in the phylogenetic tree specified by OUTPUT_TREE.  If the alignment name is followed by
# the string USE_EXISTING then we do NOT build a new alignment, but use the one already present
# in the specified file.  In this case, this file must exist and must be generated exactly 
#according to the name specifications that would be used if this program was creating the alignment.
OUTPUT_ALIGNMENT renamed_alignment.fasta
#
# The OUTPUT_TREE key specifies the name of the file to which is written the phylogenetic
# tree specifying the relationship among the sequences in OUTPUT_ALIGNMENT.  The names
# of the sequences in this tree match those used in OUTPUT_ALIGNMENT, and are detailed in
# relation to the original seauence headers by OUTPUT_SEQNAMES.  The tree is
# in Newick format (http://evolution.gs.washington.edu/phylip/newick_doc.html).  If the
# specified file already exists, the user is prompted for confirmation before the file
# is overwritten.
OUTPUT_TREE tree.newick
#
# The OUTPUT_SEQNAMES key specifies the name of a file to which we write the relationship
# between the new names given to the sequences in OUTPUT_TREE and OUTPUT_ALIGNMENT, and
# the original sequence names as given in their headers in the FASTA files specified by
# PROTSEQ and HOMOLOGOUS_SEQS.  The format is simply lines with the new name, a space,
# and then the original header in the FASTA file. 
OUTPUT_SEQNAMES seqnames.txt
#
# The PHYLIP_PROG key is used to specify information about the PHYLIP program that is
# used to generate the phylogenetic tree.  This PIPS script is guaranteed to work with 
# version 3.68 of PHYLIP for Mac OS X, downloaded from 
# http://evolution.genetics.washington.edu/phylip.html
# It will probably also work with other versions of PHYLIP on other operating
# systems, but that has not been actually tested.  There are a number of 
# strings that are present after the PHYLIP_PROG key.  These strings cannot
# contain spaces, but can have any other non-whitespace characters.  They are:
#  * The first string should give the path to which the PHYLIP executables can
#    be found. 
#  * The second string is the type of tree that is build.  Currently, the only
#    supported option is DISTANCE for methods based on pairwise distances
#    between sequences.  These trees are built using one of the 'fitch',
#    'kitsch', or 'neighbor' programs of PHYLIP.
#  * The third string specifies whether or not we are using a molecular clock.
#    It can have two values, CLOCK or NO_CLOCK depending on whether or not
#    we assume a molecular clock.  If we are not using a molecular clock,
#    then the OUTGROUP_SEQ key must be set to a value in order to provide
#    an outgroup sequence for rooting the tree.
#  * The fourth string specifies whether or not we use a neighbor joining
#    method (neighbor joining without a molecular clock, or UPGMA with
#    a molecular clock).  Neighbor joining is generally faster, but
#    considered to be less accurate.  The string should be either
#    NEIGHBOR if we are using neighbor joining / UPGMA or NO_NEIGHBOR
#    if we are not.  Neighbor joining can lead to negative branch lengths,
#    so any negative lengths are set to zero.
PHYLIP_PROG /Users/bloom/phylip-3.68/exe/ DISTANCE CLOCK NO_NEIGHBOR
#
# The ALIGNMENT_PROG key is used to specify information about the program used to make
# the multiple sequence alignment of the protein sequences.  There should be two 
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
# The HOMOLOGOUS_SEQS key specifies the homologous sequences that are aligned
# and used to build the phylogenetic tree.  The key should be followed by the
# names of one or more existing FASTA files.  All of these sequences are read.
# Any duplicate sequences (identical sequences, NOT necessarily identical headers)
# have only the first sequence retained, with subsequent duplicates discarded.  These
# sequences are then screened according to the criteria specified by HOMOLOGY_THRESHOLD,
# and any sufficiently homologous sequences are retained for the alignment and to 
# build the phylogenetic tree.  For the alignment, the sequences are renamed
# as specified by OUTPUT_SEQNAMES
HOMOLOGOUS_SEQS blast_matches.fasta
#
# The HOMOLOGY_THRESHOLD key specifies the screening of the sequences in HOMOLOGOUS_SEQS.
# It can be used to only include sequences that have a specified degree of identity with
# PROTSEQ.  You should specify both an identity cutoff and a gap cutoff.  The identity
# cutoff, which occurs after the string IDENTITY_CUTOFF, means that we only include
# sequences in which the fractions of identities among alignable residues is greater than
# or equal to this cutoff.  The gap cutoff, which comes after the string GAP_CUTOFF,
# means that we only include sequences in which the fraction of graps is less 
# than or equal to the indicated cutoff.  These values are computed
# in pairwise alignments of PROTSEQ with each sequence.  For example, if the sequence
# specified by PROTSEQ is MGCQT and the pairwise aligned homologous sequence is MG-QA,
# then the identity is 0.75 and the gap fraction is 0.2.  If you don't want to exclude any
# of the sequences, then set the identity cutoff to 0.0 and the gap cutoff 1.0, which
# will include all sequences.
HOMOLOGY_THRESHOLD IDENTITY_CUTOFF 0.55 GAP_CUTOFF 0.1
#
# The OUTGROUP_SEQ key specifies an outgroup sequence used to root the phylogenetic
# tree if the tree is being computed without a molecular clock.  If the tree is being
# computed without a molecular clock, then we need an outgroup to root the tree.  In general,
# the outgroup will be more distantly related to the homologous sequences.  In this case,
# the value should be the name of a FASTA file.  The first sequence listed in this
# FASTA file is used as the outgroup sequence.  So for example, this line would look like this:
# OUTGROUP_SEQ outgroup.fasta
# On the other hand, if we are not using a molecular clock, not outgroup sequence is 
# needed to root the tree.  In that case, the value should be NONE, as in:
OUTGROUP_SEQ NONE
#
*********

"""


import sys
import os 
import time
import traceback
import pips.version
import pips.align
import pips.fasta
import pips.tree
import pips.phylip


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
    use_existing_alignment = False # don't do this unless specified
    firstresnum = 1 # number assigned to the first residue
    logfilename = 'PIPS_BUILD_TREE_AND_ALIGNMENT.LOG' 
    infilename = 'PIPS_BUILD_TREE_AND_ALIGNMENT.IN'
    if os.path.isfile(logfilename):
        if QueryFileOverwrite(logfilename, 'log'):
            print "Deleting existing %s file." % (logfilename)
            if os.path.isfile(logfilename):
                os.remove(logfilename)
        else:
            print "Canceling program."
            sys.exit()
    logfile = open(logfilename, 'w')

    # Make sure exception information is written to logfile
    exceptionwriter = FileWriteExceptHook(logfile)
    sys.excepthook = exceptionwriter.ExceptHook

    # print introductory information
    logfile.write("Beginning execution of pips_build_tree_and_alignment.py at %s.\n" % time.asctime())
    logfile.write("Using version %s of PIPS, written by Jesse Bloom.\n" % pips.version.version)
    logfile.write("Using input/output files in the current directory, which is %s\n\n" % os.getcwd())
    logfile.write("Reading input from %s.\n" % infilename)

    # Read input data into the dictionary input_dict, which is keyed by input keywords and has
    # the input data as values
    try:
        input = open(infilename).read()
    except IOError:
        if not os.path.isfile(infilename):
            raise IOError("Cannot find pips_build_tree_and_alignnment.py input file, %s." % infilename)
        else:
            raise
    logfile.write("Contents of input file are listed below:\n\n********\n%s\n******\n\n" % input)
    input = [line for line in input.split('\n') if line and (not line.isspace()) and line[0] != '#']
    input_dict = {'PROTSEQ':None, 'OUTPUT_ALIGNMENT':None, 'OUTPUT_TREE':None, 
                  'OUTPUT_SEQNAMES':None, 'PHYLIP_PROG':None, 'ALIGNMENT_PROG':None,
                  'HOMOLOGOUS_SEQS':None, 'HOMOLOGY_THRESHOLD':None, 'OUTGROUP_SEQ':None}
    for line in input:
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
        elif entries[0] in ['OUTPUT_ALIGNMENT', 'OUTPUT_TREE', 'OUTPUT_SEQNAMES']:
            descriptions = {'OUTPUT_ALIGNMENT' : 'output alignment',
                            'OUTPUT_TREE' : 'phylogenetic tree',
                            'OUTPUT_SEQNAMES' : 'sequence names'}
            if entries[0] == 'OUTPUT_ALIGNMENT' and len(entries) == 3 and entries[2] == 'USE_EXISTING':
                logfile.write('USE_EXISTING option for OUTPUT_ALIGNMENT specifies that we use the existing protein alignment.  Note it must have been generated with all of the exact same specifications here or you will get erroneous output!\n')
                input_dict[entries[0]] = entries[1]
                logfile.write("OUTPUT_ALIGNMENT specifies that we get the alignment from the file %s.\n" % entries[1])
                if not os.path.isfile(entries[1]):
                    raise ValueError("USE_EXISTING specified for OUTPUT_ALIGNMENT, but cannot find file %s" % entries[1])
                use_existing_alignment = True
            else:
                if len(entries) != 2:
                    raise ValueError("Improper specification of %s.\nThe input file contained the following specification:\n%s" % (entries[0], line))
                input_dict[entries[0]] = entries[1]
                logfile.write("\n%s specifies that we output the %s to %s.\n" % (entries[0], descriptions[entries[0]], input_dict[entries[0]]))
                if os.path.isfile(input_dict[entries[0]]):
                    if QueryFileOverwrite(input_dict[entries[0]], descriptions[entries[0]]):
                        logfile.write("\nA file with name %s already exists.  This existing file is being deleted.\n" % input_dict[entries[0]])
                        os.remove(input_dict[entries[0]])
                    else:
                        logfile.write("\nA file with name %s already exists.  You have selected not to overwrite this file, so the script is being canceled.\n" % input_dict[entries[0]])
                        logfile.close()
                        print "Program canceled."
                        sys.exit()
        elif entries[0] == 'PHYLIP_PROG':
            if len(entries) != 5:
                raise ValueError("Improper specification of PHYLIP_PROG.\nThe correct specification has form:\nPHYLIP_PROG /Users/bloom/phylip-3.68/exe/ DISTANCE CLOCK NO_NEIGHBOR\nThe input file contained the following specification:\n%" % (line))
            logfile.write("\nReading specifications for the phylogenetic tree from PHYLIP_PROG.\n")
            dir = entries[1]
            if not os.path.isdir(dir):
                raise ValueError("Cannot find the PHYLIP directory of %s that is specified by %s" % (dir, line))
            else:
                logfile.write("The PHYLIP executables are specified to reside in the directory %s" % dir)
            if entries[2] == 'DISTANCE':
                phylip_method = entries[2]
                logfile.write("The DISTANCE keyword specifies that we will construct a tree based on pairwise sequence distances.\n")
                if entries[3] == 'CLOCK':
                    molecular_clock = True
                    logfile.write("The CLOCK keyword specifies that we will construct the tree assuming a molecular clock.\n")
                elif entries[3] == 'NO_CLOCK':
                    molecular_clock = False
                    logfile.write("The NO_CLOCK keyword specifies that we will construct the tree without assuming a molecular clock.\n")
                else:
                    raise ValueError("Third entry fails to specify CLOCK or NO_CLOCK in line:\n%s" % line)
                if entries[4] == 'NEIGHBOR':
                    neighbor_joining = True
                    logfile.write("The NEIGHBOR keyword specifies that we will build a neighbor joining or UPGMA tree.\n")
                elif entries[4] == 'NO_NEIGHBOR':
                    neighbor_joining = False
                    logfile.write("The NO_NEIGHBOR keyword specifies that we will not use neighbor joining or UPGMA.\n")
                else:
                    raise ValueError("Fourth entry fails to specify NEIGHBOR or NO_NEIGHBOR in line:\n%s\n" % line)
                input_dict['PHYLIP_PROG'] = ('DISTANCE', dir, molecular_clock, neighbor_joining)
            else:
                raise ValueError("PHYLIP_PROG specifies an invalid tree type of %s in line %s" % (entries[2], line))
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
        elif entries[0] == 'HOMOLOGOUS_SEQS':
            logfile.write("\nReading in the homologous protein sequences specified by HOMOLOGOUS_SEQS.\n")
            seqfiles = entries[1 : ]
            logfile.write("These sequences are contained in the following FASTA file(s): %s\n" % (', '.join(seqfiles)))
            homologous_seqs = []
            for seqfile in seqfiles:
                logfile.write("Reading sequences from %s...\n" % seqfile)
                if not os.path.isfile(seqfile):
                    raise IOError("Cannot find file %s." % seqfile)
                seqs = pips.fasta.Read(seqfile)
                logfile.write("Read %d sequences from %s.\n" % (len(seqs), seqfile))
                homologous_seqs += seqs
            logfile.write("Read a total of %d homologous sequences.\n" % len(homologous_seqs))
            input_dict['HOMOLOGOUS_SEQS'] = homologous_seqs
        elif entries[0] == 'HOMOLOGY_THRESHOLD':
            logfile.write("\nReading the homology thresholds specified by HOMOLOGY_THRESHOLD.\n")
            if len(entries) != 5 or entries[1] != 'IDENTITY_CUTOFF' or entries[3] != 'GAP_CUTOFF':
                raise ValueError("Invalid entry for HOMOLOGY_THRESHOLD in line:\n%s" % line)
            try:
                identity_cutoff = float(entries[2])
            except ValueError:
                raise ValueError("Failed to convert IDENTITY_CUTOFF to a number in line:\n%s" % line)
            try:
                gap_cutoff = float(entries[4])
            except ValueError:
                raise ValueError("Failed to convert GAP_CUTOFF to a number in line:\n%s" % line)
            for (cutoff, cutoff_name) in [(identity_cutoff, 'IDENTITY_CUTOFF'), (gap_cutoff, 'GAP_CUTOFF')]:
                if cutoff < 0 or cutoff > 1:
                    raise ValueError("%s of %f is not in the range from zero to one." % (cutoff_name, cutoff))
            logfile.write("Homologous sequences will be pairwise aligned with PROTSEQ, and any with less than %f fraction identities among alignable residues, or with more than %f fraction gaps, will be excluded.\n" % (identity_cutoff, gap_cutoff))
            input_dict['HOMOLOGY_THRESHOLD'] = (identity_cutoff, gap_cutoff)
        elif entries[0] == 'OUTGROUP_SEQ':
            logfile.write("\nReading the outgroup sequence specified by OUTGROUP_SEQ.  This sequence will be used to root the phylogenetic tree if we do not assume a molecular clock; it will not be used if we assume a molecular clock.\n")
            if len(entries) != 2:
                raise ValueError("Invalid entry for OUTGROUP_SEQ in line:\n%s" % line)
            if entries[1] == 'NONE':
                logfile.write("Value of NONE, meaning no outgroup sequence is specified.\n")
                input_dict['OUTGROUP_SEQ'] = 'NONE'
            else:
                if not os.path.isfile(entries[1]):
                    raise ValueError("Cannot find FASTA file specified by OUTGROUP_SEQ in line:\n%s" % line)
                logfile.write("The outgroup sequence will be the first sequence present in the FASTA file %s\n" % entries[1])
                outgroup = pips.fasta.Read(entries[1])[0]
                logfile.write("This outgroup sequence is:\n")
                pips.fasta.Write([outgroup], logfile, writable_file=True)
                input_dict['OUTGROUP_SEQ'] = outgroup
        else:
            raise IOError("Problem parsing input file keys.")
    logfile.flush()

    # We have parsed the entire input file.  Make sure we have a value for each key:
    for (inputkey, inputvalue) in input_dict.iteritems():
        if inputvalue == None:
            raise ValueError("Failed to read an input value for the required key %s." % inputkey)

    # Screen the homologous sequences for duplicates, meeting the specified homology threshold
    logfile.write("\nHOMOLOGOUS_SEQS specified %d sequences.  We will now eliminate any duplicates among these sequences...\n" % len(input_dict['HOMOLOGOUS_SEQS']))
    logfile.flush()
    seq_dict = {}
    unique_seqs = []
    for (head, seq) in input_dict['HOMOLOGOUS_SEQS']:
        seq = seq.upper()
        if seq in seq_dict:
            logfile.write("Duplicate sequence -- the sequences with the following two headers are identical; purging the latter from the sequence list:\n\t>%s\n\t>%s\n" % (seq_dict[seq], head))
        elif seq == input_dict['PROTSEQ'][1]:
            logfile.write("Duplicate sequence -- the sequence with the header shown below is identical to PROTSEQ:\n\t%s\n" % head)
        else:
            seq_dict[seq] = head
            unique_seqs.append((head, seq))
    logfile.write("After eliminating the duplicates, there are %d unique sequences.\n" % len(unique_seqs))
    (identity_cutoff, gap_cutoff) = input_dict['HOMOLOGY_THRESHOLD']
    logfile.write("\nBased on the criteria specified by HOMOLOGY_THRESHOLD, we will screen the %d unique homologous sequences...\n" % (len(unique_seqs)))
    screened_seqs = []
    if identity_cutoff <= 0.0 and gap_cutoff >= 1.0:
        screened_seqs = unique_seqs
        logfile.write("HOMOLOGY_THRESHOLD specified a IDENTITY_CUTOFF of %f and a GAP_CUTOFF of %f.  All sequences will meet this threshold, so no screening is being performed.  All %d unique sequences are being retained.\n" % (identity_cutoff, gap_cutoff, len(screened_seqs)))
    else:
        logfile.write("Each sequence will be pairwise aligned with PROTSEQ, and those with a fraction identity of at least %f among alignable residues, and with a fraction of gaps no greater than %f will be retained...\n" % (identity_cutoff, gap_cutoff))
        logfile.flush()
        iseq = 0
        for (head, seq) in unique_seqs:
            iseq += 1
            logfile.write("Performing pairwise alignment of PROTSEQ with sequence %d of %d, which has header:\n\t%s\n" % (iseq, len(unique_seqs), head))
            alignment = aligner.Align([input_dict['PROTSEQ'], (head, seq)])
            (fracidentity, fracgap) = pips.align.PairwiseStatistics(alignment)
            logfile.write("\tThis sequence aligns to PROTSEQ with %f fraction identities, and %f fraction gaps.\n" % (fracidentity, fracgap))
            if (fracidentity >= identity_cutoff) and (fracgap <= gap_cutoff):
                logfile.write("\tRetaining this sequence since it satisfies HOMOLOGY_THRESHOLD.\n")
                screened_seqs.append((head, seq))
            else:
                logfile.write("\tPurging this sequence since it does not satisfy HOMOLOGY_THRESHOLD.\n")
        logfile.write("After screening for sequences satisfying HOMOLOGY_THRESHOLD, %d unique homologous sequences remain.\n" % len(screened_seqs))

    # Now make an MSA of all PROTSEQ, all of the homologous sequences, and maybe the OUTGROUP_SEQ
    logfile.write("\nWe will now build a multiple sequence alignment (MSA) of all of the relevant sequences.\n")
    logfile.write("The first sequence in the MSA will be PROTSEQ.\n")
    logfile.write("The MSA will then include the %d screened unique homologous sequences.\n" % len(screened_seqs))
    molecular_clock = input_dict['PHYLIP_PROG'][2]
    if input_dict['OUTGROUP_SEQ'] == 'NONE' and molecular_clock:
        logfile.write("This MSA will not include an outgroup sequence, since OUTGROUP_SEQ was set to NONE, and since no outgroup is needed since we are not using a molecular clock.\n")
        outgroup = None
    elif input_dict['OUTGROUP_SEQ'] == 'NONE' and not molecular_clock:
        raise ValueError("OUTGROUP_SEQ failed to specify an outgroup sequence, yet such a sequence is needed to root the phylogenetic tree, since PHYLIP_PROG specifies that we are not using a molecular clock.\n")
    elif input_dict['OUTGROUP_SEQ'] != 'NONE' and molecular_clock:
        logfile.write("Although OUTGROUP_SEQ specifies an outgroup sequence, this sequence will be ignored since we are using a molecular clock to build the phylogenetic tree, meaning that no outgroup is necessary to root the tree.\n")
        outgroup = None
    else:
        outgroup = input_dict['OUTGROUP_SEQ']
        logfile.write("The last sequence in the MSA will be the outgroup sequence specified by OUTGROUP_SEQ, which will be used to root the phylogenetic tree.  This outgroup sequence has the following header:\n\t%s\n" % outgroup[0])
    seqnames = [('PROTSEQ', input_dict['PROTSEQ'][0])]
    sequences = [('PROTSEQ', input_dict['PROTSEQ'][1])]
    iseq = 0
    for (head, seq) in screened_seqs:
        iseq += 1
        seqname = 'SEQ%d' % iseq
        seqnames.append((seqname, head))
        sequences.append((seqname, seq))
    if outgroup:
        seqnames.append(('OUTGROUP', outgroup[0]))
        sequences.append(('OUTGROUP', outgroup[1]))
    logfile.write("The MSA will include a total of %d sequences.  These sequences have all been given short names that are acceptable to PHYLIP.  A table relating these new names to the original sequence headers has been written to the OUTPUT_SEQNAMES file of %s.\n" % (len(sequences), input_dict['OUTPUT_SEQNAMES']))
    f = open(input_dict['OUTPUT_SEQNAMES'], 'w')
    for (seqname, head) in seqnames:
        f.write("%s %s\n" % (seqname, head))
    f.close()
    logfile.write("\nHere are the contents of the OUTPUT_SEQNAMES file of %s:\n\n******\n%s\n******\n\n" % (input_dict['OUTPUT_SEQNAMES'], open(input_dict['OUTPUT_SEQNAMES']).read()))
    if use_existing_alignment:
        logfile.write("\nUsing the existing alignment in %s" % input_dict['OUTPUT_ALIGNMENT'])
        aligned_seqs = pips.fasta.Read(input_dict['OUTPUT_ALIGNMENT'])
    else:
        to_be_aligned = 'seqs_to_be_aligned.fasta'
        if os.path.isfile(to_be_aligned):
            raise IOError("File %s already exists; please delete this file." % to_be_aligned)
        logfile.write("Writing the sequences to align to 'seqs_to_be_aligned.fasta'")
        pips.fasta.Write(sequences, to_be_aligned)
        logfile.write("Now constructing the MSA (this may take a while)...\n")
        logfile.flush()
        aligned_seqs = aligner.Align(sequences)
        logfile.write("Completed aligning the sequences.  Now writing this alignment to the OUTPUT_ALIGNMENT file of %s.\n" % (input_dict['OUTPUT_ALIGNMENT']))
        logfile.flush()
        pips.fasta.Write(aligned_seqs, input_dict['OUTPUT_ALIGNMENT'])
    logfile.write("\nHere are the contents of the OUTPUT_ALIGNMENT file of %s:\n\n******\n%s\n*****\n\n" % (input_dict['OUTPUT_ALIGNMENT'], open(input_dict['OUTPUT_ALIGNMENT']).read()))
    
    # Now build the phylogenetic tree
    phylipseqfile = 'phylipsequences.temp' # temporary file to which sequences are written in PHYLIP format
    phylipdistfile = 'phylipdistances.temp' # temporary file to which sequences are written in PHYLIP format
    logfile.write("\nWe will now use the PHYLIP program specified by PHYLIP_PROG to build a phylogenetic tree of the aligned sequences.\n")
    logfile.flush()
    (treetype, phylip_path, molecular_clock, neighbor_joining) = input_dict['PHYLIP_PROG']
    if treetype == 'DISTANCE':
        logfile.write("The phylogenetic tree is a pairwise distance based tree.  First building the protein distance matrix (this may take a while)...\n")
        logfile.flush()
        try:
            pips.phylip.WritePhylipSequenceFile(aligned_seqs, phylipseqfile)
            protdistmatrix = pips.phylip.Protdist(phylipseqfile, phylip_path)
        finally:
            if os.path.isfile(phylipseqfile):
                os.remove(phylipseqfile) # remove the temporary phylip sequence file
        logfile.write("Completed constructing the protein distance matrix.  The contents of this distance matrix are as follows:\n\n******\n%s\n*******\n\n" % protdistmatrix)
        logfile.write("Now using the protein distance matrix to construct the phylogenetic tree. ")
        if molecular_clock:
            logfile.write("This tree will be constructed assuming a molecular clock. ")
        else:
            logfile.write("The tree will be constructed without assuming a molecular clock (it will be later be rooted using the outgroup sequence). ")
        if neighbor_joining:
            logfile.write("The tree will be constructed using neighbor joining or the UPGMA method.\n")
        else:
            logfile.write("No neighbor joining or the UPGMA will be used.\n")
        logfile.write("\nBeginning construction of the phylogenetic tree (this may take a while)...\n")
        logfile.flush()
        try:
            open(phylipdistfile, 'w').write(protdistmatrix)
            if outgroup:
                newick_tree = pips.phylip.DistanceTree(phylipdistfile, phylip_path, molecular_clock, neighbor_joining, outgroup=len(aligned_seqs), root_tree=True)
            else:
                newick_tree = pips.phylip.DistanceTree(phylipdistfile, phylip_path, molecular_clock, neighbor_joining, outgroup=None, root_tree=False)
        finally:
            os.remove(phylipdistfile) # remove the temporary phylip distance matrix file
        logfile.write("Completed construction of the phylogenetic tree.  This tree is given in Newick format below:\n\n******\n%s\n*****\n\n" % newick_tree)
        if outgroup:
            logfile.write("The phylogenetic tree, as well as the protein sequence alignment written to OUTPUT_ALIGNMENT contain the outgroup sequence used to root the tree.  We will now remove this outgroup sequence from the phylogenetic tree, and from the file OUTPUT_ALIGNMENT.\n")
            assert aligned_seqs[-1][0] == 'OUTGROUP'
            aligned_seqs = aligned_seqs[ : -1]
            pips.fasta.Write(aligned_seqs, input_dict['OUTPUT_ALIGNMENT'])
            newick_tree = pips.phylip.RootByRemovingOutgroup(newick_tree)
        logfile.write("The phylogenetic tree in Newick format is being written to the OUTPUT_TREE file of %s.\n" % input_dict['OUTPUT_TREE'])
        open(input_dict['OUTPUT_TREE'], 'w').write(newick_tree)
    else:
        raise ValueError("PHYLIP_PROG specifies invalid tree type of %s" % treetype)

    # Program is complete.
    logfile.write("\nCompleted execution of pips_build_tree_and_alignment.py at %s.\n" % time.asctime())
    logfile.close()


# run the script
main()
