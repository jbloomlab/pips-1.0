"""This module contains python code to run PHYLIP.

Written by Jesse Bloom, 2007."""
#------------------------------------------------------------------
import tempfile
import os
import shutil
import re
from pips import tree
from pips import align
#------------------------------------------------------------------
def WritePhylipSequenceFile(sequences, filename):
    """Writes a sequence file in Phylip format.

    'sequences' is a list of sequences and names, all of the same length.
        Each entry is the 2-tuple (name, seq), with name being the name
        and seq being the sequence.
    'filename' is a file name.  If it does not already exist, it will be
        created.  If it already exists, it will be overwritten.
    Writes 'sequences' to file in the correct format for input to Phylip
        format.  The first line contains the number of sequences and then
        the sequence length.  Each subsequent line contains a sequence, with the 
        first 10 characters consisting of the name (or the first 10 characters
        thereof), and the rest of the line being the sequence.
    """
    f = open(filename, 'w')
    length = len(sequences[0][1])
    f.write("%d %d\n" % (len(sequences), length))
    for (name, seq) in sequences:
        assert length == len(seq)
        name = name + '          '
        name = name[ : 10]
        f.write("%s%s\n" % (name, seq))
    f.close()
#------------------------------------------------------------------
def ReadPhylipSequenceFile(filename):
    """Reads a sequence file in Phylip format.

    'filename' is the name of an existing file in the Phylip format,
        as described in 'WritePhylipSequenceFile'
    The returned variable, 'sequences', is a list of 2-tuples of the
        form '(name, seq)' where 'name' is the sequence name and 'seq'
        is the sequence.
    """
    lines = open(filename).readlines()[1 : ]
    seqlength = None
    sequences = []
    for line in lines:
        assert len(line) > 10
        name = line[ : 10].strip()
        seq = line[10 : ].strip()
        if seqlength == None:
            seqlength = len(seq)
        elif seqlength != len(seq):
            raise ValueError, "Lines not all of the same length."
        sequences.append((name, seq))
    return sequences
#------------------------------------------------------------------
def RootByRemovingOutgroup(unrooted_newick):
    """This function roots a tree and removes the outgroup sequence.

    This function takes as input an unrooted tree that has been created
        using an outgroup sequence.  This tree is a string in newick
        format.  Typically, this would be the result of calling
        'DistanceTree' or 'MaximumLikelihoodTree' with 'molecular_clock'
        set to 'False' and an outgroup specified.  The outgroup sequence
        should be the last one listed in the newick string, which
        will be the case for trees returned by this function.
    The returned value of this function is a new newick string.  This
        string corresponds to the same tree that was used as input,
        except that the outgroup sequence has been removed and a
        root has been placed on the node connected to the outgroup.
    """
    assert isinstance(unrooted_newick, str)
    new_tree = unrooted_newick.split(',')
    new_tree = ('%s;' % (','.join(new_tree[ : -1])))[1 : ]
    if new_tree[0] != '(' and new_tree[-2] != ')':
        new_tree = '(%s);' % new_tree[ : -1]
    return new_tree
#------------------------------------------------------------------
def RootToOutgroup(unrooted_newick):
    """This function converts an unrooted tree to be drawn in the rooted form.

    This function takes as input an unrooted tree in newick string format.  Typically
        this would be the output of calling 'DistanceTree' or 'MaximumLikelihoodTree'
        with 'molecular_clock' set to False.  Those functions will produce
        unrooted trees.  This method roots those trees by attaching the root
        to the last node listed in the newick tree string.  Note that if either
        'DistanceTree' or 'MaximumLikelihoodTree' were called with 'outgroup'
        specified, then this last node will be the outgroup, so the tree will be
        rooted at the outgroup, as it should be.
    The returned variable is a newick string for the rooted tree.
    """
    new_tree = unrooted_newick.split(',')
    new_tree[0] = "(%s" % new_tree[0] # add leading parenthesis
    new_tree[-2] = "%s):0.0" % new_tree[-2] # add trailing parenthesis and branch length
    new_tree = ','.join(new_tree)
    return new_tree
#------------------------------------------------------------------
def EliminateNegativeBranchLengths(newick_tree):
    """This function sets any negative branch lengths to zero.

    This function takes as input a single argument specifying a
        phylogenetic tree in newick format.  It checks to see if any
        branches have negative lengths, and if they do, it sets these
        lengths to zero.  This would be of use if you were calling
        'DistanceTree' with the 'neighbor_joining' method, since
        in that case the program does not guarantee positive branch lengths.
    The returned value is a newick tree with all negative branch lengths
        set to zero.
    """
    entries = newick_tree.split(':')
    for ientry in range(len(entries)):
        entry = entries[ientry]
        if entry[0] == '-':
            # negative branch length
            if (',' in entry) and (')' in entry):
                end = min([entry.index(','), entry.index(')')])
            elif ',' in entry:
                end = entry.index(',')
            elif ')' in entry:
                end = entry.index(')')
            else:
                raise ValueError, "Can't find end to negative branchpoint in %s in %s." % (entry, newick_tree)
            entry = "0.0%s" % entry[end : ]
            entries[ientry] = entry
    return ':'.join(entries)
#------------------------------------------------------------------
def DistanceTree(infile, progpath, molecular_clock, neighbor_joining, jumble=None, outgroup=None, set_negative_branch_lengths_to_zero=True, root_tree=True):
    """This method computes a tree using the Fitch-Margoliash distance method.

    Uses either the fitch or kitsch program from the Phylip package.
    The returned variable is a string giving the inferred tree in Newick format.
    'infile' should specify the name of a file containing a distance matrix.
        For proteins, this is computed using the 'Protdist' method.
        For nucleotide sequences, this is computed using the 'DNAdist' method.
    'progpath' specifies the path to the directory where the Phylip executables
        'fitch', 'kitsch', and 'neighbor' are located.  
    'molecular_clock' specifies whether we assume a molecular clock.  If it is
        'True', we do assume a clock and use the kitsch program.  If it is 'False',
        we do not assume a clock and use the fitch program.
    'neighbor_joining' specifies whether we use a neighbor joining method.  Such
        methods are reported to be slightly less accurate, but much faster.  If
        it is 'False', then we do not use neighbor joining, and instead use the
        standard 'fitch' (no molecular clock) or 'kitcsh' (molecular clock) programs.
        If it is 'True', then we instead use the 'neighbor' program.  With a molecular
        clock, we use the UPGMA method, without a molecular clock we use the 
        neighbor joining method.
    'jumble' specifies whether we randomize the order of the input sequences.  The
        resulting tree will in general change (hopefully only slightly!) if the
        input sequences have their order jumbled.  By default, this value
        is set to 'None', which means that no jumbling is done.  If it is set to
        another value, the it should be a 2-tuple '(seed, njumbles)'.  'seed'
        must be an odd integer >= 1.  'njumbles' is the number of different jumbles
        we look at to find the best one.  It must be >= 1.  If neighbor_joining is
        True, then 'njumbles' is irrelevant since there is not an option for 
        multiple jumbles.
    'outgroup' is an option that is meaningful if and only if we are not using a 
         molecular clock.  In this case, it specifies which sequence to use as the
         outgroup.  It should be a number between 1 and the number of sequences,
         giving the line number of that sequence in the list of sequences used
         to generate the distance infile.  If we are not using an outgroup,
         it should be 'None'.
    'set_negative_branch_lengths_to_zero' is a boolean switch.  It is meaningful only if
        'neighbor_joining' is set to 'True'.  The neighbor joining algorithm does not 
        guarantee that all branch lengths are positive, so there is a possibility for having
        short negative branch lengths.  If this option is set to its default value of 'True',
        any such negative branch lengths are set to zero by calling 'EliminateNegativeBranchLengths'.
    'root_tree' is a boolean option that is meaningful if 'molecular_clock' is 'False'.  Calling the
        phylip programs with 'molecular_clock' set to 'False' produces an unrooted True.  When
        'root_tree' is at its default value of 'True', this unrooted tree is rooted by calling
        'RootToOutgroup'.  This means that if 'outgroup' is 'None' then the tree is rooted to
        a random sequence. and if 'outgroup' is set to a number, then the tree is rooted to
        this outgroup.
    """
    if not os.path.isfile(infile):
        raise IOError, "Cannot find infile %s." % infile
    assert os.path.isdir(progpath), "Cannot find directory %s." % progpath
    assert isinstance(molecular_clock, bool)
    assert isinstance(set_negative_branch_lengths_to_zero, bool)
    assert (jumble == None) or (isinstance(jumble, tuple) and len(jumble) == 2 and isinstance(jumble[0], int) and jumble[0] >= 1 and jumble[0] % 2 and isinstance(jumble[1], int) and jumble[1] >= 1)
    assert (outgroup == None) or ((not molecular_clock) and isinstance(outgroup, int) and outgroup > 0)
    assert isinstance(neighbor_joining, bool)
    if neighbor_joining:
        exe = os.path.abspath("%s/neighbor" % progpath) # the neighbor executable
    elif molecular_clock:
        exe = os.path.abspath("%s/kitsch" % progpath) # the kitsch executable
    else:
        exe = os.path.abspath("%s/fitch" % progpath) # the fitch executable
    if not os.path.isfile(exe):
        raise IOError, "Cannot find executable at %s." % exe
    currdir = os.getcwd()
    tempdir = tempfile.mkdtemp()
    try:
        # do stuff in the temporary directory
        shutil.copy(infile, "%s/infile" % tempdir)
        os.chdir(tempdir)
        (stdin, stdout_stderr) = os.popen4(exe) # run the program
        if jumble:
            stdin.write('J\n')
            stdin.write("%d\n" % jumble[0])
            if not neighbor_joining:
                stdin.write("%d\n" % jumble[1])
        if outgroup:
            stdin.write('O\n')
            stdin.write("%d\n" % outgroup)
        if neighbor_joining and molecular_clock:
            stdin.write('N\n') # set to use the UPGMA method
        stdin.write('Y\n') # accept parameter choices
        stdin.close() # all input complete
        output = stdout_stderr.read() # get the output and errors
        stdout_stderr.close()
        if not os.path.isfile('outtree'):
            raise ValueError, "Failed to create outfile, output message is %s" % output
        newick_tree = ''.join([line.strip() for line in open('outtree').readlines()])
        if not newick_tree:
            raise ValueError, "Created an empty tree, output message is %s" % output
    finally: 
        # return to the current directory
        os.chdir(currdir)
        # remove all files from the temporary directory
        for file in os.listdir(tempdir):
            os.remove("%s/%s" % (tempdir, file))
        # remove the temporary directory
        os.rmdir(tempdir)
    if set_negative_branch_lengths_to_zero and neighbor_joining: # set negative branch lengths to zero
        newick_tree = EliminateNegativeBranchLengths(newick_tree)
    if root_tree and not molecular_clock: # root the tree
        newick_tree = RootToOutgroup(newick_tree)
    return newick_tree
#------------------------------------------------------------------
def DNAdist(infile, progpath, cv=1.0):
    """This method computes a nucleotide sequence distance matrix.

    Runs the DNAdist program from the Phylip package.
    Uses the F84 model to compute distances.
    The transition / transversion ratio is assumed to be 2.0.
    Rates are assumed to be gamma distributed among positions.
    The program returns the contents of the output file containing
        the distance matrix as a string.
    'infile' specifies the text file listing the sequences in the standard Phylip
        input format.  This is the first line containing two integers (the number
        of sequences and the length of these sequences), and each subsequent line
        consisting of the ten characters specifying the sequence name followed by
        the sequence.
    'progpath' specifies the path to the directory where the Phylip executable
        'dnadist' is located.  
    'cv' specifies the coefficient of variation for unequal rates of change at
        different positions.  By default, this is set to 1.0.
    """
    if not os.path.isfile(infile):
        raise IOError, "Cannot find infile %s." % infile
    assert os.path.isdir(progpath), "Cannot find directory %s." % progpath
    exe = os.path.abspath("%s/dnadist" % progpath) # the DNAdist executable
    if not os.path.isfile(exe):
        raise IOError, "Cannot find executable at %s." % exe
    currdir = os.getcwd()
    tempdir = tempfile.mkdtemp()
    try:
        # do stuff in the temporary directory
        shutil.copy(infile, "%s/infile" % tempdir)
        os.chdir(tempdir)
        (stdin, stdout_stderr) = os.popen4(exe) # run the program
        stdin.write('G\n') # gamma distribution of rates
        stdin.write('Y\n') # accept parameter choices
        stdin.write('%f\n' % cv) # set the coefficient of variation
        stdin.close() # all input complete
        output = stdout_stderr.read() # get the output and errors
        stdout_stderr.close()
        if not os.path.isfile('outfile'):
            raise ValueError, "Failed to create outfile, output message is %s" % output
        outfile = open('outfile').read()
    finally: 
        # return to the current directory
        os.chdir(currdir)
        # remove all files from the temporary directory
        for file in os.listdir(tempdir):
            os.remove("%s/%s" % (tempdir, file))
        # remove the temporary directory
        os.rmdir(tempdir)
    return outfile
#------------------------------------------------------------------
def Protdist(infile, progpath, cv=1.0):
    """This method computes a protein sequence distance matrix.

    Runs the Protdist program from the Phylip package.
    Uses the JTT model to compute distances.
    Rates are assumed to be gamma distributed among positions.
    The program returns the contents of the output file containing
        the distance matrix as a string.
    'infile' specifies the text file listing the sequences in the standard Phylip
        input format.  This is the first line containing two integers (the number
        of sequences and the length of these sequences), and each subsequent line
        consisting of the ten characters specifying the sequence name followed by
        the sequence.
    'progpath' specifies the path to the directory where the Phylip executable
        'protdist' is located.  
    'cv' specifies the coefficient of variation for unequal rates of change at
        different amino acid positions.  By default, this is set to 1.0.
    """
    if not os.path.isfile(infile):
        raise IOError, "Cannot find infile %s." % infile
    assert os.path.isdir(progpath), "Cannot find directory %s." % progpath
    exe = os.path.abspath("%s/protdist" % progpath) # the Protdist executable
    if not os.path.isfile(exe):
        raise IOError, "Cannot find executable at %s." % exe
    currdir = os.getcwd()
    tempdir = tempfile.mkdtemp()
    try:
        # do stuff in the temporary directory
        shutil.copy(infile, "%s/infile" % tempdir)
        os.chdir(tempdir)
        (stdin, stdout_stderr) = os.popen4(exe) # run the program
        stdin.write('G\n') # gamma distribution of rates
        stdin.write('Y\n') # accept parameter choices
        stdin.write('%f\n' % cv) # set the coefficient of variation
        stdin.close() # all input complete
        output = stdout_stderr.read() # get the output and errors
        stdout_stderr.close()
        if not os.path.isfile('outfile'):
            raise ValueError, "Failed to create outfile, output message is %s" % output
        outfile = open('outfile').read()
    finally: 
        # return to the current directory
        os.chdir(currdir)
        # remove all files from the temporary directory
        for file in os.listdir(tempdir):
            os.remove("%s/%s" % (tempdir, file))
        # remove the temporary directory
        os.rmdir(tempdir)
    return outfile
#------------------------------------------------------------------
def MaximumLikelihoodTree(infile, progpath, molecular_clock, cv_ncat=(1.0, 5), run_quality=0, outgroup=None, copy_outfile=False):
    """Estimates the maximum likelihood phylogenetic tree for protein sequences.

    Uses either the proml or promlk program from Phylip.
    The returned variable is a string giving the inferred tree in Newick format.
    'infile' specifies the text file listing the sequences in the standard Phylip
        input format.  This is the first line containing two integers (the number
        of sequences and the length of these sequences), and each subsequent line
        consisting of the ten characters specifying the sequence name followed by
        the sequence.
    'progpath' specifies the path to the directory where the Phylip executables
         'proml' and 'promlk' are located. 
    'molecular_clock' specifies whether or not we assume a molecular clock.  If it 
        is 'True', then we assume a molecular clock and use the 'promlk' program.
        If it is 'False', then we do not assume a molecular clock, and use the
        'proml' program.
    'cv_ncat' specifies the coefficient of variation among the rates and the number of
        categories.  This variable should be a 2-tuple, with the first element
        a number (integer or float) and the second an integer.  The first number specifies
        the gamma distribution coefficient of variation.  The second number specifies
        the number of categories of rates.  The default value is (1.0, 5).
    'run_quality' specifies options related to the quality of the search.  The 
        default value is 0, which gives the fastest search.  Setting the value 
        to 1 means that a global search is performed, removing and re-adding each
        branch (Phylip option G).  According to the Phylip documentation, this triples 
        the run time. 
    'outgroup' gives the option of specifying an outgroup sequence if 'molecular_clock'
        is set to 'False'.  This is equivalent to the "O" option of 'proml'.  By default,
        'outgroup' is 'None', which means that no outgroup is specified.  If it is set to
        another value, it should be set to an integer >= 1 and <= the number of sequences
        in 'infile' giving the number of the sequence which is the outgroup.
    'copy_outfile' is an option specifying that we copy the Phylip output file to another file
        so that it can be examined/parsed.  By default, this option is 'False', meaning
        that the Phylip 'outfile' is simply deleted.  It can be set to another value
        giving a string representing a valid file name.  In this case, the 'outfile' from
        Phylip is copied to this file.  If an absolute path name is not given for
        'copy_outfile', it is copied to the current directory.
    """
    assert os.path.isdir(progpath), "Cannot find directory %s." % progpath
    if not os.path.isfile(infile):
        raise IOError, "Cannot find infile %s." % infile
    assert isinstance(molecular_clock, bool)
    if molecular_clock:
        exe = os.path.abspath("%s/promlk" % progpath) # the Proml executable
    else:
        exe = os.path.abspath("%s/proml" % progpath) # the Promlk executable
    if not os.path.isfile(exe):
        raise IOError, "Cannot find executable at %s." % exe
    (cv, ncat) = cv_ncat
    assert isinstance(cv, (int, float)) and cv >= 0.0
    assert isinstance(ncat, int) and ncat > 0
    assert run_quality in [0, 1]
    currdir = os.getcwd()
    tempdir = tempfile.mkdtemp()
    try:
        # do stuff in the temporary directory
        shutil.copy(infile, "%s/infile" % tempdir)
        os.chdir(tempdir)
        (stdin, stdout_stderr) = os.popen4(exe) # run the program
        if run_quality == 0:
            pass
        elif run_quality  == 1:
            stdin.write('G\n') # global rearrangements
        else:
            raise ValueError, "Invalid value of run_quality."
        if outgroup != None:
            if molecular_clock:
                raise ValueError("Cannot specify outgroup and molecular clock.")
            stdin.write('O\n')
            stdin.write('%d\n' % outgroup)
        stdin.write('1\n') # print input sequences to outfile
        stdin.write('5\n') # print node sequences to outfile
        stdin.write('R\n') # gamma distributed rates
        stdin.write('Y\n') # accept these options
        stdin.write("%f\n" % cv) # coefficient of variation
        stdin.write("%d\n" % ncat) # number of rate categories
        stdin.close() # all input complete
        output = stdout_stderr.read() # get the output and errors
        stdout_stderr.close()
        if not os.path.isfile('outtree'):
            raise ValueError, "Failed to create outtree, output message is %s" % output
        newick_tree = ''.join([line.strip() for line in open('outtree').readlines()])
        if not newick_tree:
            raise ValueError, "Created an empty tree, output message is %s" % output
        if copy_outfile:
            if not os.path.isabs(copy_outfile):
                shutil.copy('outfile', "%s/%s" % (currdir, copy_outfile))
            else:
                shutil.copy('outfile', copy_outfile)
    finally: 
        # return to the current directory
        os.chdir(currdir)
        # remove all files from the temporary directory
        for file in os.listdir(tempdir):
            os.remove("%s/%s" % (tempdir, file))
        # remove the temporary directory
        os.rmdir(tempdir)
    return newick_tree
#------------------------------------------------------------------
def BuildTreeFromPhylipOutfile(outfile, newick_tree, seqnames, rootbyremovinglast, allcaps=True):
    """This function creates a 'tree.Tree' object from a PHYLIP output file.

    This function can be used to extract information from phylogenetic trees
        produced by PHYLIP.  It is currently known to work on trees produced
        by the PHYLIP 'proml' program in which both the input and the inferred
        interior node sequences are printed to the output file.  This method
        is currently only guaranteed to work on rooted trees (or trees that
        can be rooted using the 'rootbyremovinglast' option).
    'outfile' should be the output file produced by 'proml'.  It should contain
        both the input sequences and the inferred interior nodes.
    'newick_tree' should be the Newick string giving the tree inferred by 
        PHYLIP when it produced 'outfile'.  This should be a string in
        Newick format.
    'seqnames' can be used to map the sequence names used when running
        PHYLIP to some other set of sequence names.  This would be
        desirable if names have been shortened/altered when running
        PHYLIP because of PHYLIP's requirement that names be only 
        10 characters long.  If 'seqnames' is 'None', no mapping 
        is done.  Otherwise, 'seqnames' should be a dictionary
        keyed by the sequence names used by PHYLIP, and with values
        being the desired sequence names.
    'rootbyremovinglast' is a Boolean switch.  If the tree specified
        in 'outfile' and by 'newick_tree' is already rooted, then this
        option should be set to 'False'.  Otherwise, the tree can be
        rooted by removing a specified outgroup sequence.  In this case,
        'rootbyremovinglast' should be set to the name of the string
        specifying the outgroup sequence in 'outfile' and 'newick_tree'.
        This outgroup is assumed to be the last sequence in 'newick_tree'.
        The tree is rooted by removing this outgroup sequence.
    'allcaps' is a Boolean switch specifying whether we convert all interior
        node sequences to all caps (we do if 'True').  This may be important,
        because Phylip may specify the caps state depending on the confidence
        of the reconstruction.
    The returned value is a phylogenetic tree as a pips.tree.Tree object.  The
        tip nodes have the names specified in 'seqnames', or the names in the
        PHYLIP outfile if 'seqnames' is 'None'.  The interior nodes have
        the names automatically assigned by PHYLIP.  Both tip and interior nodes
        have the sequences parsed from the PHYLIP output file.
    """
    lines = open(outfile).readlines()
    nlines = len(lines)
    iline = 0
    # Get sequences of tip nodes
    while lines[iline] != 'Name                                     Sequences\n':
        iline += 1
        if iline == nlines:
            raise IOError("Never found input sequences in %s." % outfile)
    if not (lines[iline + 1] == '----                                     ---------\n' and lines[iline + 2] == '\n'):
        raise IOError("Not sure if the input sequences were found in %s." % outfile)
    iline += 3
    tipseqs_dict = {}
    firstseq = None
    firstinterleave = True
    while not lines[iline].isspace():
        while not lines[iline].isspace():
            (name, seq) = lines[iline].split(' ', 1)
            seq = seq.strip()
            seq = seq.replace(' ', '')
            if firstseq == None:
                firstseq = name
            if firstinterleave:
                if name in tipseqs_dict:
                    raise ValueError("Multiple sequences with name %s." % name)
                tipseqs_dict[name] = seq
            else:
                tipseqs_dict[name] += seq
            iline += 1
        firstinterleave = False
        iline += 1
    tipseqs = [(firstseq, tipseqs_dict[firstseq])] + [(head, seq) for (head, seq) in tipseqs_dict.iteritems() if head != firstseq]
    tipseqs = align.RemoveDots(tipseqs)
    if rootbyremovinglast:
        n = len(tipseqs)
        tipseqs = [(head, seq) for (head, seq) in tipseqs if head != rootbyremovinglast]
        if n != len(tipseqs) + 1:
            raise ValueError("Cannot find root sequence of %s in newick tree." % rootbyremovinglast)
        newick_tree = RootByRemovingOutgroup(newick_tree)
        newick_tree = "(%s):0.0;" % newick_tree[ : -1]
    # Get sequences of interior nodes
    while lines[iline] != 'Probable sequences at interior nodes:\n':
        iline += 1
        if iline == nlines:
            raise IOError("Never found interior node sequences in %s." % outfile)
    if not (lines[iline + 1] == '\n' and lines[iline + 2] == '  node                    Reconstructed sequence (caps if > 0.95)\n' and lines[iline + 3] == '\n'):
        raise IOError("Not sure if the interior node sequences were found in %s." % outfile)
    iline += 4
    interiornameslist = []
    interiorseqs_dict = {}
    firstinterleave = True
    while iline < nlines:
        while not lines[iline].isspace():
            line = lines[iline].strip()
            (name, seq) = line.split(' ', 1)
            if name not in tipseqs_dict: # interior node
                seq = seq.strip()
                seq = seq.replace(' ', '')
                if allcaps:
                    seq = seq.upper()
                if firstinterleave:
                    interiornameslist.append(name)
                    if name in interiorseqs_dict:
                        raise ValueError("Multiple interior sequences with name %s." % name)
                    interiorseqs_dict[name] = seq
                else:
                    interiorseqs_dict[name] += seq
            iline += 1
        firstinterleave = False
        iline += 1
    assert len(interiornameslist) == len(interiorseqs_dict) == len(tipseqs) - 1
    # build the tree
    phyliptree = tree.Tree(newick_tree, tipnames_sequences=tipseqs)
    # now proceed through the tree, assigning the interior nodes and renaming the tip nodes
    root = phyliptree.GetRoot()
    assert not root.tip and not root.name and not root.seq
    root.name = interiornameslist[0]
    root.seq = interiorseqs_dict[root.name]
    # create function for recursively traversing the tree
    def _RecursivelyTraverse(node, i_interior):
        """Recursively moves through tree, assigning interior nodes and renaming tips."""
        if node.leftdescendent.tip:
            if seqnames:
                node.leftdescendent.name = seqnames[node.leftdescendent.name]
        else:
            assert not node.leftdescendent.name and not node.leftdescendent.seq
            node.leftdescendent.name = interiornameslist[i_interior]
            node.leftdescendent.seq = interiorseqs_dict[node.leftdescendent.name]
            i_interior += 1
            _RecursivelyTraverse(node.leftdescendent, i_interior)
        if node.rightdescendent.tip:
            if seqnames:
                node.rightdescendent.name = seqnames[node.rightdescendent.name]
        else:
            assert not node.rightdescendent.name and not node.rightdescendent.seq
            node.rightdescendent.name = interiornameslist[i_interior]
            node.rightdescendent.seq = interiorseqs_dict[node.rightdescendent.name]
            i_interior += 1
            _RecursivelyTraverse(node.rightdescendent, i_interior)
    # now recursively traverse the tree
    _RecursivelyTraverse(root, 1)
    # return the tree
    return phyliptree
#------------------------------------------------------------------
