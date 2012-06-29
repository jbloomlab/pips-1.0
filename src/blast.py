"""Module for running BLAST.

Written by Jesse Bloom, 2008."""


import os
import sys
import tempfile
from pips import fasta


def BlastClust(headers_seqs, identities, length_coverage, seqtype, blast_dir, ncpus=1):
    """Runs 'blastclust' program to cluster sequences.

    The method runs the 'blastclust' program to cluster sequences.
    'headers_seqs' should be a list of the sequences to be clustered.  Each
        entry is a 2-tuple '(head, seq)' where 'head' is a header for the sequence
        and 'seq' is the sequence.
    'identities' is a number between 3.0 and 100.0 specifying the percent identity
        that sequences must have to be put in the same cluster.  More precisely,
        it has the exact same meaning as the '-S' option to 'blastclust'.  This
        is the identity cutoff if it is >= 3.0; otherwise it is the
        bit score / length.  Here we assume that it is the identity, and
        >= 3.0.
    'length_coverage' is a number giving the minimum length coverage.  It
        should be between 0.0 and 1.0 (the 'blastclust' default is 0.9).
        More precisely, the is the '-L' argument to 'blastclust'.
    'seqtype' specifies whether the sequence is a DNA or protein sequence.
        It should have a value of either "DNA" or "PROT".  More precisely,
        if 'seqtype' is "DNA" then blastclust is run with '-p F'; otherwise
        it is run with '-p T'.
    'blast_dir' gives the directory where the 'blastclust' executable can
        be found.
    'ncpus' specifies the number of CPU's to use when running 'blastclust'.
        By default, it is one.  It can be set to another integer value greater
        than one if you want to use more CPU's.  This is the same as the '-a'
        option to 'blastclust'.
    The returned value is the list 'clustered_headers_seqs'.  Each entry
        in this list is another list, containing all of the '(head, seq)'
        pairs found in a specific cluster.
    """
    exe = "%s/blastclust" % blast_dir
    if not os.path.isfile(exe):
        raise ValueError("Cannot find blastclust executable of %s." % exe)
    assert 0 <= length_coverage <= 1.0
    assert 3.0 <= identities <= 100.0
    assert isinstance(ncpus, int) and ncpus > 0
    if seqtype == 'DNA':
        seqtypeswitch = '-p F'
    elif seqtype == 'PROT':
        seqtypeswitch = '-p T'
    else:
        raise ValueError("Invalid seqtype of %s." % seqtype)
    currdir = os.getcwd()
    tempdir = tempfile.mkdtemp()
    try:
        # do stuff in temporary directory
        os.chdir(tempdir)
        f = open('in.fasta', 'w')
        for iseq in range(len(headers_seqs)):
            (head, seq) = headers_seqs[iseq]
            f.write('>SEQ%d\n%s\n' % (iseq, seq))
        f.close()
        (stdin, stdout_stderr) = os.popen4('%s -i in.fasta -o out.fasta -L %f -S %f %s -a %d' % (exe, length_coverage, identities, seqtypeswitch, ncpus))
        stdin.close()
        output = stdout_stderr.read()
        stdout_stderr.close()
        if not os.path.isfile('out.fasta'):
            raise ValueError("'blastclust' failed to produce output file.  Errors are:\n%s" % output)
        clusters = open('out.fasta').readlines()
    finally:
        # return to current directory, delete temporary directory
        os.chdir(currdir)
        for file in os.listdir(tempdir):
            os.remove("%s/%s" % (tempdir, file))
        os.rmdir(tempdir)
    # Construct list giving clustered sequences
    clustered_headers_seqs = []
    seq_clustered = [False] * len(headers_seqs)
    for cluster in clusters:
        headers_seqs_in_cluster = []
        for seqid in cluster.split():
            iseq = int(seqid[3 : ])
            if seq_clustered[iseq]:
                raise ValueError("A sequence appears in multiple clusters.  Output was:\n%s" % output)
            seq_clustered[iseq] = True
            headers_seqs_in_cluster.append(headers_seqs[iseq])
        clustered_headers_seqs.append(headers_seqs_in_cluster)
    if False in seq_clustered:
        raise ValueError("At least one sequence was not clustered.  Output was:\n%s" % output)
    return clustered_headers_seqs


def GetCentralProteinSeqs(clustered_headers_seqs, AlignFunction):
    """This function is designed to find the "central" sequence in clusters of protein sequences.

    The input argument 'clustered_headers_seqs' is as would be returned by 'BlastClust'.  It is
        a list, and each item is a list containing one or more entries of the form '(head, seq').
        These entries represent the headers and protein sequences of all sequences in the cluster.
    The input argument 'AlignFunction' should be a callable function that can take as a single
        argument each sequence cluster ('clustered_headers_seqs[i]') and returns a copy of this
        argument where the sequences are all aligned to be of the same length by insertion of
        gap ('-') characters as appropriate.
    The returned value is the list 'central_headers_seqs'.  'central_headers_seqs[i]' is the 
        '(head, seq)' tuple that specifies the single central sequence in cluster
        'clustered_headers_seqs[i]'.  All sequences are returned in full uppercase.
    The central sequence is determined as follows:
        * if there is only one protein sequence in a cluster, it is obviously the central sequence
        * if there are multiple protein sequences in a cluster, all redundant sequences are first
            removed by retaining only the first of redundante sequences.  The remaining sequences
            are then aligned using 'AlignFunction'.  The Hamming distance between all pairs of
            sequences is then determined (an ambiguous nucleotide of 'X' is considered to have
            a Hamming distance of 1 even if it aligns with another 'X'), and the sequence
            with the lowest average Hamming distance is retained as the central sequence.
    """
    central_headers_seqs = []
    for cluster in clustered_headers_seqs:
        assert isinstance(cluster, list)
        cluster = fasta.PurgeDuplicates(cluster) # remove any duplicate sequences in the cluster
        if len(cluster) == 1: # only one sequence in the cluster; it is the central one
            central = cluster[0]
        elif len(cluster) > 1: # multiple sequences in the cluster
            cluster = AlignFunction(cluster)                    
            cluster = [(head, seq.upper()) for (head, seq) in cluster] # convert to all upper case
            lowest_hamming = lowest_head_seq = None
            for iseq in range(len(cluster)):
                (head, seq) = cluster[iseq]
                iavg = 0.0
                for jseq in range(len(cluster)):
                    if jseq != iseq:
                        seq2 = cluster[jseq][1]
                        assert len(seq) == len(seq2)
                        ndiffs = len([1 for (r1, r2) in zip(seq, seq2) if r1 != r2 or r1 == 'X' or r2 == 'X'])
                        iavg += ndiffs
                iavg /= (len(cluster) - 1.0)
                if lowest_hamming == None or iavg < lowest_hamming:
                    lowest_hamming = iavg
                    lowest_head_seq = (head, seq)
            central = lowest_head_seq
        else:
            raise ValueError("Cluster contains no sequences.")
        central_headers_seqs.append(central)
    return central_headers_seqs


# Test module with doctest
if __name__ == '__main__':
    import doctest
    doctest.testmod()
