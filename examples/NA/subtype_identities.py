"""This script determines the sequence identities within and between subtypes."""


import sys
import random
import pips.fasta
import pips.align
import pips.stats


def main():
    """Main body of script."""
    probconspath = '/Users/bloom/probcons/'
    refsubtype = 'N1'
    maxpersubtype = 100
    ref_seqs = pips.fasta.Read('%s.fasta' % (refsubtype))
    random.shuffle(ref_seqs)
    nref = min([maxpersubtype, len(ref_seqs)])
    subtypes = ['N%d' % i for i in range(1, 10)]
    for subtype in subtypes:
        subtype_seqs = pips.fasta.Read('%s.fasta' % (subtype))
        random.shuffle(subtype_seqs)
        identities = []
        ref_identities = []
        n = min([maxpersubtype, len(subtype_seqs)])
        for i in range(n):
            si = subtype_seqs[i]
            for j in range(i + 1, n):
                sj = subtype_seqs[j]
                alignment = pips.align.Align([si, sj], probconspath)
                identities.append(pips.align.PairwiseStatistics(alignment)[0])
            for j in range(nref):
                jref = ref_seqs[j]
                alignment = pips.align.Align([si, jref], probconspath)
                ref_identities.append(pips.align.PairwiseStatistics(alignment)[0])
        identity = pips.stats.Mean(identities)
        ref_identity = pips.stats.Mean(ref_identities)
        print "\nThe average pairwise identity of sequences within subtype %s is %.2f." % (subtype, identity)
        print "The average pairwise identity of sequences within subtype %s to those in subtype %s is %.2f." % (subtype, refsubtype, ref_identity)
        sys.stdout.flush()


main() # run the script

