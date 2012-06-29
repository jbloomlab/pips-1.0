"""Merges unique Genbank and GISAID N1 NA sequences."""

import pips.fasta
import pips.align
import pips.stats

genbank = pips.fasta.Read('N1_Genbank.fasta')
print "Read %d sequences from Genbank." % len(genbank)
gisaid = pips.fasta.Read('N1_GISAID.fasta')
print "Read %d sequences from GISAID." % len(gisaid)
d = dict([(seq, head) for (head, seq) in genbank])
print "%d of the sequences in Genbank are unique." % len(d)
for (head, seq) in gisaid:
    if seq not in d:
        d[seq] = head
print "After adding sequences from GISAID, there are now %d unique sequences." % len(d)
ca09 = pips.fasta.Read("NA_CA09.fasta")[0]
cutoff = int(0.9 * len(ca09[1]))
seqs = [(head, seq) for (seq, head) in d.iteritems()]
print "Removing all sequences that are not at least %d residues in length..." % cutoff
seqs = [(head, seq) for (head, seq) in seqs if len(seq) >= cutoff]
print "After removing these sequences, there are %d sequences remaining." % len(seqs)
identity_cutoff = 0.75
gap_cutoff = 0.2
print "Now aligning all sequences to CA09 reference sequence, and removing those without at least %.2f identity, or with more than %.2f gaps." % (identity_cutoff, gap_cutoff)
new_seqs = []
identity_list = []
gap_list = []
for (head, seq) in seqs:
    a = pips.align.Align([ca09, (head, seq)], '/Users/bloom/muscle3.8/', 'MUSCLE')
    (identities, gaps) = pips.align.PairwiseStatistics(a)
    if identities < identity_cutoff or gaps > gap_cutoff:
        print "Removing sequence with %.2f identity and %.2f gaps: %s" % (identities, gaps, head)
    else:
        new_seqs.append((head, seq))
        identity_list.append(identities)
        gap_list.append(gaps)
seqs = new_seqs
print "After removing these sequences, there are %d sequences remaining." % len(seqs)
print "Among the retained sequences, the mean and median identity is %.2f and %.2f" % (pips.stats.Mean(identity_list), pips.stats.Median(identity_list))
print "Among the retained sequences, the mean and median gap fraction is %.2f and %.2f" % (pips.stats.Mean(gap_list), pips.stats.Median(gap_list))
print "Removing sequences that are substrings of other sequences."
new_seqs = []
for (head, seq) in seqs:
    i = 0
    for (ihead, iseq) in new_seqs:
        if seq in iseq:
            break
        elif iseq in seq:
            new_seqs[i] = (head, seq)
            break
        i += 1
    else:
        new_seqs.append((head, seq))
seqs = new_seqs
print "After removing these sequences, there are %d sequences remaining." % len(seqs)
print "Writing these sequences to N1.fasta"
pips.fasta.Write(seqs, 'N1.fasta')
