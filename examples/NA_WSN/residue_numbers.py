"""Script for converting between various residue numbering schemes for influenza NA.

Jesse Bloom, 2009."""


import pips.fasta
import pips.align


def CorrespondingResidue(alignment, r):
    """Given an alignment, returns residue in second sequence corresponding to a residue in first.

    'alignment' is a pairwise sequence alignment in the form returned by pips.fasta.align.Align.
    'r' is a number corresponding to a residue of the first sequence in alignment when this 
        sequence is numbered sequentially beginning with 1 at the N-terminal residue.
    The returned value is a number corresponding to the residue in the second sequence
        that is aligned with 'r' when the second sequence is numbered sequentially beginning
        with 1 at the N-terminal residue.  If 'r' is not aligned with another residue, the
        returned value is 'None'.
    """
    i = 0 
    r1 = r2 = 0
    assert len(alignment) == 2
    (seq1, seq2) = (alignment[0][1], alignment[1][1])
    assert len(seq1) == len(seq2)
    assert isinstance(r, int) and r > 0
    for i in range(len(seq1)):
        if seq1[i] != '-':
            r1 += 1
        if seq2[i] != '-':
            r2 += 1
        if r1 == r:
            if seq2[i] == '-':
                return None
            else:
                return r2
    raise ValueError("r value of %d exceeded length of first sequence." % r)



def main():
    probconspath = '/Users/bloom/probcons/' # path to the PROBCONS alignment program executable
    """Main body of script."""
    print "This script converts between residue numbers for various influenza NA numbering schemes."
    print "These schemes are:\n"
    print "* 3BEQ: the residue numbers in the 3BEQ PDB file of the 1918 H1N1 NA.  These correspond exactly to the residue numbers in the 2HTY PDB file of the H5N1 NA.\n"
    print "* WSN: the residue numbers in the A/WSN/33 (H1N1) NA, numbering sequentially with 1 for the N-terminal Met.\n"
    print "* N2: the residue numbers in the A/Tokyo/3/67 (H3N2) N2 neuramindase.  These correspond exactly to the residue numbers in the A/Memphsi/31/98 (H3N2) N2 neuraminidase as well.\n"
    print "* SI06: the residue numbers in the A/Solomon Islands/3/2006 (H1N1) NA, numbering sequentially with 1 for the N-terminal Met.\n"
    #
    # Residue mapping between WSN and 3BEQ sequences
    wsn_to_3beq = {} # indexed by WSN residue numbers as integers, keys are 3BEQ residue numbers as strings
    seq_wsn = pips.fasta.Read('NA_WSN.fasta')[0] # sequence of WSN NA
    seq_3beq = [seq for seq in pips.fasta.Read('3BEQ_seqs.fasta') if seq[0] == 'A'] # sequence of chain A of 3BEQ PDB file
    assert len(seq_3beq) == 1
    seq_3beq = seq_3beq[0]
    alignment = pips.align.Align([seq_wsn, seq_3beq], probconspath)
    number_mapping = dict([(int(line.split()[1]), line.split()[2]) for line in open('3BEQ_numbers.txt').readlines() if line and line[0] != '#' and not line.isspace()])
    for r in range(1, len(seq_wsn[1]) + 1):
        r_3beq = CorrespondingResidue(alignment, r)
        if r_3beq != None:
            if r_3beq in number_mapping:
                r_3beq = number_mapping[r_3beq]
            else:
                r_3beq = None
        wsn_to_3beq[r] = r_3beq
    #
    # Residue mapping between N2 and WSN sequences
    n2_to_wsn = {} # indexed by N2 residue numbers as integers, keys are WSN residue numbers as integers
    seq_n2 = pips.fasta.Read('NA_N2_A_Tokyo_3_67.fasta')[0] # sequence of N2
    alignment = pips.align.Align([seq_n2, seq_wsn], probconspath)
    for r in range(1, len(seq_n2[1]) + 1):
        r_wsn = CorrespondingResidue(alignment, r)
        n2_to_wsn[r] = r_wsn
    #
    # Residue mapping between N2 and 3BEQ sequences
    n2_to_3beq = {} # indexed by N2 residue numbers as integers, keys are 3BEQ residue numbers as strings
    for (r_n2, r_wsn) in n2_to_wsn.iteritems():
        if r_wsn == None:
            continue
        r_3beq = wsn_to_3beq[r_wsn]
        n2_to_3beq[r_n2] = r_3beq
    #
    # Residue mapping between WSN and SI06 sequences
    wsn_to_si06 = {} # indexed by WSN residue numbers as integers, keys are SI06 residue numbers as integers
    seq_si06 = pips.fasta.Read('NA_SI06.fasta')[0] # sequence of SI06
    alignment = pips.align.Align([seq_wsn, seq_si06], probconspath)
    for r in range(1, len(seq_wsn[1]) + 1):
        r_si06 = CorrespondingResidue(alignment, r)
        wsn_to_si06[r] = r_si06
    #
    # Now prompt user for input
    scheme1 = raw_input("Enter the numbering scheme for the residue of known number (WSN, 3BEQ, N2, SI06): ")
    r1 = raw_input("Enter the number of the residue in scheme %s: " % scheme1)
    scheme2 = raw_input("Enter the numbering scheme for which we want to know the identity of this residue (WSN, 3BEQ, N2, SI06): ")
    if scheme1 == 'N2' and scheme2 == 'WSN':
        r1 = int(r1)
        r2 = str(n2_to_wsn[r1])
    elif scheme1 == 'WSN' and scheme2 == '3BEQ':
        r1 = int(r1)
        r2 = wsn_to_3beq[r1]
    elif scheme1 == '3BEQ' and scheme2 == 'WSN':
        map_3beq_to_wsn = dict([(y, x) for (x, y) in wsn_to_3beq.iteritems()])
        r2 = str(map_3beq_to_wsn[r1])
    elif scheme1 == 'WSN' and scheme2 == 'SI06':
        r1 = int(r1)
        r2 = str(wsn_to_si06[r1])
    elif scheme1 == 'SI06' and scheme2 == 'WSN':
        map_si06_to_wsn = dict([(y, x) for (x, y) in wsn_to_si06.iteritems()])
        r1 = int(r1)
        r2 = str(map_si06_to_wsn[r1])
    else:
        raise ValueError("Currently this program is not implemented to convert from scheme %s to %s." % (scheme1, scheme2))
    print "Residue %d in numbering scheme %s is numbered in the scheme %s as:\n%s" % (r1, scheme1, scheme2, r2)


main() # run the main program
