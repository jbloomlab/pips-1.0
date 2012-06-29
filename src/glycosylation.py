"""Module with functions for analyzing protein glycosylation.

Written by Jesse Bloom, 2008."""


def FindNGlycanSites(headers_seqs):
    """A function for finding potential N glycosylation sites in protein primary sequences.

    'headers_seqs' is a list of sequences in the form of 2-tuples '(head, seq)' where
        'head' is the sequence header and 'seq' is the sequence as one letter amino
        acid codes.
    The function finds all potential N-linked glycosylation sites, which are motifs
        of the form Asn-X-Ser or Asn-X-Thr where X can be any amino acid
        except proline.  Any gaps in the protein sequences are ignored when looking for
        the sites, so that N-DS would be a site just like NDS.
    The returned variable is a list of the same length as 'headers_seqs'.  Each entry in
        this list is another list of numbers, giving the indices (as integers ranging from
        0 to sequence length - 1) of the Asn in the N-linked glycosylation site.

    >>> FindNGlycanSites([('seq1', 'MATNWSALNQT'), ('seq2', 'matnet')])
    [[3, 8], [3]]
    
    >>> FindNGlycanSites([('seq1', 'MATNPSAN-A')])
    [[]]
    
    >>> FindNGlycanSites([('seq1', 'MA--TN-DSA')])
    [[5]]
    """
    nglycansites = []
    for (head, seq) in headers_seqs:
        seqsites = []
        seq = seq.upper()
        seqlength = len(seq)
        for ires in range(seqlength):
            if seq[ires] == 'N':
                jres = ires + 1
                while jres < seqlength and seq[jres] == '-':
                    jres += 1
                if jres == seqlength:
                    continue
                next = seq[jres]
                jres += 1
                while jres < seqlength and seq[jres] == '-':
                    jres += 1
                if jres == seqlength:
                    continue
                next2 = seq[jres]
                if next != 'P' and next2 in ['S', 'T']:
                    seqsites.append(ires)
        nglycansites.append(seqsites)
    return nglycansites


# Test with doctest
if __name__ == '__main__':
    import doctest
    doctest.testmod()
