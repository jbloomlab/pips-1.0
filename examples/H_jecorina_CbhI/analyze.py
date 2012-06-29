"""Short results analysis script.

Jesse Bloom, 2009."""


import pips.ddg_inference


def main():
    """Main body of the script."""

    # mutations to rank -- a list of "best" candidate mutations from Pete/Arvind
    mutations_to_rank = [
        'S239A', 
        'G111K',
        'N324D',
        'L343I',
        'L383I',
        'G357D',
        'A94D',
        'A123S',
        'K439R',
        'V326R',
        'V326K',
        'S64G',
        'T373L',
        'I110L',
        'W280F',
        'L299V',
        'F328I',
        'F369M',
        'K35R',
        'F435Y',
    ]

    nbest = 20 # print this many top mutations
    for prior_center in ['CUPSAT']:
        predictionsfile = 'PIPS_PREDICTIONS_PRIORS-%s.txt' % prior_center
        ddgs = pips.ddg_inference.ReadDDGs(predictionsfile)[1]
        sorted_ddgs = pips.ddg_inference.SortedDDGList(ddgs)
        i = 0
        print "\nListing the predicted %d best mutations from PIPS:" % (nbest)
        for i in range(nbest):
            (ddg, mutation) = sorted_ddgs[i]
            print "%d best is %s (predicted ddG = %.2f)" % (i + 1, mutation, ddg)
        print "\nNow listing the PIPS predicted ddG values for the mutations identified by Arvind/Pete:"
        for mutation in mutations_to_rank:
            (wt, i, mut) = (mutation[0], int(mutation[1 : -1]), mutation[-1])
            if ddgs[i][0] != wt:
                raise ValueError("Identity mismatch at residue %d" % i)
            rank = 0
            while sorted_ddgs[rank][1] != mutation:
                rank += 1
            print "Mutation %s is predicted by PIPS to have a ddG of %.2f (ranking it %d out of all %d possible mutations)." % (mutation, ddgs[i][1][mut], rank + 1, len(sorted_ddgs))


main() # run the main program
