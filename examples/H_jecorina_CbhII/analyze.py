"""Short results analysis script.

Jesse Bloom, 2009."""


import pips.ddg_inference


def main():
    """Main body of the script."""
    seqrange = (1, 36) # look at mutations to residues in this range (inclusive)
    nbest = 20 # print this many top mutations
    predictionsfile = 'PIPS_PREDICTIONS_PRIORS-REGULARIZING.txt'
    ddgs = pips.ddg_inference.ReadDDGs(predictionsfile)[1]
    sorted_ddgs = pips.ddg_inference.SortedDDGList(ddgs)
    sorted_ddgs = [(ddg, mut) for (ddg, mut) in sorted_ddgs if seqrange[0] <= int(mut[1 : -1]) <= seqrange[1]]
    i = 0
    print "\nListing the predicted %d best mutations from PIPS:" % (nbest)
    for i in range(nbest):
        (ddg, mutation) = sorted_ddgs[i]
        print "%d best is %s (predicted ddG = %.2f)" % (i + 1, mutation, ddg)


main() # run the main program
