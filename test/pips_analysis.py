"""This script performs a pips analysis of the DDG values for a protein.

Written by Jesse Bloom, 2007."""


import cProfile
import sys
import os 
import time
import re
import random
import math
import pylab
import matplotlib
import pips.cupsat
import pips.ddg_inference
import pips.fasta
import pips.align
import pips.phylip
import pips.tree
import pips.stats


def main():
    """The body of the script, performs the ddG analysis."""

    # Define the names of required input files, and other main configuration variables
    protein_w_underscores = os.getcwd().split('/')[-1]
    protein = protein_w_underscores.replace('_', ' ')
    pdbfile = 'pdb_structure.pdb' # the name of the PDB file
    pdbchain = None # chain in pdbfile -- there is only one chain, so not relevant here
    seqfile = 'protseq.txt' #  file containing the protein sequence
    ddgdatafile = 'ddG_data.txt' # file containing the literature-culled ddG values
    ddgdatafile_warning = False # warn if ddgdatafile has conflicting ddG values for a mutation
    alignment_file = "uniref_alignment-gaps_lt_0.1-identities_gt_0.5.fasta" # file with aligned sequences
    phylip_path = '/Users/bloom/phylip-3.67/exe/' # path to phylip phylogeny program

    # Define the names of files that will be created by the script if they do not already exist
    cupsatfile = 'CUPSAT_ddGs.txt' # contains the ddG values from CUPSAT
    treefile = "tree.newick" # phylogenetic tree created by phylip
    phylipsequencefile = "phylip_sequence_file" # phylip input sequence file
    phylipdistancefile = "phylip_distance_file" # phylip distance matrix
    pipsddgsfile = "pips_ddgs.txt" # pips ddgs file
    regularizingpriorpipsddgsfile = 'pips_ddgs_with_regularizing_priors.txt' # pips ddgs file calculated with regularizing priors
    hydrophobicitypriorpipsddgsfile = 'pips_ddgs_with_hydrophobicity_priors.txt' # pips ddgs file calculated with hydrophobicity priors

    # Begin execution of the program
    seq = open(seqfile).read().strip() # read in protein sequence

    # Get the ddG values from CUPSAT and store in the dictionary cupsat_ddgs.  Note that
    # in this and all subsequent ddG dictionaries, the first residue is numbered as 0.
    print "\nObtaining CUPSAT ddG values..."
    sys.stdout.flush()
    if os.path.isfile(cupsatfile): # ddG values already obtained, just read from file
        (datetime, cupsat_ddgs) = pips.ddg_inference.ReadDDGs(cupsatfile)
        print "Read the stored CUPSAT values from %s from the file %s." % (datetime, cupsatfile)
    else:  # we need to obtain the ddG values from the CUPSAT webserver
        datetime = time.asctime()
        print "Beginning to calculate and download CUPSAT ddGs at %s..." % datetime
        sys.stdout.flush()
        cupsat_ddgs = pips.cupsat.RunCUPSAT(pdbfile, seq, pdbchain)
        pips.ddg_inference.WriteDDGs(cupsat_ddgs, cupsatfile, datetime)
        print "Completed download of CUPSAT ddG values, stored in the file %s." % cupsatfile
    rescaled_cupsat_ddgs = pips.ddg_inference.RescaleDDGs(cupsat_ddgs, 10.0, '10TH_TO_90TH', recenter=5.0, min_max=(-3.0, 13.0)) 

    # Read the literature-culled ddG data from ddgdatafile and store in the dictionary ddg_data
    print "\nReading the literature-culled ddG data from %s..." % ddgdatafile
    sys.stdout.flush()
    ddgmatch = re.compile("^(?P<wt>[A-Y])(?P<r>\d+)(?P<mut>[A-Y])\s+(?P<ddg>\-{0,1}\d+\.\d+)$")
    ddg_data = {}
    for r in range(len(seq)):
        rdict = {}
        wt = seq[r]
        for aa in pips.ddg_inference.AminoAcids():
            if aa != wt:
                rdict[aa] = []
        ddg_data[r] = (wt, rdict)
    for line in open(ddgdatafile).readlines(): # loop over all lines in ddgdatafile
        if line[0] == '#':
            continue # line is a comment
        m = ddgmatch.search(line.strip()) # match the ddG value
        if not m:
            raise ValueError, "Cannot read ddG value of %s" % line
        (wt, r, mut, ddg) = (m.group('wt'), int(m.group('r')), m.group('mut'), float(m.group('ddg')))
        r -= 1 # we decrement r because we are calling the first residue 0
        if seq[r] != wt:
            raise ValueError, "Wildtype residue does not match protein sequence in %s" % line
        ddg_data[r][1][mut].append(ddg) 
    nddgs = 0
    ddgslist = []
    for (r, (wt, rddgs)) in ddg_data.iteritems():
        for mut in rddgs.iterkeys():
            if not rddgs[mut]:
                rddgs[mut] = None # no ddG value
            else:
                nddgs += 1
                ddg0 = rddgs[mut][0]
                allthesame = True
                for ddgi in rddgs[mut][1 : ]: # see if all ddG values are the same for mutation
                    if ddgi != ddg0:
                        allthesame = False
                if allthesame: # all of the ddG values are the same, take this value
                    rddgs[mut] = ddg0
                    ddgslist.append(ddg0)
                else: # ddG values differ, print warning and take the average value
                    ddg = pips.stats.Mean(rddgs[mut])
                    if ddgdatafile_warning:
                        print "WARNING: Mutation %s%d%s has multiple ddG values of" % (wt, r + 1, mut),
                        for ddgi in rddgs[mut]:
                            print "%.2f" % ddgi,
                        print "--- taking the average value of %.2f." % ddg
                        sys.stdout.flush()
                    rddgs[mut] = ddg
                    ddgslist.append(ddg)
    print "Read a total of %d different ddG values from %s.  The mean value is %.2f, the maximum value is %.2f, and the minimum value is %.2f." % (nddgs, ddgdatafile, pips.stats.Mean(ddgslist), max(ddgslist), min(ddgslist))

    # Read  the aligned sequences (into sequences), give short names for phylip
    sequences = pips.fasta.Read(alignment_file)
    nsequences = len(sequences)
    sequences = [("SEQ%d" % (i + 1), sequences[i][1]) for i in range(nsequences)] # rename 
    pips.fasta.Write(sequences, 'renamed_alignment.fasta')
    sequences = pips.align.StripGapsToFirstSequence(sequences)        
    print "\nThere are %d sequences in the alignment." % nsequences

    # Construct the phylogenetic tree
    if os.path.isfile(treefile):
        print "A phylogenetic tree has already been constructed for these sequences, and is being read from %s." % treefile
        newick_tree = open(treefile).read()
    else:
        print "Constructing a phylogenetic tree for these sequences..."
        sys.stdout.flush()
        pips.phylip.WritePhylipSequenceFile(sequences, phylipsequencefile)
        open(phylipdistancefile, 'w').write(pips.phylip.Protdist(phylipsequencefile, phylip_path))
        newick_tree = pips.phylip.DistanceTree(phylipdistancefile, phylip_path, molecular_clock=True, neighbor_joining=True)
        print "Finished constructing the phylogenetic tree, writing it to %s." % treefile
        sys.stdout.flush()
        open(treefile, 'w').write(newick_tree)

    # Perform the pips analysis
    sequences = pips.fasta.UnknownsToGaps(sequences) # replace unknown amino acids with gaps
    random.seed(1) # seed the random number generator to make output predictable
    (datetime, pips_ddgs) = pips.ddg_inference.ReadDDGs(pipsddgsfile)

    # Read things in with the new pips
    tree = pips.tree.Tree(newick_tree, tipnames_sequences=sequences) # phylogenetic tree data
    ddgset = pips.ddg_inference.DDGSet(seq, tree, ('TRANSITION_TRANSVERSION_RATIO', 0.5), ('SPECIFIED', pips_ddgs, 0, 0), ('BETA', 3, ('KYTE_DOOLITTLE_HYDROPHOBICITY', 1, 0)), 5.0, underflow=5, runtestcode=False)
    ddgset.MaximizePosterior(nrandomstarts=1, printprogress=True)
    new_pips_ddgs = ddgset.DDGDict()
    pips.ddg_inference.WriteDDGs(new_pips_ddgs, 'new_pips_ddgs.txt', time.asctime())

    # Get the consensus ddG
    consensus_ddgs = pips.ddg_inference.ConsensusDDGs(seq, sequences, pseudocounts=1)

    sys.exit()

    # Perform analysis of correlations, and make pylab plots
    print "\nAnalysis of correlations to experimental ddG values..."
    ddgtypes = ['actual', 'CUPSAT', 'consensus', '\\begin{tabular}{c} PIPS with \\\\ informative prior \end{tabular}', '\\begin{tabular}{c} PIPS with \\\\ regularizing prior \end{tabular}', '\\begin{tabular}{c} PIPS with \\\\ hydrophobicity prior \end{tabular}']
    zippedlists = pips.ddg_inference.ZippedDDGLists(ddg_data, cupsat_ddgs, consensus_ddgs, pips_ddgs, pips_ddgs_regularizing, pips_ddgs_hydrophobicity)
    mutations = zippedlists[0]
    nmutations = len(mutations)
    ddgs = dict([(ddgtypes[i], zippedlists[i + 1]) for i in range(len(ddgtypes))])
    pylab.rc('text', usetex=True)
    nplots = len(ddgtypes) - 1 # number of different plots
    invnplots = 1.0 / nplots
    (xscale, yscale) = (2.8, 2.5) # each plot covers a rectangle of this size, in inches
    bottom = 1.06
    (tmargin, bmargin, lmargin, rmargin) = (0.03, 0, 0.22, 0.03)
    fig = pylab.figure(figsize=(xscale * (1 + lmargin + rmargin), 3 * yscale * (1 + tmargin + bmargin) * bottom))
    figaxes = pylab.axes([0, 0, 1, 1])
    figaxes.axison = False
    iplot = 0
    maxticks = 5
    (xmin, xmax) = (int(round(min(ddgs['actual'])) - 1), int(round(max(ddgs['actual'])) + 1))
    xtick = 1
    while (xmax - xmin) / float(xtick) > maxticks:
        xtick += 1
    nxticks = int(math.ceil((xmax - xmin) / float(xtick)))
    xticks = [x for x in range(xmin, xmin + nxticks * xtick + 1, xtick)]
    xticklocator = matplotlib.ticker.FixedLocator(xticks)
    xtickformatter = matplotlib.ticker.FixedFormatter(["%d" % x for x in xticks])
    for ddgtype in ddgtypes[1 : ]:
        if ddgtype == ddgtypes[-1]:
            xlabel = 'experimental $\Delta\Delta G$ values'
        else:
            xlabel = ''
        (r, p, npoints) = pips.stats.PearsonCorrelation(ddgs['actual'], ddgs[ddgtype])
        axes = pylab.axes([lmargin, 1.0 - invnplots * (1 + iplot + bmargin) / bottom, 1.0 - rmargin - lmargin, invnplots * (1.0 - tmargin - bmargin) / bottom], xlabel=xlabel, ylabel=ddgtype)
        nolabels = matplotlib.ticker.NullFormatter()
        (ymin, ymax) = (int(round(min(ddgs[ddgtype])) - 1), int(round(max(ddgs[ddgtype])) + 1))
        ytick = 1
        while (ymax - ymin) / float(ytick) > maxticks:
            ytick += 1
        nyticks = int(math.ceil((ymax - ymin) / float(ytick)))
        yticks = [y for y in range(ymin, ymin + nyticks * ytick + 1, ytick)]
        yticklocator = matplotlib.ticker.FixedLocator(yticks)
        ytickformatter = matplotlib.ticker.FixedFormatter(["%d" % y for y in yticks])
        axes.xaxis.set_major_locator(xticklocator)
        axes.yaxis.set_major_locator(yticklocator)
        axes.yaxis.set_major_formatter(ytickformatter)
        if ddgtype != ddgtypes[-1]:
            axes.xaxis.set_major_formatter(nolabels)
        else:
            axes.xaxis.set_major_formatter(xtickformatter)
        iplot += 1
        pylab.text(0.64, 0.14, '$R^2 = %.2f$' % r**2, transform=axes.transAxes, ha='left', va='top', size=14)
        pylab.scatter(ddgs['actual'], ddgs[ddgtype], figure=fig, axes=axes)
    pylab.savefig("%s_vertical_plot.eps" % protein_w_underscores)

    pylab.show()


# run the script
main()
