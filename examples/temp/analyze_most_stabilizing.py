import pips.ddg_inference
import pips.stats

nlisted = 20
ddgs = pips.ddg_inference.ReadDDGs('PIPS_PREDICTIONS_PRIORS-REGULARIZING.TXT')[1]
ddgs_list = pips.ddg_inference.DDGsToList(ddgs)
sorted_ddgs = pips.ddg_inference.SortedDDGList(ddgs)
print "The mean and median for the ddG predictions are %.2f and %.2f" % (pips.stats.Mean(ddgs_list), pips.stats.Median(ddgs_list))
print "\nThe most stabilizing %d predicted mutations are:" % nlisted
for (ddg, mut) in sorted_ddgs[ : nlisted]:
    print mut, ddg
# oneper_sorted_ddgs lists only the most stabilizing predicted mutation at each site
oneper_sorted_ddgs = []
res_listed = {}
for (ddg, mut) in sorted_ddgs:
    (wt, r, m) = (mut[0], int(mut[1 : -1]), mut[-1])
    if r not in res_listed:
        oneper_sorted_ddgs.append((ddg, mut))
        res_listed[r] = True
print "\nListing only the best mutation at each residue, the most stabilizing %d predicted mutations are:" % nlisted
for (ddg, mut) in oneper_sorted_ddgs[ : nlisted]:
    print mut, ddg
# The N1 ectodomain crystal structure 3BEQ begins at residue 83 in the N2 numbering scheme, which is residue 82 of the CA09 NA.  The list ecto_oneper_sorted_ddgs lists the best mutation at each residue restricted to these ectodomain residues
ecto_oneper_sorted_ddgs = [(ddg, mut) for (ddg, mut) in oneper_sorted_ddgs if int(mut[1 : -1]) >= 82]
print "\nListing only crystallized ectodomain mutations, only the best mutation at each residue, the most stabilizing %d predicted mutations are:" % nlisted
for (ddg, mut) in ecto_oneper_sorted_ddgs[ : nlisted]:
    print mut, ddg
