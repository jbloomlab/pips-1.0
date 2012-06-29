import pips.ddg_inference
import pips.stats


regularizing = pips.ddg_inference.ReadDDGs('PIPS_PREDICTIONS_PRIORS-REGULARIZING.TXT')[1]
regularizing_list = pips.ddg_inference.DDGsToList(regularizing)
sorted_regularizing = pips.ddg_inference.SortedDDGList(regularizing)
hydrophobicity = pips.ddg_inference.ReadDDGs('PIPS_PREDICTIONS_PRIORS-HYDROPHOBICITY.TXT')[1]
hydrophobicity_list = pips.ddg_inference.DDGsToList(hydrophobicity)
sorted_hydrophobicity = pips.ddg_inference.SortedDDGList(hydrophobicity)
print "The mean and median for the regularizing priors are %.2f and %.2f" % (pips.stats.Mean(regularizing_list), pips.stats.Median(regularizing_list))
print "The mean and median for the hydrophobicity priors are %.2f and %.2f" % (pips.stats.Mean(hydrophobicity_list), pips.stats.Median(hydrophobicity_list))
print "\nhydrophobicity priors top 30:"
for (mut, ddg) in sorted_hydrophobicity[ : 30]:
    print mut, ddg
print "\nregularizing priors top 30:"
for (mut, ddg) in sorted_regularizing[ : 30]:
    print mut, ddg
