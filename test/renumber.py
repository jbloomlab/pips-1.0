import pips.ddg_inference

(datetime, ddgs) = pips.ddg_inference.ReadDDGs('CUPSAT_ddGs.txt')
renumbered_ddgs = {}
for (r, (rwt, rddgs)) in ddgs.iteritems():
    renumbered_ddgs[r + 1] = (rwt, rddgs)
ddgs = pips.ddg_inference.WriteDDGs(renumbered_ddgs, 'ddg_priors.txt', datetime)
