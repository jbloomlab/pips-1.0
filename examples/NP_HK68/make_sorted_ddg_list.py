"""Makes sorted list of ddG values."""

infile = 'PIPS_PREDICTIONS_PRIORS-REGULARIZING.TXT'
ddgs = open(infile).readlines()[2 : ]
ddgs = [(float(line.split()[1]), line.split()[0]) for line in ddgs]
ddgs.sort()
f = open('SORTED_%s' % infile, 'w')
for (ddg, m) in ddgs:
    f.write("%s\t%.3f\n" % (m, ddg))
f.close()
