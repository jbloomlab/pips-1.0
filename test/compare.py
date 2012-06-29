python = [float(x) for x in open('p_python.out').read().split(',')]
c = [float(x) for x in open('p_c.out').read().split(',')]
for (x, y) in zip(python, c):
    xform = "%.5e" % x
    yform = "%.5e" % y
    if xform != yform:
        print x, y
