#!/usr/bin/env python3

from ppp_plot import *

odir = '.'
uK = r'$\mathrm{\mu K}$'

clist = [#('300uK', -300.0, 300.0), 
         #('200uK', -200.0, 200.0),
         #('15uK',   -15.0,  15.0),
         #('10uK',   -10.0,  10.0),
         #('4uK',     -4.0,   4.0),
         #('2p5uK',   -2.5,   2.5)]
         ('1uK',   -1,   1)]

for n, vmin, vmax in clist:
    print(n)
    ofile = '{d}/colourbar_{n}.pdf'.format(d=odir, n=n)
    colourbar(vmin, vmax, ofile, units=uK)
