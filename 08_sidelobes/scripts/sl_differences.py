import os
import subprocess as sp

chains_dir = '/mn/stornext/u3/hke/BeyondPlanck/data/delivery/v6.00/chains_BP6_c1/'

first_chain = 1
last_chain = 95
burn_in = 25

pairs = [[last_chain, first_chain], [last_chain, burn_in], [burn_in, first_chain]]

freqs = [30, 44, 70]

map_editor = '/mn/stornext/d16/cmbco/bp/mathew/beyondplank/cmbco/src/f90/map_editor/map_editor' 

for pair in pairs:
    for freq in freqs:
        first_map = chains_dir + 'tod_0' + str(freq) + '_sl_c0001_k' + str(pair[0]).zfill(6) + '.fits'
        second_map = chains_dir + 'tod_0' + str(freq) + '_sl_c0001_k' + str(pair[1]).zfill(6)  + '.fits'
        out_map = './diff_map_' + str(freq) + '_' + str(pair[0]) + '_minus_' + str(pair[1]) + '.fits'
        sp.run([map_editor, 'subtract', first_map, second_map, out_map])
        sp.run(['c3pp', 'plot', out_map, '-bar', '-IQU'])
