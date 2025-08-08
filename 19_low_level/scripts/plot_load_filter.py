import matplotlib.pyplot as plt
from matplotlib import cm
from commander_tools.tod_tools.lfi import lfi
import numpy as np

for freq in lfi.freqs:

    plt.figure()

    increment = 1.0/(2*len(lfi.horns[freq]))
    color = 0

    for horn in lfi.horns[freq]:
        for hornType in lfi.hornTypes:
            realcolor = color
            ls = 'solid'
            for diode in lfi.diodeTypes[hornType]:

                filterFile = np.transpose(np.loadtxt('./load_filter_' + str(horn) + hornType + '_ref' + diode + '.dat'))

                x = filterFile[0][0:-1]
                transfer = filterFile[1][0:-1]


                plt.plot(x, transfer, label=str(horn)+hornType + '-' + diode, c=cm.tab20(realcolor), linestyle=ls)
                ls = 'dotted'
            color += increment

    col =1
    if freq == 70:
        col=2

#    plt.legend(loc='upper left', bbox_to_anchor=(1, 1.018), ncol=col)
    plt.legend(loc='best', ncol=col)
    plt.title("Load filter functions - " + str(freq) + 'GHz')
    plt.xscale('log')

    plt.xlabel('Frequency (Hz)') 
    plt.ylabel('Transfer')
    plt.tight_layout()

    plt.savefig('load_filters_' + str(freq) + '.pdf')  
 
