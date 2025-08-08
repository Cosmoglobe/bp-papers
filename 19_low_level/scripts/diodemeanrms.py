from commander_tools.tod_tools import commander_tod as cmdtod
from commander_tools.tod_tools.lfi import lfi
import matplotlib.pyplot as plt
import numpy as np

tod = cmdtod.commander_tod('/mn/stornext/d16/cmbco/bp/mathew/test/undifferenced/multiHuff', 'LFI')

for freq in [44]:

    means = {}
    stds = {}

    for od in lfi.ods:

        if(od == 1007):
            continue

        tod.init_file(freq, od)
        print(freq, od)

        for horn in lfi.horns[freq]:
            for hornType in lfi.hornTypes:
            
                for pid in tod.load_field('/common/pids'):
                   
                    print('/' + str(pid).zfill(6) +'/' + str(horn) + hornType + '/diodes') 
                    diodes = tod.load_field('/' + str(pid).zfill(6) +'/' + str(horn) + hornType + '/diodes')
                
                    i = 0
                    for diode in lfi.diodeTypes[hornType]:
                        if pid == 1:
                            means[str(horn) + hornType + '-ref' + diode] = []
                            stds[str(horn) + hornType + '-ref' + diode] = []
                            means[str(horn) + hornType + '-sky' + diode] = []
                            stds[str(horn) + hornType + '-sky' + diode] = []


                        means[str(horn) + hornType + '-ref' + diode].append(np.mean(diodes[i]))
                        stds[str(horn) + hornType + '-ref' + diode].append(np.std(diodes[i]))
                        i += 1
                        
                        means[str(horn) + hornType + '-sky' + diode].append(np.mean(diodes[i]))
                        stds[str(horn) + hornType + '-sky' + diode].append(np.std(diodes[i]))
                        i += 1

    for horn in lfi.horns[freq]:
        for hornType in lfi.hornTypes:
            for diode in lfi.diodeTypes[hornType]:
                for data, name in zip([stds, means], ['std', 'mean']):
                    np.savetxt(str(horn) + hornType + '-ref' + diode + '_' + name + '.txt', data[str(horn) + hornType + '-ref' + diode])
                    np.savetxt(str(horn) + hornType + '-sky' + diode + '_' + name + '.txt', data[str(horn) + hornType + '-sky' + diode])
    
