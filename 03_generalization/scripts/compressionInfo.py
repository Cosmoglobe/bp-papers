import h5py
import sys

indir = '/mn/stornext/d16/cmbco/bp/data/L2Data'
outdir = '/mn/stornext/d16/cmbco/bp/mathew/test/compress_tod'

in_files = {30:['LFI_030_27_L2_002_OD1000.h5', 'LFI_030_28_L2_002_OD1000.h5'], 44:['LFI_044_24_L2_002_OD1000.h5', 'LFI_044_25_L2_002_OD1000.h5', 'LFI_044_26_L2_002_OD1000.h5'], 70:['LFI_070_18_L2_002_OD1000.h5', 'LFI_070_19_L2_002_OD1000.h5', 'LFI_070_20_L2_002_OD1000.h5', 'LFI_070_21_L2_002_OD1000.h5', 'LFI_070_22_L2_002_OD1000.h5', 'LFI_070_23_L2_002_OD1000.h5']}
out_files = {30:'LFI_030_001000.h5', 44:'LFI_044_001000.h5', 70:'LFI_070_001000.h5'}

things = {'/tod':['SIGNAL'], '/pix':['PHI', 'THETA'], '/psi':['PSI'], '/scalars':['garbage'], '/flag':['FLAG'], '/time':['MJD', 'OBT', 'SCET'], '/huffsymb':['garbage'], '/hufftree':['garbage']}

freqs = [30, 44, 70]

for freq in freqs:
    outFile = h5py.File(outdir + '/' + out_files[freq])

    inFiles = []
    for in_file in in_files[freq]:
        inFiles.append(h5py.File(indir + '/' + in_file))

    for thing in things:
        outsize = 0
        #print(outFile.keys())
        for dataset in outFile.keys():
            #print(outFile[dataset].keys())
            for ds in outFile[dataset].keys():
                if('common' in dataset):
                    if 'scalars' in thing:
                        outsize += outFile[dataset + '/' + ds].size * outFile[dataset + '/' + ds].dtype.itemsize
                    continue
                for sds in outFile[dataset + '/' + ds]:
                    if thing in ('/' + sds):
                        outsize += outFile[dataset + '/' + ds + '/' + sds].size* outFile[dataset + '/' + ds + '/' + sds].dtype.itemsize

        insize = 0
        for infile in inFiles:
            for dataset in infile.keys():
                for ds in infile[dataset].keys():
                    #print(dataset, ds, things[thing], infile[dataset + '/' + ds].size, any(x in ds for x in things[thing]), (x for x in things[thing]))

                    for name in things[thing]:
                        if name in ds:
                            insize += infile[dataset + '/' + ds].size * infile[dataset + '/' + ds].dtype.itemsize

        print(freq, thing, insize/1000000, outsize/1000000)
        
