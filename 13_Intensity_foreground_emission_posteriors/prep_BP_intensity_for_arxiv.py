import os
import numpy as np

fin='13_Intensity_foreground_emission_posteriors.tex'

data=open(fin, 'r')

if (not os.path.isdir('arxiv/')):
    os.system('mkdir arxiv')
if (not os.path.isdir('arxiv/figs')):
    os.system('mkdir arxiv/figs')

outfile='arxiv/BP_XIV_Intensity_fg.tex'
if os.path.isfile(outfile) and (not outfile == fin):
    os.system('rm '+outfile)
outdata=open(outfile,'w')

for line in data:
    if '\includegraphics' in line:
        newline=False
        aline=line.split('{')
        if not len(aline[1]) < 1:
            line2=aline[1]
            aline=line2.split('}')
            if not len(aline[0]) < 1:
                fpath=aline[0]
                aline=fpath.split('/')
                f=aline[-1]
                if '.pdf' in f or '.png' in f:
                    #print(f) #debug

                    cpline='rsync -Ptr '+fpath+' arxiv/figs/.'
                    #print(cpline)
                    os.system(cpline)
                    newline=True
                    #print(templine)
        if (newline):
            templine=line.replace(fpath,'figs/'+f)
        else:
            templine=line
    else:
        templine=line

    if r'\input{../common/BP_acknowledgments.tex}' in line:
        templine=r'\input{BP_acknowledgments.tex}'
        cpline='rsync -Ptr ../common/BP_acknowledgments.tex arxiv/.'
        os.system(cpline)
    if r'\bibliography{../common/Planck_bib,../common/BP_bibliography}' in line:
        templine=r'\bibliography{Planck_bib,BP_bibliography}'
        cpline='rsync -Ptr ../common/Planck_bib.bib arxiv/.'
        os.system(cpline)
        cpline='rsync -Ptr ../common/BP_bibliography.bib arxiv/.'
        os.system(cpline)
        
    outdata.write(templine)

#copy default files
cpline='rsync -Ptr aa.cls arxiv/.'
os.system(cpline)
cpline='rsync -Ptr aa.bst arxiv/.'
os.system(cpline)
cpline='rsync -Ptr Planck.tex arxiv/.'
os.system(cpline)
cpline='rsync -Ptr BP13_authors.tex arxiv/.'
os.system(cpline)
cpline='rsync -Ptr Makefile_arxiv arxiv/Makefile'
os.system(cpline)
