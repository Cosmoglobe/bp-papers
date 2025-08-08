import os
import numpy as np

fin='13_Intensity_foreground_emission_posteriors.tex'

data=open(fin, 'r')

if (not os.path.isdir('AandA/')):
    os.system('mkdir AandA')
if (not os.path.isdir('AandA/figs')):
    os.system('mkdir AandA/figs')

outfile='AandA/BP_XIV_Intensity_fg.tex'
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

                    cpline='rsync -Ptr '+fpath+' AandA/figs/.'
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
        cpline='rsync -Ptr ../common/BP_acknowledgments.tex AandA/.'
        os.system(cpline)
    if r'\bibliography{../common/Planck_bib,../common/BP_bibliography}' in line:
        templine=r'\bibliography{Planck_bib,BP_bibliography}'
        cpline='rsync -Ptr ../common/Planck_bib.bib AandA/.'
        os.system(cpline)
        cpline='rsync -Ptr ../common/BP_bibliography.bib AandA/.'
        os.system(cpline)
        
    outdata.write(templine)

#copy default files
cpline='rsync -Ptr aa.cls AandA/.'
os.system(cpline)
cpline='rsync -Ptr aa.bst AandA/.'
os.system(cpline)
cpline='rsync -Ptr Planck.tex AandA/.'
os.system(cpline)
cpline='rsync -Ptr BP13_authors.tex AandA/.'
os.system(cpline)
cpline='rsync -Ptr Makefile_AandA AandA/Makefile'
os.system(cpline)

outdata.close()

os.system('cd AandA/')
os.system('make')
os.system('cd ../')


#convert to 1-column reviewer version
fin='AandA/BP_XIV_Intensity_fg.tex'

data=open(fin, 'r')

if (not os.path.isdir('AandA/')):
    os.system('mkdir AandA')
if (not os.path.isdir('AandA/figs')):
    os.system('mkdir AandA/figs')

outfile='AandA/BP_XIV_Intensity_fg_onecolumn.tex'
if os.path.isfile(outfile) and (not outfile == fin):
    os.system('rm '+outfile)
outdata=open(outfile,'w')

inside_fig=False
for line in data:
    #check if we are inside a figure that normally is 1-column figure
    if inside_fig:
        #check if we exit figure
        if r'\end{figure}' in line:
            inside_fig=False
            outline=line
        elif r'\includegraphics' in line: #check if \includegraphics is called
            if r'\includegraphics[width=\linewidth]' in line:
                outline=line.replace(r'\linewidth',r'0.5\linewidth')
            elif r'\linewidth' in line:
                aline=line.split(r'width=')
                temp=aline[1]
                aline2=temp.split('\linewidth')
                temp2=aline2[0]
                val=float(temp2)
                print(val,'%.3f'%(val),'%.3f'%(val/2.0))
                outline=line.replace('width='+temp2,'width='+'%.3f'%(val/2.0))
            elif r'\includegraphics[width=\textwidth]' in line:
                outline=line.replace(r'\textwidth',r'0.5\textwidth')
            elif r'\textwidth' in line:
                aline=line.split(r'width=')
                temp=aline[1]
                aline2=temp.split('\textwidth')
                temp2=aline2[0]
                val=float(temp2)
                #print(val,'%.3f'%(val),'%.3f'%(val/2.0)) #debug
                outline=line.replace('width='+temp2,'width='+'%.3f'%(val/2.0))
            else:
                outline=line
        else:
            outline=line
        #print(outline) #debug
        
    else:
        #check if we enter a 1-column figure
        if r'\begin{figure}' in line:
            inside_fig=True
            outline=line
        else:
            if r'%\modulolinenumbers[5]' in line:
                outline=r'\modulolinenumbers[5]'+'\n'
            elif r'%\linenumbers' in line:
                outline=r'\linenumbers'+'\n'
            elif r'\documentclass[twocolumn]{aa}' in line:
                outline=r'\documentclass[onecolumn]{aa}'+'\n'
            else:
                outline=line
        
    outdata.write(outline)

cpline='rsync -Ptr Makefile_AandA_onecolumn AandA/Makefile'
os.system(cpline)

outdata.close()

os.system('cd AandA/')
os.system('make')
os.system('cd ../')
