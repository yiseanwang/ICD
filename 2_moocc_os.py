import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import os

outfile = glob.glob('*.reo')
occ_init=3
occ_final=40
os.system("rm -r moocc")
########################################
tag1="# MO Occupations (alpha spin)"
tag2="# MO Occupations (beta spin)"

for ifile in outfile :
	filename=ifile.strip(".reo")
	f=open(ifile)
	lines=f.readlines()
	f.close()
	datafile = open(filename+'.moa', 'w+')
	for l in lines:
		if (tag1 in l ):
			alpha = l.strip().split()
			val1 = alpha[1:-5]
			for occ in val1 :
				datafile.write(occ +'\t')
			datafile.write('\n')
	datafile.close()
	datafile = open(filename+'.mob', 'w+')
	for l in lines:
        	if (tag2 in l ):
                	alpha = l.strip().split()
                	val1 = alpha[1:-5]
                	for occ in val1 :
                        	datafile.write(occ +'\t')
                	datafile.write('\n')
	datafile.close()
	a=np.loadtxt(filename+'.moa')
	b=np.loadtxt(filename+'.mob')
	c=a+b
	np.savetxt(filename+'.moocc', c , fmt="%f")

occfiles = glob.glob('*.moocc')
for ifile in occfiles :
	filename=ifile.strip(".moocc")
	occf=pd.read_csv(ifile ,header=None,delimiter=' ')
	f=pd.DataFrame(occf)
	f2=f.drop(columns=0)
	f_plot=f2.loc[0:,occ_init:occ_final]
	f_plot.plot(figsize=(20, 9))
	plt.title(filename,fontsize = 40)
	plt.savefig(filename+'.png',bbox_inches='tight')

os.system("mkdir moocc")
os.system("mv *.png moocc")
