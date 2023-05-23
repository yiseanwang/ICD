import sys
import math
import numpy as np
from numpy import linalg as LA
import glob
import os



np.printoptions(precision=10)
np.set_printoptions(threshold=sys.maxsize)

# os.system("python2 2_mov2asc_py2.py   H2O_scf.movecs > H2O_scf.Cmo")
# os.system("python2 2_mov2asc_py2.py H2O_noscf.movecs > H2O_noscf.Cmo")

# scf_Cmo=glob.glob("scf.Cmo")

#######################################################################################
def extract_data(filelines, tag_head, tag_end, the_num):
        val=[]
        tag1_list=[]
        tag2_list=[]
        p_head=0
        p_end=0
        for l, name in enumerate(filelines):
                if ( tag_head in name ):
                        tag1_list.append(l)
			#p_head=l
                if ( tag_end in name ):
                        tag2_list.append(l)
			#p_end=l+1
        p_head=int(tag1_list[-1])+4
        p_end=tag2_list[-1]
        print(tag1_list, p_head)
        print(tag2_list, p_end)
        p= p_head
        while p < p_end :
                S_colums=S_line[p].strip("\n").split()
                col=1
                while col < len(S_colums):
                        row=0
                        while row < the_num:
                                S_element=S_line[p+row].strip("\n").split()
                                print(S_element[col])
                                val.append(float(S_element[col]))
                                row=row+1
                        col=col+1
                p=p+(3+the_num)
        S_mat=np.array(val).reshape(the_num,the_num)
        return S_mat
########################################################################################


## 建構SCF matrix ## 
##############################################################################
def mat_scf(the_C_scf, nmo, nbf):
	scf_file=open(the_C_scf,"r")
	scf_lines=scf_file.readlines()
	scf_file.close()
	scf_list=[]
	line=0
	while line < len(scf_lines):
		element=scf_lines[line].strip('\n').split()
		i=0
		while i < len(element):
			scf_list.append(float(element[i].strip("\n")))
			i=i+1
		line=line+1
	mo_scf=np.array(scf_list,dtype=float).reshape(nmo,nbf)   # mo_scf
	mo_scf_trans=np.transpose(mo_scf)
	return mo_scf, mo_scf_trans
#############################################################################

## 建構noSCF matrix ##
#############################################################################
def mat_noscf(the_C_noscf, nmo, nbf):
	noscf_file=open(the_C_noscf,"r")
	noscf_lines=noscf_file.readlines()
	noscf_file.close()
	noscf_list=[]
	line=0
	while line < len(noscf_lines):
		element=noscf_lines[line].strip('\n').split()
		i=0
		while i < len(element):
			noscf_list.append(float(element[i].strip("\n")))
			i=i+1
		line=line+1
	mo_noscf=np.array(noscf_list,dtype=float).reshape(nmo,nbf)  # mo_noscf
	mo_noscf_trans=np.transpose(mo_noscf)
	return mo_noscf, mo_noscf_trans
#############################################################################

## 建構moocc matrix ##
##############################################################################
def moocc_mat(the_Cmat_scf, the_Cmat_scf_trans, the_Cmat_noscf, the_Cmat_noscf_trans, the_moocc, Smatrix):
	datafile = open(filename+'_noscf.occ', 'w+')
	moocc_file=open(the_moocc, 'r')
	lines=moocc_file.readlines()
	moocc_file.close()
	t=0
	while t < len(lines):
		scf_moocc0=lines[t].split()
#		print(type(scf_moocc0[3]))
		del scf_moocc0[0]
		scf_moocc=[float(x) for x in scf_moocc0]
		Pmo_scf=np.array(np.diag(scf_moocc), dtype=float)
#		print(the_Cmat_scf)
		Pao=np.dot(the_Cmat_scf_trans, np.dot(Pmo_scf, the_Cmat_scf))
		Pao_S=np.dot(Smatrix, np.dot(Pao, Smatrix))
		Pmo_noscf=np.dot(the_Cmat_noscf, np.dot(Pao_S, the_Cmat_noscf_trans))
		moocc_noscf=Pmo_noscf.diagonal()
		datafile.write(str(t) +'\t')
		for occ in moocc_noscf :
			datafile.write(str(occ) +'\t')
		datafile.write('\n')
		t=t+1
	return 
################################################################################


## THE REAL GAME !!!! ##
#tag1="Superposition of Atomic Density Guess"
tag1="global array: Temp Over"
tag2="Dispersion Parameters" #"82       0.46743     0.00000     0.00000     1.00000" #"Non-variational initial energy"
nbf=82
nmo=82
S_file=open("H2OH2O-2.895_CAM_aDZ-aDZ_noscf.nwo","r")
S_line=S_file.readlines()
S_file.close()
S_matrix=extract_data(S_line, tag1, tag2, nbf)
#################################################################################
## main ##
## 1. Pao=Cscf.Pmo.Cscf
## 2. P_mo_noscf=Cnoscf.S.Pao.S.Cnoscf
moocc_files=glob.glob("*.moocc")
# moocc_file=open("*.moocc","r")
for ifile in moocc_files :
	scf_Cmo_name="H2OH2O-2.895_CAM_aDZ-aDZ_scf.Cmo"
	noscf_Cmo_name="H2OH2O-2.895_CAM_aDZ-aDZ_noscf.Cmo"
	filename=ifile.strip(".moocc")
	moocc=open(ifile,"r")
	C_scf=scf_Cmo_name   #str(scf_Cmo_name)+"_scf.Cmo"
	C_noscf=noscf_Cmo_name #str(noscf_Cmo_name)+"_noscf.Cmo"
	Cmat_scf, Cmat_scf_trans = mat_scf(C_scf,nbf,nmo)		#construct   SCF-DM, need   SCF.Cmo
	Cmat_noscf, Cmat_noscf_trans = mat_noscf(C_noscf,nbf,nmo)	#construct noSCF-DM, need nosCF.Cmo
	Cmat_scf=np.array(Cmat_scf,dtype=float)				#
	Cmat_scf_trans=np.array(Cmat_scf_trans,dtype=float)		#
	Cmat_noscf=np.array(Cmat_noscf,dtype=float)			#
	Cmat_noscf_trans=np.array(Cmat_noscf_trans,dtype=float)		#
#	S_matrix=extract_data(S_line, tag, BF_num)
#	print(Cmat_scf)
	moocc_mat(Cmat_scf, Cmat_scf_trans, Cmat_noscf, Cmat_noscf_trans, ifile, S_matrix) 

