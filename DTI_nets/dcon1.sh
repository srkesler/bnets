#!/bin/bash

#Shelli Kesler 12/8/14 

#This script performs DTI correction and brain masking.  *Check the masks for holes before proceeding to dcon2.

#STUDY directory contains all subject folders and no other subdirectories.

#This script must be inside the STUDY directory

STUDY='/path/' #path to subject folders

SUBJS=`ls -d */`
datadir='connectome/' #name of folder containing DTI and anatomic data

for subjdir in ${SUBJS}
do 
	currentdir=${STUDY}${subjdir}${datadir}
	cd ${currentdir}
	
	#eddy current correction of raw DTI data
	eddy_correct data.nii.gz dti_corr 0

	#extract the first B0 volume
	fslroi dti_corr.nii.gz nodif 0 1
	
	#skull strip the B0 volume and create a binary mask of it
	bet nodif.nii.gz nodif_brain -f 0.5 -g 0 -m
done
